"""Prediction, consensus selection and carbon propagation stages."""

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from copy import deepcopy
from datetime import datetime
from hashlib import sha256
import json
from multiprocessing import get_context
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook, Workbook
from rdkit import Chem, RDLogger

from .atom_helpers import mapped_components
from .carbon_trace import propagate_rounds, initialize_sources
from .flow_output import FlowOutput, WorksheetRowLimitError
from .config import configure, parallel_capacity, MemoryGuard, MemoryPressureError
from .excel_io import (
    clear_formatting,
    copy_input,
    input_sheets,
    write_sheet,
    save_output,
    workbook_limit,
    WorkbookSizeLimitError,
)
from .carbon_mapping_rules import apply_rules, apply_manual_mappings, load_partials, carbon_structure_issues, carbon_presence
from .consensus import select_mappings
from .consensus import mapping_signatures, carbon_candidate
from .mapper_runner import run_mapper, validate_mapper_config
from .reactions import clean, get_smiles, prepare_reactions, parse_side, split_equation
from .smiles_completion import complete_smiles


def carbon_metabolites(metabolites, reactions):
    """Count carbon metabolite IDs participating in usable final mappings."""
    mapped_names = {
        component.metabolite
        for reaction in reactions
        if reaction.status == "mapped" and reaction.mapped_smiles
        for component in reaction.left + reaction.right
    }
    for reaction in reactions:
        if reaction.status == "partial" and reaction.mapped_smiles:
            for side in mapped_components(reaction) or ():
                for component, maps in side:
                    molecule = Chem.MolFromSmiles(component.smiles)
                    if any(number and molecule.GetAtomWithIdx(index).GetAtomicNum() == 6
                           for index, number in maps.items()):
                        mapped_names.add(component.metabolite)
    mapped_names.update(
        row[field] for reaction in reactions for row in reaction.partial_maps
        for field in ("SourceMetabolite", "ProductMetabolite")
    )
    carbon_ids, unknown_ids = defaultdict(set), set()
    for row in metabolites.drop_duplicates("MetaboliteID").itertuples(index=False):
        smiles = get_smiles(row)
        molecule = Chem.MolFromSmiles(smiles) if smiles else None
        if molecule is None:
            unknown_ids.add(clean(row.MetaboliteID))
        elif clean(row.MetaboliteName) in mapped_names and any(
            atom.GetAtomicNum() == 6 for atom in molecule.GetAtoms()
        ):
            carbon_ids[clean(row.MetaboliteName)].add(clean(row.MetaboliteID))
    return carbon_ids, len(unknown_ids)


def read_input(
    path, reaction_sheet, metabolite_sheet, structure_options=None, model=""
):
    with pd.ExcelFile(path, engine="openpyxl") as excel:
        reactions_df = pd.read_excel(excel, sheet_name=reaction_sheet)
        metabolites = pd.read_excel(excel, sheet_name=metabolite_sheet)
    required = {
        "reactions": {"ReactionID", "ReactionName", "Equation", "Flux"},
        "metabolites": {"MetaboliteID", "MetaboliteName"},
    }
    for name, table in (("reactions", reactions_df), ("metabolites", metabolites)):
        missing = required[name] - set(table.columns)
        if missing:
            raise ValueError(f"Missing {name} columns: {sorted(missing)}")
        identity = "ReactionID" if name == "reactions" else "MetaboliteID"
        if (
            table[identity].map(clean).eq("").any()
            or table[identity].duplicated().any()
        ):
            raise ValueError(f"{identity} must be nonempty and unique.")
    if (
        reactions_df.empty
        or not pd.to_numeric(reactions_df.Flux, errors="coerce")
        .map(lambda value: pd.notna(value) and abs(value) < float("inf") and value != 0)
        .all()
    ):
        raise ValueError("Reaction inputs must contain finite nonzero fluxes.")
    if "SMILES_Manual" not in metabolites:
        metabolites["SMILES_Manual"] = float("nan")
    if not {"SMILES_Model", "SMILES_Yeast9"}.intersection(metabolites.columns):
        metabolites["SMILES_Model"] = float("nan")
    if structure_options is not None:
        metabolites = complete_smiles(metabolites, model, structure_options)
    structures = [get_smiles(row) for row in metabolites.itertuples(index=False)]
    if (
        pd.DataFrame({"name": metabolites.MetaboliteName, "smiles": structures})
        .groupby("name")
        .smiles.nunique()
        .gt(1)
        .any()
    ):
        raise ValueError(
            "A compartment-qualified metabolite name has conflicting structures."
        )
    reactions = prepare_reactions(reactions_df, metabolites)
    return reactions, reactions_df, metabolites


def reuse_predictions(path, backend, settings, reactions, source, metabolites, job):
    """Reuse predictions only when inputs, direction and adapter settings match."""
    if not path.is_file():
        return False
    try:
        with pd.ExcelFile(path, engine="openpyxl") as excel:
            old_source = pd.read_excel(excel, sheet_name=job["reaction_sheet"])
            old_metabolites = pd.read_excel(excel, sheet_name=job["metabolite_sheet"])
            metadata = pd.read_excel(
                excel, sheet_name="Mapping_Info", dtype={"ReactionID": str}
            )
            info = metadata.loc[metadata.Mapper == backend]
            if info.empty:
                return False
            expected = dict(
                Mapper=backend, Adapter=settings["module"], Python=settings["python"]
            )
            if any(
                clean(info.iloc[0].get(key)) != clean(value)
                for key, value in expected.items()
            ):
                return False
            if "OptionsJSON" in info:
                if json.loads(clean(info.iloc[0].OptionsJSON)) != settings.get(
                    "options", {}
                ):
                    return False
            elif any(
                clean(info.iloc[0].get(key)) != clean(value)
                for key, value in settings.get("options", {}).items()
            ):
                return False

            pd.testing.assert_frame_equal(
                old_source[list(source.columns)],
                source,
                check_exact=False,
                rtol=1e-12,
                atol=1e-15,
            )
            if not old_metabolites.equals(metabolites):
                return False
            sheet = backend if backend in excel.sheet_names else job["reaction_sheet"]
            table = pd.read_excel(excel, sheet_name=sheet, dtype={"ReactionID": str})
            if "MappingStatus" not in table:
                details = metadata.loc[
                    (metadata.Mapper == backend) & metadata.ReactionID.notna(),
                    ["ReactionID", "MappingStatus", "MappingTemplate"],
                ]
                table = table.merge(
                    details, on="ReactionID", how="left", validate="one_to_one"
                )
                if table.MappingStatus.isna().any():
                    return False
            if table.ReactionID.tolist() != [r.reaction_id for r in reactions]:
                return False
            for reaction, row in zip(reactions, table.itertuples(index=False)):
                if clean(row.ReactionSMILES) != reaction.reaction_smiles:
                    return False
                if clean(row.MappingDirection) != (
                    "forward" if reaction.flux >= 0 else "reverse"
                ):
                    return False
                status = clean(row.MappingStatus)
                if reaction.status == "ready" and not (
                    status
                    in (
                        "mapped",
                        "Mapped fragments cannot be aligned to input metabolites",
                    )
                    or status.startswith(("Invalid mapping: ", backend + " failed: "))
                ):
                    return False
            for reaction, row in zip(reactions, table.itertuples(index=False)):
                if reaction.is_transport:
                    continue
                reaction.mapped_smiles = clean(row.MappedReactionSMILES)
                reaction.confident = (
                    "" if pd.isna(row.MappingConfident) else row.MappingConfident
                )
                if backend == "LocalMapper" and reaction.confident in (0, 1):
                    reaction.confident = bool(reaction.confident)
                reaction.template = clean(row.MappingTemplate)
                if reaction.status == "ready":
                    reaction.status = clean(row.MappingStatus)
        return True
    except (AssertionError, ValueError, KeyError, AttributeError):
        return False


def prediction_table(reactions):
    return pd.DataFrame(
        [
            {
                "ReactionID": r.reaction_id,
                "ReactionName": r.reaction_name,
                "Equation": r.equation,
                "ReactionSMILES": r.reaction_smiles,
                "MappedReactionSMILES": r.mapped_smiles,
                "MappingDirection": "forward" if r.flux >= 0 else "reverse",
                "MappingConfident": r.confident,
            }
            for r in reactions
        ]
    )


def final_mapping_table(reactions, selection):
    choices = selection.set_index("ReactionID")
    columns = [
        "ReactionID",
        "ReactionName",
        "Flux",
        "SelectedMapper",
        "ReviewRequired",
        "AIGenerated",
        "MappedReactionSMILES",
        "MappingCoverage",
        "Equation",
        "MappingDirection",
    ]
    return pd.DataFrame(
        [
            {
                "ReactionID": r.reaction_id,
                "ReactionName": r.reaction_name,
                "Flux": r.flux,
                "Equation": r.equation,
                "MappingDirection": "forward" if r.flux >= 0 else "reverse",
                "SelectedMapper": choices.loc[r.reaction_id, "SelectedMapper"],
                "MappedReactionSMILES": r.mapped_smiles,
                "ReviewRequired": bool(choices.loc[r.reaction_id, "ReviewRequired"]),
                "AIGenerated": False,
                "MappingCoverage": choices.loc[r.reaction_id].get(
                    "MappingCoverage", "Complete" if r.status == "mapped" else "Unresolved"),
            }
            for r in reactions
        ],
        columns=columns,
    )


def summarize_sources(table, sources, carbon_ids, number, previous_counts, exclusive_metabolites=None):
    """Summarize sources and exclusivity across all recorded states of each product."""
    source_ids = {name: set() for name, _, _ in sources}
    mixed_ids = set()
    production_labels = {}
    records = (table[["Production", "CarbonLabels"]].drop_duplicates().itertuples(index=False, name=None)
               if isinstance(table, pd.DataFrame) else table)
    for production, carbon_labels in records:
        labels = set(carbon_labels.split(","))
        production_labels.setdefault(production, set()).update(labels)
        matched_sources = 0
        for name, _, source_labels in sources:
            if labels & source_labels:
                source_ids[name].update(carbon_ids.get(production, ()))
                matched_sources += 1
        if matched_sources > 1:
            mixed_ids.update(carbon_ids.get(production, ()))
    total = sum(map(len, carbon_ids.values()))
    counts = {}
    if len(sources) > 1:
        counts["All sources"] = len(set().union(*source_ids.values()))
    counts.update({name: len(ids) for name, ids in source_ids.items()})
    summary = {"Round": number, "TotalCarbonMetabolites": total}
    for source, count in counts.items():
        summary[f"{source}_LabeledCarbonMetabolites"] = count
        summary[f"{source}_LabeledCarbonMetabolitesIncrement"] = (
            count - previous_counts.get(source, 0)
        )
        summary[f"{source}_LabeledCarbonPercent"] = (
            round(100 * count / total, 2) if total else None
        )
    if len(sources) > 1:
        summary["MixedSources_LabeledCarbonMetabolites"] = len(mixed_ids)
    for name, _, source_labels in sources:
        exclusive_ids = set()
        for production, labels in sorted(production_labels.items()):
            if labels and labels <= source_labels:
                ids = carbon_ids.get(production, ())
                exclusive_ids.update(ids)
                if exclusive_metabolites is not None:
                    exclusive_metabolites.extend(
                        (number, name, identifier, production) for identifier in sorted(ids))
        summary[f"{name}_ExclusiveSourceMetabolites"] = len(exclusive_ids)
    previous_counts.update(counts)
    return summary



def run_carbon_flow(jobs, rounds, output_dir, max_processes=None, memory_gb=None):
    """Read step 2 independently; protect its workbooks and limit model concurrency."""
    if type(rounds) is not int or rounds < 1:
        raise ValueError("Propagation rounds must be a positive integer.")
    if not jobs:
        raise ValueError("Configure at least one model for propagation.")
    output_dir = Path(output_dir).resolve()
    tasks = []
    inputs = set()
    for name, job in jobs.items():
        methods = job.get("analysis_methods", ["Enumeration"])
        if (not isinstance(methods, (list, tuple)) or not methods
                or len(set(methods)) != len(methods) or any(m not in ("Enumeration", "EMU") for m in methods)):
            raise ValueError("Select unique analysis methods: Enumeration, EMU")
        workbook_limit(job.get("max_excel_size_mb", 50))
        if not name or any(c in name for c in '/\\:*?"<>|') or name in (".", ".."):
            raise ValueError(f"Invalid analysis name for a file: {name}")
        sources = job.get("carbon_sources")
        if not sources or len(set(sources)) != len(sources):
            raise ValueError(f"{name} carbon sources must be nonempty and unique.")
        mappings = job.get("mapping_files", [
            job.get("mapping_file", output_dir / f"{name}_2_atom_mapping.xlsx")
        ])
        if not mappings:
            raise ValueError("Configure at least one mapping workbook for " + name)
        for mapping_file in mappings:
            mapping_file = Path(mapping_file).resolve()
            if not mapping_file.is_file():
                raise FileNotFoundError(f"Run step2_atom_mapping.py first: {mapping_file}")
            inputs.add(mapping_file)
            if "input_file" in job:
                inputs.add(Path(job["input_file"]).resolve())
            suffix = "_AIGenerated" if mapping_file.stem.endswith("_AIGenerated") else ""
            tasks.append((
                name, job, rounds, mapping_file,
                output_dir / f"{name}_3_carbon_flow{suffix}.xlsx",
            ))
    if len({task[-1] for task in tasks}) != len(tasks):
        raise ValueError("Mapping workbooks would overwrite the same carbon-flow output.")
    if any(task[-1] in inputs for task in tasks):
        raise ValueError("Output cannot overwrite an input or mapping workbook.")
    print(f"Propagation rounds: {rounds}", flush=True)
    workers, _ = parallel_capacity(
        len(tasks), max_processes, memory_gb, task_kind="carbon propagation"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = [task[-1] for task in tasks]
    previous_saves = {p: p.stat().st_mtime_ns if p.exists() else None for p in outputs}
    try:
        if workers == 1:
            saved = [trace_with_status(*task) for task in tasks]
        else:
            with get_context("spawn").Pool(processes=workers) as pool:
                saved = pool.starmap(trace_with_status, tasks)
    except Exception as error:
        changed = [p for p in outputs if p.is_file()
                   and p.stat().st_mtime_ns != previous_saves[p]]
        saved = []
        for task in tasks:
            status = dict(model=task[0], requested_rounds=rounds, completed_rounds=0,
                          stop_reason="not_completed", output_file=None)
            status["mapping_file"] = str(task[3])
            if getattr(error, "run_status", {}).get("mapping_file") == str(task[3]):
                status.update(error.run_status)
            if task[-1] in changed:
                with pd.ExcelFile(task[-1], engine="openpyxl") as excel:
                    completed = int(pd.read_excel(excel, sheet_name="Summary").Round.max()) if "Summary" in excel.sheet_names else 0
                status.update(completed_rounds=completed, output_file=str(task[-1]))
                if status["stop_reason"] == "not_completed":
                    status["stop_reason"] = "batch_interrupted"
            saved.append(status)
        write_carbon_run_report(output_dir, tasks, saved, str(error))
        try:
            collect_carbon_summaries(output_dir)
        except Exception as summary_error:
            print("Summary collection failed: {}".format(summary_error), flush=True)
        raise
    write_carbon_run_report(output_dir, tasks, saved)
    collect_carbon_summaries(output_dir)
    return saved



def collect_carbon_summaries(output_dir):
    """Copy existing Enumeration summaries without loading carbon-flow records."""
    output_dir = Path(output_dir)
    sources = sorted(output_dir.glob("*_3_carbon_flow*.xlsx"))
    sources = [path for path in sources if path.name.endswith((
        "_3_carbon_flow.xlsx",
        "_3_carbon_flow_AIGenerated.xlsx",
    )) and not path.name.startswith("~$")]
    if not sources:
        return
    output = output_dir / "carbon_flow_summary.xlsx"
    workbook = Workbook()
    workbook.remove(workbook.active)
    used_names = set()
    try:
        for source in sources:
            model, version = source.stem.rsplit("_3_carbon_flow", 1)
            model = model.replace("[", "_").replace("]", "_")
            suffix = "_AI" if version else ""
            name = model[:31 - len(suffix)] + suffix
            number = 1
            while name.lower() in used_names:
                number += 1
                tail = suffix + "_" + str(number)
                name = model[:31 - len(tail)] + tail
            used_names.add(name.lower())
            original = load_workbook(source, read_only=True)
            try:
                if "Summary" not in original.sheetnames:
                    if {"EMU_Metabolites", "EMU_BoundaryInputs"} <= set(original.sheetnames):
                        continue
                    raise KeyError("Summary is missing from " + str(source))
                sheet = workbook.create_sheet(name)
                for row in original["Summary"].iter_rows(values_only=True):
                    sheet.append(row)
            finally:
                original.close()
        if not workbook.sheetnames:
            return
        clear_formatting(workbook)
        save_output(workbook, output)
    finally:
        workbook.close()
    print("Summary workbook: {}".format(output), flush=True)


def trace_with_status(strain, job, rounds, mapping_file, output):
    """Return a small status record from each independent propagation worker."""
    methods = job.get("analysis_methods", ["Enumeration"])
    status = dict(model=strain, method="+".join(methods),
                  allow_inorganic_carbon_recycling=job.get("allow_inorganic_carbon_recycling", False),
                  requested_rounds=rounds, completed_rounds=0,
                  stop_reason="round_limit", output_file=None,
                  mapping_file=str(mapping_file), carbon_sources=job["carbon_sources"])
    try:
        digest = sha256()
        with mapping_file.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        status["mapping_sha256"] = digest.hexdigest()
        trace_mapping_file(strain, dict(job, _run_status=status), rounds, mapping_file, output)
    except Exception as error:
        status["stop_reason"] = ("memory_limit" if isinstance(error, (MemoryPressureError, MemoryError))
                                 else "error")
        error.run_status = status
        raise
    return status


def write_carbon_run_report(output_dir, tasks, statuses, error=None):
    """Keep safety stops distinguishable from convergence for automated callers."""
    output = output_dir / "carbon_flow_run.json"
    methods = sorted({row.get("method", "Enumeration") for row in statuses})
    report = dict(version=1, method=";".join(methods), finished_at_utc=datetime.utcnow().isoformat() + "Z",
                  status="error" if error else "finished", models=statuses)
    if error:
        report["error"] = error
    temporary = output.with_suffix(".tmp")
    try:
        temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        temporary.replace(output)
    finally:
        if temporary.exists():
            temporary.unlink()
    print("Run report: {}".format(output), flush=True)


def trace_mapping_file(strain, job, rounds, mapping_file, output):
    """Validate final correspondences against the structures saved in step 2."""
    RDLogger.DisableLog("rdApp.warning")
    MemoryGuard(job.get("min_available_memory_gb", 2)).check("loading " + strain, force=True)
    reactions, source, metabolites = load_final_mappings(mapping_file, job)
    methods = job.get("analysis_methods", ["Enumeration"])
    status = job.get("_run_status", {})
    extra_sheets = {}
    if "EMU" in methods:
        from .source_contributions import solve_source_contributions
        guard = MemoryGuard(job.get("min_available_memory_gb", 2))
        tables, info = solve_source_contributions(reactions, metabolites, job["carbon_sources"],
            allow_inorganic_carbon_recycling=job.get("emu_allow_inorganic_carbon_recycling", True),
            tolerance=job.get("emu_tolerance", 1e-8), memory_guard=guard)
        extra_sheets = {
            "EMU_Metabolites": tables["MetaboliteContributions"].drop(columns=["FluxResidual", "Status"]),
            "EMU_BoundaryInputs": tables["BoundaryInputs"],
        }
        info.update(mapping_version="AI" if mapping_file.stem.endswith("_AIGenerated") else "Original")
        status.update(emu=info, emu_solved=True, emu_completed=False)
        if "Enumeration" not in methods:
            guard.check("exporting EMU results", force=True)
            workbook = Workbook()
            del workbook[workbook.active.title]
            try:
                for sheet, table in extra_sheets.items():
                    write_sheet(workbook, sheet, table, guard)
                save_output(workbook, output, max_bytes=workbook_limit(job.get("max_excel_size_mb", 50)))
            finally:
                workbook.close()
            status.update(emu_completed=True, output_file=str(output))
        print("{} EMU: {} carbon positions solved; combined output {}".format(strain, info["carbon_positions"], output), flush=True)
    if "Enumeration" in methods:
        return trace_analysis(strain, dict(job, _extra_sheets=extra_sheets), reactions,
                              metabolites, rounds, Workbook(), output)
    status["stop_reason"] = "emu_solved"
    return output


def load_final_mappings(mapping_file, job):
    """Load selected mappings through the carbon-propagation input checks."""
    with pd.ExcelFile(mapping_file, engine="openpyxl") as excel:
        reaction_sheet, metabolite_sheet = input_sheets(excel)
        final = pd.read_excel(
            excel, sheet_name="Final_Mapping", dtype={"ReactionID": str}
        ).fillna("")
        partials = (pd.read_excel(excel, sheet_name="Partial_Mapping",
                     dtype={"ReactionID": str, "RuleVersion": str}).fillna("")
                    if "Partial_Mapping" in excel.sheet_names else pd.DataFrame())
    reactions, source, metabolites = read_input(
        mapping_file, reaction_sheet, metabolite_sheet
    )
    # Preserve the original-input check for callers that explicitly supply step 1.
    if "input_file" in job:
        _, original_source, original_metabolites = read_input(
            job["input_file"], job["reaction_sheet"], job["metabolite_sheet"]
        )
        pd.testing.assert_frame_equal(source, original_source, rtol=1e-12, atol=1e-15)
        pd.testing.assert_frame_equal(metabolites, original_metabolites)
    required = {"ReactionID", "MappingDirection", "Equation", "MappedReactionSMILES"}
    if not required <= set(final.columns):
        raise ValueError("Final_Mapping is missing required columns.")
    if not final.ReactionID.is_unique or not set(final.ReactionID) <= {
        r.reaction_id for r in reactions
    }:
        raise ValueError(
            "Final_Mapping contains duplicated or unknown ReactionID values."
        )
    lookup = final.set_index("ReactionID")
    structure_issues = carbon_structure_issues(metabolites)
    carbon_by_name = {clean(row.MetaboliteName): carbon_presence(row)
                      for row in metabolites.itertuples(index=False)}
    for reaction in reactions:
        reaction.status = "Not selected"
        if reaction.reaction_id not in lookup.index:
            continue
        row = lookup.loc[reaction.reaction_id]
        if (
            row.MappingDirection != ("forward" if reaction.flux >= 0 else "reverse")
            or row.Equation != reaction.equation
        ):
            raise ValueError(
                "Final_Mapping input or direction changed; rerun the mapping stage: {} ({!r} vs {!r})".format(reaction.reaction_id, row.Equation, reaction.equation)
            )
        coverage = row.get("MappingCoverage", "Complete" if row.MappedReactionSMILES else "Unresolved")
        if coverage not in ("Complete", "Partial", "Unresolved", "Boundary", "NoCarbon"):
            raise ValueError("Unknown MappingCoverage: " + str(coverage))
        reaction.mapped_smiles = row.MappedReactionSMILES
        inconsistent = [c.metabolite for c in reaction.left + reaction.right
                        if c.metabolite in structure_issues]
        if coverage in ("Complete", "Partial") and reaction.mapped_smiles and inconsistent:
            raise ValueError("Final mapping uses carbon-inconsistent structures: "
                             + ", ".join(sorted(set(inconsistent)))
                             + "; review structures and rerun step 2 selection")
        if coverage != "Complete":
            if coverage == "NoCarbon":
                names = [name for side in split_equation(reaction.equation) for _, name in parse_side(side)]
                if not names or any(carbon_by_name.get(name) is not False for name in names):
                    raise ValueError("NoCarbon coverage conflicts with saved inputs: " + reaction.reaction_id)
            if coverage == "Boundary" and all(split_equation(reaction.equation)):
                raise ValueError("Boundary coverage requires an open reaction side: " + reaction.reaction_id)
            if coverage in ("Unresolved", "Boundary", "NoCarbon"):
                # Retained failed candidates are audit data, never usable atom transfers.
                reaction.mapped_smiles = ""
                continue
            if coverage == "Partial" and reaction.reaction_id in set(partials.get("ReactionID", [])):
                continue  # Reduced rule subreactions are checked against their exact records below.
            if coverage == "Partial" and reaction.mapped_smiles:
                checked = carbon_candidate(reaction)
                if ((checked.status != "partial" and not reaction.omitted_components)
                        or checked.mapped_smiles != reaction.mapped_smiles):
                    raise ValueError("Saved partial mapper correspondence must be a validated carbon-only view")
                reaction.status = "partial"
            elif reaction.mapped_smiles:
                raise ValueError("Unresolved or boundary coverage cannot contain a mapped reaction")
            continue
        if reaction.omitted_components:
            raise ValueError("Complete coverage cannot omit original reaction participants: " + reaction.reaction_id)
        if not reaction.mapped_smiles:
            raise ValueError("Complete coverage requires a full mapped reaction")
        mapping_signatures(reaction)
        reaction.status = "mapped"
    load_partials(reactions, final, partials, metabolites)
    return reactions, source, metabolites


def run_analyses(
    jobs,
    rounds,
    batch_size,
    max_reaction_atoms,
    output_dir,
    backends,
    mapper_configs,
    refresh=False,
    max_mapper_processes=None,
    mapper_memory_gb=6,
    max_carbon_processes=None,
    carbon_memory_gb=None,
    stage="all",
    structure_options=None,
    manual_mapping_file=None,
):
    """Predict independently, select one mapping per reaction and propagate."""
    if stage == "trace":
        return run_carbon_flow(
            jobs, rounds, output_dir, max_carbon_processes, carbon_memory_gb
        )
    if not backends or len(set(backends)) != len(backends):
        raise ValueError("Mapper list must be nonempty and unique.")
    if stage not in ("all", "map", "trace"):
        raise ValueError("Stage must be all, map or trace.")
    for backend in backends:
        if (
            backend
            in (
                "Mapping_Selection",
                "Final_Mapping",
                "Partial_Mapping",
                "Summary",
                "Mapping_Info",
                "Identity",
            )
            or len(backend) > 31
            or any(c in backend for c in "[]")
            or backend.startswith("Carbon_Flow_")
        ):
            raise ValueError(
                "Mapper names must not conflict with output sheets and must fit Excel limits."
            )
        if stage != "trace":
            validate_mapper_config(backend, mapper_configs[backend])
    configure(batch_size, max_reaction_atoms)
    if stage != "map" and (type(rounds) is not int or rounds < 1):
        raise ValueError("Propagation rounds must be a positive integer.")
    for name, job in jobs.items():
        if not name or any(c in name for c in '/\\:*?"<>|') or name in (".", ".."):
            raise ValueError(f"Invalid analysis name for a file: {name}")
        if not job["input_file"].is_file():
            raise FileNotFoundError(job["input_file"])
        if stage != "map" and (
            not job["carbon_sources"]
            or len(set(job["carbon_sources"])) != len(job["carbon_sources"])
        ):
            raise ValueError(f"{name} carbon sources must be nonempty and unique.")
        if any(
            backend in (job["reaction_sheet"], job["metabolite_sheet"])
            for backend in backends
        ):
            raise ValueError("Mapper names cannot match input sheet names.")
    input_paths = {job["input_file"].resolve() for job in jobs.values()}
    if any(
        (output_dir / filename).resolve() in input_paths
        for name in jobs
        for filename in (f"{name}_2_atom_mapping.xlsx", f"{name}_3_carbon_flow.xlsx")
    ):
        raise ValueError("Output cannot overwrite any input workbook.")
    if stage != "map":
        print(f"Propagation rounds: {rounds}", flush=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    RDLogger.DisableLog("rdApp.warning")
    prepared, pending = {}, []
    for strain, job in jobs.items():
        output = output_dir / f"{strain}_2_atom_mapping.xlsx"
        reactions, source, metabolites = read_input(
            job["input_file"],
            job["reaction_sheet"],
            job["metabolite_sheet"],
            structure_options,
            strain,
        )
        predictions, metadata = {}, []
        for backend in backends:
            settings = mapper_configs[backend]
            candidates = deepcopy(reactions)
            reused = False
            if not refresh:
                reused = reuse_predictions(
                    output, backend, settings, candidates, source, metabolites, job
                )
                if reused:
                    print(
                        f"{strain} {backend}: inputs and configuration match; reusing predictions.",
                        flush=True,
                    )
            if not reused and any(r.status == "ready" for r in candidates):
                pending.append((strain, backend, candidates, settings))
            predictions[backend] = candidates
            info = dict(
                Mapper=backend,
                Adapter=settings["module"],
                Python=settings["python"],
                ReusedPredictions=reused,
                Rounds=rounds,
                BatchSize=batch_size,
                MaxReactionAtoms=max_reaction_atoms,
                OptionsJSON=json.dumps(
                    settings.get("options", {}), ensure_ascii=False, sort_keys=True
                ),
                SelectionPolicy="Select carbon agreement supported by at least two tools; otherwise require review; do not compare scores across tools.",
            )
            info.update(settings.get("options", {}))
            info["SourceInputSHA256"] = job.get("input_sha256", "")
            metadata.append(info)
        prepared[strain] = metabolites, predictions, metadata
    workers, cpu_threads = parallel_capacity(
        len(pending), max_mapper_processes, mapper_memory_gb
    )
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                run_mapper, strain, candidates, backend, settings, cpu_threads
            ): (strain, backend)
            for strain, backend, candidates, settings in pending
        }
        for future in as_completed(futures):
            strain, backend = futures[future]
            future.result()
            print(f"{strain} {backend}: predictions complete.", flush=True)
    workers, _ = parallel_capacity(
        len(jobs),
        max_carbon_processes,
        carbon_memory_gb,
        task_kind="carbon propagation",
    )
    tasks = [
        (strain, job, *prepared[strain], rounds, backends, output_dir, stage, manual_mapping_file)
        for strain, job in jobs.items()
    ]
    if workers == 1:
        for task in tasks:
            finish_analysis(*task)
    else:
        with get_context("spawn").Pool(processes=workers) as pool:
            pool.starmap(finish_analysis, tasks)


def finish_analysis(
    strain,
    job,
    metabolites,
    predictions,
    metadata,
    rounds,
    backends,
    output_dir,
    stage="all",
    manual_mapping_file=None,
):
    """Write one strain workbook without sharing its state catalog."""
    RDLogger.DisableLog("rdApp.warning")
    print(f"{strain}: writing selected atom mappings.", flush=True)
    output = output_dir / f"{strain}_2_atom_mapping.xlsx"
    reactions, selection = select_mappings(predictions, carbon_structure_issues(metabolites))
    apply_manual_mappings(strain, reactions, selection, metabolites, manual_mapping_file, job.get("source_model"))
    partials = apply_rules(reactions, selection, metabolites)
    if "StructureReviewRequired" in metabolites:
        review_names = set(
            metabolites.loc[
                metabolites.StructureReviewRequired.fillna(False).eq(True),
                "MetaboliteName",
            ]
        )
        for index, reaction in enumerate(reactions):
            if any(
                component.metabolite in review_names
                for component in reaction.left + reaction.right
            ):
                selection.at[index, "ReviewRequired"] = True
                selection.at[index, "SelectionReason"] += (
                    "; Automatically selected structure requires review"
                )

    temporary = output.with_suffix(".tmp.xlsx")
    copy_input(job["input_file"], temporary)
    workbook = load_workbook(temporary)
    clear_formatting(workbook)
    write_sheet(workbook, job["metabolite_sheet"], metabolites)
    for backend in backends:
        write_sheet(workbook, backend, prediction_table(predictions[backend]))
    selection_columns = ["ReactionID", "ReactionName", "Flux", "Equation", "MappingDirection",
                         "SelectedMapper", "ReviewRequired", "SelectionReason", "SelectedMappedReactionSMILES"]
    write_sheet(workbook, "Mapping_Selection", selection[selection_columns])
    write_sheet(workbook, "Final_Mapping", final_mapping_table(reactions, selection))
    if not partials.empty:
        write_sheet(workbook, "Partial_Mapping", partials)
    metadata.extend(
        dict(
            Mapper=backend,
            ReactionID=r.reaction_id,
            MappingStatus=r.status,
            MappingTemplate=r.template,
            MappingInputEquation=r.mapping_equation,
            OmittedComponentsJSON=json.dumps(r.omitted_components, sort_keys=True),
        )
        for backend in backends
        for r in predictions[backend]
    )
    write_sheet(workbook, "Mapping_Info", pd.DataFrame(metadata))
    workbook.move_sheet(
        "Mapping_Selection", offset=-workbook.sheetnames.index("Mapping_Selection")
    )
    workbook.move_sheet(
        "Final_Mapping", offset=1 - workbook.sheetnames.index("Final_Mapping")
    )
    save_output(workbook, output)
    chosen = selection.loc[selection.SelectedMapper != ""]
    print(
        f"{strain} consensus: {int(selection.MappingCoverage.eq('Complete').sum())}/{len(reactions)} complete, "
        f"{int(selection.MappingCoverage.eq('Partial').sum())} partial carbon mappings; "
        f"{int(selection.ReviewRequired.sum())} require review; "
        f"selected tools {chosen.SelectedMapper.value_counts().to_dict()}",
        flush=True,
    )
    if stage == "map":
        workbook.close()
        return
    trace_analysis(
        strain,
        job,
        reactions,
        metabolites,
        rounds,
        workbook,
        output_dir / f"{strain}_3_carbon_flow.xlsx",
    )


def trace_analysis(strain, job, reactions, metabolites, rounds, workbook, output):
    """Reuse final atom correspondences; never invoke a prediction tool here."""
    workbook.close()
    summary = []
    previous_counts = {}
    carbon_ids, unknown_count = carbon_metabolites(metabolites, reactions)
    carbon_total = sum(len(ids) for ids in carbon_ids.values())
    print(
        f"{strain}: {carbon_total} carbon metabolites in the mapped network; {unknown_count} missing/invalid structures in Metabolites.",
        flush=True,
    )
    smiles = {
        clean(row.MetaboliteName): get_smiles(row)
        for row in metabolites.itertuples(index=False)
    }
    sources = initialize_sources(job["carbon_sources"], smiles)
    memory_guard = MemoryGuard(job.get("min_available_memory_gb", 2))
    completed_round = 0
    max_bytes = workbook_limit(job.get("max_excel_size_mb", 50))
    status = job.get("_run_status", {})
    batch_rows = job.get("output_batch_rows", 50000)
    status.update(output_batch_rows=batch_rows, max_excel_bytes=max_bytes)
    detail_rounds = job.get("export_detail_rounds", 3)
    status["export_detail_rounds"] = detail_rounds
    store = FlowOutput(output, max_bytes, batch_rows, summary, memory_guard,
                       detail_rounds=detail_rounds, extra_sheets=job.get("_extra_sheets"))
    iterator = propagate_rounds(
        strain, reactions, metabolites, rounds, job["carbon_sources"], sources,
        allow_inorganic_carbon_recycling=job.get("allow_inorganic_carbon_recycling", False),
        memory_guard=memory_guard, row_store=store,
    )
    try:
        for number, _, added, converged in iterator:
            exclusive_metabolites = [] if job.get("export_exclusive_source_metabolites", False) else None
            row = summarize_sources(store.source_records(), sources, carbon_ids, number,
                                    previous_counts, exclusive_metabolites)
            # Each published round is provisional final; a rejected round leaves the old file intact.
            if summary:
                for key in list(summary[-1]):
                    if key.endswith("_ExclusiveSourceMetabolites"):
                        del summary[-1][key]
            store.exclusive_metabolites = exclusive_metabolites
            row.update(Method="Enumeration", InorganicCarbonRecycling=job.get("allow_inorganic_carbon_recycling", False))
            summary.append(row)
            store.complete_round()
            if job.get("_extra_sheets"):
                status["emu_completed"] = True
            completed_round = number
            status["stop_reason"] = "converged" if converged else "round_limit"
            print(f"{strain} R{number} completed (detail export: {detail_rounds is None or number <= detail_rounds}): {store.row_count} rows, {added} new states -> {output}", flush=True)
    except (WorkbookSizeLimitError, WorksheetRowLimitError) as error:
        status["stop_reason"] = ("file_size_limit" if isinstance(error, WorkbookSizeLimitError) else "row_limit")
        print("{}: {}; {} complete round(s) retained. This stop is not convergence."
              .format(strain, error, completed_round), flush=True)
    except (MemoryPressureError, MemoryError) as error:
        print(
            "{}: memory safety stop; {} round(s) saved in this run. "
            "The last saved workbook is unchanged by the incomplete round. "
            "If no rounds were saved, any existing output belongs to an earlier run. "
            "Rerunning starts from round 1; this is not a resumable checkpoint. {}"
            .format(strain, completed_round, error), flush=True,
        )
        raise
    finally:
        status.update(completed_rounds=completed_round,
                      output_file=str(output) if completed_round or status.get("emu_completed") else None)
        iterator.close()
        store.close()
    return output if completed_round or status.get("emu_completed") else None
