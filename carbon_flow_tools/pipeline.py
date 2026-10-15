"""Prediction, consensus selection and carbon propagation stages."""

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from copy import deepcopy
import json
from multiprocessing import get_context
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook, Workbook
from rdkit import Chem, RDLogger

from .carbon_trace import propagate_rounds, initialize_sources
from .config import configure, parallel_capacity
from .excel_io import clear_formatting, copy_input, input_sheets, write_sheet
from .consensus import select_mappings
from .consensus import mapping_signatures
from .mapper_runner import run_mapper, validate_mapper_config
from .reactions import clean, get_smiles, prepare_reactions
from .smiles_completion import complete_smiles


def carbon_metabolites(metabolites, reactions):
    """Count carbon metabolite IDs participating in usable final mappings."""
    mapped_names = {
        component.metabolite
        for reaction in reactions
        if reaction.status == "mapped" and reaction.mapped_smiles
        for component in reaction.left + reaction.right
    }
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
            }
            for r in reactions
            if r.status == "mapped"
        ],
        columns=columns,
    )


def summarize_sources(table, sources, carbon_ids, number, previous_counts):
    """Return one row per round, with separate columns for each source."""
    source_ids = {name: set() for name, _, _ in sources}
    for row in (
        table[["Production", "CarbonLabels"]].drop_duplicates().itertuples(index=False)
    ):
        labels = set(row.CarbonLabels.split(","))
        for name, _, source_labels in sources:
            if labels & source_labels:
                source_ids[name].update(carbon_ids.get(row.Production, ()))
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
    previous_counts.update(counts)
    return summary


def write_flow(workbook, number, table):
    # Paginate instead of truncating at Excel's row limit.
    for page, start in enumerate(range(0, max(1, len(table)), 1048575), 1):
        suffix = "" if page == 1 else f"_p{page}"
        write_sheet(
            workbook,
            f"Carbon_Flow_R{number}" + suffix,
            table.iloc[start : start + 1048575],
        )


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
        if not name or any(c in name for c in '/\\:*?"<>|') or name in (".", ".."):
            raise ValueError(f"Invalid analysis name for a file: {name}")
        sources = job.get("carbon_sources")
        if not sources or len(set(sources)) != len(sources):
            raise ValueError(f"{name} carbon sources must be nonempty and unique.")
        mapping_file = Path(
            job.get("mapping_file", output_dir / f"{name}_2_atom_mapping.xlsx")
        ).resolve()
        if not mapping_file.is_file():
            raise FileNotFoundError(f"Run step2_atom_mapping.py first: {mapping_file}")
        inputs.add(mapping_file)
        if "input_file" in job:
            inputs.add(Path(job["input_file"]).resolve())
        suffix = "_AIGenerated" if mapping_file.stem.endswith("_AIGenerated") else ""
        tasks.append(
            (
                name,
                job,
                rounds,
                mapping_file,
                output_dir / f"{name}_3_carbon_flow{suffix}.xlsx",
            )
        )
    if any(task[-1] in inputs for task in tasks):
        raise ValueError("Output cannot overwrite an input or mapping workbook.")
    print(f"Propagation rounds: {rounds}", flush=True)
    workers, _ = parallel_capacity(
        len(tasks), max_processes, memory_gb, task_kind="carbon propagation"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    if workers == 1:
        for task in tasks:
            trace_mapping_file(*task)
    else:
        with get_context("spawn").Pool(processes=workers) as pool:
            pool.starmap(trace_mapping_file, tasks)


def trace_mapping_file(strain, job, rounds, mapping_file, output):
    """Validate final correspondences against the structures saved in step 2."""
    RDLogger.DisableLog("rdApp.warning")
    with pd.ExcelFile(mapping_file, engine="openpyxl") as excel:
        reaction_sheet, metabolite_sheet = input_sheets(excel)
        final = pd.read_excel(
            excel, sheet_name="Final_Mapping", dtype={"ReactionID": str}
        ).fillna("")
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
                "Final_Mapping input or direction changed; rerun the mapping stage."
            )
        reaction.mapped_smiles = row.MappedReactionSMILES
        mapping_signatures(reaction)
        reaction.status = "mapped"
    trace_analysis(strain, job, reactions, metabolites, rounds, Workbook(), output)


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
        raise ValueError("Output cannot overwrite any MATLAB input.")
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
        (strain, job, *prepared[strain], rounds, backends, output_dir, stage)
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
):
    """Write one strain workbook without sharing its state catalog."""
    RDLogger.DisableLog("rdApp.warning")
    print(f"{strain}: writing selected atom mappings.", flush=True)
    output = output_dir / f"{strain}_2_atom_mapping.xlsx"
    reactions, selection = select_mappings(predictions)
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
    write_sheet(workbook, "Mapping_Selection", selection)
    write_sheet(workbook, "Final_Mapping", final_mapping_table(reactions, selection))
    metadata.extend(
        dict(
            Mapper=backend,
            ReactionID=r.reaction_id,
            MappingStatus=r.status,
            MappingTemplate=r.template,
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
        f"{strain} consensus: {len(chosen)}/{len(reactions)} reactions have complete carbon mappings; "
        f"{int(chosen.ReviewRequired.sum())} require review; "
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


def save_output(workbook, output):
    """Keep the previous completed workbook until its replacement is saved."""
    # Use an explicit font rather than the locale-dependent theme font.
    for font in workbook._fonts:
        font.name = "Calibri"
        font.scheme = None
        font.charset = None
    workbook._fonts._dict = {font: i for i, font in enumerate(workbook._fonts)}
    temporary = output.with_suffix(".tmp.xlsx")
    workbook.save(temporary)
    try:
        temporary.replace(output)
    except PermissionError:
        raise PermissionError(
            f"Close {output.name} in Excel; completed output is in {temporary.name}."
        )


def trace_analysis(strain, job, reactions, metabolites, rounds, workbook, output):
    """Reuse final atom correspondences; never invoke a prediction tool here."""
    workbook.close()
    workbook = Workbook()
    del workbook[workbook.active.title]
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
    for number, table, added, converged in propagate_rounds(
        strain,
        reactions,
        metabolites,
        rounds,
        job["carbon_sources"],
        sources,
        allow_inorganic_carbon_recycling=job.get(
            "allow_inorganic_carbon_recycling", False
        ),
    ):
        write_flow(workbook, number, table)
        source_summary = summarize_sources(
            table, sources, carbon_ids, number, previous_counts
        )
        summary.append(source_summary)
        write_sheet(workbook, "Summary", pd.DataFrame(summary))
        workbook.move_sheet("Summary", offset=-workbook.sheetnames.index("Summary"))
        save_output(workbook, output)
        print(
            f"{strain} R{number} saved: {len(table)} rows, {added} new states -> {output}",
            flush=True,
        )
    workbook.close()
