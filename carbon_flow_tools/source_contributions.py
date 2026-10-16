"""Fixed-flux, single-carbon origin balances; not a full isotopomer simulator."""

from collections import defaultdict, deque
import math
import re
import warnings

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import MatrixRankWarning, spsolve

from .atom_helpers import mapped_components
from .consensus import mapping_signatures
from .reactions import clean, get_smiles, parse_side, split_equation

TOL = 1e-8
UNKNOWN = "Unknown"


def solve_source_contributions(reactions, metabolites, source_names):
    """Balance carbon positions at metabolic and origin steady state.

    Fractions are unconditional: Unknown is retained, never redistributed.
    Only external inflows receive prescribed origins; internal source pools mix.
    Use the validated selected mapping as the positional representative.
    Graph symmetry is recorded, not converted into additional Unknown.
    """
    if not source_names or len(set(source_names)) != len(source_names):
        raise ValueError("Configure unique, nonempty carbon sources.")
    if (metabolites.MetaboliteName.map(clean).eq("").any()
            or metabolites.MetaboliteName.map(clean).duplicated().any()
            or metabolites.MetaboliteID.map(clean).eq("").any()
            or metabolites.MetaboliteID.map(clean).duplicated().any()):
        raise ValueError("Contribution analysis requires unique metabolite IDs and names.")
    if len({r.reaction_id for r in reactions}) != len(reactions):
        raise ValueError("Duplicate reaction IDs.")
    nodes, positions, ambiguous, structure_issues = [], {}, set(), {}
    rows_by_name = {clean(r.MetaboliteName): r for r in metabolites.itertuples(index=False)}
    for name, row in rows_by_name.items():
        structure = get_smiles(row)
        molecule = Chem.MolFromSmiles(structure) if structure else None
        if molecule is None:
            structure_issues[name] = "Missing or invalid structure; carbon positions unavailable"
            continue
        formula = clean(getattr(row, "MolecularFormula", ""))
        if formula and re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", formula):
            counts = re.findall(r"([A-Z][a-z]?)(\d*)", formula)
            expected = sum(int(count or 1) for element, count in counts if element == "C")
            actual = sum(a.GetAtomicNum() == 6 for a in molecule.GetAtoms())
            if actual != expected:
                structure_issues[name] = "Structure/formula carbon count mismatch"
                continue
        positions[name] = []
        for part, fragment in enumerate(structure.split(".")):
            mol = Chem.MolFromSmiles(fragment)
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
            ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=True))
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() != 6:
                    continue
                node = (name, part, atom.GetIdx())
                positions[name].append(node)
                nodes.append(node)
                if ranks.count(ranks[atom.GetIdx()]) > 1:
                    ambiguous.add(node)
    for name in source_names:
        if not positions.get(name):
            raise ValueError("Carbon source has no usable carbon structure: " + name)
    if not nodes:
        raise ValueError("No carbon positions available.")

    incoming, outgoing = defaultdict(float), defaultdict(float)
    directed, imports = [], []
    for reaction in reactions:
        if not math.isfinite(reaction.flux) or reaction.flux == 0:
            raise ValueError("Fluxes must be finite and nonzero: " + reaction.reaction_id)
        left, right = split_equation(reaction.equation)
        if reaction.flux < 0:
            left, right = right, left
        left, right = parse_side(left), parse_side(right)
        rate = abs(reaction.flux)
        for side, totals in ((left, outgoing), (right, incoming)):
            for coefficient, name in side:
                if name not in rows_by_name or coefficient <= 0 or not math.isfinite(coefficient):
                    raise ValueError("Invalid or unknown reaction term: " + name)
                totals[name] += coefficient * rate
        if not left:
            for coefficient, name in right:
                if positions.get(name):
                    origin = name if name in source_names else "External:" + name
                    imports.append((reaction.reaction_id, name, coefficient * rate, origin))
        directed.append((reaction, left, right, rate))
    origins = list(source_names) + sorted({x[3] for x in imports} - set(source_names))
    if UNKNOWN in origins or len(set(origins + [UNKNOWN])) != len(origins) + 1:
        raise ValueError("Source names collide with reserved origin categories.")
    reserved = {"MetaboliteID", "MetaboliteName", "CarbonCount", "Inflow", "Outflow",
                "FluxResidual", "Status", "CarbonPosition", "FragmentIndex", "AtomIndex", "AmbiguousPosition"}
    if reserved.intersection(origins):
        raise ValueError("Source names collide with result metadata columns.")
    origins.append(UNKNOWN)
    index = {node: i for i, node in enumerate(nodes)}
    origin_index = {name: i for i, name in enumerate(origins)}
    n = len(nodes)
    rhs = np.zeros((n, len(origins)))
    edge_rows, edge_cols, edge_rates = [], [], []
    diagnostics, symmetry_mapped = [], 0
    for reaction, left, right, rate in directed:
        status, unknown_rate = "sink", 0.0
        if not left:
            status = "external_input"
            for coefficient, name in right:
                origin = name if name in source_names else "External:" + name
                for node in positions.get(name, []):
                    rhs[index[node], origin_index[origin]] += coefficient * rate
        elif right:
            status = "mapped"
            transfers = None
            if reaction.status == "mapped":
                mapping_signatures(reaction)
                components = mapped_components(reaction)
                source_by_map = {}
                for component, atom_maps in components[0]:
                    mol = Chem.MolFromSmiles(component.smiles)
                    for atom in mol.GetAtoms():
                        if atom.GetAtomicNum() == 6:
                            source_by_map[atom_maps[atom.GetIdx()]] = (
                                component.metabolite, component.component_index, atom.GetIdx())
                transfers = []
                for component, atom_maps in components[1]:
                    mol = Chem.MolFromSmiles(component.smiles)
                    for atom in mol.GetAtoms():
                        if atom.GetAtomicNum() == 6:
                            target = (component.metabolite, component.component_index, atom.GetIdx())
                            if target in index:
                                transfers.append((target, source_by_map[atom_maps[atom.GetIdx()]]))
                            else:
                                status = "unavailable_product_structure"
                expected = sum(c * len(positions.get(name, [])) for c, name in right)
                if len(transfers) != expected:
                    raise ValueError("Mapping/stoichiometry mismatch: " + reaction.reaction_id)
            if transfers is None:
                status = "missing_mapping"
                for coefficient, name in right:
                    for node in positions.get(name, []):
                        rhs[index[node], -1] += coefficient * rate
                        unknown_rate += coefficient * rate
            else:
                for target, source in transfers:
                    if source not in index:
                        rhs[index[target], -1] += rate
                        unknown_rate += rate
                        status = "unavailable_precursor_structure"
                    else:
                        # The selected map defines correspondence, including symmetry representatives.
                        edge_rows.append(index[target])
                        edge_cols.append(index[source])
                        edge_rates.append(rate)
                        if target in ambiguous or source in ambiguous:
                            symmetry_mapped += 1
        diagnostics.append(dict(ReactionID=reaction.reaction_id, Flux=reaction.flux,
                                Status=status, UnknownCarbonInputRate=unknown_rate))

    q = np.array([incoming[node[0]] for node in nodes])
    balance = {name: incoming[name] - outgoing[name] for name in rows_by_name}
    invalid = {name for name in rows_by_name if abs(balance[name]) >
               1e-6 + 1e-7 * max(incoming[name], outgoing[name])}

    local_rhs = rhs.copy()
    weights = coo_matrix((edge_rates, (edge_rows, edge_cols)), shape=(n, n)).tocsr()
    fixed = {i for i, node in enumerate(nodes) if q[i] == 0 or node[0] in invalid}
    # An unanchored carbon cycle cannot acquire an origin from steady-state fluxes.
    reverse = weights.transpose().tocsr()
    reachable = set(np.flatnonzero(local_rhs.sum(axis=1) > 0)) | fixed
    queue = deque(reachable)
    while queue:
        i = queue.popleft()
        for target in reverse.indices[reverse.indptr[i]:reverse.indptr[i + 1]]:
            if int(target) not in reachable:
                reachable.add(int(target))
                queue.append(int(target))
    unanchored = set(range(n)) - reachable
    fixed |= unanchored
    matrix = (diags(q) - weights).tolil()
    for i in fixed:
        matrix.rows[i], matrix.data[i] = [i], [1.0]
        local_rhs[i] = 0
        local_rhs[i, -1] = 1
    matrix = matrix.tocsr()
    scales = np.asarray(abs(matrix).sum(axis=1)).ravel()
    scaled = diags(1 / scales).dot(matrix).tocsc()
    scaled_rhs = local_rhs / scales[:, None]
    with warnings.catch_warnings():
        warnings.simplefilter("error", MatrixRankWarning)
        fractions = np.asarray(spsolve(scaled, scaled_rhs)).reshape(n, len(origins))
    residual = float(np.max(np.abs(scaled.dot(fractions) - scaled_rhs)))
    if (not np.isfinite(fractions).all() or residual > TOL
            or fractions.min() < -TOL or fractions.max() > 1 + TOL
            or np.max(np.abs(fractions.sum(axis=1) - 1)) > TOL):
        raise ValueError("Origin balance solution failed probability/residual checks.")
    fractions = np.clip(fractions, 0, 1)

    atom_rows, metabolite_rows = [], []
    for name, row in rows_by_name.items():
        atom_ids = [index[node] for node in positions.get(name, [])]
        if atom_ids:
            means = fractions[atom_ids].mean(axis=0)
            status = ("invalid_flux_balance" if name in invalid else
                      "no_inflow" if incoming[name] == 0 else
                      "unanchored_cycle" if any(i in unanchored for i in atom_ids) else
                      "complete" if means[-1] <= TOL else "partial_unknown")
        else:
            means, status = [np.nan] * len(origins), structure_issues.get(name, "no_carbon")
        common = dict(MetaboliteID=clean(row.MetaboliteID), MetaboliteName=name)
        metabolite_rows.append(dict(common, CarbonCount=len(atom_ids),
            Inflow=incoming[name], Outflow=outgoing[name], FluxResidual=balance[name],
            Status=status, **{o: float(v) for o, v in zip(origins, means)}))
        for ordinal, node in enumerate(positions.get(name, []), 1):
            i = index[node]
            atom_rows.append(dict(common, CarbonPosition=ordinal, FragmentIndex=node[1],
                AtomIndex=node[2], AmbiguousPosition=node in ambiguous,
                **{o: float(v) for o, v in zip(origins, fractions[i])}))
    atoms = pd.DataFrame(atom_rows)
    compounds = pd.DataFrame(metabolite_rows)
    info = dict(method="fixed_flux_single_carbon_origin_balance", max_scaled_residual=residual,
                carbon_positions=n, ambiguous_carbon_positions=len(ambiguous),
                symmetry_mapped_transfers=symmetry_mapped,
                symmetry_policy="Follow validated Final_Mapping representatives; graph symmetry does not add Unknown or imply random scrambling.",
                unanchored_positions=len(unanchored), invalid_balance_metabolites=sorted(invalid),
                unavailable_structure_metabolites=len(structure_issues),
                unknown_interpretation="Known fractions are contributions resolved by this conservative model. Unknown may also belong to any source; known fractions are not renormalized.",
                configured_sources_without_import=sorted(set(source_names) - {x[3] for x in imports}),
                assumptions="Metabolic and origin steady state; well-mixed pools; net-direction fluxes; no isotope effects; CO2 refixation included where mapped; Validated selected atom correspondences are positional representatives, including graph-equivalent atoms; no extra symmetry scrambling is assumed.",
                scope="Mean carbon-origin fractions, not molecular source combinations, experimental measurements or source essentiality.",
                position_convention="CarbonPosition is 1-based carbon order in saved SMILES; FragmentIndex and AtomIndex are 0-based, not biochemical numbering.")
    return dict(AtomContributions=atoms, MetaboliteContributions=compounds,
                ReactionDiagnostics=pd.DataFrame(diagnostics),
                BoundaryInputs=pd.DataFrame(imports, columns=["ReactionID", "MetaboliteName", "MoleculeInputRate", "Source"])), info


def run_source_contributions(jobs, output_dir, analysis_plan=None):
    """Export each mapping version independently of Enumeration."""
    from datetime import datetime
    from hashlib import sha256
    import json
    from pathlib import Path
    from openpyxl import Workbook
    from rdkit import RDLogger
    from .cli import project_path
    from .config import MemoryGuard
    from .excel_io import input_sheets, save_output, write_sheet, workbook_limit
    from .pipeline import load_final_mappings

    RDLogger.DisableLog("rdApp.warning")
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "source_contributions_run.json"
    report = dict(version=1, method="EMU", status="running", models=[])
    if analysis_plan:
        plan = project_path(analysis_plan)
        report.update(analysis_plan=str(plan), plan_sha256=sha256(plan.read_bytes()).hexdigest())
    try:
        for name, job in jobs.items():
            if not name or any(c in name for c in '/\\:*?"<>|') or name in (".", ".."):
                raise ValueError("Invalid analysis name for a file: " + str(name))
            for mapping in job.get("mapping_files", [job["mapping_file"]]):
                mapping = Path(mapping).resolve()
                guard = MemoryGuard(job.get("min_available_memory_gb", 2))
                guard.check("loading EMU inputs for " + name, force=True)
                if "source_input_file" not in job:
                    raise ValueError("EMU requires the original step-1 input from an analysis plan.")
                source_file = Path(job["source_input_file"])
                digest = sha256(source_file.read_bytes()).hexdigest()
                if digest != job.get("input_sha256"):
                    raise ValueError("Step-1 input hash changed: " + str(source_file))
                suffix = "_AIGenerated" if mapping.stem.endswith("_AIGenerated") else ""
                output = output_dir / (name + "_3_carbon_flow_EMU" + suffix + ".xlsx")
                if output.resolve() in (mapping.resolve(), source_file.resolve()):
                    raise ValueError("Output cannot overwrite an input.")
                with pd.ExcelFile(source_file, engine="openpyxl") as excel:
                    reaction_sheet, _ = input_sheets(excel)
                    original_reactions = pd.read_excel(excel, sheet_name=reaction_sheet)
                with pd.ExcelFile(mapping, engine="openpyxl") as excel:
                    metadata = pd.read_excel(excel, sheet_name="Mapping_Info")
                hashes = set(metadata.SourceInputSHA256.dropna()) if "SourceInputSHA256" in metadata else set()
                if hashes != {digest}:
                    raise ValueError("Mapping provenance does not match step 1: " + str(mapping))
                reactions, source, metabolites = load_final_mappings(mapping, {})
                pd.testing.assert_frame_equal(source, original_reactions, rtol=1e-12, atol=1e-15)
                guard.check("solving EMU source balances for " + name, force=True)
                tables, info = solve_source_contributions(reactions, metabolites, job["carbon_sources"])
                info.update(analysis_method="EMU", emu_scope="Size-1 EMU-inspired mean origin balances; no multi-carbon isotopomer simulation",
                            model=name, mapping_version="AI" if suffix else "Original", mapping_file=str(mapping), mapping_sha256=sha256(mapping.read_bytes()).hexdigest(),
                            source_input_sha256=digest, carbon_sources=job["carbon_sources"],
                            allow_inorganic_carbon_recycling=True, output_file=str(output),
                            generated_at_utc=datetime.utcnow().isoformat() + "Z",
                            execution="Deterministic numerical code; no LLM or mapper calls; each mapping version solved separately",
                            mapping_quality="Uses selected mappings, including unresolved review flags; no new biochemical review performed.",
                            metabolite_status_counts=tables["MetaboliteContributions"].Status.value_counts().to_dict())
                guard.check("exporting EMU source balances for " + name, force=True)
                workbook = Workbook()
                del workbook[workbook.active.title]
                try:
                    for sheet in ("MetaboliteContributions", "BoundaryInputs"):
                        table = tables[sheet]
                        if sheet == "MetaboliteContributions":
                            table = table.drop(columns=["FluxResidual", "Status"])
                        write_sheet(workbook, sheet, table)
                    save_output(workbook, output, max_bytes=workbook_limit(job.get("max_excel_size_mb", 50)))
                finally:
                    workbook.close()
                report["models"].append(info)
                print("{} EMU: {} carbon positions; {} complete metabolites; {}".format(
                    name, info["carbon_positions"], int(tables["MetaboliteContributions"].Status.eq("complete").sum()),
                    output.name), flush=True)
        report["status"] = "finished"
    except Exception as error:
        report.update(status="error", error=str(error))
        raise
    finally:
        report["finished_at_utc"] = datetime.utcnow().isoformat() + "Z"
        temporary = report_path.with_suffix(".tmp")
        temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        temporary.replace(report_path)
    return report["models"]
