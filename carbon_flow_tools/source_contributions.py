"""Fixed-flux, single-carbon origin balances; not a full isotopomer simulator."""

from collections import defaultdict, deque
import math
import warnings

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import MatrixRankWarning, spsolve

from .atom_helpers import mapped_components
from .consensus import mapping_signatures
from .carbon_mapping_rules import structure_kind, STRUCTURE_LIMITS, carbon_structure_issues
from .reactions import clean, get_smiles, parse_side, split_equation

def inorganic_carbon_names(smiles):
    """Recognize CO2 and carbonate protonation forms independently of model IDs."""
    structures = {"O=C=O", "O=C([O-])O", "O=C(O)O", "O=C([O-])[O-]"}
    names = set()
    for name, structure in smiles.items():
        molecule = Chem.MolFromSmiles(structure) if structure else None
        if molecule is not None:
            for atom in molecule.GetAtoms():
                atom.SetAtomMapNum(0)
            if Chem.MolToSmiles(molecule) in structures:
                names.add(name)
    return names


TOL = 1e-8
UNKNOWN = "Unknown"


def solve_source_contributions(reactions, metabolites, source_names,
                               allow_inorganic_carbon_recycling=True, tolerance=TOL,
                               memory_guard=None):
    """Balance carbon positions at metabolic and origin steady state.

    Fractions are unconditional: Unknown is retained, never redistributed.
    Only external inflows receive prescribed origins; internal source pools mix.
    Use the validated selected mapping as the positional representative.
    Graph symmetry is recorded, not converted into additional Unknown.
    """
    if type(allow_inorganic_carbon_recycling) is not bool:
        raise ValueError("EMU inorganic-carbon recycling must be True or False")
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("EMU numerical tolerance must be finite and positive")
    if memory_guard:
        memory_guard.check("preparing EMU carbon positions", force=True)
    terminal = set() if allow_inorganic_carbon_recycling else inorganic_carbon_names(
        {clean(r.MetaboliteName): get_smiles(r) for r in metabolites.itertuples(index=False)})
    if set(source_names) & terminal:
        raise ValueError("Enable EMU inorganic-carbon recycling for an inorganic carbon source")
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
    carbon_issues = carbon_structure_issues(metabolites)
    for name, row in rows_by_name.items():
        structure = get_smiles(row)
        kind = structure_kind(structure, clean(getattr(row, "MolecularFormula", "")))
        if kind == "polymer":
            structure_issues[name] = STRUCTURE_LIMITS[kind]
            continue
        molecule = Chem.MolFromSmiles(structure) if structure else None
        if molecule is None:
            structure_issues[name] = "Missing or invalid structure; carbon positions unavailable"
            continue
        if name in carbon_issues:
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
        row = rows_by_name.get(name)
        if row is not None and structure_kind(get_smiles(row), clean(getattr(row, "MolecularFormula", ""))) != "complete":
            raise ValueError("Carbon source requires a complete, explicit structure: " + name)
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
    if memory_guard:
        memory_guard.check("allocating EMU origin arrays", force=True, extra_bytes=n * len(origins) * 8 * 5)
    rhs = np.zeros((n, len(origins)))
    edge_rows, edge_cols, edge_rates = [], [], []
    diagnostics, symmetry_mapped = [], 0
    for reaction_number, (reaction, left, right, rate) in enumerate(directed):
        if memory_guard and reaction_number % 25 == 0:
            memory_guard.check("building EMU transfers")
        status, unknown_rate = "sink", 0.0
        if not left:
            status = "external_input"
            for coefficient, name in right:
                origin = name if name in source_names else "External:" + name
                for node in positions.get(name, []):
                    rhs[index[node], origin_index[origin]] += coefficient * rate
        elif right:
            status = "missing_mapping"
            transfers = []
            if reaction.status == "partial" and reaction.partial_maps:
                import json
                status = "partial_rule_mapping"
                for record in reaction.partial_maps:
                    amount = float(record["TransferCoefficient"])
                    for a, i, b, j in json.loads(record["CarbonAtomMap"]):
                        transfers.append(((record["ProductMetabolite"], b, j),
                                          (record["SourceMetabolite"], a, i), amount))
            elif reaction.status in ("mapped", "partial") and reaction.mapped_smiles:
                if reaction.status == "mapped":
                    mapping_signatures(reaction)
                components = mapped_components(reaction)
                if components is None:
                    raise ValueError("EMU mapping does not align with saved structures: " + reaction.reaction_id)
                status = "mapped" if reaction.status == "mapped" else "partial_mapper_mapping"
                source_by_map = {}
                for component, atom_maps in components[0]:
                    mol = Chem.MolFromSmiles(component.smiles)
                    for atom in mol.GetAtoms():
                        number = atom_maps[atom.GetIdx()]
                        if atom.GetAtomicNum() == 6 and number:
                            if number in source_by_map:
                                raise ValueError("Duplicated EMU carbon map ID: " + reaction.reaction_id)
                            source_by_map[number] = (component.metabolite, component.component_index, atom.GetIdx())
                for component, atom_maps in components[1]:
                    mol = Chem.MolFromSmiles(component.smiles)
                    for atom in mol.GetAtoms():
                        number = atom_maps[atom.GetIdx()]
                        if atom.GetAtomicNum() == 6 and number in source_by_map:
                            transfers.append(((component.metabolite, component.component_index, atom.GetIdx()),
                                              source_by_map[number], 1.0))
            # Keep full-network turnover; unresolved product positions receive Unknown.
            supplied, used = defaultdict(float), defaultdict(float)
            expected = defaultdict(float)
            available = defaultdict(float)
            for coefficient, name in right:
                for node in positions.get(name, []):
                    expected[node] += coefficient
            for coefficient, name in left:
                for node in positions.get(name, []):
                    available[node] += coefficient
            for target, source, amount in transfers:
                if target not in index:
                    continue
                supplied[target] += amount
                used[source] += amount
                if supplied[target] > expected[target] + tolerance:
                    raise ValueError("EMU product carbon is over-assigned: " + reaction.reaction_id)
                if source in index and used[source] > available[source] + tolerance:
                    raise ValueError("EMU precursor carbon is over-assigned: " + reaction.reaction_id)
                if source not in index or source[0] in terminal:
                    rhs[index[target], -1] += amount * rate
                    unknown_rate += amount * rate
                    status = "unavailable_precursor_or_stopped_inorganic_carbon"
                else:
                    edge_rows.append(index[target])
                    edge_cols.append(index[source])
                    edge_rates.append(amount * rate)
                    if target in ambiguous or source in ambiguous:
                        symmetry_mapped += 1
            for target, coefficient in expected.items():
                missing = max(0.0, coefficient - supplied[target])
                rhs[index[target], -1] += missing * rate
                unknown_rate += missing * rate
        diagnostics.append(dict(ReactionID=reaction.reaction_id, Flux=reaction.flux,
                                Status=status, UnknownCarbonInputRate=unknown_rate))

    q = np.array([incoming[node[0]] for node in nodes])
    balance = {name: incoming[name] - outgoing[name] for name in rows_by_name}
    invalid = {name for name in rows_by_name if abs(balance[name]) >
               1e-6 + 1e-7 * max(incoming[name], outgoing[name])}

    if memory_guard:
        memory_guard.check("building EMU sparse matrices", force=True,
            extra_bytes=len(edge_rates) * 256 + n * len(origins) * 8 * 5)
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
    if memory_guard:
        # Sparse LU fill-in can be much denser than the transfer matrix.
        memory_guard.check("solving EMU balances", force=True, extra_bytes=n * n * 16)
    with warnings.catch_warnings():
        warnings.simplefilter("error", MatrixRankWarning)
        fractions = np.asarray(spsolve(scaled, scaled_rhs)).reshape(n, len(origins))
    residual = float(np.max(np.abs(scaled.dot(fractions) - scaled_rhs)))
    if (not np.isfinite(fractions).all() or residual > tolerance
            or fractions.min() < -tolerance or fractions.max() > 1 + tolerance
            or np.max(np.abs(fractions.sum(axis=1) - 1)) > tolerance):
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
                      "complete" if means[-1] <= tolerance else "partial_unknown")
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
    if memory_guard:
        memory_guard.check("exporting EMU tables", force=True,
            extra_bytes=(len(atom_rows) + len(metabolite_rows)) * (len(origins) + 12) * 128)
    atoms = pd.DataFrame(atom_rows)
    compounds = pd.DataFrame(metabolite_rows)
    kinds = {name: structure_kind(get_smiles(row), clean(getattr(row, "MolecularFormula", "")))
             for name, row in rows_by_name.items()}
    if any(kind in ("partial", "polymer") for kind in kinds.values()):
        compounds["StructureType"] = compounds.MetaboliteName.map(kinds)
        compounds["KnownCarbonCount"] = compounds.CarbonCount
        compounds.loc[compounds.StructureType.eq("polymer"), "KnownCarbonCount"] = np.nan
        compounds["ContributionScope"] = compounds.StructureType.map(
            lambda kind: "explicit carbon positions only; unknown region not quantified" if kind == "partial"
            else "unavailable; repeat count unspecified" if kind == "polymer" else "whole explicit structure")
        compounds.loc[compounds.StructureType.isin(["partial", "polymer"]), "CarbonCount"] = np.nan
    info = dict(method="fixed_flux_single_carbon_origin_balance", max_scaled_residual=residual,
                carbon_positions=n, ambiguous_carbon_positions=len(ambiguous),
                symmetry_mapped_transfers=symmetry_mapped,
                symmetry_policy="Follow validated Final_Mapping representatives; graph symmetry does not add Unknown or imply random scrambling.",
                unanchored_positions=len(unanchored), invalid_balance_metabolites=sorted(invalid),
                unavailable_structure_metabolites=len(structure_issues),
                unknown_interpretation="Known fractions are contributions resolved by this conservative model. Unknown may also belong to any source; known fractions are not renormalized.",
                configured_sources_without_import=sorted(set(source_names) - {x[3] for x in imports}),
                assumptions="Metabolic and origin steady state; well-mixed pools; net-direction fluxes; no isotope effects; CO2 refixation follows the EMU setting; Validated selected atom correspondences are positional representatives, including graph-equivalent atoms; no extra symmetry scrambling is assumed.",
                allow_inorganic_carbon_recycling=allow_inorganic_carbon_recycling,
                numerical_tolerance=tolerance,
                scope="Mean carbon-origin fractions, not molecular source combinations, experimental measurements or source essentiality.",
                position_convention="CarbonPosition is 1-based carbon order in saved SMILES; FragmentIndex and AtomIndex are 0-based, not biochemical numbering.")
    return dict(AtomContributions=atoms, MetaboliteContributions=compounds,
                ReactionDiagnostics=pd.DataFrame(diagnostics),
                BoundaryInputs=pd.DataFrame(imports, columns=["ReactionID", "MetaboliteName", "MoleculeInputRate", "Source"])), info


