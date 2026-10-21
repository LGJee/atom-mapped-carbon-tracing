"""Propagate source labels for concrete substrate states."""

from collections import defaultdict
from itertools import product

import pandas as pd
from rdkit import Chem

from .atom_helpers import (
    glucose_labels,
    isoprenol_labels,
    labeled_smiles,
    mapped_components,
    assignment_from_labeled_smiles,
)
from .carbon_mapping_rules import partial_reactions
from .reactions import clean, get_smiles, parse_side, split_equation

COLUMNS = [
    "Step",
    "ReactionID",
    "Source",
    "SourceIDs",
    "Production",
    "ProductionID",
    "CarbonLabels",
    "LabeledSMILES",
]


def group_components(terms, components):
    """Separate molecular copies according to integer stoichiometry."""
    groups, cursor = [], 0
    for term_index, (coefficient, name) in enumerate(terms):
        count = (
            max(c.component_index for c, _ in components if c.metabolite == name) + 1
        )
        for _ in range(int(round(coefficient))):
            group = components[cursor : cursor + count]
            assert len(group) == count and all(c.metabolite == name for c, _ in group)
            groups.append((term_index, name, group))
            cursor += count
    assert cursor == len(components)
    return groups


def map_one_reaction(
    reaction,
    components,
    states,
    step,
    catalog,
    next_id,
    smiles_by_name,
    first_steps=None,
    required_labels=None,
    unordered_sources=False,
    memory_guard=None,
    record_sink=None,
):
    if memory_guard:
        memory_guard.check(reaction.reaction_id)
    left, right = components
    left_groups = group_components(parse_side(reaction.source_text), left)
    sides = split_equation(reaction.mapping_equation or reaction.equation)
    right_groups = group_components(
        parse_side(sides[1 if reaction.mapping_equation or reaction.flux >= 0 else 0]), right
    )
    choices = [
        sorted(states.get(name, {}).items()) or [None] for _, name, _ in left_groups
    ]
    rows, produced = [], []
    for index, selected in enumerate(product(*choices)):
        if memory_guard and index % 256 == 0:
            memory_guard.check(reaction.reaction_id)
        if not any(choice is not None for choice in selected):
            continue
        if required_labels is not None and not any(
            label in required_labels
            for choice in selected
            if choice
            for _, _, label in choice[1]
        ):
            continue
        label_by_map, by_term = {}, defaultdict(list)
        used_steps = []
        for (term, name, group), choice in zip(left_groups, selected):
            by_term[term].append(0 if choice is None else choice[0])
            if choice is None:
                continue
            if first_steps is not None:
                used_steps.append(first_steps[(name, choice[0])])
            assignment = {(part, atom): label for part, atom, label in choice[1]}
            for component, atom_maps in group:
                for atom, number in atom_maps.items():
                    label = assignment.get((component.component_index, atom))
                    if label and number > 0:
                        label_by_map[number] = label

        # Copy order affects asymmetric products; enumerate it before deduplication.
        source_ids = " + ".join(
            "&".join(map(str, sorted(by_term[i]) if unordered_sources else by_term[i]))
            for i in sorted(by_term)
        )
        actual_step = (
            max(used_steps, default=0) + 1 if first_steps is not None else step
        )
        for _, name, group in right_groups:
            assignment = tuple(
                sorted(
                    (component.component_index, atom, label_by_map[number])
                    for component, atom_maps in group
                    for atom, number in atom_maps.items()
                    if number > 0 and number in label_by_map
                )
            )
            if not assignment:
                continue
            marked = labeled_smiles(smiles_by_name[name], assignment)
            key = (name, marked)
            if key not in catalog:
                catalog[key] = next_id
                next_id += 1
            production_id = catalog[key]

            row = dict(zip(COLUMNS, [
                actual_step, reaction.parent_reaction_id or reaction.reaction_id, reaction.source_text, source_ids,
                name, production_id,
                ",".join(sorted((a[2] for a in assignment), key=lambda x: int(x[1:]))),
                marked,
            ]))
            if record_sink is None:
                restored = assignment_from_labeled_smiles(smiles_by_name[name], marked)
                produced.append((name, production_id, restored))
                rows.append(row)
            else:
                record_sink(row)
    return produced, rows, next_id


def initialize_sources(source_names, smiles):
    """Allocate non-overlapping origin labels in source-list order."""
    if not source_names or len(set(source_names)) != len(source_names):
        raise ValueError(
            "Carbon sources must be nonempty and compartment-qualified names unique."
        )
    sources, offset = [], 0
    rules = {"D-glucose": glucose_labels, "isoprenol": isoprenol_labels}
    for name in source_names:
        structure = smiles.get(name, "")
        molecule = Chem.MolFromSmiles(structure) if structure else None
        if molecule is None:
            raise ValueError(f"Carbon source has no valid SMILES: {name}")
        rule = rules.get(name.rsplit("[", 1)[0])
        if rule is not None:
            local = rule(structure)
            atoms = sorted(local, key=lambda atom: int(local[atom][1:]))
        else:
            atoms = [
                atom.GetIdx()
                for atom in molecule.GetAtoms()
                if atom.GetAtomicNum() == 6
            ]
            print(
                f"{name}: origin labels follow carbon atom order in the input SMILES.",
                flush=True,
            )
        if not atoms:
            raise ValueError(f"Selected source contains no carbon: {name}")
        assignment = tuple(
            (0, atom, f"C{offset + i}") for i, atom in enumerate(atoms, 1)
        )

        if "." in structure:
            raise ValueError(f"Provide a single-component carbon source SMILES: {name}")
        labels = {label for _, _, label in assignment}
        sources.append((name, assignment, labels))
        offset += len(atoms)
        print(f"Carbon source {name}: C{offset - len(atoms) + 1}閳ユ弲{offset}", flush=True)
    return sources


def prepare_mappings(reactions):
    prepared = []
    used_ids = {reaction.reaction_id for reaction in reactions}
    for reaction in reactions:
        if reaction.status == "partial" and reaction.partial_maps:
            for unit, components in partial_reactions(reaction):
                while unit.reaction_id in used_ids:
                    unit.reaction_id += ":"
                used_ids.add(unit.reaction_id)
                prepared.append((unit, components))
            continue
        if reaction.status not in ("mapped", "partial"):
            continue
        components = mapped_components(reaction)
        if components is None or any(
            not maps for side in components for _, maps in side
        ):
            reaction.status = "Mapped fragments cannot be aligned to input metabolites"
        else:
            prepared.append((reaction, components))
    return prepared


def propagate_once(
    prepared,
    states,
    first_steps,
    catalog,
    next_id,
    smiles,
    required_labels=None,
    unordered_sources=False,
    progress=None,
    memory_guard=None,
    record_sink=None,
):
    """Read one state snapshot; return products for the caller to merge."""
    produced, rows = [], []
    for index, (reaction, components) in enumerate(prepared, 1):
        new_states, new_rows, next_id = map_one_reaction(
            reaction,
            components,
            states,
            0,
            catalog,
            next_id,
            smiles,
            first_steps,
            required_labels,
            unordered_sources,
            memory_guard,
            record_sink=record_sink,
        )
        produced.extend(new_states)
        rows.extend(new_rows)
        del new_states, new_rows
        if progress and index % 25 == 0:
            print(
                f"{progress}: {index}/{len(prepared)} reactions; "
                + ("records staged on disk" if record_sink is not None else f"{len(rows)} records"),
                flush=True,
            )
    if record_sink is not None:
        return None, None, next_id
    if memory_guard:
        memory_guard.check("building propagation table", force=True,
                           extra_bytes=len(rows) * (len(COLUMNS) * 8 + 64) * 3)
    table = pd.DataFrame(rows, columns=COLUMNS)
    del rows
    if memory_guard:
        memory_guard.check("deduplicating propagation table", force=True,
                           extra_bytes=int(table.memory_usage(index=True).sum()) * 3)
    table = table.drop_duplicates(ignore_index=True)
    if memory_guard:
        memory_guard.check("sorting propagation table", force=True,
                           extra_bytes=int(table.memory_usage(index=True).sum()) * 3)
    table = table.sort_values(
        ["Step", "ReactionID", "SourceIDs", "ProductionID"], ignore_index=True
    )
    return produced, table, next_id


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


def merge_states(
    produced, table, states, first_steps, terminal_names=(), memory_guard=None
):
    for index, (name, number, assignment) in enumerate(produced):
        if memory_guard and index % 1024 == 0:
            memory_guard.check("merging propagation states")
        # Keep terminal products in the output/catalog, not as labelled substrates.
        if name not in terminal_names:
            states[name].setdefault(number, assignment)
    for index, row in enumerate(table.itertuples(index=False)):
        if memory_guard and index % 1024 == 0:
            memory_guard.check("merging propagation steps")
        key = (row.Production, row.ProductionID)
        first_steps[key] = min(first_steps.get(key, row.Step), row.Step)


def propagate_rounds(
    strain,
    reactions,
    metabolites,
    rounds,
    source_names,
    sources=None,
    allow_inorganic_carbon_recycling=False,
    memory_guard=None,
    row_store=None,
):
    """Propagate by source and layer in round one, then by round snapshots."""
    if memory_guard:
        memory_guard.check("starting propagation", force=True)
    smiles = {
        clean(r.MetaboliteName): get_smiles(r)
        for r in metabolites.itertuples(index=False)
    }
    if type(allow_inorganic_carbon_recycling) is not bool:
        raise ValueError("allow_inorganic_carbon_recycling must be True or False.")
    terminal_names = (
        set() if allow_inorganic_carbon_recycling else inorganic_carbon_names(smiles)
    )
    sources = initialize_sources(source_names, smiles) if sources is None else sources
    if any(name in terminal_names for name, _, _ in sources):
        raise ValueError(
            "For an inorganic carbon source, set ALLOW_INORGANIC_CARBON_RECYCLING = True "
            "at the top of step3_carbon_flow.py."
        )
    print(
        "Inorganic-carbon label recycling: "
        + (
            "enabled."
            if allow_inorganic_carbon_recycling
            else "disabled; CO2/bicarbonate labels terminate."
        ),
        flush=True,
    )
    prepared = prepare_mappings(reactions)
    states, catalog, first_steps = defaultdict(dict), {}, {}
    next_id, phase_tables = 1, []
    if row_store is not None:
        row_store.start_round(1)
    for phase, (seed, assignment, labels) in enumerate(sources):
        states[seed][0] = assignment
        first_steps[(seed, 0)] = 0
        pending, tables, layer = list(prepared), [], 0
        required = None if phase == 0 else labels
        while pending:
            if memory_guard:
                memory_guard.check("preparing first-round layer", force=True)
            ready = []
            for item in pending:
                available = [
                    a
                    for _, name in parse_side(item[0].source_text)
                    for a in states.get(name, {}).values()
                ]
                if available and (
                    required is None
                    or any(label in required for a in available for _, _, label in a)
                ):
                    ready.append(item)
            if not ready:
                break
            ready_ids = {r.reaction_id for r, _ in ready}
            pending = [item for item in pending if item[0].reaction_id not in ready_ids]
            if row_store is not None:
                row_store.start_layer(phase)
            produced, table, next_id = propagate_once(
                ready, states, first_steps, catalog, next_id, smiles, required,
                memory_guard=memory_guard,
                record_sink=row_store.append if row_store is not None else None,
            )
            if row_store is None:
                merge_states(produced, table, states, first_steps, terminal_names, memory_guard)
                tables.append(table)
                layer_records = len(table)
            else:
                row_store.merge_states(states, first_steps, terminal_names, smiles)
                layer_records = row_store.db.execute(
                    "SELECT COUNT(*) FROM records WHERE round_no=1 AND layer=?",
                    (row_store.layer,),
                ).fetchone()[0]
            del produced
            layer += 1
            print(
                f"{strain} {seed} layer {layer}: {len(ready)} reactions,{layer_records} records",
                flush=True,
            )
        if row_store is None:
            if memory_guard:
                memory_guard.check("combining first-round layers", force=True)
            phase_table = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(columns=COLUMNS)
            phase_tables.append(phase_table.drop_duplicates(ignore_index=True).sort_values(
                ["Step", "ReactionID", "SourceIDs", "ProductionID"], ignore_index=True))
    if row_store is None:
        if memory_guard:
            memory_guard.check("combining first-round sources", force=True)
        first = pd.concat(phase_tables, ignore_index=True)
        del phase_tables, tables, phase_table
    else:
        first = None
    yield 1, first, len(catalog), False
    del first
    table = None
    for number in range(2, rounds + 1):
        if memory_guard:
            memory_guard.check("starting round {}".format(number), force=True)
        if row_store is not None:
            row_store.start_round(number)
            row_store.start_layer()
        before, old_steps = len(catalog), dict(first_steps)
        produced, table, next_id = propagate_once(
            prepared,
            states,
            first_steps,
            catalog,
            next_id,
            smiles,
            unordered_sources=True,
            progress=f"{strain} R{number}",
            memory_guard=memory_guard,
            record_sink=row_store.append if row_store is not None else None,
        )
        if row_store is None:
            merge_states(produced, table, states, first_steps, terminal_names, memory_guard)
        else:
            row_store.merge_states(states, first_steps, terminal_names, smiles)
        converged = len(catalog) == before and first_steps == old_steps
        del produced, old_steps
        yield number, table, len(catalog) - before, converged
        del table
        if converged:
            print(
                f"{strain} R{number}: no new states or earlier steps; propagation converged.",
                flush=True,
            )
            break
