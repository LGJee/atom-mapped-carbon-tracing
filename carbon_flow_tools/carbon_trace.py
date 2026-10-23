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
from .console import Progress, detail
from .enumeration import candidate_rows
from .propagation_state import ordered_product
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


def reaction_transfer(reaction, components):
    """Compile fixed molecular-copy/atom transfers, never growing label states."""
    signature = (reaction.source_text, reaction.equation, reaction.mapping_equation,
        reaction.flux >= 0, tuple(tuple((c.metabolite, c.component_index, c.smiles,
            tuple(maps.items())) for c, maps in side) for side in components))
    if getattr(reaction, "_transfer_signature", None) != signature:
        sides = split_equation(reaction.mapping_equation or reaction.equation)
        groups = (group_components(parse_side(reaction.source_text), components[0]),
                  group_components(parse_side(sides[1 if reaction.mapping_equation or reaction.flux >= 0 else 0]), components[1]))
        reaction._transfer_plan = tuple(tuple((term, name, tuple(
            (c.component_index, atom, number) for c, maps in group
            for atom, number in maps.items() if number > 0))
            for term, name, group in side) for side in groups)
        reaction._transfer_signature = signature
    return reaction._transfer_plan


def projected_choices(left_groups, states, first_steps, required_labels):
    """Project each substrate state once per molecular copy, before combinations."""
    views, choices = {}, []
    for _, name, transfers in left_groups:
        if name not in views:
            views[name] = [(number, {(part, atom): label for part, atom, label in assignment},
                first_steps[(name, number)] if first_steps is not None else None,
                required_labels is None or any(label in required_labels for _, _, label in assignment))
                for number, assignment in sorted(states.get(name, {}).items())]
        choices.append([(number, tuple((mapped, assignment[(part, atom)])
            for part, atom, mapped in transfers if assignment.get((part, atom))),
            first, required) for number, assignment, first, required in views[name]] or [None])
    return choices


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
    progress=None,
    progress_label="",
    candidate_pool=None,
    compact_sink=None,
):
    if memory_guard:
        memory_guard.check(reaction.reaction_id)
    left_groups, right_groups = reaction_transfer(reaction, components)
    disk_states = hasattr(states, "choices")
    choices = (states.choices(left_groups, required_labels) if disk_states else
               projected_choices(left_groups, states, first_steps, required_labels))
    combinations = 1
    for options in choices:
        combinations *= len(options)
    if progress is not None:
        progress.begin(progress_label + " | current reaction " + reaction.reaction_id,
                       combinations)
    parallel = candidate_pool is not None and candidate_pool.eligible(combinations)
    def checked_selections():
        selections = ordered_product(choices) if disk_states else product(*choices)
        for index, selected in enumerate(selections):
            if index % 256 == 0:
                if progress is not None and not parallel:
                    progress.update(index)
                if memory_guard:
                    memory_guard.check(reaction.reaction_id)
            yield selected

    args = (left_groups, right_groups, checked_selections(), smiles_by_name,
            step, first_steps is not None, required_labels is not None, unordered_sources)
    candidates = (candidate_pool.rows(*args, memory_guard=memory_guard, progress=progress,
                                     compact=compact_sink is not None)
                  if parallel else candidate_rows(*args, compact=compact_sink is not None))
    rows, produced = [], []
    try:
        for actual_step, source_ids, name, carbon_labels, marked in candidates:
            key = (name, marked)
            if key not in catalog:
                catalog[key] = next_id
                next_id += 1
            production_id = catalog[key]
            if compact_sink is not None:
                compact_sink(name, production_id, marked, carbon_labels, actual_step)
                continue
            row = dict(zip(COLUMNS, [actual_step,
                reaction.parent_reaction_id or reaction.reaction_id, reaction.source_text,
                source_ids, name, production_id, carbon_labels, marked]))
            if record_sink is None:
                restored = assignment_from_labeled_smiles(smiles_by_name[name], marked)
                produced.append((name, production_id, restored))
                rows.append(row)
            else:
                record_sink(row)
    finally:
        candidates.close()
    if progress is not None:
        progress.update(combinations)
    return produced, rows, next_id


def initialize_sources(source_names, smiles):
    """Allocate non-overlapping origin labels in source-list order."""
    if not source_names or len(set(source_names)) != len(source_names):
        raise ValueError(
            "Carbon sources must be nonempty and compartment-qualified names unique."
        )
    from .carbon_mapping_rules import structure_kind
    sources, offset = [], 0
    rules = {"D-glucose": glucose_labels, "isoprenol": isoprenol_labels}
    for name in source_names:
        structure = smiles.get(name, "")
        if structure_kind(structure) != "complete":
            raise ValueError("Carbon source requires a complete, explicit structure: " + name)
        molecule = Chem.MolFromSmiles(structure)
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
            detail(
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
        detail(f"Carbon source {name}: C{offset - len(atoms) + 1}-C{offset}", flush=True)
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
    candidate_pool=None,
    compact_sink=None,
):
    """Read one state snapshot; return products for the caller to merge."""
    tracker = progress if isinstance(progress, Progress) else (
        Progress(progress, unit="combinations") if progress else None)
    context = tracker.label if tracker is not None else ""
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
            progress=tracker,
            progress_label="{} | reaction {}/{}".format(context, index, len(prepared)),
            candidate_pool=candidate_pool,
            compact_sink=compact_sink,
        )
        produced.extend(new_states)
        rows.extend(new_rows)
        del new_states, new_rows
        if progress and index % 25 == 0:
            detail(
                f"{context}: {index}/{len(prepared)} reactions; "
                + ("records staged on disk" if record_sink is not None else f"{len(rows)} records"),
                flush=True,
            )
    if record_sink is not None or compact_sink is not None:
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
    candidate_pool=None,
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
    from .carbon_mapping_rules import structure_kind
    for row in metabolites.itertuples(index=False):
        if clean(row.MetaboliteName) in source_names and structure_kind(
                get_smiles(row), clean(getattr(row, "MolecularFormula", ""))) != "complete":
            raise ValueError("Carbon source requires a complete, explicit structure: " + clean(row.MetaboliteName))
    sources = initialize_sources(source_names, smiles) if sources is None else sources
    if any(name in terminal_names for name, _, _ in sources):
        raise ValueError(
            "For an inorganic carbon source, set ALLOW_INORGANIC_CARBON_RECYCLING = True "
            "at the top of step3_carbon_flow.py."
        )
    detail(
        "Inorganic-carbon label recycling: "
        + (
            "enabled."
            if allow_inorganic_carbon_recycling
            else "disabled; CO2/bicarbonate labels terminate."
        ),
        flush=True,
    )
    prepared = prepare_mappings(reactions)
    states, catalog, first_steps = ((row_store.states, row_store.catalog, row_store.first_steps)
                                   if row_store is not None else (defaultdict(dict), {}, {}))
    next_id, phase_tables = 1, []
    tracker = Progress(strain, unit="combinations")
    print("{}: calculating round 1 of {}...".format(strain, rounds), flush=True)
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
                available = (a for _, name in parse_side(item[0].source_text)
                             for a in states.get(name, {}).values())
                if any(required is None or any(label in required for _, _, label in a)
                       for a in available):
                    ready.append(item)
            if not ready:
                break
            ready_ids = {r.reaction_id for r, _ in ready}
            pending = [item for item in pending if item[0].reaction_id not in ready_ids]
            if row_store is not None:
                row_store.start_layer(phase)
            tracker.label = "{} | R1/{} | {} | layer {}".format(strain, rounds, seed, layer + 1)
            produced, table, next_id = propagate_once(
                ready, states, first_steps, catalog, next_id, smiles, required,
                progress=tracker, memory_guard=memory_guard,
                record_sink=row_store.append if row_store is not None else None,
                candidate_pool=candidate_pool,
                compact_sink=row_store.append_state if row_store is not None and row_store.compact else None,
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
            detail(
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
        tracker.label = "{} | R{}/{}".format(strain, number, rounds)
        print("{}: calculating round {} of {}..."
              .format(strain, number, rounds), flush=True)
        if memory_guard:
            memory_guard.check("starting round {}".format(number), force=True)
        if row_store is not None:
            row_store.start_round(number)
            row_store.start_layer()
        before = len(catalog)
        # Steps only appear or decrease; a revision detects changes without copying the table.
        old_steps = first_steps.revision if row_store is not None else dict(first_steps)
        produced, table, next_id = propagate_once(
            prepared,
            states,
            first_steps,
            catalog,
            next_id,
            smiles,
            unordered_sources=True,
            progress=tracker,
            memory_guard=memory_guard,
            record_sink=row_store.append if row_store is not None else None,
            candidate_pool=candidate_pool,
            compact_sink=row_store.append_state if row_store is not None and row_store.compact else None,
        )
        if row_store is None:
            merge_states(produced, table, states, first_steps, terminal_names, memory_guard)
        else:
            row_store.merge_states(states, first_steps, terminal_names, smiles)
        converged = len(catalog) == before and (
            first_steps.revision == old_steps if row_store is not None else first_steps == old_steps)
        del produced, old_steps
        yield number, table, len(catalog) - before, converged
        del table
        if converged:
            print(
                f"{strain}: no further changes after round {number}; calculation complete.",
                flush=True,
            )
            break
