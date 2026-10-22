"""Select mapper correspondences; always retain a nonempty generated candidate."""

from collections import Counter, defaultdict
from copy import deepcopy

import pandas as pd
from rdkit import Chem

from .atom_helpers import mapped_components
from .mapper_runner import mapping_error
from .reactions import clean, parse_side, split_equation


def mapping_signatures(reaction):
    """Compare structurally aligned correspondences, ignoring IDs and equivalent symmetry."""
    from .carbon_mapping_rules import structure_kind
    error = mapping_error(reaction.mapped_smiles)
    if error:
        raise ValueError(error)
    sides = mapped_components(reaction)
    if sides is None:
        raise ValueError(
            "Mapped structures do not match input identity, charge or stereochemistry"
        )
    molecule = Chem.Mol()
    positions, carbon_maps, colors = [], [], []
    for side_index, side in enumerate(sides):
        side_maps, carbons = {}, set()
        for component, atom_maps in side:
            if structure_kind(component.smiles) != "complete":
                raise ValueError("Partial structure: wildcard groups have no established atom correspondence")
            part = Chem.MolFromSmiles(component.smiles)
            if len(atom_maps) != part.GetNumAtoms():
                raise ValueError("Input atoms were not completely aligned")
            offset = molecule.GetNumAtoms()
            for atom in part.GetAtoms():
                number = atom_maps[atom.GetIdx()]
                if atom.GetAtomicNum() == 6:
                    if not number:
                        raise ValueError("Unmapped carbon atom")
                    carbons.add(number)
                if number:
                    side_maps[number] = (offset + atom.GetIdx(), atom.GetAtomicNum())
                colors.append(
                    (
                        side_index,
                        component.metabolite,
                        atom.GetIsotope(),
                        atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else "",
                    )
                )
                atom.SetAtomMapNum(0)

                atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            molecule = Chem.CombineMols(molecule, part)
        positions.append(side_maps)
        carbon_maps.append(carbons)
    if carbon_maps[0] != carbon_maps[1]:
        raise ValueError("Carbon correspondence is not complete across the two sides")
    color_ids = {color: i for i, color in enumerate(sorted(set(colors)), 1)}
    for atom, color in zip(molecule.GetAtoms(), colors):
        atom.SetIsotope(color_ids[color])

    def signature(carbon_only):
        graph = Chem.RWMol(molecule)
        for number in positions[0].keys() & positions[1].keys():
            left, element = positions[0][number]
            right, _ = positions[1][number]
            if element == 6 or (not carbon_only and element > 1):
                graph.AddBond(left, right, Chem.BondType.ZERO)
        return Chem.MolToSmiles(graph.GetMol(), canonical=True, isomericSmiles=True)

    return signature(True), signature(False)


def score_type(backend, reaction):
    if reaction.is_transport:
        return "Identity transport; prediction tools were not called"
    if backend == "LocalMapper":
        return "Boolean template verification flag; not a correctness probability"
    if backend == "RXNMapper":
        return "Attention confidence; not calibrated on this dataset"
    if backend == "EnzymeMap":
        return (
            "Quality: relative rule frequency; not a correctness probability"
            if reaction.template.startswith("EnzymeMap CSV")
            else "Rule or identity mapping; no numeric confidence score"
        )
    return "Uncalibrated native tool score"


def preferred_backend(backends, predictions, index):
    """Use explicit tool priority, not a cross-tool numeric score."""

    def priority(backend):
        return {"EnzymeMap": 0, "RXNMapper": 1, "LocalMapper": 2}.get(backend, 3)

    return min(
        backends,
        key=lambda backend: (priority(backend), list(predictions).index(backend)),
    )


def carbon_candidate(reaction):
    """Keep only unique C-to-C correspondences in an exactly aligned mapper result."""
    view = deepcopy(reaction)
    sides = mapped_components(view)
    if sides is None:
        raise ValueError("Mapper structures do not match input identity, charge or stereochemistry")
    counts = []
    for side in sides:
        numbers = []
        for component, maps in side:
            molecule = Chem.MolFromSmiles(component.smiles)
            numbers.extend(number for index, number in maps.items()
                           if number and molecule.GetAtomWithIdx(index).GetAtomicNum() == 6)
        counts.append(Counter(numbers))
    shared = {n for n in counts[0] if counts[0][n] == counts[1][n] == 1}
    if not shared:
        raise ValueError("Mapper has no unique carbon-to-carbon correspondence")
    marked = []
    for text in reaction.mapped_smiles.split(">>"):
        mol = Chem.MolFromSmiles(text)
        if mol is None:
            raise ValueError("Invalid mapper SMILES")
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6 or atom.GetAtomMapNum() not in shared:
                atom.SetAtomMapNum(0)
        marked.append(Chem.MolToSmiles(mol))
    view.mapped_smiles = ">>".join(marked)
    view.status = "partial"
    try:
        mapping_signatures(view)
        view.status = "mapped"
    except ValueError:
        pass
    return view


def select_mappings(predictions, invalid_structures=None):
    """Prefer usable carbon consensus; retain failed candidates for review, never blank them."""
    from .carbon_mapping_rules import structure_kind
    if not predictions:
        raise ValueError("At least one mapper result is required")
    reference = next(iter(predictions.values()))
    ids = [reaction.reaction_id for reaction in reference]
    if len(set(ids)) != len(ids):
        raise ValueError("ReactionID must be unique")
    for backend, reactions in predictions.items():
        if len(reactions) != len(reference) or any(
            (r.reaction_id, r.reaction_smiles, r.equation, r.flux)
            != (base.reaction_id, base.reaction_smiles, base.equation, base.flux)
            for r, base in zip(reactions, reference)
        ):
            raise ValueError(backend + " has inconsistent reaction inputs or direction")
    selected, records = [], []
    for index, base in enumerate(reference):
        try:
            invalid = {name: invalid_structures[name] for side in split_equation(base.equation)
                       for _, name in parse_side(side) if name in (invalid_structures or {})}
        except ValueError:
            invalid = {}  # Preparation already records malformed equations for review.
        input_error = ("Input carbon formula/structure mismatch: " + "; ".join(
                       name + " (" + detail + ")" for name, detail in sorted(invalid.items()))) if invalid else ""
        partial_input = any(structure_kind(component.smiles) == "partial"
                            for component in base.left + base.right)
        groups, heavy_keys, valid = defaultdict(list), {}, []
        candidates, partial_candidates, carbon_views = {}, {}, set()
        record = {
            "ReactionID": base.reaction_id,
            "ReactionName": base.reaction_name,
            "Flux": base.flux,
            "Equation": base.equation,
            "MappingDirection": "forward" if base.flux >= 0 else "reverse",
        }
        for backend, reactions in predictions.items():
            reaction = reactions[index]
            error = ""
            if invalid:
                error = input_error
            elif reaction.is_transport and not partial_input:
                pass
            elif reaction.status == "mapped":
                try:
                    carbon_key, heavy_keys[backend] = mapping_signatures(reaction)
                    groups[carbon_key].append(backend)
                    valid.append(backend)
                except ValueError as exc:
                    error = str(exc)
            else:
                error = reaction.status
            if not invalid and (not reaction.is_transport or partial_input) and reaction.mapped_smiles:
                if not error:
                    candidates[backend] = reaction
                else:
                    try:
                        view = carbon_candidate(reaction)
                        if view.status == "mapped":
                            carbon_key, _ = mapping_signatures(view)
                            groups[carbon_key].append(backend)
                            valid.append(backend)
                            candidates[backend] = view
                            carbon_views.add(backend)
                        else:
                            partial_candidates[backend] = view
                        error += "; usable " + view.status + " carbon-only correspondence retained"
                    except ValueError:
                        pass
            reaction.comparison_error = error
            record[backend + "_Score"] = reaction.confident or None
            if reaction.confident is False or reaction.confident == 0:
                record[backend + "_Score"] = reaction.confident
            record[backend + "_ScoreType"] = score_type(backend, reaction)
            record[backend + "_Validation"] = (
                "Identity transport; prediction tools were not called"
                if reaction.is_transport and not invalid and not partial_input
                else error or "Structure and carbon mapping checks passed"
            )
        final = deepcopy(base)
        final.mapped_smiles, final.confident, final.template = "", "", ""
        winner, supporters, review, reason = (
            "",
            [],
            False,
            "No mapping available for carbon propagation",
        )
        if base.is_transport and not invalid and not partial_input:
            final = deepcopy(base)
            winner = "Identity"
            reason = "Pure transport: unchanged identity and stoichiometry; use identity correspondence"
        elif valid:
            support = max(map(len, groups.values()))
            tied = [
                backend
                for group in groups.values()
                if len(group) == support
                for backend in group
            ]
            winner = preferred_backend(tied, predictions, index)
            supporters = next(group for group in groups.values() if winner in group)
            review = support < 2
            if support >= 2:
                reason = "{} tools agree on carbon mapping; selected directly".format(
                    support
                )
            elif len(groups) > 1:
                reason = "Carbon mapping conflict; selected {} by EnzymeMap > RXNMapper > LocalMapper priority; review required".format(winner)
            else:
                reason = "Only {} provides a complete usable carbon mapping; review required".format(winner)
            final = deepcopy(candidates[winner])
            final.status = "mapped"
            if winner in carbon_views:
                review = True
                reason += "; selected a validated carbon-only view; review required"
        elif partial_candidates:
            winner = preferred_backend(partial_candidates, predictions, index)
            final = deepcopy(partial_candidates[winner])
            review = True
            reason = "No complete candidate; retain partial carbon correspondences from {}; EnzymeMap > RXNMapper > LocalMapper priority; review required".format(winner)
        else:
            generated = [name for name, rows in predictions.items() if rows[index].mapped_smiles]
            review = True
            if generated:
                winner = preferred_backend(generated, predictions, index)
                final = deepcopy(predictions[winner][index])
                final.status = "candidate"
                reason = (
                    "No usable carbon correspondence passed validation; retain {} nonempty prediction by "
                    "EnzymeMap > RXNMapper > LocalMapper priority; review required; "
                    "excluded from propagation: {}"
                ).format(winner, predictions[winner][index].comparison_error or input_error or "Unverified mapper result")
            else:
                reason = input_error or reason
                final.status = "Not selected: " + reason
        if final.omitted_components:
            review = True
            if final.status == "mapped" and final.mapped_smiles:
                try:
                    final = carbon_candidate(final)
                    final.status = "partial"
                except ValueError:
                    final.status = "candidate"  # No usable carbon correspondence in the retained view.
            reason += (
                "; known-structure mapper input omits missing-SMILES participants: "
                + "; ".join("{Side}: {Coefficient:g} {MetaboliteName}".format(**item)
                            for item in final.omitted_components)
                + "; partial original-reaction coverage; optional AI review must assess applicability")
        if partial_input:
            review = True
            final.is_transport = False
            reason += "; explicit carbon positions only; wildcard groups are unresolved, not paired by their symbols"
        record.update(
            ValidMappers="; ".join(valid),
            SupportingMappers="; ".join(supporters),
            CarbonMappingGroups=len(groups),
            HeavyAtomAgreement=(
                "identity"
                if base.is_transport and not invalid and not partial_input
                else "not comparable"
                if carbon_views
                else "agreement"
                if len(set(heavy_keys.values())) == 1 and len(valid) > 1
                else "disagreement"
                if len(set(heavy_keys.values())) > 1
                else "not comparable"
            ),
            SelectedMapper=winner,
            ReviewRequired=review,
            SelectionReason=reason,
            SelectedMappedReactionSMILES=final.mapped_smiles,
        )
        selected.append(final)
        records.append(record)
    table = pd.DataFrame(records)
    first = [
        "ReactionID",
        "ReactionName",
        "Flux",
        "SelectedMapper",
        "ReviewRequired",
        "SelectionReason",
        "Equation",
        "MappingDirection",
        "ValidMappers",
        "SupportingMappers",
        "CarbonMappingGroups",
        "HeavyAtomAgreement",
    ]
    return selected, table[first + [column for column in table if column not in first]]
