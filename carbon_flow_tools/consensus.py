"""Compare candidate correspondences and select one complete original mapping."""

from collections import defaultdict
from copy import deepcopy

import pandas as pd
from rdkit import Chem

from .atom_helpers import mapped_components
from .mapper_runner import mapping_error
from .reactions import clean


def mapping_signatures(reaction):
    """Compare structurally aligned correspondences, ignoring IDs and equivalent symmetry."""
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
        score = predictions[backend][index].confident
        if backend == "EnzymeMap":
            return 0
        if backend == "LocalMapper" and clean(score).lower() in ("true", "1", "1.0"):
            return 1
        if backend == "RXNMapper":
            return 2
        if backend == "LocalMapper":
            return 3
        return 4

    return min(
        backends,
        key=lambda backend: (priority(backend), list(predictions).index(backend)),
    )


def select_mappings(predictions):
    """Select a carbon consensus; fewer than two supporters require review."""
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
        groups, heavy_keys, valid = defaultdict(list), {}, []
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
            if reaction.is_transport:
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
            reaction.comparison_error = error
            record[backend + "_Score"] = reaction.confident or None
            if reaction.confident is False or reaction.confident == 0:
                record[backend + "_Score"] = reaction.confident
            record[backend + "_ScoreType"] = score_type(backend, reaction)
            record[backend + "_Validation"] = (
                "Identity transport; prediction tools were not called"
                if reaction.is_transport
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
        if base.is_transport:
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
                reason = "Carbon mapping conflict; use documented tool priority; review required"
            else:
                reason = "Only one tool provides a usable mapping; review required"
            final = deepcopy(predictions[winner][index])
            final.status = "mapped"
        else:
            review = True
            final.status = "Not selected: " + reason
        record.update(
            ValidMappers="; ".join(valid),
            SupportingMappers="; ".join(supporters),
            CarbonMappingGroups=len(groups),
            HeavyAtomAgreement=(
                "identity"
                if base.is_transport
                else "agreement"
                if len(set(heavy_keys.values())) == 1 and len(valid) > 1
                else "disagreement"
                if len(set(heavy_keys.values())) > 1
                else "not comparable"
            ),
            SelectedMapper=winner,
            ReviewRequired=review,
            AIGenerated=False,
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
        "AIGenerated",
        "SelectionReason",
        "Equation",
        "MappingDirection",
        "ValidMappers",
        "SupportingMappers",
        "CarbonMappingGroups",
        "HeavyAtomAgreement",
    ]
    return selected, table[first + [column for column in table if column not in first]]
