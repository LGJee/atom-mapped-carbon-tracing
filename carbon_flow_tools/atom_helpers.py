from typing import Dict

from rdkit import Chem

from .reactions import PreparedReaction










def component_atom_maps(base_smiles: str, mapped_smiles: str) -> Dict[int, int]:
    """Align input atom indices with reaction-local map identifiers."""
    base = Chem.MolFromSmiles(base_smiles)
    mapped = Chem.MolFromSmiles(mapped_smiles)
    if base is None or mapped is None or base.GetNumAtoms() != mapped.GetNumAtoms():
        return {}
    match = mapped.GetSubstructMatch(base)
    if not match:
        match = base.GetSubstructMatch(mapped)
        if match:
            return {
                base_index: mapped.GetAtomWithIdx(mapped_index).GetAtomMapNum()
                for mapped_index, base_index in enumerate(match)
            }
        return {}
    return {
        base_index: mapped.GetAtomWithIdx(mapped_index).GetAtomMapNum()
        for base_index, mapped_index in enumerate(match)
    }


def mapped_components(reaction: PreparedReaction):
    # Cache the fixed correspondence, not the growing substrate combinations.
    signature = (
        reaction.mapped_smiles,
        tuple(
            (c.metabolite, c.component_index, c.smiles)
            for c in reaction.left + reaction.right
        ),
    )
    if getattr(reaction, "_mapping_signature", None) == signature:
        return reaction._mapped_components
    try:
        left_text, right_text = reaction.mapped_smiles.split(">>", 1)
        left_smiles = left_text.split(".") if left_text else []
        right_smiles = right_text.split(".") if right_text else []
        if len(left_smiles) != len(reaction.left) or len(right_smiles) != len(
            reaction.right
        ):
            return None

        # Mappers may reorder molecules; preserve original stoichiometric order.
        def align(components, mapped_parts):
            remaining = list(mapped_parts)
            result = []
            for component in components:
                base = Chem.MolToSmiles(Chem.MolFromSmiles(component.smiles))
                for index, mapped in enumerate(remaining):
                    molecule = Chem.MolFromSmiles(mapped)
                    if molecule is None:
                        continue
                    for atom in molecule.GetAtoms():
                        atom.SetAtomMapNum(0)
                    if Chem.MolToSmiles(molecule) == base:
                        result.append(
                            (component, component_atom_maps(component.smiles, mapped))
                        )
                        remaining.pop(index)
                        break
                else:
                    return None
            return result

        left = align(reaction.left, left_smiles)
        right = align(reaction.right, right_smiles)
        if left is None or right is None:
            return None
        reaction._mapping_signature = signature
        reaction._mapped_components = (left, right)
        return reaction._mapped_components
    except Exception:
        return None
