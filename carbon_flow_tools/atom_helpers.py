from typing import Dict, Iterable, Tuple

from rdkit import Chem

from .reactions import PreparedReaction


def labeled_smiles(smiles: str, assignment: Iterable[Tuple[int, int, str]]) -> str:
    """Write source-origin carbon labels into the product SMILES atom-map fields."""
    labels = {
        (int(part), int(atom)): int(str(label)[1:]) for part, atom, label in assignment
    }
    output = []
    for part, fragment in enumerate(smiles.split(".")):
        molecule = Chem.MolFromSmiles(fragment)
        if molecule is None:
            raise ValueError(f"Invalid SMILES: {fragment}")
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(labels.get((part, atom.GetIdx()), 0))
        output.append(Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True))
    return ".".join(output)


def assignment_from_labeled_smiles(
    base_smiles: str, marked_smiles: str
) -> Tuple[Tuple[int, int, str], ...]:
    """Recover origin labels at exact positions in the input structure."""
    base_parts = base_smiles.split(".")
    marked_parts = marked_smiles.split(".")
    if len(base_parts) != len(marked_parts):
        return tuple()
    assignment = []
    for part, (base, marked) in enumerate(zip(base_parts, marked_parts)):
        for atom, number in component_atom_maps(base, marked).items():
            if number > 0:
                assignment.append((part, atom, f"C{number}"))
    return tuple(sorted(assignment))


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


def glucose_labels(smiles: str) -> Dict[int, str]:
    """Identify biochemical C1-C6 in cyclic glucose."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid D-glucose SMILES")
    six_rings = [ring for ring in mol.GetRingInfo().AtomRings() if len(ring) == 6]
    for ring in six_rings:
        ring_set = set(ring)
        oxygens = [i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "O"]
        if len(oxygens) != 1:
            continue
        ring_o = oxygens[0]
        ring_c = [i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "C"]
        c5 = None
        c6 = None
        for carbon in ring_c:
            outside_c = [
                n.GetIdx()
                for n in mol.GetAtomWithIdx(carbon).GetNeighbors()
                if n.GetIdx() not in ring_set and n.GetSymbol() == "C"
            ]
            if outside_c:
                c5, c6 = carbon, outside_c[0]
                break
        if c5 is None:
            continue
        o_neighbors = [
            n.GetIdx()
            for n in mol.GetAtomWithIdx(ring_o).GetNeighbors()
            if n.GetIdx() in ring_set
        ]
        c1_candidates = [i for i in o_neighbors if i != c5]
        if len(c1_candidates) != 1:
            continue
        chain = [c1_candidates[0]]
        previous = ring_o
        while chain[-1] != c5:
            current = chain[-1]
            next_c = [
                n.GetIdx()
                for n in mol.GetAtomWithIdx(current).GetNeighbors()
                if n.GetIdx() in ring_set
                and n.GetSymbol() == "C"
                and n.GetIdx() != previous
            ]
            if len(next_c) != 1:
                break
            previous, current = current, next_c[0]
            chain.append(current)
        if len(chain) == 5 and chain[-1] == c5:
            labels = {atom: f"C{i}" for i, atom in enumerate(chain, 1)}
            labels[c6] = "C6"
            return labels
    raise ValueError("Cannot identify glucose C1-C6 from the structure")


def isoprenol_labels(smiles: str) -> Dict[int, str]:
    """Identify isoprenol from the hydroxyl end (legacy local labels C7-C11)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid isoprenol SMILES")
    alcohol_oxygen = next(
        (atom for atom in mol.GetAtoms() if atom.GetSymbol() == "O"), None
    )
    if alcohol_oxygen is None:
        raise ValueError("No oxygen atom found in isoprenol")
    c1 = next(
        (n.GetIdx() for n in alcohol_oxygen.GetNeighbors() if n.GetSymbol() == "C"),
        None,
    )
    if c1 is None:
        raise ValueError("Cannot identify isoprenol C1")
    c2 = next(
        n.GetIdx()
        for n in mol.GetAtomWithIdx(c1).GetNeighbors()
        if n.GetSymbol() == "C"
    )
    c3 = next(
        n.GetIdx()
        for n in mol.GetAtomWithIdx(c2).GetNeighbors()
        if n.GetSymbol() == "C" and n.GetIdx() != c1
    )
    branches = [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(c3).GetNeighbors()
        if n.GetSymbol() == "C" and n.GetIdx() != c2
    ]
    c4 = next(
        i
        for i in branches
        if mol.GetBondBetweenAtoms(c3, i).GetBondType() == Chem.BondType.DOUBLE
    )
    c5 = next(i for i in branches if i != c4)
    return {c1: "C7", c2: "C8", c3: "C9", c4: "C10", c5: "C11"}
