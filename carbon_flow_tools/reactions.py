import math
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

import pandas as pd
from rdkit import Chem

from . import config


@dataclass
class Component:
    metabolite: str
    component_index: int
    smiles: str


@dataclass
class PreparedReaction:
    reaction_id: str
    reaction_name: str
    flux: float
    equation: str
    reaction_smiles: str
    left: List[Component]
    right: List[Component]
    source_text: str = ""
    status: str = "ready"
    mapped_smiles: str = ""
    template: str = ""
    confident: Any = ""
    is_transport: bool = False


def clean(value) -> str:
    return "" if pd.isna(value) else str(value).strip()


def get_smiles(row):
    """Manual/input structures precede automatically retrieved candidates."""
    return next(
        (
            clean(getattr(row, field, ""))
            for field in (
                "SMILES_Manual",
                "SMILES_Model",
                "SMILES_Yeast9",
                "SMILES_Retrieved",
            )
            if clean(getattr(row, field, ""))
        ),
        "",
    )


def split_equation(equation: str) -> Tuple[str, str]:
    for arrow in ("<=>", "=>", "->"):
        if arrow in equation:
            return tuple(part.strip() for part in equation.split(arrow, 1))
    raise ValueError("Unrecognized reaction arrow")


def parse_side(text: str) -> List[Tuple[float, str]]:
    if not text:
        return []
    terms = []
    for raw in text.split(" + "):
        match = re.match(r"^\s*(?:(\d+(?:\.\d+)?)\s+)?(.+?)\s*$", raw)
        if not match:
            raise ValueError(f"Cannot parse: {raw}")
        terms.append((float(match.group(1) or 1), match.group(2)))
    return terms


def build_side(
    terms: List[Tuple[float, str]], smiles_by_name: Dict[str, str]
) -> Tuple[List[Component], str]:
    components = []  # type: List[Component]
    for coefficient, name in terms:
        if not math.isclose(coefficient, round(coefficient), abs_tol=1e-9):
            return [], f"Noninteger stoichiometric coefficient: {coefficient:g} {name}"
        smiles = smiles_by_name.get(name, "")
        if not smiles:
            return [], f"Missing SMILES: {name}"
        fragments = smiles.split(".")
        if any(Chem.MolFromSmiles(fragment) is None for fragment in fragments):
            return [], f"Invalid SMILES: {name}"
        for _ in range(int(round(coefficient))):
            for component_index, fragment in enumerate(fragments):
                components.append(Component(name, component_index, fragment))
    return components, ""


def transport_mapping(reaction):
    """Preserve atoms for identity- and stoichiometry-preserving compartment transport."""
    if not reaction.left or len(reaction.left) != len(reaction.right):
        return ""
    remaining = list(enumerate(reaction.right))
    left_smiles, right_smiles = [], [None] * len(reaction.right)
    number, changed_compartment = 1, False
    for source in reaction.left:
        source_name = re.fullmatch(r"(.+)\[([^\[\]]+)\]", source.metabolite)
        if source_name is None:
            return ""
        molecule = Chem.MolFromSmiles(source.smiles)
        structure = Chem.MolToSmiles(molecule)
        for index, (position, target) in enumerate(remaining):
            target_name = re.fullmatch(r"(.+)\[([^\[\]]+)\]", target.metabolite)
            if target_name is None or target_name.group(1) != source_name.group(1):
                continue
            product = Chem.MolFromSmiles(target.smiles)
            if Chem.MolToSmiles(product) != structure:
                continue
            match = product.GetSubstructMatch(molecule, useChirality=True)
            if len(match) != molecule.GetNumAtoms():
                continue
            for atom, target_index in zip(molecule.GetAtoms(), match):
                atom.SetAtomMapNum(number)
                product.GetAtomWithIdx(target_index).SetAtomMapNum(number)
                number += 1
            left_smiles.append(Chem.MolToSmiles(molecule, canonical=False))
            right_smiles[position] = Chem.MolToSmiles(product, canonical=False)
            changed_compartment |= source_name.group(2) != target_name.group(2)
            remaining.pop(index)
            break
        else:
            return ""
    return (
        ".".join(left_smiles) + ">>" + ".".join(right_smiles)
        if changed_compartment
        else ""
    )


def prepare_reactions(
    reactions: pd.DataFrame, metabolites: pd.DataFrame
) -> List[PreparedReaction]:
    smiles_by_name = {}
    for row in metabolites.itertuples(index=False):
        smiles_by_name[clean(row.MetaboliteName)] = get_smiles(row)

    prepared = []
    for row in reactions.itertuples(index=False):
        reaction = PreparedReaction(
            clean(row.ReactionID),
            clean(row.ReactionName),
            float(row.Flux),
            clean(row.Equation),
            "",
            [],
            [],
        )
        try:
            left_text, right_text = split_equation(reaction.equation)
            if reaction.flux < 0:
                left_text, right_text = right_text, left_text
            reaction.source_text = left_text
            left_terms, right_terms = parse_side(left_text), parse_side(right_text)
            if not left_terms or not right_terms:
                reaction.status = "Boundary reaction; no atom mapping required"
            else:
                reaction.left, error = build_side(left_terms, smiles_by_name)
                if not error:
                    reaction.right, error = build_side(right_terms, smiles_by_name)
                if error:
                    reaction.status = error
                else:
                    reaction.reaction_smiles = (
                        ".".join(x.smiles for x in reaction.left)
                        + ">>"
                        + ".".join(x.smiles for x in reaction.right)
                    )
                    reaction.mapped_smiles = transport_mapping(reaction)
                    if reaction.mapped_smiles:
                        reaction.is_transport = True
                        reaction.status = "mapped"
                        reaction.template = "Identity transport mapping"
                    else:
                        atom_count = sum(
                            Chem.MolFromSmiles(component.smiles).GetNumAtoms()
                            for component in reaction.left + reaction.right
                        )
                        if atom_count > config.MAX_REACTION_ATOMS:
                            reaction.status = f"Reaction exceeds the configured size limit: {atom_count} atoms"
        except Exception as exc:
            reaction.status = f"Equation parsing failed: {exc}"
        prepared.append(reaction)
    return prepared
