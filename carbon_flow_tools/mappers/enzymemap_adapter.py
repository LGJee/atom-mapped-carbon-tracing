"""Match EnzymeMap records, then use official single-step rules."""

import csv
from functools import lru_cache
from pathlib import Path

from rdkit import Chem, RDLogger


@lru_cache(maxsize=50000)
def structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES: " + smiles)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)

    return Chem.MolToSmiles(mol, isomericSmiles=False)


def key(reaction):
    return tuple(
        tuple(sorted(structure(s) for s in side.split(".")))
        for side in reaction.split(">>")
    )


def transfer(input_reaction, mapped_reaction, normalizer=None):
    """Transfer database map IDs while retaining input charge and stereochemistry."""
    output = []
    for given, stored in zip(input_reaction.split(">>"), mapped_reaction.split(">>")):
        available = stored.split(".")
        parts = []
        for smiles in given.split("."):
            base = Chem.MolFromSmiles(smiles)
            query = normalizer(Chem.Mol(base)) if normalizer else base
            if query.GetNumAtoms() != base.GetNumAtoms():
                raise ValueError(
                    "Normalization changed the atom count; cannot restore the input."
                )
            for index, candidate in enumerate(available):
                if structure(Chem.MolToSmiles(query)) != structure(candidate):
                    continue
                mapped = Chem.MolFromSmiles(candidate)
                match = mapped.GetSubstructMatch(query, useChirality=True)
                if len(match) != base.GetNumAtoms():
                    continue
                for atom in base.GetAtoms():
                    atom.SetAtomMapNum(
                        mapped.GetAtomWithIdx(match[atom.GetIdx()]).GetAtomMapNum()
                    )
                parts.append(Chem.MolToSmiles(base))
                available.pop(index)
                break
            else:
                if (
                    normalizer
                    and base.GetNumAtoms() == 1
                    and base.GetAtomWithIdx(0).GetAtomicNum() == 1
                ):
                    parts.append(Chem.MolToSmiles(base))
                    continue
                raise ValueError(
                    "Database molecules do not match input structures or specified stereochemistry."
                )
        if available and not (
            normalizer
            and all(
                Chem.MolFromSmiles(s).GetNumAtoms() == 1
                and Chem.MolFromSmiles(s).GetAtomWithIdx(0).GetAtomicNum() == 1
                for s in available
            )
        ):
            raise ValueError("Database stoichiometry does not match input.")
        output.append(".".join(parts))
    if len(output) != 2:
        raise ValueError("Reaction must contain >>.")
    return ">>".join(output)


def identity(reaction):
    """Map identical structures without claiming a database prediction."""
    left, right = key(reaction)
    if left != right:
        return None
    mapped, number = [], 1
    for smiles in reaction.split(">>")[0].split("."):
        mol = Chem.MolFromSmiles(smiles)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(number)
            number += 1
        mapped.append(Chem.MolToSmiles(mol))
    side = ".".join(mapped)
    try:
        return transfer(reaction, side + ">>" + side)
    except ValueError:
        return None


def load_model(options):
    path = Path(options["dataset"])
    if not path.is_file():
        raise FileNotFoundError("EnzymeMap dataset not found: " + str(path))
    requested = {key(r) for r in options["reactions"]}
    matches = {}
    RDLogger.DisableLog("rdApp.warning")
    with path.open(encoding="utf-8-sig", newline="") as source:
        rows = csv.DictReader(source)
        if not {"mapped", "unmapped"} <= set(rows.fieldnames or []):
            raise ValueError("EnzymeMap CSV requires mapped and unmapped columns.")
        for number, row in enumerate(rows, 2):
            if number % 10000 == 0:
                print(
                    "[EnzymeMap lookup] scanned {} records; {} reactions matched.".format(
                        number - 1, len(matches)
                    ),
                    flush=True,
                )
            try:
                reaction_key = key(row["unmapped"])
                candidates = [
                    (reaction_key, row["mapped"], "forward"),
                    (
                        reaction_key[::-1],
                        ">>".join(row["mapped"].split(">>")[::-1]),
                        "reverse",
                    ),
                ]
                for candidate, mapped, direction in candidates:
                    if candidate not in requested:
                        continue
                    quality = float(row.get("quality") or 0)
                    old = matches.get(candidate)
                    if old is None or quality > old["quality"]:
                        matches[candidate] = {
                            "mapped": mapped,
                            "quality": quality,
                            "template": "EnzymeMap CSV line {}; {}; rxn_idx={}; rule_id={}; {}".format(
                                number,
                                direction,
                                row.get("rxn_idx", ""),
                                row.get("rule_id", ""),
                                row.get("source", ""),
                            ),
                        }
            except (ValueError, TypeError):
                continue
    print(
        "[EnzymeMap lookup] {} requested reactions; {} database matches.".format(
            len(requested), len(matches)
        ),
        flush=True,
    )
    predictor = None
    if options.get("rules"):
        from .enzymemap_rules import RuleMapper

        predictor = RuleMapper(options)
    return {"matches": matches, "predictor": predictor}


def map_batch(model, reaction_smiles):
    results = []
    for reaction in reaction_smiles:
        same = identity(reaction)
        if same:
            results.append(
                {
                    "mapped_rxn": same,
                    "template": "Identity transport; not an EnzymeMap prediction",
                }
            )
            continue
        match = model["matches"].get(key(reaction))
        error = (
            "No EnzymeMap record with matching structures, charges and stoichiometry."
        )
        if match is not None:
            try:
                results.append(
                    {
                        "mapped_rxn": transfer(reaction, match["mapped"]),
                        "confidence": match["quality"],
                        "template": match["template"],
                    }
                )
                continue
            except ValueError as exc:
                error = str(exc)
        if model["predictor"] is not None:
            results.append(model["predictor"].map(reaction))
        else:
            results.append({"error": error})
    return results
