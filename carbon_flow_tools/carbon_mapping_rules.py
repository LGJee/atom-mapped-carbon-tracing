"""User-maintained carbon correspondences and automatic mapper fallbacks.

Rules, in order:
* Classify carbon from explicit structures/formulas; missing evidence is unknown.
* Open-boundary reactions need no internal carbon mapping.
* Preserve mapper predictions, including failed candidates retained for review.
* Recover only structurally recognized, uniquely paired carbon skeletons.
* Save a reduced unit subreaction plus exact decimal amounts in hidden Mapping_Info.
* Unknown pools and transfers remain unresolved; never average or invent atoms.

Cofactor skeleton retention does not establish attached-group transfer. Generic
biomass/protein pools have no unique positional structure. Graph similarity is
not biochemical proof; every rule prediction requires review. See
``docs/carbon-mapping-rules.md`` for class-specific assumptions and references.
"""

from collections import Counter
from functools import lru_cache
from pathlib import Path
import csv
import json
import math
import re

import pandas as pd
from rdkit import Chem

from .reactions import Component, PreparedReaction, clean, get_smiles, parse_side, split_equation

RULE_VERSION = "1"
COLUMNS = ["ReactionID", "MappingDirection", "SourceMetabolite", "ProductMetabolite",
           "TransferCoefficient", "SourceSMILES", "ProductSMILES", "CarbonAtomMap",
           "Rule", "RuleVersion"]
# Full connectivity references distinguish ADP from positional isomers such as PAP.
_REFERENCES = {
    "ATP": "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
    "ADP": "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
    "AMP": "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)O)O)O)N",
    "NAD": "C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N",
    "NADP": "C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)([O-])OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)([O-])[O-])O)O)O)C(=O)N",
}
_REFERENCES.update({
    "GTP": "C1=NC2=C(N1C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N=C(NC2=O)N",
    "CTP": "C1=CN(C(=O)N=C1N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O",
    "UTP": "C1=CN(C(=O)NC1=O)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O",
    "CoA": "CC(C)(COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)C(C(=O)NCCC(=O)NCCS)O",
})
NUCLEOTIDE_FAMILIES = {"adenylate", "guanylate", "cytidylate", "uridylate"}
_FAMILY_BY_NAME = {}
for _prefix, _family_name in (("A", "adenylate"), ("G", "guanylate"),
                             ("C", "cytidylate"), ("U", "uridylate")):
    for _count, _suffix in ((1, "MP"), (2, "DP"), (3, "TP")):
        _name = _prefix + _suffix
        _REFERENCES[_name] = _REFERENCES[_prefix + "TP"].replace(
            "OP(=O)(O)OP(=O)(O)OP(=O)(O)O", "OP(=O)(O)" * _count + "O")
        _FAMILY_BY_NAME[_name] = _family_name


# A formula is metadata, never a substitute for a positional molecular graph.
_FORMULA = r"(?:[A-Z][a-z]?\d*)+"
_REPEAT_FORMULA = r"(?:[A-Z][a-z]?\d*|\((?:[A-Z][a-z]?\d*)+\)[nm]\d*|[nm])+"
STRUCTURE_LIMITS = {
    "partial": "Partial structure: explicit carbon positions only; wildcard groups remain unresolved",
    "polymer": "Unspecified polymer repeat count; no unique positional structure (n is not set to 1)",
    "invalid": "Invalid SMILES syntax or unsupported valence",
    "missing": "Missing SMILES",
}


def repeat_formula(formula):
    """Recognize variable repeats only in formula syntax, not aromatic SMILES n."""
    text = re.sub(r"[+-]\d*$", "", clean(formula))
    return bool(re.fullmatch(_REPEAT_FORMULA, text) and re.search(
        r"[nm]", re.sub(r"[A-Z][a-z]?\d*", "", text)))


@lru_cache(maxsize=4096)
def structure_kind(smiles, formula=""):
    """Classify without editing the original SMILES, charges or stereochemistry.

    RDKit dummy atoms (atomic number zero) are legal attachment placeholders,
    not carbon atoms. Matching two dummies does not establish carrier identity.
    A variable-repeat formula cannot define the omitted carbon positions.
    """
    text = clean(smiles)
    if repeat_formula(formula):
        return "polymer"
    if not text:
        return "missing"
    molecule = Chem.MolFromSmiles(text)
    if molecule is None or not molecule.GetNumAtoms():
        return "polymer" if repeat_formula(text) else "invalid"
    return "partial" if any(a.GetAtomicNum() == 0 for a in molecule.GetAtoms()) else "complete"


def formula_counts(formula):
    """Count explicit elements; variable groups and repeat units stay unresolved."""
    text = clean(formula)
    if not re.fullmatch(_FORMULA, text):
        return None
    totals = Counter()
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", text):
        if element in ("R", "X"):
            return None
        totals[element] += int(count or 1)
    return {element: count for element, count in totals.items() if count}


def formulas_match(first, second):
    """Compare elemental composition, accepting explicit 1 and reordered tokens.

    Charge is checked separately by callers; unsupported formulas require exact
    text equality, never equality between two failed parses.
    """
    left, right = formula_counts(first), formula_counts(second)
    return left == right if left is not None and right is not None else clean(first) == clean(second)


def protonation_compatible(first_formula, first_charge, second_formula, second_charge):
    """Allow a reviewed H+/charge mismatch, never a generic charge bypass.

    Non-hydrogen element counts must match and delta H must equal delta charge.
    This is compatibility evidence, not proof of compound identity or pH state;
    missing/variable formulas and redox-like differences remain unresolved.
    """
    first, second = formula_counts(first_formula), formula_counts(second_formula)
    if first is None or second is None or not clean(first_charge) or not clean(second_charge):
        return False
    if {e: n for e, n in first.items() if e != "H"} != {e: n for e, n in second.items() if e != "H"}:
        return False
    try:
        charges = [float(first_charge), float(second_charge)]
    except (TypeError, ValueError):
        return False
    if not all(math.isfinite(q) and q.is_integer() for q in charges):
        return False
    hydrogen_difference = first.get("H", 0) - second.get("H", 0)
    return hydrogen_difference != 0 and hydrogen_difference == charges[0] - charges[1]


def fixed_formula_carbons(formula):
    """Return a fixed total only for an explicit formula without variable groups."""
    text = re.sub(r"[+-]\d*$", "", clean(formula))
    counts = formula_counts(text)
    return None if counts is None else counts.get("C", 0)


def partial_transport_pairs(left, right, smiles):
    """Use explicit identity-preserving transport, never wildcard identity alone.

    Keep only carbon pairs invariant over all exact graph alignments. Unknown
    groups and symmetry-ambiguous positions are not assigned. Chemical changes
    of ACP/lipid variable groups require independent reaction-specific evidence.
    """
    if len(left) != 1 or len(right) != 1:
        return None
    source, target = left[0][1], right[0][1]
    a = re.fullmatch(r"(.+)\[([^\[\]]+)\]", source)
    b = re.fullmatch(r"(.+)\[([^\[\]]+)\]", target)
    if not a or not b or a.group(1) != b.group(1) or a.group(2) == b.group(2):
        return None
    if structure_kind(smiles.get(source, "")) != "partial" or structure_kind(smiles.get(target, "")) != "partial":
        return None
    first, second = _structure(smiles[source]), _structure(smiles[target])
    if first is None or second is None or Chem.MolToSmiles(first) != Chem.MolToSmiles(second):
        return None
    matches = second.GetSubstructMatches(first, uniquify=False, useChirality=True, maxMatches=257)
    if not matches or len(matches) >= 257:
        return None
    assignments = [{(atom.GetIdx(), match[atom.GetIdx()]) for atom in first.GetAtoms()
                    if atom.GetAtomicNum() == 6 and second.GetAtomWithIdx(match[atom.GetIdx()]).GetAtomicNum() == 6}
                   for match in matches]
    pairs = tuple(sorted(set.intersection(*assignments)))
    return ("partial_identity_transport", pairs) if pairs else None


def _graph(molecule, adenosine=False):
    """Normalize only an internal matching view; saved structures stay unchanged."""
    graph = Chem.RWMol(molecule)
    for atom in graph.GetAtoms():
        atom.SetIntProp("original_index", atom.GetIdx())
        atom.SetAtomMapNum(0)
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
        atom.SetIsAromatic(False)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for bond in graph.GetBonds():
        bond.SetIsAromatic(False)
        bond.SetBondType(Chem.BondType.SINGLE)
    if adenosine:
        removed = [a.GetIdx() for a in graph.GetAtoms()
                   if a.GetAtomicNum() == 15 or (a.GetAtomicNum() == 8 and
                       not any(n.GetAtomicNum() == 6 for n in a.GetNeighbors()))]
        for index in reversed(removed):
            graph.RemoveAtom(index)
    result = graph.GetMol()
    result.UpdatePropertyCache(strict=False)
    return result


@lru_cache(maxsize=1024)
def _structure(smiles):
    molecule = Chem.MolFromSmiles(smiles) if smiles and "." not in smiles else None
    return molecule


@lru_cache(maxsize=1024)
def _family(smiles):
    molecule = _structure(smiles)
    if molecule is None:
        return ""
    key = Chem.MolToSmiles(_graph(molecule), isomericSmiles=True)
    for name, reference in _REFERENCES.items():
        if key == Chem.MolToSmiles(_graph(_structure(reference)), isomericSmiles=True):
            return _FAMILY_BY_NAME.get(name, name)
    return ""


@lru_cache(maxsize=1024)
def _coa_core(smiles):
    molecule = _structure(smiles)
    if molecule is None:
        return None
    graph = _graph(molecule)
    reference = Chem.MolToSmiles(_graph(_structure(_REFERENCES["CoA"])))
    if Chem.MolToSmiles(graph) == reference:
        return graph
    candidates = []
    for bond in molecule.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetAtomicNum() == 6 and b.GetAtomicNum() == 16:
            a, b = b, a
        if a.GetAtomicNum() != 16 or b.GetAtomicNum() != 6:
            continue
        if not any(n.GetAtomicNum() == 8 and molecule.GetBondBetweenAtoms(
                b.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.DOUBLE
                   for n in b.GetNeighbors()):
            continue
        cut = Chem.RWMol(graph)
        cut.RemoveBond(a.GetIdx(), b.GetIdx())
        for fragment in Chem.GetMolFrags(cut.GetMol(), asMols=True, sanitizeFrags=False):
            if Chem.MolToSmiles(fragment) == reference:
                candidates.append(fragment)
    return candidates[0] if len(candidates) == 1 else None


@lru_cache(maxsize=2048)
def carbon_pairs(source, product):
    """Return only uniquely aligned, structurally recognized carbon skeletons."""
    left, right = _structure(source), _structure(product)
    if (left is None or right is None or structure_kind(source) != "complete"
            or structure_kind(product) != "complete"):
        return None
    family = _family(source)
    target_family = _family(product)
    identical = Chem.MolToSmiles(left) == Chem.MolToSmiles(right)
    a, b = _coa_core(source), _coa_core(product)
    coa = a is not None and b is not None
    if coa:
        rule = "CoA_core"
    else:
        if not identical and (not family or family != target_family):
            return None
        rule = "identity_skeleton" if identical else family + "_skeleton"
        a = _graph(left, family in NUCLEOTIDE_FAMILIES)
        b = _graph(right, family in NUCLEOTIDE_FAMILIES)
    if a.GetNumAtoms() != b.GetNumAtoms():
        return None
    matches = b.GetSubstructMatches(a, uniquify=False, maxMatches=257)
    if not matches or len(matches) >= 257:
        return None
    assignments = set()
    for match in matches:
        pairs = tuple(sorted((atom.GetIntProp("original_index"),
                              b.GetAtomWithIdx(match[atom.GetIdx()]).GetIntProp("original_index"))
                             for atom in a.GetAtoms() if atom.GetAtomicNum() == 6))
        if any(left.GetAtomWithIdx(i).HasProp("_CIPCode") and
               right.GetAtomWithIdx(j).HasProp("_CIPCode") and
               left.GetAtomWithIdx(i).GetProp("_CIPCode") !=
               right.GetAtomWithIdx(j).GetProp("_CIPCode") for i, j in pairs):
            continue
        assignments.add(pairs)
    if not assignments:
        return None
    if coa:
        # Equivalent pantothenate methyl positions stay unresolved, not arbitrarily paired.
        pairs = tuple(sorted(set.intersection(*(set(p) for p in assignments))))
    elif len(assignments) == 1:
        pairs = next(iter(assignments))
    else:
        return None
    if not pairs:
        return None
    return rule, pairs


def carbon_presence(row):
    """Return True/False/None; formula-only carbon evidence does not define positions."""
    formula = clean(getattr(row, "MolecularFormula", ""))
    count = fixed_formula_carbons(formula)
    text = get_smiles(row)
    kind = structure_kind(text, formula)
    molecule = Chem.MolFromSmiles(text) if kind in ("complete", "partial") else None
    if molecule is not None:
        actual = sum(a.GetAtomicNum() == 6 for a in molecule.GetAtoms())
        if count is not None and (actual > count or (kind == "complete" and bool(actual) != bool(count))):
            return None
        return (True if actual or count else None) if kind == "partial" else bool(actual)
    if count is not None:
        return bool(count)
    # Carbon in a repeat formula proves presence, not a count or atom ordering.
    return True if re.search(r"C(?=\d|[A-Z(]|$)", formula) else None


def partial_mapped_smiles(records):
    """Serialize reliable unit subreactions; actual amounts remain in the records."""
    sides = [[], []]
    number = 1
    for row in records:
        molecules = [_structure(row[k]) for k in ("SourceSMILES", "ProductSMILES")]
        molecules = [Chem.Mol(m) for m in molecules]
        for molecule in molecules:
            for atom in molecule.GetAtoms():
                atom.SetAtomMapNum(0)
        for part_a, i, part_b, j in json.loads(row["CarbonAtomMap"]):
            if part_a or part_b:
                raise ValueError("Rule unit structures must be single fragments")
            molecules[0].GetAtomWithIdx(i).SetAtomMapNum(number)
            molecules[1].GetAtomWithIdx(j).SetAtomMapNum(number)
            number += 1
        for side, molecule in zip(sides, molecules):
            side.append(Chem.MolToSmiles(molecule))
    return ">>".join(".".join(side) for side in sides) if records else ""


def carbon_structure_issues(metabolites):
    """Detect a known carbon-count mismatch without guessing a corrected structure."""
    issues = {}
    for row in metabolites.itertuples(index=False):
        formula = clean(getattr(row, "MolecularFormula", ""))
        text = get_smiles(row)
        kind = structure_kind(text, formula)
        if kind == "polymer":
            continue
        expected = fixed_formula_carbons(formula)
        molecule = Chem.MolFromSmiles(text) if kind in ("complete", "partial") else None
        if expected is None or molecule is None:
            continue
        actual = sum(a.GetAtomicNum() == 6 for a in molecule.GetAtoms())
        # A partial graph can contain fewer carbons than the whole compound.
        if actual > expected or (kind == "complete" and expected != actual):
            issues[clean(row.MetaboliteName)] = "formula C{} versus SMILES C{}".format(expected, actual)
    return issues


def rule_smiles(metabolites):
    """Exclude unavailable positional graphs without blocking other participants."""
    issues = carbon_structure_issues(metabolites)
    return {clean(row.MetaboliteName): get_smiles(row)
            for row in metabolites.itertuples(index=False)
            if clean(row.MetaboliteName) not in issues and structure_kind(
                get_smiles(row), clean(getattr(row, "MolecularFormula", ""))) in ("complete", "partial")}


def generic_pool(name):
    """A generic pool name is not the identity of a discrete small molecule."""
    name = re.sub(r"\[[^\[\]]+\]$", "", name).strip().lower()
    return (name in {"biomass", "protein", "rna", "dna", "lipid", "lipids",
                     "carbohydrate", "carbohydrates", "cofactor", "cofactors", "ion", "ions"}
            or "trna" in name)


def rule_transfers(reaction, smiles):
    """Recover only unique, structurally recognized carbon-skeleton transfers.

    Work in the flux direction and require matching original stoichiometric
    amounts. Generic biomass/polymer pools do not define positional structures.
    Unpaired or ambiguous carbon stays unresolved; this is not a full-reaction map.
    """
    left, right = split_equation(reaction.equation)
    if reaction.flux < 0:
        left, right = right, left
    left, right = parse_side(left), parse_side(right)
    candidates = []
    for amount, source in left:
        if generic_pool(source):
            continue
        if not math.isfinite(amount) or amount <= 0:
            continue
        for target_amount, product in right:
            if generic_pool(product):
                continue
            if (not math.isfinite(target_amount) or target_amount <= 0 or
                    not math.isclose(amount, target_amount, rel_tol=1e-12, abs_tol=1e-12)):
                continue
            mapping = partial_transport_pairs(left, right, smiles)
            if mapping is None:
                mapping = carbon_pairs(smiles.get(source, ""), smiles.get(product, ""))
            if mapping:
                candidates.append((source, product, amount, mapping))
    # Reject competing partners instead of choosing a plausible but unsupported pair.
    donors = Counter(item[0] for item in candidates)
    products = Counter(item[1] for item in candidates)
    return [dict(zip(COLUMNS, [reaction.reaction_id,
                "forward" if reaction.flux >= 0 else "reverse", source, product,
                amount, smiles[source], smiles[product],
                json.dumps([[0, i, 0, j] for i, j in pairs]), rule, RULE_VERSION]))
            for source, product, amount, (rule, pairs) in candidates
            if donors[source] == 1 and products[product] == 1]


def validate_manual_records(reaction, records, smiles):
    """Bind human-provided carbon positions to exact saved structures and amounts."""
    sides = list(split_equation(reaction.equation))
    if reaction.flux < 0:
        sides.reverse()
    amounts = [Counter() for _ in sides]
    for totals, side in zip(amounts, sides):
        for amount, name in parse_side(side):
            totals[name] += amount
    used, supplied = Counter(), Counter()
    expected_direction = "forward" if reaction.flux >= 0 else "reverse"
    for row in records:
        if (row["MappingDirection"] != expected_direction or row["Rule"] != "ManualVerified"
                or str(row["RuleVersion"]) != RULE_VERSION or not clean(row.get("Evidence", ""))):
            raise ValueError("Manual mapping requires matching direction, version and evidence: " + reaction.reaction_id)
        source, target = row["SourceMetabolite"], row["ProductMetabolite"]
        molecules = []
        for name, field in ((source, "SourceSMILES"), (target, "ProductSMILES")):
            text = row[field]
            if not text or text != smiles.get(name) or "." in text:
                raise ValueError("Manual atom indices require the exact single-component input SMILES: " + name)
            molecule = _structure(text)
            if molecule is None:
                raise ValueError("Invalid manual structure: " + name)
            molecules.append(molecule)
        amount = float(row["TransferCoefficient"])
        if (not math.isfinite(amount) or amount <= 0 or amount > amounts[0][source] + 1e-9
                or amount > amounts[1][target] + 1e-9):
            raise ValueError("Manual transfer exceeds original directional stoichiometry: " + reaction.reaction_id)
        pairs = json.loads(row["CarbonAtomMap"])
        if not isinstance(pairs, list) or not pairs:
            raise ValueError("Manual mapping requires carbon-position pairs")
        for pair in pairs:
            if (not isinstance(pair, list) or len(pair) != 4 or
                    any(type(number) is not int for number in pair) or pair[0] != 0 or pair[2] != 0):
                raise ValueError("Manual atom pairs require single-component integer positions")
            for molecule, position in zip(molecules, (pair[1], pair[3])):
                if position < 0 or position >= molecule.GetNumAtoms() or molecule.GetAtomWithIdx(position).GetAtomicNum() != 6:
                    raise ValueError("Manual pair does not identify carbon atoms")
            used[(source, pair[1])] += amount
            supplied[(target, pair[3])] += amount
        if any(value > amounts[0][name] + 1e-9 for (name, _), value in used.items()) or any(
                value > amounts[1][name] + 1e-9 for (name, _), value in supplied.items()):
            raise ValueError("Manual mapping assigns an original carbon more than once")
    return records


def manual_mapping_records(reaction, mapped, smiles, evidence):
    """Derive positional records and amounts from a supplied carbon subreaction."""
    from .atom_helpers import component_atom_maps
    from .mapper_runner import mapping_error
    error = mapping_error(mapped)
    if error:
        raise ValueError("Invalid manual mapping: " + error)
    texts = list(split_equation(reaction.equation))
    if reaction.flux < 0:
        texts.reverse()
    sides = []
    for text, marked_side in zip(texts, mapped.split(">>")):
        amounts = Counter()
        for amount, name in parse_side(text):
            amounts[name] += amount
        aligned = []
        for marked in marked_side.split("."):
            molecule = Chem.MolFromSmiles(marked)
            for atom in molecule.GetAtoms():
                atom.SetAtomMapNum(0)
            canonical = Chem.MolToSmiles(molecule)
            matches = [name for name in amounts if smiles.get(name) and "." not in smiles[name]
                       and _structure(smiles[name]) is not None
                       and Chem.MolToSmiles(_structure(smiles[name])) == canonical]
            if len(matches) != 1:
                raise ValueError("Manual mapping component must identify one original metabolite: " + reaction.reaction_id)
            name = matches[0]
            aligned.append((name, component_atom_maps(smiles[name], marked)))
        copies = Counter(name for name, _ in aligned)
        positions = {}
        for copy, (name, maps) in enumerate(aligned):
            base = _structure(smiles[name])
            for atom, number in maps.items():
                if number and base.GetAtomWithIdx(atom).GetAtomicNum() == 6:
                    positions[number] = (copy, name, atom, amounts[name] / copies[name])
        sides.append(positions)
    shared = sorted(sides[0].keys() & sides[1].keys())
    if not shared:
        raise ValueError("Manual mapping has no carbon-to-carbon correspondence: " + reaction.reaction_id)
    groups = {}
    for number in shared:
        source, target = sides[0][number], sides[1][number]
        if not math.isclose(source[3], target[3], rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError("Manual mapping molecular copies conflict with original stoichiometry: " + reaction.reaction_id)
        key = (source[0], target[0], source[1], target[1], source[3])
        groups.setdefault(key, []).append([0, source[2], 0, target[2]])
    records = {}
    for (_, _, source, target, amount), pairs in groups.items():
        key = (source, target, json.dumps(sorted(pairs)))
        if key in records:
            records[key]["TransferCoefficient"] += amount
        else:
            records[key] = dict(zip(COLUMNS, [reaction.reaction_id,
                "forward" if reaction.flux >= 0 else "reverse", source, target, amount,
                smiles[source], smiles[target], key[2], "ManualVerified", RULE_VERSION]),
                Evidence=evidence, Notes="Automatically derived from the user-maintained mapped reaction")
    return validate_manual_records(reaction, list(records.values()), smiles)


def apply_manual_mappings(model, reactions, selection, metabolites, filename=None, source=None):
    """Match species/reaction IDs; supplied mappings precede software and fallbacks."""
    from .smiles_completion import species_name
    from .consensus import mapping_signatures
    from .mapper_runner import mapping_error
    path = Path(filename) if filename is not None else Path(__file__).resolve().parents[1] / "data/manual_atom_mappings.csv"
    if not path.is_file():
        return
    fields = ["Species", "ReactionID", "MappedReactionSMILES"]
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if set(reader.fieldnames or []) != set(fields) or len(reader.fieldnames) != 3:
            raise ValueError("Manual atom-mapping CSV requires exactly: " + ",".join(fields))
        rows = {}
        for entry in reader:
            row = {key: clean(entry.get(key)) for key in fields}
            key = (row["Species"], row["ReactionID"])
            if not all(row.values()) or key in rows or mapping_error(row["MappedReactionSMILES"]):
                raise ValueError("Invalid, empty or duplicated manual mapping: " + str(key))
            rows[key] = row["MappedReactionSMILES"]
    species = species_name(model, source)
    smiles = rule_smiles(metabolites)
    issues = carbon_structure_issues(metabolites)
    for index, reaction in enumerate(reactions):
        mapped = rows.get((species, reaction.reaction_id))
        if mapped is None:
            continue
        if reaction.flux < 0:
            mapped = ">>".join(reversed(mapped.split(">>")))
        reaction.mapped_smiles = mapped
        reaction.partial_maps = []
        try:
            mapping_signatures(reaction)
            complete = (not reaction.omitted_components and not any(
                component.metabolite in issues for component in reaction.left + reaction.right))
        except ValueError:
            complete = False
        if not complete:
            reaction.partial_maps = manual_mapping_records(reaction, mapped, smiles,
                "User-maintained mapping: " + path.name + "; Species=" + species)
            if any(structure_kind(row[field]) == "partial" for row in reaction.partial_maps
                   for field in ("SourceSMILES", "ProductSMILES")):
                # Preserve the CSV; the usable view must not assert dummy-atom identity.
                mapped = partial_mapped_smiles(reaction.partial_maps)
                reaction.mapped_smiles = mapped
        reaction.status = "mapped" if complete else "partial"
        selection.at[index, "SelectedMapper"] = "Manual"
        selection.at[index, "ReviewRequired"] = not complete
        selection.at[index, "SelectionReason"] = (
            "User-maintained mapped reaction takes priority; " + ("complete" if complete else "partial")
            + " original-reaction coverage; source: " + path.name + "; Species=" + species)
        selection.at[index, "SelectedMappedReactionSMILES"] = mapped


def apply_rules(reactions, selection, metabolites):
    """Preserve mapper candidates; use rules only when no prediction is available.

    A rule supplies only local pairs, never invents the missing remainder,
    and always requires review. Exact amounts are persisted in Mapping_Info;
    unit mapped SMILES alone cannot represent fractional original coefficients.
    """
    smiles = {clean(r.MetaboliteName): get_smiles(r)
              for r in metabolites.itertuples(index=False)}
    issues = carbon_structure_issues(metabolites)
    usable_smiles = rule_smiles(metabolites)
    kinds = {clean(row.MetaboliteName): structure_kind(
                get_smiles(row), clean(getattr(row, "MolecularFormula", "")))
             for row in metabolites.itertuples(index=False)}
    internal_consumers = set()
    for item in reactions:
        try:
            a, b = split_equation(item.equation)
            if item.flux < 0:
                a, b = b, a
            if a and b:
                internal_consumers.update(name for _, name in parse_side(a))
        except ValueError:
            pass  # The affected reaction retains its own parse diagnostic below.
    records = []
    for index, reaction in enumerate(reactions):
        try:
            sides = split_equation(reaction.equation)
            terms = [parse_side(side) for side in sides]
        except ValueError as error:
            selection.at[index, "MappingCoverage"] = "Unresolved"
            selection.at[index, "ReviewRequired"] = True
            selection.at[index, "SelectionReason"] += "; equation cannot be parsed: " + str(error)
            continue
        if reaction.status == "candidate" and reaction.mapped_smiles:
            coverage = "Unresolved"
        elif reaction.status == "mapped" and reaction.mapped_smiles:
            coverage = "Complete"
        elif reaction.status == "partial" and reaction.mapped_smiles:
            coverage = "Partial"
            records.extend(reaction.partial_maps)
        elif not all(sides):
            coverage = "Boundary"
            selection.at[index, "SelectedMapper"] = "Boundary"
            selection.at[index, "ReviewRequired"] = False
            selection.at[index, "SelectionReason"] = "Open input/output boundary; no internal atom correspondence required"
        else:
            reaction.partial_maps = rule_transfers(reaction, usable_smiles)
            coverage = "Partial" if reaction.partial_maps else "Unresolved"
            if reaction.partial_maps:
                reaction.status = "partial"
                reaction.mapped_smiles = partial_mapped_smiles(reaction.partial_maps)
                records.extend(reaction.partial_maps)
                selection.at[index, "SelectedMapper"] = "CarbonRules"
                selection.at[index, "ReviewRequired"] = True
                selection.at[index, "SelectionReason"] = (
                    "All mapper predictions are empty; CarbonRules fallback: "
                    + "; ".join(row["Rule"] + ": " + row["SourceMetabolite"] + " -> "
                                + row["ProductMetabolite"] for row in reaction.partial_maps)
                    + "; reduced unit carbon-mapping subreaction; original amounts retained in Mapping_Info"
                    + "; remaining carbon unresolved; review required")
            involved_issues = sorted({name for side in terms
                                      for _, name in side if name in issues})
            if involved_issues and reaction.partial_maps:
                selection.at[index, "SelectionReason"] += (
                    "; input carbon formula/structure mismatch: " + "; ".join(
                        name + " (" + issues[name] + ")" for name in involved_issues))
            pools = sorted({name for side in terms
                            for _, name in side if generic_pool(name) and not smiles.get(name)})
            if pools:
                selection.at[index, "SelectionReason"] += (
                    "; generic pools have no positional structure: " + ", ".join(pools))
                product_names = {name for _, name in terms[1 if reaction.flux >= 0 else 0]}
                terminal = [name for name in pools if name in product_names and name not in internal_consumers]
                if terminal:
                    selection.at[index, "SelectionReason"] += (
                        "; terminal unresolved pools in the full nonzero network: " + ", ".join(terminal))
        limitations = sorted({name + ": " + STRUCTURE_LIMITS[kind]
                              for side in terms for _, name in side
                              for kind in [kinds.get(name, "missing")]
                              if kind in ("partial", "polymer", "invalid")})
        if limitations:
            selection.at[index, "SelectionReason"] += "; " + "; ".join(limitations)
            if coverage != "Boundary":
                selection.at[index, "ReviewRequired"] = True
        if not selection.at[index, "SelectedMapper"] and coverage == "Unresolved":
            selection.at[index, "SelectedMapper"] = "Unresolved"
        selection.at[index, "SelectedMappedReactionSMILES"] = reaction.mapped_smiles
        selection.at[index, "MappingCoverage"] = coverage
    fields = COLUMNS + (["Evidence", "Notes"] if any(row["Rule"] == "ManualVerified" for row in records) else [])
    return pd.DataFrame(records, columns=fields)


def load_partials(reactions, final, records, metabolites):
    """Validate saved structures, direction, amounts and rule version before use.

    Metadata is audit context, not permission to propagate: only a Partial final
    row activates its saved pairs. An Unresolved/rejected final row disables them.
    """
    if records.empty:
        if ("MappingCoverage" in final and (final.MappingCoverage.eq("Partial") &
                (final.MappedReactionSMILES.eq("") | final.get("SelectedMapper", pd.Series("", index=final.index)).eq("CarbonRules"))).any()):
            raise ValueError("Partial Final_Mapping rows require saved partial correspondence records")
        return
    if not set(COLUMNS) <= set(records.columns):
        raise ValueError("Saved partial correspondence records are missing required columns")
    if records[COLUMNS].duplicated().any():
        raise ValueError("Duplicated partial carbon correspondence")
    lookup = {r.reaction_id: r for r in reactions}
    if not set(records.ReactionID) <= set(lookup):
        raise ValueError("Unknown partial correspondence ReactionID")
    smiles = rule_smiles(metabolites)
    choices = final.set_index("ReactionID")
    for key, group in records.groupby("ReactionID", sort=False):
        if key not in choices.index or choices.loc[key].get("MappingCoverage", "") != "Partial":
            continue
        reaction = lookup[key]
        manual = group.Rule.eq("ManualVerified").any()
        fields = COLUMNS + (["Evidence", "Notes"] if manual else [])
        if not set(fields) <= set(group.columns):
            raise ValueError("Manual partial records require saved evidence")
        saved = group[fields].to_dict("records")
        expected = validate_manual_records(reaction, saved, smiles) if manual else rule_transfers(reaction, smiles)
        # Compare numeric values as numbers and JSON as positions, not cell spelling.
        def normalized(row):
            result = dict(row)
            result["RuleVersion"] = str(row["RuleVersion"])
            result["TransferCoefficient"] = float(row["TransferCoefficient"])
            result["CarbonAtomMap"] = json.loads(row["CarbonAtomMap"])
            return result
        if [normalized(r) for r in saved] != [normalized(r) for r in expected]:
            raise ValueError("Saved partial rules do not match input structures, amounts or rule version: " + key)
        reduced = partial_mapped_smiles(saved)
        if reaction.mapped_smiles and reaction.mapped_smiles != reduced:
            if not manual:
                raise ValueError("Saved reduced rule mapping does not match saved partial correspondences: " + key)
            derived = manual_mapping_records(reaction, reaction.mapped_smiles, smiles, "Saved manual mapping")
            def transfers(rows):
                return sorted(json.dumps({field: normalized(row)[field] for field in COLUMNS}, sort_keys=True)
                              for row in rows)
            if transfers(derived) != transfers(saved):
                raise ValueError("Saved manual mapping does not match saved partial correspondences: " + key)
        reaction.mapped_smiles = reaction.mapped_smiles if manual and reaction.mapped_smiles else reduced
        reaction.partial_maps = saved
        reaction.status = "partial"
    for key, row in choices.iterrows():
        if (row.get("MappingCoverage", "") == "Partial" and not row.MappedReactionSMILES
                and not lookup[key].partial_maps):
            raise ValueError("Missing partial correspondence for " + key)


def partial_reactions(reaction):
    """Use a unit transfer for reachability; never round fractional stoichiometry."""
    for index, row in enumerate(reaction.partial_maps):
        source, target = row["SourceMetabolite"], row["ProductMetabolite"]
        left = _structure(row["SourceSMILES"])
        right = _structure(row["ProductSMILES"])
        maps = [{a.GetIdx(): 0 for a in mol.GetAtoms()} for mol in (left, right)]
        for number, (_, i, _, j) in enumerate(json.loads(row["CarbonAtomMap"]), 1):
            maps[0][i] = maps[1][j] = number
        unit = PreparedReaction(reaction.reaction_id + "::carbon_rule:" + str(index),
                    reaction.reaction_name, 1, source + " => " + target, "",
                    [Component(source, 0, row["SourceSMILES"])],
                    [Component(target, 0, row["ProductSMILES"])], source_text=source,
                    status="partial", parent_reaction_id=reaction.reaction_id)
        yield unit, ([(unit.left[0], maps[0])], [(unit.right[0], maps[1])])
