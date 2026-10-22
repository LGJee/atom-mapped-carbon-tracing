"""Discover model structures within the source model's own directory."""

from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re

from rdkit import Chem

from .console import detail


def read_model_smiles(path):
    """Read scalar/first-list SMILES from COBRA/RAVEN ordered-map YAML exports."""
    lines = Path(path).read_text(encoding="utf-8-sig").splitlines()
    try:
        start = lines.index("- metabolites:") + 1
        end = lines.index("- reactions:", start)
    except ValueError:
        raise ValueError(
            "Expected model metabolite/reaction YAML sections: " + str(path)
        )
    structures, current = {}, None
    for position in range(start, end):
        match = re.match(
            r"^\s+- (id|name|compartment|formula|smiles):\s*(.*)$", lines[position]
        )
        if not match:
            continue
        field, value = match.groups()
        if field == "smiles" and not value and position + 1 < end:
            first = re.match(r"^\s+- (.+)$", lines[position + 1])
            value = first.group(1) if first else ""
        if value.startswith('"'):
            value = json.loads(value)
        if field == "id":
            if not value or value in structures:
                raise ValueError("Empty or duplicate model metabolite ID: " + value)
            current = value
            structures[current] = {}
        elif current:
            structures[current][field] = value
    if not structures:
        raise ValueError("No metabolite IDs in the configured model YAML: " + str(path))
    return structures


def first_text(value):
    if isinstance(value, (list, tuple)):
        value = value[0] if value else ""
    return value.strip() if isinstance(value, str) else ""


def compartment(key):
    match = re.search(r"\[([^\]]+)\]$|_([^_]+)$", key)
    return next((part for part in match.groups() if part), "") if match else ""


@lru_cache(maxsize=8)
def read_annotations(filename, modified, size, variable_name):
    """Cache small extracted records, not full numerical models, across non-IUPD/IUPD."""
    path = Path(filename)
    if path.suffix.lower() in (".yml", ".yaml"):
        return read_model_smiles(path)
    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text(encoding="utf-8-sig"))
        if not isinstance(data, dict) or not isinstance(data.get("metabolites"), list):
            return {}  # Source/license records are not model annotations.
        records = {}
        for row in data["metabolites"]:
            key = row["id"]
            if key in records:
                raise ValueError("Duplicate model metabolite ID: " + key)
            annotations = row.get("annotation") or {}
            if not isinstance(annotations, dict):
                annotations = {}
            inchi_key = annotations.get("inchi_key") or annotations.get("inchikey")
            if isinstance(inchi_key, (list, tuple)) and len(inchi_key) != 1:
                inchi_key = ""
            records[key] = dict(
                name=row.get("name", ""),
                compartment=row.get("compartment") or compartment(key),
                formula=row.get("formula", ""),
                smiles=first_text(row.get("smiles"))
                or first_text(annotations.get("smiles")),
                InChIKey=first_text(inchi_key),
            )
        return records
    if path.suffix.lower() == ".mat":
        import numpy as np
        from scipy.io import loadmat

        data = loadmat(str(path), squeeze_me=True, struct_as_record=False)
        if variable_name:
            model = data[variable_name]
        else:
            models = [
                value
                for value in data.values()
                if hasattr(value, "mets") and hasattr(value, "S")
            ]
            if len(models) != 1:
                raise ValueError("Expected one source model in " + str(path))
            model = models[0]
        ids = np.atleast_1d(model.mets).reshape(-1)

        def values(*fields):
            for field in fields:
                if hasattr(model, field):
                    result = np.atleast_1d(getattr(model, field)).reshape(-1)
                    if len(result) == len(ids):
                        return result
            return [""] * len(ids)

        names = values("metNames")
        formulas = values("metFormulas")
        smiles = values("metSmiles", "metSMILES", "smiles", "SMILES")
        keys = values("metInChIKey", "metInChIKeys")
        compartments = values("metComps")
        codes = np.atleast_1d(getattr(model, "comps", [])).reshape(-1)
        records = {}
        for key, name, formula, structure, inchi_key, code in zip(
            ids, names, formulas, smiles, keys, compartments
        ):
            key = first_text(key)
            if not key or key in records:
                raise ValueError("Empty or duplicate model metabolite ID: " + key)
            comp = compartment(key)
            if (
                not comp
                and isinstance(code, (int, float, np.number))
                and 1 <= code <= len(codes)
            ):
                comp = first_text(codes[int(code) - 1])
            records[key] = dict(
                name=first_text(name),
                compartment=comp,
                formula=first_text(formula),
                smiles=first_text(structure),
                InChIKey=first_text(inchi_key),
            )
        return records
    return {}


def compatible(annotation, row, require_name=False):
    from .carbon_mapping_rules import formulas_match
    name = annotation.get("name", "")
    if name == "D-Glucose":
        name = "D-glucose"
    comp = annotation.get("compartment", "")
    if name and comp and not name.endswith("[" + comp + "]"):
        name += "[" + comp + "]"
    formula = getattr(row, "MolecularFormula", "")
    formula = formula if isinstance(formula, str) else ""
    return (
        (bool(name) or not require_name)
        and (not name or name == row.MetaboliteName)
        and (
            not formula
            or not annotation.get("formula")
            or formulas_match(formula, annotation["formula"])
        )
    )


def discover_annotations(source, metabolites):
    """Prefer the recorded source; flag disagreements between compatible structures."""
    path = Path(source["path"])
    if hashlib.sha256(path.read_bytes()).hexdigest() != source["sha256"]:
        raise ValueError(
            "Source model changed after step 1; rerun step 1: " + str(path)
        )
    files = [path] + sorted(
        (
            item
            for item in path.parent.iterdir()
            if item.is_file()
            and item != path
            and item.suffix.lower() in (".yml", ".yaml", ".json")
        ),
        key=lambda item: (item.stem != path.stem, item.name),
    )
    candidates, identifier_candidates = {}, {}
    for item in files:
        try:
            stat = item.stat()
            records = read_annotations(
                str(item),
                stat.st_mtime_ns,
                stat.st_size,
                source.get("variable_name") if item == path else None,
            )
        except (ValueError, KeyError, TypeError, OSError) as error:
            print(
                "Skipping model annotations {}: {}".format(item.name, error), flush=True
            )
            continue
        matched = 0
        for row in metabolites.itertuples(index=False):
            annotation = records.get(row.MetaboliteID, {})
            if not annotation or not compatible(
                annotation, row, require_name=item != path
            ):
                continue
            if annotation.get("InChIKey"):
                identifier_candidates.setdefault(row.MetaboliteID, set()).add(
                    annotation["InChIKey"]
                )
            smiles = annotation.get("smiles", "")
            molecule = Chem.MolFromSmiles(smiles) if smiles else None
            from .carbon_mapping_rules import structure_kind
            if structure_kind(smiles) not in ("complete", "partial"):
                continue
            candidates.setdefault(row.MetaboliteID, []).append(
                dict(
                    smiles=smiles,
                    reference=item.name,
                    canonical=Chem.MolToSmiles(molecule),
                )
            )
            matched += 1
        if matched:
            detail(
                "Model structures: {} matched metabolites from {}".format(
                    matched, item.name
                ),
                flush=True,
            )
    chosen = {}
    for key, records in candidates.items():
        chosen[key] = dict(
            records[0],
            conflict=len({row["canonical"] for row in records}) > 1,
            candidates=json.dumps(
                [dict(SMILES=row["smiles"], File=row["reference"]) for row in records]
            ),
        )
    # Conflicting identifiers must not silently choose a different database compound.
    identifiers = {
        key: next(iter(values))
        for key, values in identifier_candidates.items()
        if len(values) == 1
    }
    return chosen, identifiers
