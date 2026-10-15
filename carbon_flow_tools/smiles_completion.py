"""Prepare structures from manual overrides, scoped model YAML and PubChem."""

import csv
import json
from pathlib import Path
import re
import time
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

import pandas as pd
from rdkit import Chem

from .reactions import clean, get_smiles

API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
PROPERTIES = "SMILES,MolecularFormula,Charge,InChIKey"


def load_manual(path):
    """Require model-scoped identities so local IDs cannot collide across models."""
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"Model", "MetaboliteID", "MetaboliteName", "SMILES"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError(
                "Manual structure CSV requires: " + ", ".join(sorted(required))
            )
        rows = list(reader)
    for row in rows:
        if not row["Model"] or not (row["MetaboliteID"] or row["MetaboliteName"]):
            raise ValueError(
                "Manual structures require a model and an ID or qualified name."
            )
        molecule = Chem.MolFromSmiles(row["SMILES"])
        if (
            molecule is None
            or not molecule.GetNumAtoms()
            or any(atom.GetAtomicNum() == 0 for atom in molecule.GetAtoms())
        ):
            raise ValueError("Invalid or unspecified manual structure: " + str(row))
    return rows


def read_model_smiles(path):
    """Read scalar/first-list SMILES from supported COBRA and RAVEN YAML exports."""
    lines = Path(path).read_text(encoding="utf-8-sig").splitlines()
    try:
        start = lines.index("- metabolites:") + 1
        end = lines.index("- reactions:", start)
    except ValueError:
        raise ValueError("Expected model metabolite/reaction YAML sections: " + str(path))
    structures, current = {}, None
    for position in range(start, end):
        match = re.match(r"^\s+- (id|name|compartment|formula|smiles):\s*(.*)$", lines[position])
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


def query_key(row):
    """Prefer optional exact identifiers; names are an explicitly reviewed fallback."""
    for field, namespace in (("PubChemCID", "cid"), ("InChIKey", "inchikey")):
        value = clean(getattr(row, field, ""))
        if value:
            if namespace == "cid":
                if value.endswith(".0"):
                    value = value[:-2]
                if not value.isdigit():
                    raise ValueError("PubChemCID must be a positive integer: " + value)
                if int(value) < 1:
                    raise ValueError("PubChemCID must be positive.")
            return namespace, value
    return "name", re.sub(r"\[[^\[\]]+\]$", "", clean(row.MetaboliteName)).strip()


def fetch_properties(namespace, identifier, timeout):
    url = "{}/{}/{}/property/{}/JSON".format(
        API, namespace, quote(identifier, safe=""), PROPERTIES
    )
    request = Request(url, headers={"User-Agent": "atom-mapped-carbon-tracing"})
    try:
        with urlopen(request, timeout=timeout) as response:
            payload = json.load(response)
    except HTTPError as error:
        if error.code == 404:
            return []
        raise
    return payload["PropertyTable"]["Properties"]


def select_candidate(properties, row):
    """Use available formula/charge constraints; resolve remaining name ties by CID."""
    formula = re.sub(r"[+-]\d*$", "", clean(getattr(row, "MolecularFormula", "")))
    charge = getattr(row, "Charge", "")
    candidates = []
    for record in properties:
        smiles = record.get("SMILES") or record.get("IsomericSMILES", "")
        molecule = Chem.MolFromSmiles(smiles) if smiles else None
        if molecule is None or not molecule.GetNumAtoms():
            continue
        if any(atom.GetAtomicNum() == 0 for atom in molecule.GetAtoms()):
            continue
        database_formula = re.sub(
            r"[+-]\d*$", "", clean(record.get("MolecularFormula", ""))
        )
        if formula and formula != database_formula:
            continue
        if clean(charge) and float(charge) != record.get("Charge"):
            continue
        candidates.append(dict(record, SMILES=smiles))
    candidates.sort(key=lambda record: int(record["CID"]))
    return candidates


def complete_smiles(metabolites, model, options):
    """Return an enriched copy; never write into the input workbook."""
    timeout = options["timeout"]
    interval = options["min_interval"]
    if timeout <= 0 or interval < 0:
        raise ValueError(
            "Structure timeout must be positive; query interval cannot be negative."
        )
    manual = load_manual(options["manual_file"])
    cache_file = Path(options["cache_file"])
    cache = (
        json.loads(cache_file.read_text(encoding="utf-8"))
        if options["lookup"] and cache_file.is_file()
        else {}
    )
    if not isinstance(cache, dict):
        raise ValueError("Structure cache must be a JSON object: " + str(cache_file))
    yaml_file = options.get("model_yaml_files", {}).get(model)
    annotations = read_model_smiles(yaml_file) if yaml_file else {}
    table = metabolites.copy()
    if annotations:
        if "SMILES_Model" not in table:
            table["SMILES_Model"] = float("nan")
        table["SMILES_Model"] = table.SMILES_Model.astype(object)
    if "SMILES_Manual" not in table:
        table["SMILES_Manual"] = float("nan")
    table["SMILES_Manual"] = table.SMILES_Manual.astype(object)
    for field in (
        "SMILES_Retrieved",
        "StructureSource",
        "StructureReference",
        "StructureCandidates",
        "StructureStatus",
    ):
        if field not in table:
            table[field] = pd.Series(float("nan"), index=table.index, dtype=object)
    if "StructureReviewRequired" not in table:
        table["StructureReviewRequired"] = False
    last_request = 0.0
    fetched, hits, filled = 0, 0, 0
    for position, (index, row) in enumerate(
        zip(table.index, table.itertuples(index=False)), 1
    ):
        if position == 1 or position % 25 == 0:
            print(
                "{} structures {}/{}: {} additions, {} requests, {} cache hits".format(
                    model, position, len(table), filled, fetched, hits
                ),
                flush=True,
            )
        matches = [
            entry
            for entry in manual
            if entry["Model"] == model
            and (
                not entry["MetaboliteID"]
                or entry["MetaboliteID"] == clean(row.MetaboliteID)
            )
            and (
                not entry["MetaboliteName"]
                or entry["MetaboliteName"] == clean(row.MetaboliteName)
            )
        ]
        if len({entry["SMILES"] for entry in matches}) > 1:
            raise ValueError(
                "Conflicting manual structures for {} {}".format(
                    model, row.MetaboliteID
                )
            )
        if matches:
            table.at[index, "SMILES_Manual"] = matches[0]["SMILES"]
            table.at[index, "StructureSource"] = "manual table"
            table.at[index, "StructureReference"] = Path(options["manual_file"]).name
            table.at[index, "StructureReviewRequired"] = False
            table.at[index, "StructureStatus"] = "Verified override"
            continue
        existing = get_smiles(row)
        annotation = annotations.get(clean(row.MetaboliteID), {})
        # Scoped IDs alone do not establish identity for companion annotations.
        expected_name = annotation.get("name", "")
        if annotation.get("compartment"):
            expected_name += "[" + annotation["compartment"] + "]"
        compatible = (
            (not expected_name or expected_name == clean(row.MetaboliteName))
            and (
                not annotation.get("formula")
                or not clean(getattr(row, "MolecularFormula", ""))
                or annotation["formula"] == clean(row.MolecularFormula)
            )
        )
        if not existing and compatible and annotation.get("smiles"):
            existing = annotation["smiles"]
            table.at[index, "SMILES_Model"] = existing
            table.at[index, "StructureSource"] = "Model YAML"
            table.at[index, "StructureReference"] = Path(yaml_file).name
        if existing:
            molecule = Chem.MolFromSmiles(existing)
            table.at[index, "StructureStatus"] = (
                "Existing structure" if molecule else "Invalid existing structure"
            )
            if not clean(table.at[index, "StructureSource"]):
                table.at[index, "StructureSource"] = "input workbook"
            continue
        table.at[index, "StructureStatus"] = "Missing structure; lookup disabled"
        if not options["lookup"]:
            continue
        try:
            namespace, identifier = query_key(row)
            key = namespace + ":" + identifier
            if key in cache:
                properties = cache[key]
                hits += 1
            else:
                time.sleep(max(0, interval - (time.monotonic() - last_request)))
                last_request = time.monotonic()
                fetched += 1
                properties = fetch_properties(namespace, identifier, timeout)
                cache[key] = properties
            candidates = select_candidate(properties, row)
            if not candidates:
                table.at[index, "StructureStatus"] = "No compatible PubChem structure"
                continue
            chosen = candidates[0]
            table.at[index, "SMILES_Retrieved"] = chosen["SMILES"]
            table.at[index, "StructureSource"] = "PubChem " + namespace
            table.at[index, "StructureReference"] = "PubChem:{}".format(chosen["CID"])
            table.at[index, "StructureCandidates"] = "; ".join(
                "PubChem:{}".format(item["CID"]) for item in candidates
            )
            table.at[index, "StructureReviewRequired"] = (
                namespace == "name" or len(candidates) > 1
            )
            table.at[index, "StructureStatus"] = (
                "Automatically selected; review required"
                if table.at[index, "StructureReviewRequired"]
                else "Identifier matched"
            )
            filled += 1
        except (
            HTTPError,
            URLError,
            TimeoutError,
            OSError,
            ValueError,
            KeyError,
        ) as error:
            table.at[index, "StructureStatus"] = "Lookup failed: " + str(error)
    if fetched and cache:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        temporary = cache_file.with_suffix(".tmp.json")
        temporary.write_text(json.dumps(cache, sort_keys=True), encoding="utf-8")
        temporary.replace(cache_file)
    print(
        "{} structures: {} database additions; {} requests, {} cache hits".format(
            model, filled, fetched, hits
        ),
        flush=True,
    )
    return table.infer_objects()
