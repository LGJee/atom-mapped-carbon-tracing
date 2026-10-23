"""Prepare structures from manual overrides, source models and PubChem."""

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
from rdkit.Chem import rdMolDescriptors

from .carbon_mapping_rules import generic_pool, structure_kind, STRUCTURE_LIMITS, formulas_match, protonation_compatible
from .console import Progress, detail
from .config import MemoryGuard
from .reactions import clean, get_smiles
from .model_structures import compatible, discover_annotations, read_model_smiles

API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
PROPERTIES = "SMILES,MolecularFormula,Charge,InChIKey"


def metabolite_name(value):
    """Use the compound name without a terminal compartment suffix."""
    return re.sub(r"\[[^\[\]]+\]$", "", clean(value)).strip()


def species_name(analysis, source=None):
    """Prefer a full source-directory species name, otherwise its analysis prefix."""
    names = [Path(source["path"]).parent.name] if source else []
    names.append(clean(analysis))
    matches = [re.match(r"^([A-Z][a-z]+_[a-z][a-z0-9]+)(?:_|$)", name) for name in names]
    found = [match.group(1) for match in matches if match]
    if len(set(found)) > 1:
        raise ValueError("Analysis species conflicts with the source model directory")
    return found[0] if found else clean(analysis)


def load_manual(path, with_source=False, species=None, names=None):
    """Read species/name structures; local metabolite IDs are not shared identifiers."""
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"Species", "MetaboliteName", "SMILES"} | ({"Source"} if with_source else set())
        if set(reader.fieldnames or []) != required or len(reader.fieldnames) != len(required):
            raise ValueError(
                "Structure CSV requires exactly: Species,MetaboliteName,SMILES"
                + (",Source" if with_source else "")
                + "\nFile: " + str(Path(path).resolve())
                + "\nFound columns: " + ", ".join(reader.fieldnames or [])
                + "\nAction: use the required CSV header; keep one verified structure per species/metabolite."
            )
        rows, guard = [], MemoryGuard()
        for index, row in enumerate(reader):
            if index % 1024 == 0:
                guard.check("reading structure CSV")
            entry = {field: clean(row.get(field)) for field in required}
            # Filter while streaming, rather than retaining unrelated species first.
            if ((species is None or entry["Species"] == species)
                    and (names is None or metabolite_name(entry["MetaboliteName"]) in names)):
                rows.append((reader.line_num, entry))
    structures = {}
    for line, row in rows:
        if not row["Species"] or not row["MetaboliteName"]:
            raise ValueError(
                "Manual structures require a full species name and metabolite name.\n"
                "File: {}\nLine: {}\nAction: fill Species and MetaboliteName."
                .format(Path(path).resolve(), line)
            )
        molecule = Chem.MolFromSmiles(row["SMILES"])
        kind = structure_kind(row["SMILES"])
        if kind in ("invalid", "missing"):
            raise ValueError(
                "Invalid manual structure.\nFile: {}\nLine: {}\n"
                "Species: {}\nMetabolite: {}\nSMILES: {}\n"
                "Action: correct the SMILES syntax; legal wildcard structures are accepted as partial."
                .format(Path(path).resolve(), line, row["Species"], row["MetaboliteName"], row["SMILES"])
            )
        key = (row["Species"], metabolite_name(row["MetaboliteName"]))
        canonical = Chem.MolToSmiles(molecule) if molecule is not None else row["SMILES"]
        if key in structures and structures[key][0] != canonical:
            if kind == "polymer" or structure_kind(structures[key][2]["SMILES"]) == "polymer":
                raise ValueError("Conflicting polymer/structure entries: " + str(key))
            _, previous_line, previous = structures[key]
            previous_molecule = Chem.MolFromSmiles(previous["SMILES"])
            raise ValueError(
                "Conflicting SMILES entries\n"
                "File: {}\nSpecies: {}\nMetabolite: {}\n\n"
                "Line {}: {}\n  Formula: {}; formal charge: {}\n"
                "Line {}: {}\n  Formula: {}; formal charge: {}\n\n"
                "These entries represent different structures, not just different SMILES spellings.\n"
                "Check connectivity, stereochemistry and protonation; keep one verified structure "
                "for this species/metabolite. No conflicting structure was selected; "
                "the source CSV was not modified."
                .format(Path(path).resolve(), key[0], key[1],
                        previous_line, previous["SMILES"],
                        rdMolDescriptors.CalcMolFormula(previous_molecule),
                        Chem.GetFormalCharge(previous_molecule),
                        line, row["SMILES"], rdMolDescriptors.CalcMolFormula(molecule),
                        Chem.GetFormalCharge(molecule))
            )
        if key not in structures:
            structures[key] = (canonical, line, row)
    return [row for _, row in rows]


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
    return "name", metabolite_name(row.MetaboliteName)


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
        if structure_kind(smiles, formula) != "complete":
            continue
        database_formula = re.sub(
            r"[+-]\d*$", "", clean(record.get("MolecularFormula", ""))
        )
        if formula and not formulas_match(formula, database_formula):
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
    source = options.get("source_models", {}).get(model)
    species = species_name(model, source)
    names = {metabolite_name(name) for name in metabolites.MetaboliteName}
    manual = load_manual(options["manual_file"], species=species, names=names)
    ai_file = Path(options["manual_file"]).with_name("ai_smiles.csv")
    ai_structures = load_manual(ai_file, with_source=True, species=species, names=names) if ai_file.is_file() else []
    cache_file = Path(options["cache_file"])
    guard = MemoryGuard()
    if options["lookup"] and cache_file.is_file():
        guard.check_file(cache_file, "loading structure cache")
    cache = (
        json.loads(cache_file.read_text(encoding="utf-8"))
        if options["lookup"] and cache_file.is_file()
        else {}
    )
    if not isinstance(cache, dict):
        raise ValueError("Structure cache must be a JSON object: " + str(cache_file))
    yaml_file = options.get("model_yaml_files", {}).get(model)
    annotations = read_model_smiles(yaml_file) if yaml_file else {}
    discovered, identifiers = (
        discover_annotations(source, metabolites) if source else ({}, {})
    )
    table = metabolites.copy()
    if identifiers and "InChIKey" not in table:
        table["InChIKey"] = ""
    if identifiers:
        table["InChIKey"] = table.InChIKey.astype(object)
    if annotations or discovered:
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
    # Cache hits and remote queries have unrelated costs; a shared ETA is misleading.
    progress = Progress(model + " | structures", len(table), "metabolites", estimate=False)
    last_request = 0.0
    fetched, hits, filled = 0, 0, 0
    for position, (index, row) in enumerate(
        zip(table.index, table.itertuples(index=False)), 1
    ):
        guard.check("completing compound structures")
        progress.update(position - 1)
        if position == 1 or position % 25 == 0:
            detail(
                "{} structures {}/{}: {} additions, {} requests, {} cache hits".format(
                    model, position, len(table), filled, fetched, hits
                ),
                flush=True,
            )
        matches = [entry for entry in manual if entry["Species"] == species
                   and metabolite_name(entry["MetaboliteName"]) == metabolite_name(row.MetaboliteName)]
        ai_match = False
        if not matches and not get_smiles(row):
            matches = [entry for entry in ai_structures if entry["Species"] == species
                       and metabolite_name(entry["MetaboliteName"]) == metabolite_name(row.MetaboliteName)]
            ai_match = bool(matches)
        if matches:
            molecule = Chem.MolFromSmiles(matches[0]["SMILES"])
            formula = re.sub(r"[+-]\d*$", "", clean(getattr(row, "MolecularFormula", "")))
            kind = structure_kind(matches[0]["SMILES"], formula)
            observed = re.sub(r"[+-]\d*$", "", rdMolDescriptors.CalcMolFormula(molecule)) if molecule is not None else ""
            charge = clean(getattr(row, "Charge", ""))
            observed_charge = Chem.GetFormalCharge(molecule) if molecule is not None else None
            protonation_note = ""
            if kind == "complete" and protonation_compatible(formula, charge, observed, observed_charge):
                protonation_note = (
                    "Protonation-compatible mismatch: model {} (charge {}); SMILES {} (charge {}); "
                    "original structure retained; review required"
                ).format(formula, charge, observed, observed_charge)
                print("Warning: {} / {}: {}".format(model, clean(row.MetaboliteName), protonation_note), flush=True)
            if (kind == "complete" and not protonation_note and formula and re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", formula)
                    and not re.search(r"[RX]", formula) and not formulas_match(formula, observed)):
                raise ValueError(
                    "Manual structure/formula conflict.\nFile: {}\nSpecies: {}\nMetabolite: {}\n"
                    "Input formula: {}; SMILES formula: {}\n"
                    "Action: verify compound identity and chemical form; do not bypass the formula check."
                    .format(ai_file if ai_match else options["manual_file"], species,
                            clean(row.MetaboliteName), formula, observed)
                )
            if kind == "complete" and not protonation_note and charge and float(charge) != observed_charge:
                raise ValueError(
                    "Manual structure/charge conflict.\nFile: {}\nSpecies: {}\nMetabolite: {}\n"
                    "Input charge: {}; SMILES formal charge: {}\n"
                    "Action: verify the protonation/charge state before rerunning."
                    .format(ai_file if ai_match else options["manual_file"], species,
                            clean(row.MetaboliteName), charge, Chem.GetFormalCharge(molecule))
                )
            table.at[index, "SMILES_Retrieved" if ai_match else "SMILES_Manual"] = matches[0]["SMILES"]
            table.at[index, "StructureSource"] = "AI structure table" if ai_match else "manual table"
            table.at[index, "StructureReference"] = matches[0]["Source"] if ai_match else Path(options["manual_file"]).name
            table.at[index, "StructureReviewRequired"] = ai_match or bool(protonation_note)
            table.at[index, "StructureStatus"] = protonation_note or (
                "AI-assisted structure; review required" if ai_match else "Verified override")
            if ai_match:
                table.at[index, "AIStructureCompleted"] = True
            continue
        existing = get_smiles(row)
        annotation = annotations.get(clean(row.MetaboliteID), {})
        if not existing and compatible(annotation, row) and annotation.get("smiles"):
            existing = annotation["smiles"]
            table.at[index, "SMILES_Model"] = existing
            table.at[index, "StructureSource"] = "Model YAML"
            table.at[index, "StructureReference"] = Path(yaml_file).name
        found = discovered.get(clean(row.MetaboliteID))
        structure_conflict = bool(found and found["conflict"])
        if found:
            table.at[index, "StructureCandidates"] = found["candidates"]
            if not existing:
                existing = found["smiles"]
                table.at[index, "SMILES_Model"] = existing
                table.at[index, "StructureSource"] = "Source model directory"
                table.at[index, "StructureReference"] = found["reference"]
            else:
                molecule = Chem.MolFromSmiles(existing)
                structure_conflict = structure_conflict or bool(
                    molecule and Chem.MolToSmiles(molecule) != found["canonical"]
                )
            if structure_conflict:
                table.at[index, "StructureReviewRequired"] = True
        inchi_key = identifiers.get(clean(row.MetaboliteID))
        if inchi_key and not clean(getattr(row, "InChIKey", "")):
            table.at[index, "InChIKey"] = inchi_key
        if existing:
            molecule = Chem.MolFromSmiles(existing)
            table.at[index, "StructureStatus"] = (
                "Conflicting model structures; review required"
                if structure_conflict
                else "Existing structure"
                if molecule
                else "Invalid existing structure"
            )
            if not clean(table.at[index, "StructureSource"]):
                table.at[index, "StructureSource"] = "input workbook"
            if generic_pool(clean(row.MetaboliteName)):
                table.at[index, "StructureReviewRequired"] = True
                table.at[index, "StructureStatus"] = "Generic pool annotation; verify explicit macromolecular identity"
            continue
        table.at[index, "StructureStatus"] = "Missing structure; lookup disabled"
        if not options["lookup"]:
            continue
        try:
            # Include identifiers supplied by the source model directory.
            row = table.loc[index]
            namespace, identifier = query_key(row)
            if namespace == "name" and generic_pool(clean(row.MetaboliteName)):
                table.at[index, "StructureStatus"] = "Generic biomass/macromolecule pool; no unique molecular structure"
                table.at[index, "StructureReviewRequired"] = True
                continue
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
    progress.update(len(table))
    if fetched and cache:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        temporary = cache_file.with_suffix(".tmp.json")
        temporary.write_text(json.dumps(cache, sort_keys=True), encoding="utf-8")
        temporary.replace(cache_file)
    detail(
        "{} structures: {} database additions; {} requests, {} cache hits".format(
            model, filled, fetched, hits
        ),
        flush=True,
    )
    for index, row in zip(table.index, table.itertuples(index=False)):
        kind = structure_kind(get_smiles(row), clean(getattr(row, "MolecularFormula", "")))
        table.at[index, "StructureType"] = kind
        if kind in STRUCTURE_LIMITS and kind != "missing":
            table.at[index, "StructureStatus"] = STRUCTURE_LIMITS[kind]
            table.at[index, "StructureReviewRequired"] = True
    return table.infer_objects()
