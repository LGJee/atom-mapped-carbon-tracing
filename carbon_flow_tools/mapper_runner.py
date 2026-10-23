"""Shared isolated-process runner for mapper adapters."""

import json
import os
import re
import subprocess
import tempfile
from pathlib import Path
from rdkit import Chem

from . import config
from .reactions import clean


def mapping_error(mapped):
    """Check mapping syntax, ID uniqueness and element consistency, not mechanism."""
    try:
        sides = mapped.split(">>")
        if len(sides) != 2 or not all(sides):
            return "Nonempty reactants and products separated by >> are required"
        atom_tables = []
        for side in sides:
            atoms = {}
            for smiles in side.split("."):
                molecule = Chem.MolFromSmiles(smiles)
                if molecule is None:
                    return "Invalid mapped SMILES"
                for atom in molecule.GetAtoms():
                    number = atom.GetAtomMapNum()
                    if not number:
                        continue
                    if number in atoms:
                        return "Repeated map ID {} on one reaction side".format(number)
                    atoms[number] = atom.GetSymbol()
            atom_tables.append(atoms)
        for number in sorted(atom_tables[0].keys() & atom_tables[1].keys()):
            left, right = atom_tables[0][number], atom_tables[1][number]
            if left != right:
                return "Element mismatch for map ID {}: {} -> {}".format(
                    number, left, right
                )
    except Exception as error:
        return str(error)
    return ""


def validate_mapper_config(name, settings):
    if not name or any(c in name for c in '/\\:*?"<>|') or name in (".", ".."):
        raise ValueError("Invalid mapper name for an output directory: " + str(name))
    module = settings.get("module", "")
    if not re.fullmatch(r"[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*", module):
        raise ValueError("Configure a valid Python module name for the mapper adapter.")
    python = Path(settings.get("python", ""))
    if not python.is_file():
        raise FileNotFoundError(
            "Python interpreter for {} not found: {}".format(name, python)
        )
    if not isinstance(settings.get("options", {}), dict):
        raise ValueError("Mapper options must be a dictionary.")


def group_mapper_jobs(pending):
    """Share one database scan across pending EnzymeMap jobs with identical settings."""
    tasks, groups = [], {}
    for strain, name, candidates, settings in pending:
        if settings["module"] != "carbon_flow_tools.mappers.enzymemap_adapter":
            tasks.append((strain, candidates, name, settings))
            continue
        signature = (name, str(settings["python"]),
                     json.dumps(settings.get("options", {}), sort_keys=True))
        if signature not in groups:
            shared = dict(settings, log_file=str(Path(settings["log_file"]).with_name("EnzymeMap.log")))
            task = ("All models", [], name, shared)
            groups[signature] = task
            tasks.append(task)
        groups[signature][1].extend(candidates)
    return tasks


def run_mapper(strain, reactions, name, settings, cpu_threads=1, memory_gb=6):
    """Run an isolated mapper process without biochemical corrections."""
    validate_mapper_config(name, settings)
    from .carbon_mapping_rules import structure_kind, STRUCTURE_LIMITS
    for reaction in reactions:
        if reaction.status != "ready":
            continue
        for side in reaction.reaction_smiles.split(">>"):
            kind = structure_kind(side)
            if kind not in ("complete", "partial"):
                reaction.status = "Mapper input unavailable: " + STRUCTURE_LIMITS[kind]
                break
    ready = [r for r in reactions if r.status == "ready"]
    if not ready:
        return
    # Exact strings retain stoichiometry, direction, charges, stereo and atom order.
    # Reuse predictions only inside this call, whose adapter/options are identical.
    unique = list(dict.fromkeys(r.reaction_smiles for r in ready))
    if len(unique) != len(ready):
        print("{} / {}: {} identical reaction inputs will share predictions."
              .format(strain, name, len(ready) - len(unique)), flush=True)
    worker = Path(__file__).with_name("mapper_worker.py")
    with tempfile.TemporaryDirectory() as folder:
        request, response = Path(folder) / "input.json", Path(folder) / "output.json"
        request.write_text(
            json.dumps(
                {
                    "reactions": unique,
                    "cpu_threads": cpu_threads,
                    "memory_gb": memory_gb,
                    "options": settings.get("options", {}),
                    "log_file": str(settings.get("log_file", Path(folder) / "mapper.log")),
                }
            ),
            encoding="utf-8",
        )
        environment = os.environ.copy()
        for key in (
            "OMP_NUM_THREADS",
            "MKL_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ):
            environment[key] = str(1 if settings["module"] ==
                "carbon_flow_tools.mappers.enzymemap_adapter" else cpu_threads)
        result = subprocess.run(
            [
                str(settings["python"]),
                str(worker),
                settings["module"],
                str(request),
                str(response),
                str(config.BATCH_SIZE),
                strain,
                name,
            ],
            check=False,
            env=environment,
        )
        if result.returncode:
            raise RuntimeError(
                "Atom-mapping worker failed (exit {}).\nModel: {}\nMapper: {}\nPython: {}\n"
                "Details: {}\nAction: check the log, or rerun step2_atom_mapping.py --verbose."
                .format(result.returncode, strain, name, settings["python"],
                        settings.get("log_file", "rerun with --verbose for details"))
            )
        results = json.loads(response.read_text(encoding="utf-8"))
    if not isinstance(results, list) or len(results) != len(unique):
        raise RuntimeError(name + " returned a result count inconsistent with input.")
    lookup = dict(zip(unique, results))
    for reaction in ready:
        result = lookup[reaction.reaction_smiles]
        if not isinstance(result, dict):
            raise RuntimeError(name + " must return a dictionary per reaction.")
        reaction.mapped_smiles = clean(result.get("mapped_rxn"))
        reaction.confident = result.get("confidence", "")
        reaction.template = clean(result.get("template"))
        error = mapping_error(reaction.mapped_smiles) if reaction.mapped_smiles else ""
        if reaction.mapped_smiles and not error:
            error = clean(result.get("error"))
        reaction.status = (
            "Invalid mapping: " + error
            if error
            else "mapped"
            if reaction.mapped_smiles
            else name + " failed: " + str(result.get("error", "Empty result"))
        )
