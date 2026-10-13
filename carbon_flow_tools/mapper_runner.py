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


def run_mapper(strain, reactions, name, settings, cpu_threads=1):
    """Run an isolated mapper process without biochemical corrections."""
    validate_mapper_config(name, settings)
    ready = [r for r in reactions if r.status == "ready"]
    if not ready:
        return
    worker = Path(__file__).with_name("mapper_worker.py")
    with tempfile.TemporaryDirectory() as folder:
        request, response = Path(folder) / "input.json", Path(folder) / "output.json"
        request.write_text(
            json.dumps(
                {
                    "reactions": [r.reaction_smiles for r in ready],
                    "options": settings.get("options", {}),
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
            environment[key] = str(cpu_threads)
        subprocess.run(
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
            check=True,
            env=environment,
        )
        results = json.loads(response.read_text(encoding="utf-8"))
    if not isinstance(results, list) or len(results) != len(ready):
        raise RuntimeError(name + " returned a result count inconsistent with input.")
    for reaction, result in zip(ready, results):
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
