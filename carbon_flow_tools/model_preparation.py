"""COBRApy input loading, numerical checks and plain Excel export."""

import hashlib
import json
import logging
from pathlib import Path
import re

import cobra.io
from cobra.util.array import create_stoichiometric_matrix
from cobra.util.solver import linear_reaction_coefficients
import numpy as np
from openpyxl import Workbook
import pandas as pd

from .excel_io import save_output, write_sheet
from .console import detail, verbose


class _MatCompatibilityFilter(logging.Filter):
    """Hide known legacy-format notices, never unrelated warnings or errors."""

    def __init__(self):
        super().__init__()
        self.hidden = False

    def filter(self, record):
        routine = record.getMessage().startswith((
            "This model seems to have metCharge instead of metCharges field.",
            "No defined compartments in model ",
            "Using regular expression found the following compartments:",
        ))
        if routine and record.levelno < logging.ERROR:
            self.hidden = True
            return False
        return True


def read_model(path, solver="glpk", variable_name=None):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            "Model not found: {}. Supply a local model file; see README.md.".format(path)
        )
    suffix = path.suffix.lower()
    if suffix == ".mat":
        logger = logging.getLogger("cobra.io.mat")
        compatibility = _MatCompatibilityFilter()
        if not verbose():
            logger.addFilter(compatibility)
        try:
            model = cobra.io.load_matlab_model(str(path), variable_name=variable_name)
        finally:
            logger.removeFilter(compatibility)
        if compatibility.hidden:
            print("  Source model: legacy MAT fields/compartments adapted (--verbose for details).",
                  flush=True)
        # Retain source annotations that the generic MAT reader does not import.
        from scipy.io import loadmat

        raw = loadmat(path, simplify_cells=True)
        record = (
            raw[variable_name]
            if variable_name
            else next(
                value
                for value in raw.values()
                if isinstance(value, dict) and "S" in value
            )
        )
        if "version" in record:
            model.annotation["source_version"] = str(record["version"])
        if "inchis" in record:
            for key, inchi in zip(
                np.atleast_1d(record["mets"]), np.atleast_1d(record["inchis"])
            ):
                if isinstance(inchi, str) and inchi:
                    model.metabolites.get_by_id(str(key)).annotation["inchi"] = inchi
    elif suffix == ".json":
        model = cobra.io.load_json_model(str(path))
    elif suffix in (".xml", ".sbml"):
        model = cobra.io.read_sbml_model(str(path))
    elif suffix in (".yml", ".yaml"):
        model = cobra.io.load_yaml_model(str(path))
    else:
        raise ValueError("Use a COBRApy-readable MAT, JSON, SBML or YAML model.")
    model.solver = solver
    # Record provenance automatically; the user only supplies the model path.
    model._tracing_source = dict(
        path=str(path.resolve()),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        variable_name=variable_name,
    )
    model._tracing_original_reactions = {reaction.id for reaction in model.reactions}
    record_path = path.parent / "source.json"
    if record_path.is_file():
        model._tracing_source["source_record"] = dict(
            path=str(record_path.resolve()),
            sha256=hashlib.sha256(record_path.read_bytes()).hexdigest(),
        )
    return model


def check_solution(model, solution):
    values = solution.fluxes.reindex([r.id for r in model.reactions]).to_numpy()
    if solution.status != "optimal" or not np.isfinite(values).all():
        raise ValueError("Optimization did not return a finite optimal solution.")
    residual = create_stoichiometric_matrix(model, array_type="lil").tocsr().dot(values)
    if np.max(np.abs(residual), initial=0) >= 1e-6:
        raise ValueError("Steady-state residual exceeds tolerance.")
    if any(
        value < reaction.lower_bound - 1e-6 or value > reaction.upper_bound + 1e-6
        for reaction, value in zip(model.reactions, values)
    ):
        raise ValueError("A reaction flux violates its bounds.")
    objective = sum(
        coefficient * solution.fluxes[reaction.id]
        for reaction, coefficient in linear_reaction_coefficients(model).items()
    )
    if not np.isfinite(objective) or objective <= 0:
        raise ValueError("The configured biological objective must have positive flux.")
    return float(objective)


def glucose_bound(model, reaction_id):
    value = model.reactions.get_by_id(reaction_id).lower_bound
    if not np.isfinite(value) or value >= 0:
        raise ValueError("The source model must allow glucose uptake: " + reaction_id)
    return value


def close_other_carbon_uptake(model, glucose_exchange):
    """Keep one exchange or a collection of simultaneous carbon-source exchanges."""
    allowed = {glucose_exchange} if isinstance(glucose_exchange, str) else set(glucose_exchange)
    for name in allowed:
        model.reactions.get_by_id(name)
    for reaction in model.reactions:
        if reaction.id in allowed or len(reaction.metabolites) != 1:
            continue
        metabolite, coefficient = next(iter(reaction.metabolites.items()))
        if not re.search(r"C(?![a-z])", metabolite.formula or ""):
            continue
        if coefficient < 0 and reaction.lower_bound < 0:
            reaction.lower_bound = 0
        elif coefficient > 0 and reaction.upper_bound > 0:
            reaction.upper_bound = 0


def metabolite_labels(model):
    labels = {}
    for metabolite in model.metabolites:
        match = re.search(r"\[([^\]]+)\]$", metabolite.id)
        compartment = match.group(1) if match else metabolite.compartment
        if not compartment:
            match = re.search(r"_([^_]+)$", metabolite.id)
            compartment = match.group(1) if match else None
        if not compartment:
            raise ValueError("No identifiable compartment: " + metabolite.id)
        name = metabolite.name or metabolite.id
        if name == "D-Glucose":
            name = "D-glucose"
        labels[metabolite.id] = "{}[{}]".format(name, compartment)
    return labels


def export_model(
    model, solution, name, sources, project, output_dir="results", tolerance=1e-8
):
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name):
        raise ValueError(
            "Use an analysis name starting with a letter, followed by letters, digits or underscores."
        )
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("Nonzero tolerance must be finite and positive.")
    if (
        not isinstance(sources, list)
        or not sources
        or any(not isinstance(s, str) or not s.strip() for s in sources)
        or len(set(sources)) != len(sources)
    ):
        raise ValueError("Supply a nonempty ordered carbon-source list.")
    project = Path(project).resolve()
    directory = Path(str(output_dir).replace("\\", "/"))
    if directory.is_absolute() or re.match(r"^[A-Za-z]:", str(directory)):
        raise ValueError("Use a project-relative output directory.")
    directory = (project / directory).resolve()
    directory.relative_to(project)
    objective = check_solution(model, solution)
    active = [r for r in model.reactions if abs(solution.fluxes[r.id]) > tolerance]
    participating = {m.id for reaction in active for m in reaction.metabolites}
    labels = metabolite_labels(model)
    if any(source not in set(labels.values()) for source in sources):
        raise ValueError(
            "Every carbon source must be present in the model: " + name
        )
    # Keep an unused supplied source so later summaries can report zero reachability.
    source_ids = {key for key, label in labels.items() if label in sources}
    included = participating | source_ids
    indices = {m.id: i for i, m in enumerate(model.metabolites)}
    rows = []
    for reaction in active:
        sides = [[], []]
        for metabolite, coefficient in sorted(
            reaction.metabolites.items(), key=lambda item: indices[item[0].id]
        ):
            term = labels[metabolite.id]
            if abs(coefficient) != 1:
                term = "{:.15g} ".format(abs(coefficient)) + term
            sides[coefficient > 0].append(term)
        arrow = " <=> " if reaction.lower_bound < 0 < reaction.upper_bound else " => "
        equation = " + ".join(sides[0]) + arrow + " + ".join(sides[1])
        rows.append(
            [
                reaction.id,
                reaction.name or reaction.id,
                float(solution.fluxes[reaction.id]),
                equation,
            ]
        )
    reactions = pd.DataFrame(
        rows, columns=["ReactionID", "ReactionName", "Flux", "Equation"]
    )
    metabolites = []
    for metabolite in model.metabolites:
        if metabolite.id not in included:
            continue
        row = dict(
            MetaboliteID=metabolite.id,
            MetaboliteName=labels[metabolite.id],
            MolecularFormula=metabolite.formula or "",
            Charge=metabolite.charge,
        )
        for key, column in (("inchi_key", "InChIKey"), ("inchi", "InChI")):
            value = metabolite.annotation.get(key)
            if isinstance(value, str) and value:
                row[column] = value
        metabolites.append(row)
    workbook = Workbook()
    del workbook[workbook.active.title]
    write_sheet(workbook, "Nonzero", reactions)
    write_sheet(workbook, "Metabolites", pd.DataFrame(metabolites))
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / (name + "_1_pFBA_nonzero_reactions.xlsx")
    try:
        save_output(workbook, path)
    finally:
        workbook.close()
    objectives = linear_reaction_coefficients(model)
    direction = {"max": "maximize", "min": "minimize"}.get(
        model.objective_direction, model.objective_direction)
    print("  Model: {}\n  Solution: {}".format(name, solution.status), flush=True)
    for reaction, coefficient in objectives.items():
        label = "{} [{}]".format(reaction.name, reaction.id) if reaction.name else reaction.id
        print("  Target: {} ({})\n  Target flux: {:.9g} (model flux units)"
              .format(label, direction, float(solution.fluxes[reaction.id])), flush=True)
        if len(objectives) > 1 or coefficient != 1:
            print("  Objective coefficient: {:.9g}".format(coefficient), flush=True)
    if len(objectives) != 1 or next(iter(objectives.values())) != 1:
        print("  Weighted objective value: {:.9g}".format(objective), flush=True)
    print("  Carbon sources: {}\n  Exported network: {} reactions, {} metabolites\n  File: {}"
          .format(", ".join(sources), len(reactions), len(metabolites), path), flush=True)
    detail("  Total absolute flux (diagnostic): {:.9g}"
           .format(float(solution.fluxes.abs().sum())), flush=True)
    entry = dict(
        input_file=path.relative_to(project).as_posix(),
        carbon_sources=sources,
        input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    )
    exchanges = {
        reaction.id: dict(
            lower_bound=float(reaction.lower_bound),
            upper_bound=float(reaction.upper_bound),
            flux=float(solution.fluxes[reaction.id]),
        )
        for reaction in model.exchanges
    }
    source_uptake = {}
    for source_name in sources:
        ids = {key for key in source_ids if labels[key] == source_name}
        imported = {}
        for reaction_id in exchanges:
            reaction = model.reactions.get_by_id(reaction_id)
            net_import = sum(coefficient * solution.fluxes[reaction_id]
                             for metabolite, coefficient in reaction.metabolites.items()
                             if metabolite.id in ids)
            if any(metabolite.id in ids for metabolite in reaction.metabolites):
                imported[reaction_id] = float(net_import)
        net_import = sum(imported.values()) if imported else None
        source_uptake[source_name] = dict(
            metabolite_ids=sorted(ids), exchange_imports=imported,
            net_import=net_import, present_in_nonzero_network=bool(ids & participating),
        )
        if net_import is not None and net_import <= tolerance:
            print("{}: {} has no net uptake in this pFBA solution; uptake is not forced."
                  .format(name, source_name), flush=True)
    entry["preparation"] = dict(
        model_id=model.id, solver=model.solver.interface.__name__,
        objective_direction=model.objective_direction,
        objective_reactions={r.id: float(c) for r, c in linear_reaction_coefficients(model).items()},
        objective_value=objective, total_absolute_flux=float(solution.fluxes.abs().sum()),
        nonzero_tolerance=tolerance, carbon_supply_mode="simultaneous",
        source_uptake=source_uptake, exchange_conditions=exchanges,
        knocked_out_genes=sorted(g.id for g in model.genes if not g.functional),
    )
    if hasattr(model, "_tracing_original_reactions"):
        entry["preparation"]["added_reactions"] = sorted(
            {r.id for r in model.reactions} - model._tracing_original_reactions
        )
    if hasattr(model, "_tracing_source"):
        entry["source_model"] = dict(model._tracing_source)
        source = Path(entry["source_model"]["path"])
        try:
            entry["source_model"]["path"] = source.relative_to(project).as_posix()
        except ValueError:
            entry["source_model"]["path"] = source.as_posix()
        if "source_record" in entry["source_model"]:
            record = dict(entry["source_model"]["source_record"])
            try:
                record["path"] = Path(record["path"]).relative_to(project).as_posix()
            except ValueError:
                pass
            entry["source_model"]["source_record"] = record
    return entry


def plan_path(project, plan_file="data/analysis_plan.json"):
    """Keep generated plans inside their project, separate from model data."""
    project = Path(project).resolve()
    target = (project / plan_file).resolve()
    target.relative_to(project)
    if target.suffix.lower() != ".json":
        raise ValueError("The analysis plan must be a JSON file.")
    return target


def invalidate_plan(project, plan_file="data/analysis_plan.json"):
    """Invalidate this case only; never delete an unrelated model JSON file."""
    target = plan_path(project, plan_file)
    if target.exists():
        previous = json.loads(target.read_text(encoding="utf-8"))
        if not isinstance(previous, dict) or previous.get("version") not in (2, 3):
            raise ValueError("Refusing to remove an unrelated JSON file: " + str(target))
        target.unlink()
    return target


def save_plan(project, models, output_dir, plan_file="data/analysis_plan.json", script_file=None):
    if not models:
        raise ValueError("Export at least one model before saving the analysis plan.")
    target = plan_path(project, plan_file)
    if script_file is not None:
        script = Path(script_file).resolve()
        try:
            relative_script = script.relative_to(Path(project).resolve()).as_posix()
        except ValueError:
            relative_script = script.as_posix()
        provenance = dict(path=relative_script, sha256=hashlib.sha256(script.read_bytes()).hexdigest())
        models = {name: dict(entry, preparation_script=provenance) for name, entry in models.items()}
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(".tmp")
    try:
        temporary.write_text(
            json.dumps(
                dict(
                    version=3,
                    models=models,
                    output_dir=str(output_dir).replace("\\", "/"),
                ),
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        temporary.replace(target)
    finally:
        if temporary.exists():
            temporary.unlink()
    print("Analysis plan saved: " + str(target), flush=True)
    return target
