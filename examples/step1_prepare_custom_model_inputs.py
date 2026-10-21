"""Independent example: glucose versus simultaneous glucose/xylose growth pFBA."""

from pathlib import Path
import importlib.util
import argparse
import os
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = next((parent for parent in SCRIPT_DIR.parents
             if (parent / "carbon_flow_tools/model_preparation.py").is_file()), None)
if ROOT is None:
    raise SystemExit("Keep this script inside the project, e.g. results/results_by_AI/my_case/.")
sys.path.insert(0, str(ROOT))

# Output locations follow the case folder; copying the script creates a new case.
CASE_DIR = (SCRIPT_DIR if SCRIPT_DIR != ROOT / "examples"
            else ROOT / "results/results_by_AI/ecoli_glucose_xylose")
PLAN_FILE = (CASE_DIR / "analysis_plan.json").relative_to(ROOT).as_posix()
SOLVER = "glpk"
NONZERO_TOLERANCE = 1e-8


def prepare_models():
    from cobra.flux_analysis import pfba

    from carbon_flow_tools.model_preparation import (
        close_other_carbon_uptake,
        read_model,
        metabolite_labels,
    )

    models = {}

    # 1. Escherichia_coli_glucose: illustrative single-carbon-source reference.
    model_file = "data/Escherichia_coli/iML1515.json"
    model = read_model(ROOT / model_file, SOLVER)
    labels = metabolite_labels(model)
    carbon_sources = [labels["glc__D_e"]]  # D-glucose, extracellular.
    close_other_carbon_uptake(model, "EX_glc__D_e")
    model.reactions.get_by_id("EX_glc__D_e").lower_bound = -10  # Glucose exchange.
    # Inherit oxygen and other inorganic nutrients from this model.
    model.objective = "BIOMASS_Ec_iML1515_core_75p37M"  # Core biomass reaction.
    model.objective_direction = "max"
    solution = pfba(model, fraction_of_optimum=1.0)
    models["Escherichia_coli_glucose"] = (model, solution, carbon_sources)

    # 2. Escherichia_coli_glucose_xylose: simultaneous supply, not sequential feeding.
    model_file = "data/Escherichia_coli/iML1515.json"
    model = read_model(ROOT / model_file, SOLVER)
    labels = metabolite_labels(model)
    carbon_sources = [labels["glc__D_e"], labels["xyl__D_e"]]
    close_other_carbon_uptake(model, ["EX_glc__D_e", "EX_xyl__D_e"])
    # Equal allowed carbon input to the reference: 5*6 + 6*5 = 10*6.
    # These uptake capacities are examples, not measured fluxes or forced uptake.
    model.reactions.get_by_id("EX_glc__D_e").lower_bound = -5  # Glucose exchange.
    model.reactions.get_by_id("EX_xyl__D_e").lower_bound = -6  # D-xylose exchange.
    model.objective = "BIOMASS_Ec_iML1515_core_75p37M"  # Core biomass reaction.
    model.objective_direction = "max"
    solution = pfba(model, fraction_of_optimum=1.0)
    models["Escherichia_coli_glucose_xylose"] = (model, solution, carbon_sources)
    return models


def main(argv=None):
    arguments = sys.argv[1:] if argv is None else list(argv)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", default=PLAN_FILE, help="Generated case plan path")
    parser.add_argument("--output-dir", help="Project-relative output folder")
    args = parser.parse_args(arguments)
    # Share the project's preparation environment; do not create one per case.
    from step1_prepare_carbon_tracing_inputs import PREPARATION_PYTHON
    python = ROOT / PREPARATION_PYTHON
    if python.is_file() and python.resolve() != Path(sys.executable).resolve():
        os.execv(str(python), [str(python), str(Path(__file__).resolve())] + arguments)
    if importlib.util.find_spec("cobra") is None:
        raise SystemExit("COBRApy is unavailable. Run bash scripts/setup_environments.sh first.")
    from carbon_flow_tools.model_preparation import export_model, save_plan, invalidate_plan

    target = invalidate_plan(ROOT, args.plan)
    output_dir = args.output_dir or target.parent.relative_to(ROOT).as_posix()
    entries = {}
    for name, (model, solution, sources) in prepare_models().items():
        entries[name] = export_model(model, solution, name, sources, ROOT,
                                     output_dir, NONZERO_TOLERANCE)
    save_plan(ROOT, entries, output_dir, target, script_file=__file__)
    selected_plan = target.relative_to(ROOT).as_posix()
    print("Next: python3 step2_atom_mapping.py --plan {}\n"
          "Then: python3 step3_carbon_flow.py --plan {}".format(selected_plan, selected_plan))


if __name__ == "__main__":
    main()
