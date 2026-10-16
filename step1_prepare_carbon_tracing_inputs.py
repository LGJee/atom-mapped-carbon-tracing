"""Prepare six reference networks with COBRApy, then export the analysis plan."""

import importlib.util
import argparse
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent

# User settings: model/input preparation only.
OUTPUT_DIR = "results"
SOLVER = "glpk"
PREPARATION_PYTHON = ".envs/preparation/bin/python"
NONZERO_TOLERANCE = 1e-8


def prepare_models():
    from cobra import Metabolite, Reaction
    from cobra.flux_analysis import pfba
    from carbon_flow_tools.model_preparation import (
        check_solution,
        close_other_carbon_uptake,
        read_model,
    )

    models = {}
    # Stoichiometry: negative coefficients consume substrates; positive ones form products.

    # ================================================================
    # 1. Saccharomyces_cerevisiae_non_IUPD
    # ================================================================
    # Load the model.
    model_file = "data/yeast9/yeast-GEM.mat"
    model_non_iupd = read_model(ROOT / model_file, SOLVER, "model")

    # Carbon sources (label introduction order).
    carbon_sources = ["D-glucose[e]"]

    # Culture conditions and model modifications.
    glucose_non_iupd = model_non_iupd.reactions.get_by_id(
        "r_1714"
    ).lower_bound  # D-glucose exchange.
    # Inherit native glucose, oxygen and remaining medium for growth optimization.

    # Growth constraint: calculate non-IUPD maximum, then fix 50%.
    model_non_iupd.objective = "r_2111"  # Biomass growth.
    growth_non_iupd = model_non_iupd.optimize(objective_sense="maximize")
    check_solution(model_non_iupd, growth_non_iupd)
    max_growth_saccharomyces = growth_non_iupd.objective_value
    fixed_growth_saccharomyces = 0.5 * max_growth_saccharomyces
    biomass = model_non_iupd.reactions.get_by_id("r_2111")
    biomass.lower_bound = fixed_growth_saccharomyces
    biomass.upper_bound = fixed_growth_saccharomyces
    model_non_iupd.reactions.get_by_id("r_1714").lower_bound = 2 * glucose_non_iupd

    # Production objective and pFBA.
    model_non_iupd.objective = "r_1012"  # Squalene synthase.
    # pFBA maximizes squalene first, then minimizes total absolute flux.
    model_non_iupd.objective_direction = "max"
    solution_non_iupd = pfba(model_non_iupd, fraction_of_optimum=1.0)

    # Register this analysis for export.
    models["Saccharomyces_cerevisiae_non_IUPD"] = (
        model_non_iupd,
        solution_non_iupd,
        carbon_sources,
    )

    # ================================================================
    # 2. Saccharomyces_cerevisiae_IUPD
    # ================================================================
    # Load the model.
    model_file = "data/yeast9/yeast-GEM.mat"
    # Load independently of the non-IUPD production constraints.
    model_iupd = read_model(ROOT / model_file, SOLVER, "model")

    # Carbon sources (label introduction order).
    carbon_sources = ["D-glucose[e]", "isoprenol[e]"]

    # Culture conditions and model modifications.
    glucose_iupd = model_iupd.reactions.get_by_id(
        "r_1714"
    ).lower_bound  # D-glucose exchange.
    model_iupd.genes.get_by_id("YML126C").knock_out()  # ERG13, HMG-CoA synthase.
    # Add the IU-pathway metabolites.
    isoprenol_e = Metabolite(
        "isoprenol[e]", name="isoprenol", formula="C5H10O", compartment="e"
    )
    isoprenol_c = Metabolite(
        "isoprenol[c]", name="isoprenol", formula="C5H10O", compartment="c"
    )
    isopentenyl_p = Metabolite(
        "isopentenyl_phosphate[c]",
        name="isopentenyl phosphate",
        formula="C5H11O4P",
        compartment="c",
    )
    model_iupd.add_metabolites([isoprenol_e, isoprenol_c, isopentenyl_p])
    isoprenol_exchange = Reaction(
        "Ex_isopentenol_e", lower_bound=2.09 * glucose_iupd, upper_bound=1000
    )
    # Isoprenol exchange; negative flux is uptake.
    isoprenol_exchange.add_metabolites({isoprenol_e: -1})
    isoprenol_transport = Reaction(
        "Ex_isopentenol_e_t", lower_bound=0, upper_bound=1000
    )
    # Isoprenol transport: extracellular -> cytosol.
    isoprenol_transport.add_metabolites({isoprenol_e: -1, isoprenol_c: 1})
    # s_0434: ATP; s_0394: ADP; s_0943: isopentenyl diphosphate (IPP), all cytosolic.
    atp = model_iupd.metabolites.get_by_id("s_0434")
    adp = model_iupd.metabolites.get_by_id("s_0394")
    ipp = model_iupd.metabolites.get_by_id("s_0943")
    # Isoprenol + ATP -> isopentenyl phosphate + ADP.
    isoprenol_kinase = Reaction("isopentenolK", lower_bound=0, upper_bound=1000)
    isoprenol_kinase.add_metabolites(
        {isoprenol_c: -1, atp: -1, isopentenyl_p: 1, adp: 1}
    )
    isoprenol_kinase.gene_reaction_rule = "YLR133W"
    # Isopentenyl phosphate + ATP -> IPP + ADP.
    isopentenyl_phosphate_kinase = Reaction(
        "PisopentenolK", lower_bound=0, upper_bound=1000
    )
    isopentenyl_phosphate_kinase.add_metabolites(
        {isopentenyl_p: -1, atp: -1, ipp: 1, adp: 1}
    )
    isopentenyl_phosphate_kinase.gene_reaction_rule = "AtIPK"
    model_iupd.add_reactions(
        [
            isoprenol_exchange,
            isoprenol_transport,
            isoprenol_kinase,
            isopentenyl_phosphate_kinase,
        ]
    )

    # Growth constraint: use the paired non-IUPD half-maximum rate.
    biomass = model_iupd.reactions.get_by_id("r_2111")
    biomass.lower_bound = fixed_growth_saccharomyces
    biomass.upper_bound = fixed_growth_saccharomyces

    # Production objective and pFBA.
    model_iupd.objective = "r_1012"  # Squalene synthase.
    # pFBA maximizes squalene first, then minimizes total absolute flux.
    model_iupd.objective_direction = "max"
    solution_iupd = pfba(model_iupd, fraction_of_optimum=1.0)

    # Register this analysis for export.
    models["Saccharomyces_cerevisiae_IUPD"] = (
        model_iupd,
        solution_iupd,
        carbon_sources,
    )

    # ================================================================
    # 3. Escherichia_coli_non_IUPD
    # ================================================================
    # Load the model.
    model_file = "data/ecoli/iML1515.mat"
    model_non_iupd = read_model(ROOT / model_file, SOLVER, "iML1515")

    # Carbon sources (label introduction order).
    carbon_sources = ["D-glucose[e]"]

    # Culture conditions and model modifications.
    glucose_non_iupd = model_non_iupd.reactions.get_by_id(
        "EX_glc__D_e"
    ).lower_bound  # D-glucose exchange.
    # Both comparators require functional SQS activity and a squalene output.
    squalene = Metabolite(
        "squalene[c]", name="squalene", formula="C30H50", charge=0, compartment="c"
    )
    sqs = Reaction(
        "SQS", name="Heterologous squalene synthase", lower_bound=0, upper_bound=1000
    )
    # RHEA:32295: 2 FPP + NADPH + H+ -> squalene + 2 PPi + NADP+.
    fpp = model_non_iupd.metabolites.get_by_id("frdp_c")  # Farnesyl diphosphate.
    nadph = model_non_iupd.metabolites.get_by_id("nadph_c")
    proton = model_non_iupd.metabolites.get_by_id("h_c")
    pyrophosphate = model_non_iupd.metabolites.get_by_id("ppi_c")
    nadp = model_non_iupd.metabolites.get_by_id("nadp_c")
    sqs.add_metabolites(
        {fpp: -2, nadph: -1, proton: -1, pyrophosphate: 2, nadp: 1, squalene: 1}
    )
    demand = Reaction(
        "DM_squalene", name="Squalene output", lower_bound=0, upper_bound=1000
    )
    demand.add_metabolites({squalene: -1})
    model_non_iupd.add_reactions([sqs, demand])
    close_other_carbon_uptake(model_non_iupd, "EX_glc__D_e")
    # Inherit native oxygen, inorganic medium and maintenance bounds.

    # Growth constraint: calculate non-IUPD maximum, then fix 50%.
    model_non_iupd.objective = "BIOMASS_Ec_iML1515_core_75p37M"  # Biomass growth.
    growth_non_iupd = model_non_iupd.optimize(objective_sense="maximize")
    check_solution(model_non_iupd, growth_non_iupd)
    max_growth_escherichia = growth_non_iupd.objective_value
    fixed_growth_escherichia = 0.5 * max_growth_escherichia
    biomass = model_non_iupd.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M")
    biomass.lower_bound = fixed_growth_escherichia
    biomass.upper_bound = fixed_growth_escherichia
    model_non_iupd.reactions.get_by_id("EX_glc__D_e").lower_bound = 2 * glucose_non_iupd

    # Production objective and pFBA.
    model_non_iupd.objective = "SQS"  # Introduced squalene synthase.
    # pFBA maximizes squalene first, then minimizes total absolute flux.
    model_non_iupd.objective_direction = "max"
    solution_non_iupd = pfba(model_non_iupd, fraction_of_optimum=1.0)

    # Register this analysis for export.
    models["Escherichia_coli_non_IUPD"] = (
        model_non_iupd,
        solution_non_iupd,
        carbon_sources,
    )

    # ================================================================
    # 4. Escherichia_coli_IUPD
    # ================================================================
    # Load the model.
    model_file = "data/ecoli/iML1515.mat"
    # Load independently of the non-IUPD production constraints.
    model_iupd = read_model(ROOT / model_file, SOLVER, "iML1515")

    # Carbon sources (label introduction order).
    carbon_sources = ["D-glucose[e]", "isoprenol[e]"]

    # Culture conditions and model modifications.
    glucose_iupd = model_iupd.reactions.get_by_id(
        "EX_glc__D_e"
    ).lower_bound  # D-glucose exchange.
    squalene = Metabolite(
        "squalene[c]", name="squalene", formula="C30H50", charge=0, compartment="c"
    )
    sqs = Reaction(
        "SQS", name="Heterologous squalene synthase", lower_bound=0, upper_bound=1000
    )
    # RHEA:32295: 2 FPP + NADPH + H+ -> squalene + 2 PPi + NADP+.
    fpp = model_iupd.metabolites.get_by_id("frdp_c")  # Farnesyl diphosphate.
    nadph = model_iupd.metabolites.get_by_id("nadph_c")
    proton = model_iupd.metabolites.get_by_id("h_c")
    pyrophosphate = model_iupd.metabolites.get_by_id("ppi_c")
    nadp = model_iupd.metabolites.get_by_id("nadp_c")
    sqs.add_metabolites(
        {fpp: -2, nadph: -1, proton: -1, pyrophosphate: 2, nadp: 1, squalene: 1}
    )
    demand = Reaction(
        "DM_squalene", name="Squalene output", lower_bound=0, upper_bound=1000
    )
    demand.add_metabolites({squalene: -1})
    model_iupd.add_reactions([sqs, demand])
    close_other_carbon_uptake(model_iupd, "EX_glc__D_e")
    model_iupd.genes.get_by_id("b2515").knock_out()  # ispG, MECDP reductase.
    # Add the IU-pathway metabolites.
    isoprenol_e = Metabolite(
        "isoprenol[e]", name="isoprenol", formula="C5H10O", charge=0, compartment="e"
    )
    isoprenol_c = Metabolite(
        "isoprenol[c]", name="isoprenol", formula="C5H10O", charge=0, compartment="c"
    )
    isopentenyl_p = Metabolite(
        "isopentenyl_phosphate[c]",
        name="isopentenyl phosphate",
        formula="C5H11O4P",
        charge=0,
        compartment="c",
    )
    model_iupd.add_metabolites([isoprenol_e, isoprenol_c, isopentenyl_p])
    isoprenol_exchange = Reaction(
        "Ex_isopentenol_e", lower_bound=2.09 * glucose_iupd, upper_bound=1000
    )
    isoprenol_exchange.add_metabolites({isoprenol_e: -1})  # Isoprenol exchange.
    isoprenol_transport = Reaction(
        "Ex_isopentenol_e_t", lower_bound=0, upper_bound=1000
    )
    # Isoprenol transport: extracellular -> cytosol.
    isoprenol_transport.add_metabolites({isoprenol_e: -1, isoprenol_c: 1})
    # atp_c: ATP; adp_c: ADP; ipdp_c: IPP, all cytosolic.
    atp = model_iupd.metabolites.get_by_id("atp_c")
    adp = model_iupd.metabolites.get_by_id("adp_c")
    ipp = model_iupd.metabolites.get_by_id("ipdp_c")
    # Isoprenol + ATP -> isopentenyl phosphate + ADP.
    isoprenol_kinase = Reaction("isopentenolK", lower_bound=0, upper_bound=1000)
    isoprenol_kinase.add_metabolites(
        {isoprenol_c: -1, atp: -1, isopentenyl_p: 1, adp: 1}
    )
    isoprenol_kinase.gene_reaction_rule = "YLR133W"
    # Isopentenyl phosphate + ATP -> IPP + ADP.
    isopentenyl_phosphate_kinase = Reaction(
        "PisopentenolK", lower_bound=0, upper_bound=1000
    )
    isopentenyl_phosphate_kinase.add_metabolites(
        {isopentenyl_p: -1, atp: -1, ipp: 1, adp: 1}
    )
    isopentenyl_phosphate_kinase.gene_reaction_rule = "AtIPK"
    model_iupd.add_reactions(
        [
            isoprenol_exchange,
            isoprenol_transport,
            isoprenol_kinase,
            isopentenyl_phosphate_kinase,
        ]
    )

    # Growth constraint: use the paired non-IUPD half-maximum rate.
    biomass = model_iupd.reactions.get_by_id("BIOMASS_Ec_iML1515_core_75p37M")
    biomass.lower_bound = fixed_growth_escherichia
    biomass.upper_bound = fixed_growth_escherichia

    # Production objective and pFBA.
    model_iupd.objective = "SQS"  # Introduced squalene synthase.
    # pFBA maximizes squalene first, then minimizes total absolute flux.
    model_iupd.objective_direction = "max"
    solution_iupd = pfba(model_iupd, fraction_of_optimum=1.0)

    # Register this analysis for export.
    models["Escherichia_coli_IUPD"] = (model_iupd, solution_iupd, carbon_sources)

    # ================================================================
    # 5. Yarrowia_lipolytica_non_IUPD
    # ================================================================
    # Load the model.
    model_file = "data/yarrowia/iYali5.mat"
    model_non_iupd = read_model(ROOT / model_file, SOLVER, "model")

    # Carbon sources (label introduction order).
    carbon_sources = ["D-glucose[e]"]

    # Culture conditions and model modifications.
    glucose_uptake = 1.0  # mmol/gDW/h; uptake is closed in the source iYali5.
    model_non_iupd.reactions.get_by_id("1714").bounds = (
        -glucose_uptake,
        1000,
    )  # D-glucose exchange.
    glucose_non_iupd = model_non_iupd.reactions.get_by_id("1714").lower_bound
    close_other_carbon_uptake(model_non_iupd, "1714")

    # Growth constraint: calculate non-IUPD maximum, then fix 50%.
    model_non_iupd.objective = "xBIOMASS"  # Biomass growth.
    growth_non_iupd = model_non_iupd.optimize(objective_sense="maximize")
    check_solution(model_non_iupd, growth_non_iupd)
    max_growth_yarrowia = growth_non_iupd.objective_value
    fixed_growth_yarrowia = 0.5 * max_growth_yarrowia
    biomass = model_non_iupd.reactions.get_by_id("xBIOMASS")
    biomass.lower_bound = fixed_growth_yarrowia
    biomass.upper_bound = fixed_growth_yarrowia
    model_non_iupd.reactions.get_by_id("1714").lower_bound = 2 * glucose_non_iupd

    # Production objective and pFBA.
    model_non_iupd.objective = "1012"  # Squalene synthase.
    # pFBA maximizes squalene first, then minimizes total absolute flux.
    model_non_iupd.objective_direction = "max"
    solution_non_iupd = pfba(model_non_iupd, fraction_of_optimum=1.0)

    # Register this analysis for export.
    models["Yarrowia_lipolytica_non_IUPD"] = (
        model_non_iupd,
        solution_non_iupd,
        carbon_sources,
    )

    # ================================================================
    # 6. Yarrowia_lipolytica_IUPD
    # ================================================================
    # Load the model.
    model_file = "data/yarrowia/iYali5.mat"
    # Load independently of the non-IUPD production constraints.
    model_iupd = read_model(ROOT / model_file, SOLVER, "model")

    # Carbon sources (label introduction order).
    carbon_sources = ["D-glucose[e]", "isoprenol[e]"]

    # Culture conditions and model modifications.
    glucose_uptake = 1.0  # mmol/gDW/h; uptake is closed in the source iYali5.
    model_iupd.reactions.get_by_id("1714").bounds = (
        -glucose_uptake,
        1000,
    )  # D-glucose exchange.
    glucose_iupd = model_iupd.reactions.get_by_id("1714").lower_bound
    close_other_carbon_uptake(model_iupd, "1714")
    # YALI0F30481g: ERG13, HMG-CoA synthase (reactions 559 and 560).
    model_iupd.genes.get_by_id("YALI0F30481g").knock_out()
    # Add the IU-pathway metabolites.
    isoprenol_e = Metabolite(
        "isoprenol[e]", name="isoprenol", formula="C5H10O", compartment="e"
    )
    isoprenol_c = Metabolite(
        "isoprenol[c]", name="isoprenol", formula="C5H10O", compartment="c"
    )
    isopentenyl_p = Metabolite(
        "isopentenyl_phosphate[c]",
        name="isopentenyl phosphate",
        formula="C5H11O4P",
        compartment="c",
    )
    model_iupd.add_metabolites([isoprenol_e, isoprenol_c, isopentenyl_p])
    isoprenol_exchange = Reaction(
        "Ex_isopentenol_e", lower_bound=2.09 * glucose_iupd, upper_bound=1000
    )
    isoprenol_exchange.add_metabolites({isoprenol_e: -1})  # Isoprenol exchange.
    isoprenol_transport = Reaction(
        "Ex_isopentenol_e_t", lower_bound=0, upper_bound=1000
    )
    # Isoprenol transport: extracellular -> cytosol.
    isoprenol_transport.add_metabolites({isoprenol_e: -1, isoprenol_c: 1})
    # m141: ATP; m143: ADP; m491: IPP, all cytosolic.
    atp = model_iupd.metabolites.get_by_id("m141")
    adp = model_iupd.metabolites.get_by_id("m143")
    ipp = model_iupd.metabolites.get_by_id("m491")
    # Isoprenol + ATP -> isopentenyl phosphate + ADP.
    isoprenol_kinase = Reaction("isopentenolK", lower_bound=0, upper_bound=1000)
    isoprenol_kinase.add_metabolites(
        {isoprenol_c: -1, atp: -1, isopentenyl_p: 1, adp: 1}
    )
    isoprenol_kinase.gene_reaction_rule = "YLR133W"
    # Isopentenyl phosphate + ATP -> IPP + ADP.
    isopentenyl_phosphate_kinase = Reaction(
        "PisopentenolK", lower_bound=0, upper_bound=1000
    )
    isopentenyl_phosphate_kinase.add_metabolites(
        {isopentenyl_p: -1, atp: -1, ipp: 1, adp: 1}
    )
    isopentenyl_phosphate_kinase.gene_reaction_rule = "AtIPK"
    model_iupd.add_reactions(
        [
            isoprenol_exchange,
            isoprenol_transport,
            isoprenol_kinase,
            isopentenyl_phosphate_kinase,
        ]
    )

    # Growth constraint: use the paired non-IUPD half-maximum rate.
    biomass = model_iupd.reactions.get_by_id("xBIOMASS")
    biomass.lower_bound = fixed_growth_yarrowia
    biomass.upper_bound = fixed_growth_yarrowia

    # Production objective and pFBA.
    model_iupd.objective = "1012"  # Squalene synthase.
    # pFBA maximizes squalene first, then minimizes total absolute flux.
    model_iupd.objective_direction = "max"
    solution_iupd = pfba(model_iupd, fraction_of_optimum=1.0)

    # Register this analysis for export.
    models["Yarrowia_lipolytica_IUPD"] = (model_iupd, solution_iupd, carbon_sources)
    return models


# Automatic environment selection and export; model edits belong above.
def main(argv=None):
    arguments = sys.argv[1:] if argv is None else list(argv)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", default="data/analysis_plan.json",
                        help="Generated plan path inside this project")
    parser.add_argument("--output-dir", help="Project-relative output folder")
    args = parser.parse_args(arguments)
    python = ROOT / PREPARATION_PYTHON
    if python.is_file() and python.resolve() != Path(sys.executable).resolve():
        os.execv(
            str(python), [str(python), str(Path(__file__).resolve())] + arguments
        )
    if importlib.util.find_spec("cobra") is None:
        raise SystemExit(
            "COBRApy is unavailable. Run bash scripts/setup_environments.sh first."
        )
    from carbon_flow_tools.model_preparation import export_model, save_plan, invalidate_plan

    # Invalidate the old plan before computing so a partial run cannot reuse it.
    target = invalidate_plan(ROOT, args.plan)
    output_dir = args.output_dir or (
        OUTPUT_DIR if target == (ROOT / "data/analysis_plan.json").resolve()
        else target.parent.relative_to(ROOT.resolve()).as_posix()
    )
    models = prepare_models()
    plan = {}
    for name, (model, solution, carbon_sources) in models.items():
        plan[name] = export_model(
            model,
            solution,
            name,
            carbon_sources,
            ROOT,
            output_dir,
            NONZERO_TOLERANCE,
        )
    save_plan(ROOT, plan, output_dir, target, script_file=__file__)
    relative_plan = target.relative_to(ROOT.resolve()).as_posix()
    print("Next: python3 step2_atom_mapping.py --plan {}\n"
          "Then: python3 step3_carbon_flow.py --plan {}"
          .format(relative_plan, relative_plan))


if __name__ == "__main__":
    main()
