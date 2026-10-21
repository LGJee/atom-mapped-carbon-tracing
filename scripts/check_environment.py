"""Check installed versions and run a real pFBA or mapping smoke test."""

import argparse
import importlib
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "backend", choices=("preparation", "localmapper", "rxnmapper", "enzymemap")
    )
    args = parser.parse_args()

    try:
        from importlib.metadata import version
    except ImportError:
        from pkg_resources import get_distribution

        version = lambda package: get_distribution(package).version

    if args.backend == "preparation":
        from cobra import Model, Metabolite, Reaction
        from cobra.flux_analysis import pfba

        model = Model("preparation_check")
        model.solver = "glpk"
        substrate = Metabolite("A_c", compartment="c")
        supply = Reaction("supply", lower_bound=0, upper_bound=10)
        supply.add_metabolites({substrate: 1})
        growth = Reaction("growth", lower_bound=0, upper_bound=1000)
        growth.add_metabolites({substrate: -1})
        model.add_reactions([supply, growth])
        model.objective = growth
        result = pfba(model)
        if result.status != "optimal" or abs(result.fluxes["growth"] - 10) > 1e-8:
            raise RuntimeError("COBRApy/GLPK pFBA smoke test failed.")
        for package in ("cobra", "swiglpk", "numpy", "scipy", "pandas", "openpyxl"):
            print("{}={}".format(package, version(package)))
        print("preparation: real COBRApy/GLPK pFBA test passed")
        return
    from rdkit import Chem

    if args.backend == "localmapper":
        for package in ("pandas", "openpyxl", "scipy.sparse.linalg"):
            importlib.import_module(package)
        from openpyxl import LXML

        if not LXML:
            raise RuntimeError(
                "Native Excel backend unavailable: install lxml==4.9.4 in this "
                "environment and ensure OPENPYXL_LXML is not disabled."
            )
        from carbon_flow_tools.mappers import localmapper_adapter as adapter

        result = adapter.map_batch(adapter.load_model(), ["CC=O>>CCO"])[0]
        versions = (
            "localmapper",
            "torch",
            "dgl",
            "dgllife",
            "rdkit-pypi",
            "pandas",
            "openpyxl",
            "lxml",
            "scipy",
        )
    elif args.backend == "rxnmapper":
        from carbon_flow_tools.mappers import rxnmapper_adapter as adapter

        result = adapter.map_batch(adapter.load_model(), ["CC=O>>CCO"])[0]
        versions = ("rxnmapper", "torch", "transformers", "rdkit")
    else:
        helpers_map = importlib.import_module("enzymemap.helpers_map")
        helpers_rdkit = importlib.import_module("enzymemap.helpers_rdkit")
        if not callable(helpers_map.map) or not callable(
            helpers_rdkit.neutralize_atoms
        ):
            raise RuntimeError("Official EnzymeMap interfaces are unavailable.")
        from rdchiral.main import rdchiralRunText

        results = rdchiralRunText("[C:1]=[O:2]>>[C:1][O:2]", "CC=O")
        if not results:
            raise RuntimeError(
                "EnzymeMap rule dependency smoke test produced no product."
            )
        if any(Chem.MolFromSmiles(product) is None for product in results):
            raise RuntimeError("Rule dependency produced an invalid product.")
        result = None
        versions = ("enzymemap", "pandas", "wrapt-timeout-decorator")
    if result is not None:
        sides = result.get("mapped_rxn", "").split(">>")
        if len(sides) != 2 or any(Chem.MolFromSmiles(side) is None for side in sides):
            raise RuntimeError("Smoke test returned invalid mapped reaction SMILES.")
    for package in versions:
        print("{}={}".format(package, version(package)))
    print(
        "{}: {} passed".format(
            args.backend,
            "real rule dependency test" if result is None else "real mapper test",
        )
    )


if __name__ == "__main__":
    main()
