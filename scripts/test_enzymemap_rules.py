"""Exercise the official rule fallback without using any database match."""

from pathlib import Path

from carbon_flow_tools.config import configure
from carbon_flow_tools.consensus import mapping_signatures
from carbon_flow_tools.pipeline import read_input
from carbon_flow_tools.mappers.enzymemap_adapter import map_batch
from carbon_flow_tools.mappers.enzymemap_rules import RuleMapper
from carbon_flow_tools.cli import build_mapper_configs
from step2_atom_mapping import settings


def main():
    root = Path(__file__).resolve().parents[1]
    configs = build_mapper_configs(settings())
    configure(1, 800)
    reactions, _, _ = read_input(
        root / "examples/Example_1_input.xlsx", "Nonzero", "Metabolites"
    )
    reaction = reactions[1]
    model = {"matches": {}, "predictor": RuleMapper(configs["EnzymeMap"]["options"])}
    result = map_batch(model, [reaction.reaction_smiles])[0]
    assert result.get("mapped_rxn") and not result.get("error"), result
    assert result["template"].startswith("EnzymeMap rule prediction")
    reaction.mapped_smiles = result["mapped_rxn"]
    mapping_signatures(reaction)
    print(
        "PASS: native EnzymeMap single-step rule, input restoration and complete carbon correspondence"
    )


if __name__ == "__main__":
    main()
