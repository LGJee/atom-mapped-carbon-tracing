"""Regression checks for preserved scientific behavior."""

import sys
from unittest.mock import patch

import pandas as pd

from carbon_flow_tools import config
from carbon_flow_tools.carbon_trace import prepare_mappings, map_one_reaction
from carbon_flow_tools.mapper_runner import (
    run_mapper,
    validate_mapper_config,
    mapping_error,
)
from carbon_flow_tools.reactions import (
    PreparedReaction,
    build_side,
    parse_side,
    prepare_reactions,
)
from carbon_flow_tools.atom_helpers import isoprenol_labels


def load_model(options=None):
    if options is not None:
        assert options["probe"] == 7 and len(options["reactions"]) == 4
    return {
        "C>>C": "[CH4:101]>>[CH4:101]",
        "CC>>CC": "[CH3:40][CH3:80]>>[CH3:80][CH3:40]",
    }


def map_batch(model, reaction_smiles):
    if len(reaction_smiles) > 1:
        raise RuntimeError("test batch failure")
    if reaction_smiles == ["O>>O"]:
        return [
            {
                "mapped_rxn": "[OH2:1]>>[OH2:1]",
                "error": "test: prediction cannot restore input",
            }
        ]
    return [{"mapped_rxn": model[reaction_smiles[0]], "confidence": 0.8}]


def main():
    assert not mapping_error("[CH4:1]>>[CH4:1]")
    assert "Element mismatch" in mapping_error("[CH4:1]>>[OH2:1]")
    assert "Repeated" in mapping_error("[CH4:1].[CH4:1]>>[CH4:1]")
    assert mapping_error("C>C")
    settings = {
        "module": "scripts.test_mapper_backends",
        "python": sys.executable,
        "options": {"probe": 7},
    }
    reactions = [
        PreparedReaction(str(i), "", 1, "", smiles, [], [])
        for i, smiles in enumerate(("C>>C", "CC>>CC", "N>>N", "O>>O"))
    ]
    skipped = PreparedReaction("skip", "", 1, "", "", [], [], status="Missing SMILES")
    reactions.append(skipped)
    config.configure(3, 800)
    run_mapper("Test", reactions, "ExampleMapper", settings)
    assert [r.status for r in reactions[:2]] == ["mapped", "mapped"]
    assert "ExampleMapper failed" in reactions[2].status
    assert skipped.status == "Missing SMILES"
    assert reactions[1].mapped_smiles == load_model()["CC>>CC"]
    assert reactions[3].mapped_smiles == "[OH2:1]>>[OH2:1]"
    assert reactions[3].status.startswith("Invalid mapping: test:")
    try:
        validate_mapper_config("../unsafe", settings)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid output name was accepted")
    try:
        validate_mapper_config("Example", dict(settings, python="/no/such/python"))
    except FileNotFoundError:
        pass
    else:
        raise AssertionError("missing environment was accepted")

    structures = {
        "isoprenol[e]": "CC(=C)CCO",
        "isoprenol[c]": "OCCC(=C)C",
        "A[e]": "CCO",
        "A[c]": "OCC",
        "B[e]": "CCO",
        "B[c]": "CCO",
        "X[e]": "[NH4+]",
        "X[c]": "N",
        "S[e]": "C[C@H](O)F",
        "S[c]": "C[C@@H](O)F",
        "ATP[c]": "CC",
        "ADP[c]": "C",
        "phosphate[c]": "OP(=O)(O)O",
    }
    metabolites = pd.DataFrame(
        [
            dict(MetaboliteName=name, SMILES_Yeast9=smiles, SMILES_Manual="")
            for name, smiles in structures.items()
        ]
    )
    equations = [
        ("transport", "isoprenol[e] => isoprenol[c]", 1),
        ("reverse", "isoprenol[e] => isoprenol[c]", -1),
        ("copies", "2 A[e] => 2 A[c]", 1),
        ("antiport", "A[e] + B[c] => A[c] + B[e]", 1),
        ("chemical", "A[c] => B[c]", 1),
        ("charge", "X[e] => X[c]", 1),
        ("stereo", "S[e] => S[c]", 1),
        ("coupled", "ATP[c] + A[e] => ADP[c] + phosphate[c] + A[c]", 1),
        ("exchange", "isoprenol[e] <=>", -1),
    ]
    table = pd.DataFrame(
        [
            dict(ReactionID=name, ReactionName=name, Equation=equation, Flux=flux)
            for name, equation, flux in equations
        ]
    )
    prepared = prepare_reactions(table, metabolites)
    assert [r.is_transport for r in prepared] == [True] * 4 + [False] * 5
    assert all(not mapping_error(r.mapped_smiles) for r in prepared[:4])
    assert (
        prepared[-1].status == "Boundary reaction; no atom mapping required"
        and not prepared[-1].mapped_smiles
    )
    assert all(r.status == "ready" for r in prepared[4:8])
    with patch("carbon_flow_tools.mapper_runner.subprocess.run") as worker:
        run_mapper("Test", prepared[:4] + prepared[-1:], "ExampleMapper", settings)
        worker.assert_not_called()
    labels = tuple(
        (0, atom, label)
        for atom, label in isoprenol_labels(structures["isoprenol[e]"]).items()
    )
    mapping = prepare_mappings(prepared[:1])[0][1]
    _, rows, _ = map_one_reaction(
        prepared[0], mapping, {"isoprenol[e]": {1: labels}}, 1, {}, 2, structures
    )
    assert rows[0]["Production"] == "isoprenol[c]" and rows[0]["Step"] == 1
    assert rows[0]["LabeledSMILES"] == "O[CH2:7][CH2:8][C:9](=[CH2:10])[CH3:11]"
    config.configure(1, 1)
    assert prepare_reactions(table.iloc[:1], metabolites)[0].is_transport
    config.configure(3, 800)

    structures = {"D-glucose[c]": "CCO", "D-glucose 6-phosphate[c]": "CCO"}
    left, error = build_side(parse_side("D-glucose[c]"), structures)
    assert not error
    right, error = build_side(parse_side("D-glucose 6-phosphate[c]"), structures)
    assert not error
    raw = "[CH3:40][CH2:80][OH:90]>>[CH3:80][CH2:40][OH:90]"
    reaction = PreparedReaction(
        "r_0534",
        "",
        1,
        "D-glucose[c] -> D-glucose 6-phosphate[c]",
        "",
        left,
        right,
        source_text="D-glucose[c]",
        status="mapped",
        mapped_smiles=raw,
    )
    prepared = prepare_mappings([reaction])
    assert len(prepared) == 1 and reaction.mapped_smiles == raw
    _, rows, _ = map_one_reaction(
        reaction,
        prepared[0][1],
        {"D-glucose[c]": {1: ((0, 0, "C1"),)}},
        1,
        {},
        2,
        structures,
    )
    assert rows[0]["LabeledSMILES"] == "C[CH2:1]O"
    assert rows[0]["CarbonLabels"] == "C1"
    print(
        "PASS: new adapter, batch fallback, isolated failure, raw mapping and source labels"
    )


if __name__ == "__main__":
    main()
