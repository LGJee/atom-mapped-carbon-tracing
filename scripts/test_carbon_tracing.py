"""Regression checks for preserved scientific behavior."""

from carbon_flow_tools.carbon_trace import (
    map_one_reaction,
    initialize_sources,
    propagate_once,
    merge_states,
    inorganic_carbon_names,
    propagate_rounds,
)
from carbon_flow_tools.atom_helpers import mapped_components
from carbon_flow_tools.reactions import PreparedReaction, build_side, parse_side


def read_carbon_round(path, number=1, **kwargs):
    """Compare one exported round with the unchanged propagation records."""
    import pandas as pd
    kwargs.setdefault("engine", "openpyxl")
    table = pd.read_excel(path, sheet_name="Carbon_Flow", **kwargs)
    return table.loc[table.Carbon_Flow_Round.eq(number)].drop(columns="Carbon_Flow_Round").reset_index(drop=True)


def reaction(left, right, smiles, mapped):
    l, error = build_side(parse_side(left), smiles)
    assert not error
    r, error = build_side(parse_side(right), smiles)
    assert not error
    rxn = PreparedReaction(
        "test",
        "",
        1,
        left + " -> " + right,
        "",
        l,
        r,
        source_text=left,
        mapped_smiles=mapped,
    )
    return rxn, mapped_components(rxn)


def check_smiles_reuse():
    from carbon_flow_tools.atom_helpers import (
        _labeled_smiles,
        labeled_smiles,
        assignment_from_labeled_smiles,
    )

    _labeled_smiles.cache_clear()
    assignment_from_labeled_smiles.cache_clear()
    cases = [
        ("CCO", ((0, 0, "C1"),)),
        ("CCO", ((0, 1, "C1"),)),
        ("OCC", ((0, 0, "C1"),)),
        ("CCC", ((0, 0, "C1"), (0, 2, "C2"))),
        ("CCC", ((0, 0, "C2"), (0, 2, "C1"))),
        ("CC", ((0, 0, "C7"), (0, 1, "C7"))),
        ("C.CO", ((0, 0, "C1"), (1, 0, "C2"))),
    ]
    for base, assignment in cases:
        expected = _labeled_smiles.__wrapped__(base, assignment)
        assert labeled_smiles(base, list(assignment)) == expected
        assert labeled_smiles(base, iter(assignment)) == expected
        assert assignment_from_labeled_smiles(base, expected) == (
            assignment_from_labeled_smiles.__wrapped__(base, expected)
        )
    assert labeled_smiles(*cases[0]) != labeled_smiles(*cases[1])
    assert labeled_smiles(*cases[3]) == labeled_smiles(*cases[4])
    assert [
        label
        for _, _, label in assignment_from_labeled_smiles(
            "CC", labeled_smiles(*cases[5])
        )
    ] == ["C7", "C7"]
    for number in range(1, 8194):
        marked = labeled_smiles("C", ((0, 0, "C" + str(number)),))
        assignment_from_labeled_smiles("C", marked)
    assert _labeled_smiles.cache_info().currsize == 8192
    assert assignment_from_labeled_smiles.cache_info().currsize == 8192
    assert _labeled_smiles.cache_info().hits > 0
    _labeled_smiles.cache_clear()
    assignment_from_labeled_smiles.cache_clear()


def check_inorganic_carbon_termination():
    import pandas as pd

    structures = {
        "formate[c]": "O=CO",
        "co2[c]": "O=C=O",
        "co2[m]": "[O:7]=[C:8]=[O:9]",
        "bicarbonate[m]": "O=C(O)[O-]",
        "carbonic_acid[e]": "OC(O)=O",
        "carbonate[e]": "[O-]C([O-])=O",
        "pyruvate[m]": "CC(=O)C(=O)O",
        "oxaloacetate[m]": "O=C(O)CC(=O)C(=O)O",
        "oxaloacetate[c]": "O=C(O)CC(=O)C(=O)O",
        "water[m]": "O",
        "hydroxide[m]": "[OH-]",
    }
    assert inorganic_carbon_names(structures) == {
        "co2[c]",
        "co2[m]",
        "bicarbonate[m]",
        "carbonic_acid[e]",
        "carbonate[e]",
    }
    # Use unnumbered input structures; only the reaction mappings carry map IDs.
    structures["co2[m]"] = "O=C=O"
    pyr = "[CH3:1][C:2](=[O:3])[C:4](=[O:5])[OH:6]"
    oaa = "[O:7]=[C:8]([OH:9])[CH2:1][C:2](=[O:3])[C:4](=[O:5])[OH:6]"
    specs = [
        ("formate[c]", "co2[c]", "[CH:1](=[O:2])[OH:3]>>[O:2]=[C:1]=[O:3]"),
        ("co2[c]", "co2[m]", "[O:1]=[C:2]=[O:3]>>[O:1]=[C:2]=[O:3]"),
        (
            "co2[m] + water[m]",
            "bicarbonate[m]",
            "[O:1]=[C:2]=[O:3].[OH2:4]>>[O:1]=[C:2]([OH:3])[O-:4]",
        ),
        (
            "formate[c] + water[m]",
            "bicarbonate[m]",
            "[CH:1](=[O:2])[OH:3].[OH2:4]>>[O:2]=[C:1]([OH:3])[O-:4]",
        ),
        (
            "pyruvate[m] + bicarbonate[m]",
            "oxaloacetate[m] + hydroxide[m]",
            pyr + ".[O:7]=[C:8]([OH:9])[O-:10]>>" + oaa + ".[OH-:10]",
        ),
        ("pyruvate[m] + co2[m]", "oxaloacetate[m]", pyr + ".[O:7]=[C:8]=[O:9]>>" + oaa),
        ("oxaloacetate[m]", "oxaloacetate[c]", oaa + ">>" + oaa),
    ]
    reactions = []
    for index, (left, right, mapped) in enumerate(specs):
        item, components = reaction(left, right, structures, mapped)
        assert components is not None, index
        item.reaction_id = "r" + str(index)
        item.status = "mapped"
        reactions.append(item)
    metabolites = pd.DataFrame(
        [
            dict(MetaboliteID=name, MetaboliteName=name, SMILES_Model=smiles)
            for name, smiles in structures.items()
        ]
    )
    sources = [
        ("formate[c]", ((0, 1, "C11"),), {"C11"}),
        ("pyruvate[m]", ((0, 0, "C1"), (0, 1, "C2"), (0, 3, "C3")), {"C1", "C2", "C3"}),
    ]
    for ordered_sources in (sources, list(reversed(sources))):
        for enabled in (False, True):
            results = list(
                propagate_rounds(
                    "Test",
                    reactions,
                    metabolites,
                    3,
                    [s[0] for s in ordered_sources],
                    sources=ordered_sources,
                    allow_inorganic_carbon_recycling=enabled,
                )
            )
            assert len(results) >= 2
            for _, table, _, _ in results:
                assert (
                    table.loc[table.Production == "co2[c]", "CarbonLabels"]
                    .eq("C11")
                    .all()
                )
                assert "co2[c]" in set(table.Production)
                assert "bicarbonate[m]" in set(table.Production)
                assert ("co2[m]" in set(table.Production)) == enabled
                fixed = table.loc[table.Production.str.startswith("oxaloacetate[")]
                assert {"oxaloacetate[m]", "oxaloacetate[c]"} <= set(fixed.Production)
                if enabled:
                    assert fixed.CarbonLabels.eq("C1,C2,C3,C11").any()
                    assert set(fixed.CarbonLabels) <= {
                        "C11",
                        "C1,C2,C3",
                        "C1,C2,C3,C11",
                    }
                else:
                    assert fixed.CarbonLabels.eq("C1,C2,C3").all()
                    assert not fixed.SourceIDs.str.contains("&").any()
                    # Both carboxylases still transfer the pyruvate labels.
                    assert {"r4", "r5"} <= set(fixed.ReactionID)
                    assert (
                        table.loc[table.ReactionID.isin(["r4", "r5"]), "SourceIDs"]
                        .str.endswith(" + 0")
                        .all()
                    )
    try:
        list(propagate_rounds("Test", reactions, metabolites, 2, ["co2[c]"]))
    except ValueError as error:
        assert "ALLOW_INORGANIC_CARBON_RECYCLING = True" in str(error)
    else:
        raise AssertionError("An inorganic-carbon source was silently blocked")
    enabled = list(
        propagate_rounds(
            "Test",
            reactions,
            metabolites,
            2,
            ["co2[c]"],
            allow_inorganic_carbon_recycling=True,
        )
    )
    assert "oxaloacetate[c]" in set(enabled[0][1].Production)
    try:
        list(
            propagate_rounds(
                "Test",
                reactions,
                metabolites,
                2,
                ["formate[c]"],
                allow_inorganic_carbon_recycling="False",
            )
        )
    except ValueError as error:
        assert "must be True or False" in str(error)
    else:
        raise AssertionError("A string was accepted as a boolean setting")


def main():
    check_smiles_reuse()
    check_inorganic_carbon_termination()
    structures = {"A[c]": "CC", "water[c]": "O", "B[c]": "CC"}
    rxn, parts = reaction(
        "A[c] + water[c]",
        "B[c] + water[c]",
        structures,
        "[OH2:3].[CH3:1][CH3:2]>>[OH2:3].[CH3:1][CH3:2]",
    )
    assert parts is not None
    assert parts[0][0][1] == {0: 1, 1: 2}
    assert parts[0][1][1] == {0: 3}
    _, rows, _ = map_one_reaction(
        rxn, parts, {"A[c]": {1: ((0, 0, "C1"),)}}, 1, {}, 2, structures
    )
    assert len(rows) == 1 and rows[0]["Production"] == "B[c]"
    assert rows[0]["CarbonLabels"] == "C1"
    import pandas as pd
    from carbon_flow_tools.pipeline import carbon_metabolites, summarize_sources

    sources = [("glucose", (), {"C1"}), ("isoprenol", (), {"C7"})]
    ids = {"A[c]": {"a"}, "B[m]": {"b"}, "C[c]": {"c"}}
    flow = pd.DataFrame(
        [["A[c]", "C1"], ["A[c]", "C1"], ["B[m]", "C7"], ["C[c]", "C1,C7"]],
        columns=["Production", "CarbonLabels"],
    )
    previous = {}
    summary = summarize_sources(flow, sources, ids, 1, previous)
    assert summary["All sources_LabeledCarbonMetabolites"] == 3
    assert summary["glucose_LabeledCarbonMetabolites"] == 2
    assert summary["isoprenol_LabeledCarbonMetabolites"] == 2
    assert summary["MixedSources_LabeledCarbonMetabolites"] == 1
    summary = summarize_sources(flow.iloc[:3], sources, ids, 2, previous)
    assert summary["MixedSources_LabeledCarbonMetabolites"] == 0
    assert summary["glucose_LabeledCarbonMetabolitesIncrement"] == -1
    assert summary["isoprenol_LabeledCarbonMetabolitesIncrement"] == -1
    single = summarize_sources(flow, sources[:1], ids, 1, {})
    assert single["glucose_LabeledCarbonMetabolites"] == 2
    assert not any(key.startswith("All sources_") for key in single)
    assert len(pd.DataFrame([single, single])) == 2
    assert "MixedSources_LabeledCarbonMetabolites" not in single
    assert single == {
        "Round": 1, "TotalCarbonMetabolites": 3,
        "glucose_LabeledCarbonMetabolites": 2,
        "glucose_LabeledCarbonMetabolitesIncrement": 2,
        "glucose_LabeledCarbonPercent": 66.67,
        "glucose_ExclusiveSourceMetabolites": 1,
    }
    mixed_flow = pd.DataFrame(
        [
            ["A[c]", "C1"], ["A[c]", "C7"],
            ["B[c]", "C1,C1"], ["B[c]", "C1,C99"],
            ["C[c]", "C1,C7"], ["C[c]", "C1,C7"],
            ["C[c]", "C1,C8"], ["C[m]", "C7,C12"],
            ["unmapped[c]", "C1,C7"],
        ], columns=["Production", "CarbonLabels"],
    )
    mixed_sources = [
        ("glucose", (), {"C1"}), ("xylose", (), {"C7", "C8"}),
        ("acetate", (), {"C12"}),
    ]
    mixed_carbon_ids = {"A[c]": {"a"}, "B[c]": {"b"}, "C[c]": {"c"}, "C[m]": {"cm"}}
    summary = summarize_sources(mixed_flow, mixed_sources, mixed_carbon_ids, 1, {})
    assert summary["TotalCarbonMetabolites"] == 4
    assert summary["All sources_LabeledCarbonMetabolites"] == 4
    assert summary["MixedSources_LabeledCarbonMetabolites"] == 2
    assert summarize_sources(
        mixed_flow.iloc[:4], mixed_sources, mixed_carbon_ids, 2, {}
    )["MixedSources_LabeledCarbonMetabolites"] == 0
    assert summarize_sources(
        mixed_flow.iloc[:0], mixed_sources, mixed_carbon_ids, 3, {}
    )["MixedSources_LabeledCarbonMetabolites"] == 0

    structures = {
        name: "C" for name in ("A[c]", "A[m]", "B[c]", "D[c]", "Z[c]", "E[e]", "U[c]")
    }
    structures["water[c]"] = "O"
    mapped = [
        reaction(
            "A[c] + water[c]",
            "B[c] + water[c]",
            structures,
            "[CH4:1].[OH2:2]>>[CH4:1].[OH2:2]",
        )[0],
        reaction("2 A[c]", "2 A[m]", structures, "[CH4:1].[CH4:2]>>[CH4:1].[CH4:2]")[0],
        reaction("D[c]", "D[c]", structures, "[CH4:1]>>[CH4:1]")[0],
    ]
    for item in mapped:
        item.status = "mapped"
    mapped[1].is_transport = True
    failed = reaction("Z[c]", "B[c]", structures, "[CH4:1]>>[CH4:1]")[0]
    failed.status = "mapping failed"
    boundary = reaction("E[e]", "", structures, "")[0]
    boundary.status = "Boundary reaction; no atom mapping required"
    empty = reaction("U[c]", "B[c]", structures, "")[0]
    empty.status = "mapped"
    metabolites = pd.DataFrame(
        [
            dict(
                MetaboliteID=name,
                MetaboliteName=name,
                SMILES_Yeast9=smiles,
                SMILES_Manual="",
            )
            for name, smiles in structures.items()
        ]
        + [
            dict(
                MetaboliteID="missing",
                MetaboliteName="missing[c]",
                SMILES_Yeast9="",
                SMILES_Manual="",
            )
        ]
    )
    metabolites = pd.concat([metabolites, metabolites.iloc[:1]], ignore_index=True)
    ids, unknown = carbon_metabolites(metabolites, mapped + [failed, boundary, empty])
    assert (
        set(ids) == {"A[c]", "A[m]", "B[c]", "D[c]"}
        and sum(map(len, ids.values())) == 4
    )
    assert unknown == 1
    flow = pd.DataFrame(
        [["B[c]", "C1"], ["B[c]", "C1"]], columns=["Production", "CarbonLabels"]
    )
    summary = summarize_sources(flow, [("glucose", (), {"C1"})], ids, 1, {})
    assert summary["TotalCarbonMetabolites"] == 4
    assert summary["glucose_LabeledCarbonMetabolites"] == 1
    assert summary["glucose_LabeledCarbonPercent"] == 25
    assert (
        summarize_sources(flow, [("glucose", (), {"C1"})], {}, 1, {})[
            "glucose_LabeledCarbonPercent"
        ]
        is None
    )
    structures = {"D-glucose[e]": "OC1OC(CO)C(O)C(O)C1O", "isoprenol[e]": "CC(=C)CCO"}
    sources = initialize_sources(["D-glucose[e]", "isoprenol[e]"], structures)
    assert sources[0][2] == {"C" + str(i) for i in range(1, 7)}
    assert sources[1][2] == {"C" + str(i) for i in range(7, 12)}
    reverse = initialize_sources(["isoprenol[e]", "D-glucose[e]"], structures)
    assert reverse[0][2] == {"C" + str(i) for i in range(1, 6)}
    assert reverse[1][2] == {"C" + str(i) for i in range(6, 12)}

    from collections import defaultdict

    structures = {"A[c]": "C", "B[c]": "C", "C[c]": "C"}
    ab = reaction("A[c]", "B[c]", structures, "[CH4:1]>>[CH4:1]")
    bc = reaction("B[c]", "C[c]", structures, "[CH4:1]>>[CH4:1]")
    states = defaultdict(dict, {"A[c]": {0: ((0, 0, "C1"),)}})
    steps, catalog = {("A[c]", 0): 0}, {}
    produced, table, next_id = propagate_once(
        [ab, bc], states, steps, catalog, 1, structures
    )
    assert set(table.Production) == {"B[c]"} and "B[c]" not in states
    merge_states(produced, table, states, steps)
    _, table, _ = propagate_once([bc], states, steps, catalog, next_id, structures)
    assert set(table.Production) == {"C[c]"} and table.iloc[0].Step == 2
    smiles = {"A[c]": "CCC", "B[c]": "CCC"}
    rxn, parts = reaction(
        "A[c]", "B[c]", smiles, "[CH3:1][CH2:2][CH3:3]>>[CH3:1][CH2:2][CH3:3]"
    )
    states = {
        "A[c]": {36: ((0, 0, "C1"), (0, 1, "C6"), (0, 2, "C5")), 39: ((0, 0, "C1"),)}
    }
    _, rows, _ = map_one_reaction(rxn, parts, states, 8, {}, 1, smiles)
    assert len(rows) == 2
    assert {r["SourceIDs"]: r["CarbonLabels"] for r in rows} == {
        "36": "C1,C5,C6",
        "39": "C1",
    }
    smiles = {"A[c]": "C", "B[c]": "C", "P[c]": "CC"}
    rxn, parts = reaction(
        "A[c] + B[c]", "P[c]", smiles, "[CH4:1].[CH4:2]>>[CH3:1][CH3:2]"
    )
    _, rows, _ = map_one_reaction(
        rxn,
        parts,
        {"A[c]": {1: ((0, 0, "C7"),)}, "B[c]": {2: ((0, 0, "C7"),)}},
        2,
        {},
        3,
        smiles,
    )
    assert rows[0]["SourceIDs"] == "1 + 2"
    assert rows[0]["CarbonLabels"] == "C7,C7"
    rxn, parts = reaction("P[c]", "2 A[c]", smiles, "[CH3:1][CH3:2]>>[CH4:1].[CH4:2]")
    _, rows, _ = map_one_reaction(
        rxn, parts, {"P[c]": {3: ((0, 0, "C1"), (0, 1, "C2"))}}, 3, {}, 4, smiles
    )
    assert {r["CarbonLabels"] for r in rows} == {"C1", "C2"}
    assert len(rows) == 2

    smiles = {"A[c]": "C", "B[c]": "C", "C[c]": "C", "P[c]": "CCCC"}
    rxn, parts = reaction(
        "A[c] + 2 B[c] + C[c]",
        "P[c]",
        smiles,
        "[CH4:1].[CH4:2].[CH4:3].[CH4:4]>>[CH3:1][CH2:2][CH2:3][CH3:4]",
    )
    _, rows, _ = map_one_reaction(
        rxn,
        parts,
        {
            "A[c]": {1: ((0, 0, "C1"),)},
            "B[c]": {2: ((0, 0, "C2"),), 3: ((0, 0, "C3"),)},
            "C[c]": {4: ((0, 0, "C4"),), 5: ((0, 0, "C5"),)},
        },
        2,
        {},
        6,
        smiles,
        unordered_sources=True,
    )
    assert {r["SourceIDs"] for r in rows} == {
        "1 + 2&2 + 4",
        "1 + 2&3 + 4",
        "1 + 3&3 + 4",
        "1 + 2&2 + 5",
        "1 + 2&3 + 5",
        "1 + 3&3 + 5",
    }
    assert (
        len({r["LabeledSMILES"] for r in rows if r["SourceIDs"] == "1 + 2&3 + 4"}) == 2
    )
    print("PASS: separate states, repeated labels, separate product copies")


if __name__ == "__main__":
    main()
