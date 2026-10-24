"""Check structure priority, reviewed name matches, caching and input protection."""


import csv
import hashlib
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch
from types import SimpleNamespace
from urllib.error import URLError

import pandas as pd
from openpyxl import load_workbook

from carbon_flow_tools.pipeline import read_input, run_analyses, run_carbon_flow
from carbon_flow_tools.reactions import get_smiles
from carbon_flow_tools.smiles_completion import (
    complete_smiles,
    query_key,
    read_model_smiles,
    select_candidate,
)
from scripts.test_parallel_carbon import make_input


def record(cid, smiles="CCO", formula="C2H6O", charge=0):
    return dict(CID=cid, SMILES=smiles, MolecularFormula=formula, Charge=charge)


def main():
    assert get_smiles(SimpleNamespace(SMILES_Model="CCO", SMILES_Yeast9="CC")) == "CCO"
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        manual = root / "manual.csv"
        with manual.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["Species", "MetaboliteName", "SMILES"])
            writer.writerow(["Checked", "isoprenol", "CC(=C)CCO"])
        options = dict(
            manual_file=manual,
            cache_file=root / "cache.json",
            lookup=True,
            timeout=10,
            min_interval=0,
        )
        yaml_file = root / "yeast-GEM.yml"
        yaml_file.write_text(
            "- metabolites:\n  - !!omap\n    - id: s_0\n"
            "      - smiles:\n        - CC\n"
            "  - !!omap\n    - id: s_1\n"
            '      - smiles:\n        - "CCO"\n        - C\n'
            "  - !!omap\n    - id: s_2\n"
            "      - smiles:\n        - CCO\n"
            "- reactions:\n    - id: reaction_only\n",
            encoding="utf-8",
        )
        assert read_model_smiles(yaml_file) == {
            "s_0": {"smiles": "CC"},
            "s_1": {"smiles": "CCO"},
            "s_2": {"smiles": "CCO"},
        }
        yaml_table = pd.DataFrame(
            [
                ["s_0", "isoprenol[c]", ""],
                ["s_1", "ethanol[c]", ""],
                ["s_2", "known[c]", "CC"],
            ],
            columns=["MetaboliteID", "MetaboliteName", "SMILES_Model"],
        )
        yaml_before = yaml_table.copy()
        yaml_options = dict(
            options, lookup=False, model_yaml_files={"Checked": yaml_file}
        )
        with patch(
            "carbon_flow_tools.smiles_completion.fetch_properties",
            side_effect=AssertionError("YAML lookup must work offline"),
        ):
            loaded = complete_smiles(yaml_table, "Checked", yaml_options)
            unrelated = complete_smiles(yaml_table, "AnotherSpecies", yaml_options)
        pd.testing.assert_frame_equal(yaml_table, yaml_before)
        assert [get_smiles(row) for row in loaded.itertuples(index=False)] == [
            "CC(=C)CCO",
            "CCO",
            "CC",
        ]
        assert loaded.iloc[1].StructureSource == "Model YAML"
        assert loaded.iloc[1].StructureReference == "yeast-GEM.yml"
        assert unrelated.SMILES_Model.eq(yaml_before.SMILES_Model).all()
        invalid_yaml = root / "invalid.yml"
        invalid_yaml.write_text("metabolites: []\n", encoding="utf-8")
        try:
            read_model_smiles(invalid_yaml)
        except ValueError as error:
            assert "model metabolite/reaction" in str(error)
        else:
            raise AssertionError("Unsupported YAML layout accepted")

        raven_yaml = root / "companion.yml"
        raven_yaml.write_text(
            '- metabolites:\n    - !!omap\n      - id: "m1"\n'
            '      - name: "ethanol"\n      - compartment: "c"\n'
            '      - formula: "C2H6O"\n      - smiles: "CCO"\n'
            '- reactions:\n    - !!omap\n      - id: "r1"\n',
            encoding="utf-8",
        )
        native = pd.DataFrame(
            [["m1", "ethanol[c]", "C2H6O"]],
            columns=["MetaboliteID", "MetaboliteName", "MolecularFormula"],
        )
        native["SMILES_Model"] = float("nan")  # read_input creates this empty numeric column.
        native_options = dict(options, lookup=False, model_yaml_files={"Checked": raven_yaml})
        with patch("carbon_flow_tools.smiles_completion.fetch_properties",
                   side_effect=AssertionError("Native annotations must work offline")):
            loaded = complete_smiles(native, "Checked", native_options)
            assert loaded.iloc[0].SMILES_Model == "CCO"
            assert loaded.iloc[0].StructureReference == "companion.yml"
            for field, value in [("MetaboliteName", "methanol[c]"),
                                 ("MetaboliteName", "ethanol[e]"),
                                 ("MolecularFormula", "CH4O")]:
                wrong = native.copy()
                wrong[field] = value
                rejected = complete_smiles(wrong, "Checked", native_options)
                assert not get_smiles(next(rejected.itertuples(index=False)))
            override = native.copy()
            override["SMILES_Model"] = "CC"
            preserved = complete_smiles(override, "Checked", native_options)
            assert preserved.iloc[0].SMILES_Model == "CC"
        duplicate = root / "duplicate.yml"
        duplicate.write_text('- metabolites:\n    - id: m1\n    - id: m1\n- reactions:\n')
        try:
            read_model_smiles(duplicate)
        except ValueError as error:
            assert "duplicate" in str(error)
        else:
            raise AssertionError("Duplicate annotation ID accepted")

        table = pd.DataFrame(
            [
                ["s_0", "isoprenol[c]", "CC", ""],
                ["s_1", "ethanol[c]", "", ""],
                ["s_2", "ethanol[e]", "", ""],
                ["s_3", "known[c]", "CC", ""],
            ],
            columns=[
                "MetaboliteID",
                "MetaboliteName",
                "SMILES_Yeast9",
                "SMILES_Manual",
            ],
        )
        before = table.copy()
        with patch(
            "carbon_flow_tools.smiles_completion.fetch_properties",
            return_value=[record(7)],
        ) as fetch:
            completed = complete_smiles(table, "Checked", options)
            assert fetch.call_count == 1  # Compartments share one chemical lookup.
        pd.testing.assert_frame_equal(table, before)
        assert get_smiles(next(completed.itertuples(index=False))) == "CC(=C)CCO"
        assert completed.loc[0, "StructureSource"] == "manual table"
        assert completed.loc[1, "StructureReviewRequired"]
        assert completed.loc[2, "StructureReference"] == "PubChem:7"
        assert completed.loc[3, "SMILES_Yeast9"] == "CC"
        with patch(
            "carbon_flow_tools.smiles_completion.fetch_properties",
            side_effect=AssertionError("Cache bypassed"),
        ):
            again = complete_smiles(table, "OtherModel", options)
        assert (
            get_smiles(next(again.itertuples(index=False))) == "CC"
        )  # Species-scoped names.

        exact = table.iloc[[1]].copy()
        exact["PubChemCID"] = 7
        with patch(
            "carbon_flow_tools.smiles_completion.fetch_properties",
            return_value=[record(7)],
        ):
            matched = complete_smiles(exact, "OtherModel", options)
        assert not matched.iloc[0].StructureReviewRequired
        assert matched.iloc[0].StructureSource == "PubChem cid"
        annotated = exact.copy()
        annotated["MolecularFormula"] = "C2H6O"
        annotated["Charge"] = 0
        row = next(annotated.itertuples(index=False))
        assert [
            item["CID"]
            for item in select_candidate(
                [
                    record(8),
                    record(7),
                    record(6, formula="C3H8O"),
                    record(9, charge=1),
                    record(10, "*"),
                ],
                row,
            )
        ] == [7, 8]
        charged = SimpleNamespace(MolecularFormula="C3H5O6P", Charge=-2)
        assert select_candidate(
            [record(4643300, "OCC(=O)COP(=O)([O-])[O-]", "C3H5O6P-2", -2)], charged
        )
        assert not select_candidate(
            [record(4643300, "OCC(=O)COP(=O)([O-])[O-]", "C3H5O6P-2", -1)], charged
        )
        key = exact.copy()
        key["PubChemCID"] = ""
        key["InChIKey"] = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        assert query_key(next(key.itertuples(index=False)))[0] == "inchikey"

        failed = table.iloc[[1]].copy()
        failed.MetaboliteName = "not_a_compound[c]"
        with patch(
            "carbon_flow_tools.smiles_completion.fetch_properties",
            side_effect=URLError("offline"),
        ):
            unavailable = complete_smiles(failed, "OtherModel", options)
        assert "Lookup failed" in unavailable.iloc[0].StructureStatus
        assert not get_smiles(next(unavailable.itertuples(index=False)))
        with patch(
            "carbon_flow_tools.smiles_completion.fetch_properties",
            return_value=[record(17)],
        ) as fetch:
            complete_smiles(failed, "OtherModel", options)
            assert fetch.call_count == 1  # A transient outage is not cached as missing.

        source = root / "Test_1_input.xlsx"
        make_input(source, "B[c]")
        workbook = load_workbook(source)
        workbook["Metabolites"].delete_cols(3, 2)  # No SMILES columns at all.
        workbook.save(source)
        workbook.close()
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        jobs = {
            "Test": dict(
                input_file=source,
                reaction_sheet="Nonzero",
                metabolite_sheet="Metabolites",
            )
        }
        configs = {
            name: dict(
                module="scripts.test_parallel_carbon",
                python=sys.executable,
                options={"test": True},
            )
            for name in ("FirstMapper", "SecondMapper")
        }
        with patch(
            "carbon_flow_tools.smiles_completion.fetch_properties",
            return_value=[record(7)],
        ):
            run_analyses(
                jobs,
                1,
                800,
                root,
                list(configs),
                configs,
                stage="map",
                max_mapper_processes=1,
                structure_options=options,
            )
        mapping = root / "Test_2_atom_mapping.xlsx"
        selected = pd.read_excel(
            mapping, sheet_name="Mapping_Selection", engine="openpyxl"
        )
        assert selected.ReviewRequired.eq(True).all()  # Even two agreeing tools.
        assert selected.SelectionReason.str.contains("structure requires review").all()
        assert (
            pd.read_excel(mapping, sheet_name="Final_Mapping", engine="openpyxl")
            .ReviewRequired.eq(True)
            .all()
        )
        _, _, enriched = read_input(mapping, "Nonzero", "Metabolites")
        assert enriched.SMILES_Retrieved.eq("CCO").all()
        assert hashlib.sha256(source.read_bytes()).hexdigest() == digest
        mapping_digest = hashlib.sha256(mapping.read_bytes()).hexdigest()
        run_carbon_flow(
            {"Test": dict(mapping_file=mapping, carbon_sources=["A[e]"])},
            root,
            max_processes=1,
        )
        assert hashlib.sha256(mapping.read_bytes()).hexdigest() == mapping_digest
        assert "EMU_Metabolites" in pd.ExcelFile(root / "Test_3_carbon_flow.xlsx", engine="openpyxl").sheet_names
        with patch(
            "carbon_flow_tools.pipeline.run_mapper",
            side_effect=AssertionError("Reuse failed"),
        ), patch(
            "carbon_flow_tools.smiles_completion.fetch_properties",
            side_effect=AssertionError("Cache failed"),
        ):
            run_analyses(
                jobs,
                1,
                800,
                root,
                list(configs),
                configs,
                stage="map",
                max_mapper_processes=1,
                structure_options=options,
            )
    print(
        "PASS: manual priority, species/name scoping, PubChem ID/name review, cache/outages, full stages and input protection"
    )


if __name__ == "__main__":
    main()
