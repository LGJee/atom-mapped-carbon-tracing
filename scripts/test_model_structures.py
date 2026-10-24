"""Check automatic source discovery without network access or mapper calls."""

import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.io import savemat

from carbon_flow_tools.model_structures import discover_annotations
from carbon_flow_tools.smiles_completion import complete_smiles


def source(path, variable=None):
    return dict(
        path=str(path),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        variable_name=variable,
    )


def yaml_model(path, smiles="CCO", name="ethanol", compartment="c", formula="C2H6O"):
    path.write_text(
        '- metabolites:\n  - !!omap\n    - id: "m1"\n    - name: "{}"\n'
        '    - compartment: "{}"\n    - formula: "{}"\n    - smiles: "{}"\n'
        "- reactions:\n".format(name, compartment, formula, smiles),
        encoding="utf-8",
    )


class DiscoveryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.folder = self.root / "species"
        self.folder.mkdir()
        self.path = self.folder / "network.json"
        self.path.write_text(
            json.dumps(
                dict(
                    metabolites=[
                        dict(id="m1", name="ethanol", compartment="c", formula="C2H6O")
                    ]
                )
            ),
            encoding="utf-8",
        )
        self.table = pd.DataFrame(
            [
                dict(
                    MetaboliteID="m1",
                    MetaboliteName="ethanol[c]",
                    MolecularFormula="C2H6O",
                )
            ]
        )
        self.manual = self.root / "manual.csv"
        self.manual.write_text(
            "Species,MetaboliteName,SMILES\n", encoding="utf-8"
        )

    def tearDown(self):
        self.temporary.cleanup()

    def options(self, lookup=False):
        return dict(
            manual_file=self.manual,
            cache_file=self.root / "cache.json",
            timeout=1,
            min_interval=0,
            lookup=lookup,
            source_models={"Example": source(self.path)},
        )

    def test_different_stem_yaml_metadata_json_and_directory_isolation(self):
        yaml_model(self.folder / "companion-different-name.yml")
        (self.folder / "source.json").write_text('{"source":"metadata only"}')
        yaml_model(self.root / "unrelated.yml", smiles="CO", name="methanol")
        before = self.table.copy()
        with patch(
            "carbon_flow_tools.smiles_completion.fetch_properties",
            side_effect=AssertionError("Offline discovery requested network access"),
        ):
            loaded = complete_smiles(self.table, "Example", self.options())
        self.assertEqual(loaded.iloc[0].SMILES_Model, "CCO")
        self.assertEqual(
            loaded.iloc[0].StructureReference, "companion-different-name.yml"
        )
        self.assertFalse(loaded.iloc[0].StructureReviewRequired)
        pd.testing.assert_frame_equal(self.table, before)

    def test_companion_identity_must_match(self):
        path = self.folder / "companion.yml"
        for field in (
            dict(name="methanol"),
            dict(compartment="e"),
            dict(formula="CH4O"),
        ):
            yaml_model(path, **field)
            found, _ = discover_annotations(source(self.path), self.table)
            self.assertEqual(found, {})

    def test_native_json_priority_conflict_review_and_manual_override(self):
        self.path.write_text(
            json.dumps(
                dict(
                    metabolites=[
                        dict(
                            id="m1",
                            name="ethanol",
                            compartment="c",
                            formula="C2H6O",
                            annotation=dict(smiles=["CCO"]),
                        )
                    ]
                )
            )
        )
        yaml_model(self.folder / "companion.yml", smiles="COC")
        loaded = complete_smiles(self.table, "Example", self.options())
        self.assertEqual(loaded.iloc[0].SMILES_Model, "CCO")
        self.assertTrue(loaded.iloc[0].StructureReviewRequired)
        self.assertIn("Conflicting", loaded.iloc[0].StructureStatus)
        self.assertEqual(len(json.loads(loaded.iloc[0].StructureCandidates)), 2)
        self.manual.write_text(
            "Species,MetaboliteName,SMILES\nExample,ethanol,OCC\n"
        )
        loaded = complete_smiles(self.table, "Example", self.options())
        self.assertEqual(loaded.iloc[0].SMILES_Manual, "OCC")
        self.assertFalse(loaded.iloc[0].StructureReviewRequired)

    def test_equivalent_smiles_and_existing_input_priority(self):
        yaml_model(self.folder / "one.yml", smiles="CCO")
        yaml_model(self.folder / "two.yml", smiles="OCC")
        table = self.table.copy()
        table["SMILES_Model"] = "OCC"
        loaded = complete_smiles(table, "Example", self.options())
        self.assertEqual(loaded.iloc[0].SMILES_Model, "OCC")
        self.assertFalse(loaded.iloc[0].StructureReviewRequired)

    def test_mat_native_smiles_and_compartments(self):
        path = self.folder / "network.mat"
        savemat(
            str(path),
            dict(
                model=dict(
                    mets=np.array(["m1"], dtype=object),
                    metNames=np.array(["ethanol"], dtype=object),
                    metFormulas=np.array(["C2H6O"], dtype=object),
                    metSmiles=np.array(["CCO"], dtype=object),
                    metComps=np.array([1]),
                    comps=np.array(["c"], dtype=object),
                    S=np.array([[1]]),
                )
            ),
        )
        found, _ = discover_annotations(source(path, "model"), self.table)
        self.assertEqual(found["m1"]["smiles"], "CCO")
        self.assertEqual(found["m1"]["reference"], "network.mat")

    def test_missing_structures_use_source_identifiers_for_existing_lookup(self):
        self.path.write_text(
            json.dumps(
                dict(
                    metabolites=[
                        dict(
                            id="m1",
                            name="ethanol",
                            compartment="c",
                            annotation=dict(inchi_key=["LFQSCWFLJHTTHZ-UHFFFAOYSA-N"]),
                        )
                    ]
                )
            )
        )
        with patch(
            "carbon_flow_tools.smiles_completion.fetch_properties",
            return_value=[
                dict(CID=702, SMILES="CCO", MolecularFormula="C2H6O", Charge=0)
            ],
        ) as fetch:
            loaded = complete_smiles(self.table, "Example", self.options(lookup=True))
        self.assertEqual(
            fetch.call_args[0][:2], ("inchikey", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
        )
        self.assertEqual(loaded.iloc[0].SMILES_Retrieved, "CCO")
        self.assertFalse(loaded.iloc[0].StructureReviewRequired)

    def test_source_hash_prevents_stale_model_association(self):
        recorded = source(self.path)
        self.path.write_text('{"metabolites": []}')
        with self.assertRaisesRegex(ValueError, "Source model changed"):
            discover_annotations(recorded, self.table)

    def test_ambiguous_source_identifiers_fall_back_to_reviewed_name_lookup(self):
        for keys in (["first", "second"], ["first"]):
            self.path.write_text(
                json.dumps(
                    dict(
                        metabolites=[
                            dict(
                                id="m1",
                                name="ethanol",
                                compartment="c",
                                annotation=dict(inchi_key=keys),
                            )
                        ]
                    )
                )
            )
            if len(keys) == 1:
                (self.folder / "companion.json").write_text(
                    json.dumps(
                        dict(
                            metabolites=[
                                dict(
                                    id="m1",
                                    name="ethanol",
                                    compartment="c",
                                    annotation=dict(inchi_key="second"),
                                )
                            ]
                        )
                    )
                )
            with patch(
                "carbon_flow_tools.smiles_completion.fetch_properties",
                return_value=[
                    dict(CID=702, SMILES="CCO", MolecularFormula="C2H6O", Charge=0)
                ],
            ) as fetch:
                loaded = complete_smiles(
                    self.table, "Example", self.options(lookup=True)
                )
            if fetch.called:
                self.assertEqual(fetch.call_args[0][0], "name")
            self.assertTrue(loaded.iloc[0].StructureReviewRequired)

    def test_structure_conflict_marks_final_mapping_for_review(self):
        from carbon_flow_tools.pipeline import run_analyses
        from scripts.test_parallel_carbon import make_input

        workbook = self.root / "Example_1_input.xlsx"
        make_input(workbook, "B[c]")
        original = workbook.read_bytes()
        self.path.write_text(
            json.dumps(
                dict(
                    metabolites=[
                        dict(
                            id="s_1",
                            name="A",
                            compartment="c",
                            annotation=dict(smiles="CCO"),
                        )
                    ]
                )
            )
        )
        companion = self.folder / "companion.yml"
        yaml_model(companion, smiles="COC", name="A")
        companion.write_text(companion.read_text().replace('id: "m1"', 'id: "s_1"'))
        run_analyses(
            {
                "Example": dict(
                    input_file=workbook,
                    reaction_sheet="Nonzero",
                    metabolite_sheet="Metabolites",
                )
            },
            1,
            800,
            self.root,
            ["LocalMapper", "RXNMapper"],
            {
                "LocalMapper": dict(
                    python=sys.executable,
                    module="scripts.test_parallel_carbon",
                    options={"test": True},
                ),
                "RXNMapper": dict(
                    python=sys.executable,
                    module="scripts.test_parallel_carbon",
                    options={"test": True},
                ),
            },
            stage="map",
            max_mapper_processes=1,
            structure_options=self.options(),
        )
        result = self.root / "Example_2_atom_mapping.xlsx"
        final = pd.read_excel(result, sheet_name="Final_Mapping", engine="openpyxl")
        self.assertTrue(final.ReviewRequired.all())
        selection = pd.read_excel(
            result, sheet_name="Mapping_Selection", engine="openpyxl"
        )
        self.assertTrue(
            selection.SelectionReason.str.contains(
                "Automatically selected structure"
            ).all()
        )
        self.assertEqual(workbook.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
