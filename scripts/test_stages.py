"""Mapping output reuse and an independent, mapper-free EMU stage."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from carbon_flow_tools.pipeline import run_analyses, run_carbon_flow
from scripts.test_parallel_carbon import prepare_fixture


class StageTests(unittest.TestCase):
    def test_saved_mapping_reuse_and_independent_trace(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mapping = prepare_fixture(root)
            configs = {"ExampleMapper": dict(python=sys.executable, module="scripts.test_parallel_carbon", options={"test":True})}
            job = dict(input_file=root / "Test_1_input.xlsx", reaction_sheet="Nonzero", metabolite_sheet="Metabolites",
                       carbon_sources=["A[e]"])
            with patch("carbon_flow_tools.pipeline.run_mapper", side_effect=AssertionError("Prediction reuse failed")):
                run_analyses({"Test":job}, 1, 800, root, ["ExampleMapper"], configs, stage="map")
            saved = mapping.read_bytes()
            with patch("carbon_flow_tools.pipeline.run_mapper", side_effect=AssertionError("Mapper called by EMU")):
                statuses = run_carbon_flow({"Test":dict(mapping_file=mapping, carbon_sources=["A[e]"])},
                                          root / "output", max_processes=1)
            self.assertTrue(statuses[0]["emu_completed"])
            self.assertEqual(mapping.read_bytes(), saved)


def main():
    unittest.main(module=__name__)


if __name__ == "__main__":
    main()
