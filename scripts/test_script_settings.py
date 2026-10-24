"""Entry settings contain only active stage-specific options."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from carbon_flow_tools.cli import discover_inputs, main
import step2_atom_mapping as mapping
import step3_carbon_flow as tracing


class ScriptSettingsTests(unittest.TestCase):
    def test_discovery_and_mapper_free_trace(self):
        settings = tracing.settings()
        self.assertNotIn("rounds", settings)
        self.assertNotIn("analysis_methods", settings)
        self.assertNotIn("batch_size", settings)
        self.assertTrue(settings["emu_allow_inorganic_carbon_recycling"])
        self.assertNotIn("rounds", mapping.settings())
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            normal = root / "Culture_2_atom_mapping.xlsx"
            ai = root / "Culture_2_atom_mapping_AIGenerated.xlsx"
            normal.touch()
            ai.touch()
            (root / "Ignore_2_atom_mapping_backup.xlsx").touch()
            self.assertEqual(set(discover_inputs(root, "trace")), {"Culture"})
            settings.pop("analysis_plan")
            settings.update(input_dir=str(root), output_dir=str(root), carbon_sources={"Culture":["A[e]"]},
                            tracing_python=sys.executable)
            with patch("carbon_flow_tools.pipeline.run_carbon_flow") as run:
                self.assertEqual(main("trace", settings, ["Culture"]), 0)
                job = run.call_args[0][0]["Culture"]
                self.assertEqual(job["mapping_files"], [normal, ai])
                self.assertEqual(job["carbon_sources"], ["A[e]"])
                self.assertEqual(run.call_args[0][1], root)
            with self.assertRaises(SystemExit):
                main("trace", settings, ["--rounds", "3"])
            settings["carbon_sources"] = {}
            with self.assertRaises(SystemExit):
                main("trace", settings, ["Culture"])


def run():
    unittest.main(module=__name__)


if __name__ == "__main__":
    run()
