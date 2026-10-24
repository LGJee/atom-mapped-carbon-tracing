"""Validate generated plans and an isolated map-to-EMU handoff."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import pandas as pd
from carbon_flow_tools.cli import load_analysis_plan, main
from scripts.test_parallel_carbon import make_input
import step2_atom_mapping as mapping
import step3_carbon_flow as tracing


class AnalysisPlanTests(unittest.TestCase):
    def test_plan_hash_validation_case_isolation_and_environment_handoff(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.xlsx"
            make_input(source, "B[c]")
            from openpyxl import load_workbook
            wb = load_workbook(source)
            wb["Nonzero"].append(["in", "Import", " => A[e]", 1])
            wb["Nonzero"].append(["out", "Drain", "B[c] => ", 1])
            wb.save(source)
            wb.close()
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            (root / "missing.csv").write_text("Species,MetaboliteName,SMILES\n")
            case = root / "results/results_by_AI/case"
            plan = root / "analysis_plan.json"
            content = dict(version=3, output_dir=str(case), models={"Case": dict(input_file=str(source),
                           input_sha256=digest, carbon_sources=["A[e]"])})
            plan.write_text(json.dumps(content))
            ms = mapping.settings()
            ms.update(analysis_plan=str(plan), mappers=["LocalMapper"], max_mapper_processes=1,
                      structure_options=dict(lookup=False, manual_file=str(root / "missing.csv"),
                                             cache_file=str(root / "cache.json"), timeout=1, min_interval=0))
            ts = tracing.settings()
            ts.update(analysis_plan=str(plan), tracing_python=sys.executable, max_carbon_processes=1)
            configs = {"LocalMapper": dict(python=sys.executable, module="scripts.test_parallel_carbon", options={"test": True})}
            with patch("carbon_flow_tools.cli.build_mapper_configs", return_value=configs):
                self.assertEqual(main("map", ms, []), 0)
            with patch("carbon_flow_tools.cli.build_mapper_configs", side_effect=AssertionError("Mapper requested by step 3")):
                self.assertEqual(main("trace", ts, []), 0)
            out = case / "Case_3_carbon_flow.xlsx"
            table = pd.read_excel(out, sheet_name="EMU_Metabolites", engine="openpyxl").set_index("MetaboliteName")
            self.assertAlmostEqual(table.loc["B[c]", "A[e]"], 1)
            report = json.loads((case / "carbon_flow_run.json").read_text())
            self.assertTrue(report["models"][0]["emu_completed"])
            content["models"]["Case"]["input_sha256"] = "0" * 64
            plan.write_text(json.dumps(content))
            with self.assertRaisesRegex(ValueError, "rerun step2"):
                main("trace", ts, [])
            self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), digest)
            content["version"] = 2
            plan.write_text(json.dumps(content))
            with self.assertRaisesRegex(ValueError, "rerun step1"):
                load_analysis_plan(plan, "trace", ts)
            plan.unlink()
            with self.assertRaises(FileNotFoundError):
                load_analysis_plan(plan, "trace", ts)


def run():
    unittest.main(module=__name__)


if __name__ == "__main__":
    run()
