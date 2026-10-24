"""EMU export remains atomic under file-size, memory and disk limits."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from openpyxl import Workbook
from carbon_flow_tools.excel_io import save_output, workbook_limit, WorkbookSizeLimitError
from carbon_flow_tools.config import MemoryPressureError, DiskSpaceLimitError
from carbon_flow_tools.pipeline import run_carbon_flow
from scripts.test_parallel_carbon import prepare_fixture
import step3_carbon_flow as step3


class OutputLimitTests(unittest.TestCase):
    def test_limit_and_shared_settings(self):
        self.assertEqual(workbook_limit(50), 50000000)
        self.assertEqual(step3.settings()["max_excel_size_mb"], 50)
        self.assertNotIn("rounds", step3.settings())
        for value in (None, True, 0, -.1, float("inf"), float("nan"), .0000001):
            with self.assertRaises(ValueError):
                workbook_limit(value)

    def test_rejected_save_preserves_exact_previous_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.xlsx"
            wb = Workbook()
            save_output(wb, path)
            before = path.read_bytes()
            with self.assertRaises(WorkbookSizeLimitError):
                save_output(wb, path, max_bytes=1)
            self.assertEqual(path.read_bytes(), before)
            self.assertFalse(path.with_suffix(".tmp.xlsx").exists())
            wb.close()

    def test_resource_failure_keeps_history_and_reports_no_new_completion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mapping = prepare_fixture(root)
            output = root / "Test_3_carbon_flow.xlsx"
            jobs = {"Test": dict(mapping_file=mapping, carbon_sources=["A[e]"], min_available_memory_gb=.01)}
            run_carbon_flow(jobs, root, max_processes=1)
            previous = output.read_bytes()
            failures = [("carbon_flow_tools.source_contributions.solve_source_contributions",
                         MemoryPressureError("test memory pressure"), "memory_limit"),
                        ("carbon_flow_tools.config.DiskGuard.check", DiskSpaceLimitError("test disk reserve"), "disk_space_limit")]
            for target, error, reason in failures:
                with patch(target, side_effect=error), self.assertRaises(type(error)):
                    run_carbon_flow(jobs, root, max_processes=1)
                status = json.loads((root / "carbon_flow_run.json").read_text())["models"][0]
                self.assertEqual(status["stop_reason"], reason)
                self.assertFalse(status["emu_completed"])
                self.assertIsNone(status["output_file"])
                self.assertEqual(output.read_bytes(), previous)
            jobs["Test"]["max_excel_size_mb"] = .001
            with self.assertRaises(WorkbookSizeLimitError):
                run_carbon_flow(jobs, root, max_processes=1)
            status = json.loads((root / "carbon_flow_run.json").read_text())["models"][0]
            self.assertEqual(status["stop_reason"], "file_size_limit")
            self.assertTrue(status["emu_solved"])
            self.assertFalse(status["emu_completed"])
            self.assertEqual(output.read_bytes(), previous)

    def test_duplicate_output_and_input_overwrite_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            p = root / "Test_3_carbon_flow.xlsx"
            p.touch()
            with self.assertRaisesRegex(ValueError, "overwrite an input"):
                run_carbon_flow({"Test": dict(mapping_file=p, carbon_sources=["A"])}, root)
            with self.assertRaisesRegex(ValueError, "same carbon-flow"):
                run_carbon_flow({"Test": dict(mapping_files=[p, p], carbon_sources=["A"])}, root)


if __name__ == "__main__":
    unittest.main()
