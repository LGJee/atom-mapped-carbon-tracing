"""Check memory safety stops without exhausting real machine memory."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from openpyxl import Workbook

from carbon_flow_tools.config import MemoryGuard, MemoryPressureError, available_memory
from carbon_flow_tools.carbon_trace import map_one_reaction, propagate_rounds
from carbon_flow_tools.pipeline import trace_analysis
from scripts.test_carbon_tracing import reaction
import step3_carbon_flow as step3


class MemoryGuardTests(unittest.TestCase):
    def test_validation_and_reserve(self):
        for value in (0, -1, True, None, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                MemoryGuard(value)
        with patch("carbon_flow_tools.config.available_memory", return_value=(3, 40)):
            with self.assertRaises(MemoryPressureError):
                MemoryGuard(2).check(force=True)
        with patch("carbon_flow_tools.config.available_memory", return_value=(5, 40)):
            guard = MemoryGuard(2)
            guard.check(force=True)
            with self.assertRaises(MemoryPressureError):
                guard.check(force=True, extra_bytes=2 * 1024 ** 3)

    def test_throttling_and_unreadable_memory(self):
        with patch("carbon_flow_tools.config.available_memory", return_value=(8, 16)) as read:
            with patch("carbon_flow_tools.config.time.monotonic", side_effect=[0, .1, .6]):
                guard = MemoryGuard()
                guard.check()
                guard.check()
                guard.check()
            self.assertEqual(read.call_count, 2)
        with patch("carbon_flow_tools.config.available_memory", side_effect=OSError("unreadable")):
            with self.assertRaisesRegex(MemoryPressureError, "Cannot read"):
                MemoryGuard().check()

    def test_linux_and_nested_cgroup_limits(self):
        files = {
            "/proc/meminfo": "MemTotal: 33554432 kB\nMemAvailable: 16777216 kB\n",
            "/proc/self/cgroup": "0::/parent/child\n",
            "/sys/fs/cgroup/parent/child/memory.max": "max",
            "/sys/fs/cgroup/parent/memory.max": str(8 * 1024 ** 3),
            "/sys/fs/cgroup/parent/memory.current": str(7 * 1024 ** 3),
        }
        def read(path):
            if str(path) not in files:
                raise FileNotFoundError(str(path))
            return files[str(path)]
        with patch.object(Path, "read_text", read):
            self.assertEqual(available_memory(), (1, 8))
            del files["/sys/fs/cgroup/parent/memory.max"]
            self.assertEqual(available_memory(), (16, 32))

    def test_inner_combination_loop_stops(self):
        structures = {"A[c]": "C", "P[c]": "CC"}
        rxn, parts = reaction(
            "2 A[c]", "P[c]", structures,
            "[CH4:1].[CH4:2]>>[CH3:1][CH3:2]",
        )
        states = {"A[c]": {i: ((0, 0, "C" + str(i)),) for i in range(1, 18)}}
        with patch("carbon_flow_tools.config.available_memory",
                   side_effect=[(8, 16), (8, 16), (.5, 16)]) as read:
            with patch("carbon_flow_tools.config.time.monotonic", side_effect=[0, 1, 2]):
                with self.assertRaises(MemoryPressureError):
                    map_one_reaction(rxn, parts, states, 1, {}, 1, structures,
                                     memory_guard=MemoryGuard())
        self.assertEqual(read.call_count, 3)

    def test_guard_preserves_results_and_empty_round(self):
        structures = {"A[c]": "C", "B[c]": "C", "C[c]": "C"}
        reactions = []
        for i, (left, right) in enumerate((("A[c]", "B[c]"), ("B[c]", "C[c]"))):
            item, _ = reaction(left, right, structures, "[CH4:1]>>[CH4:1]")
            item.reaction_id = "r" + str(i)
            item.status = "mapped"
            reactions.append(item)
        metabolites = pd.DataFrame([
            dict(MetaboliteID=name, MetaboliteName=name, SMILES_Model=value)
            for name, value in structures.items()
        ])
        with patch("carbon_flow_tools.config.available_memory", return_value=(12, 16)):
            expected = list(propagate_rounds("Test", reactions, metabolites, 3, ["A[c]"]))
            actual = list(propagate_rounds("Test", reactions, metabolites, 3,
                                         ["A[c]"], memory_guard=MemoryGuard()))
            self.assertEqual(len(expected), len(actual))
            for a, b in zip(expected, actual):
                self.assertEqual((a[0], a[2], a[3]), (b[0], b[2], b[3]))
                pd.testing.assert_frame_equal(a[1], b[1])
            empty = list(propagate_rounds("Test", [], metabolites, 2,
                                         ["A[c]"], memory_guard=MemoryGuard()))
            self.assertTrue(empty[0][1].empty)

    def test_saved_round_survives_stop_and_existing_file_is_preserved(self):
        structures = {"A[c]": "C", "B[c]": "C"}
        rxn, _ = reaction("A[c]", "B[c]", structures, "[CH4:1]>>[CH4:1]")
        rxn.status = "mapped"
        metabolites = pd.DataFrame([
            dict(MetaboliteID=name, MetaboliteName=name, SMILES_Model=value)
            for name, value in structures.items()
        ])
        original = propagate_rounds
        def stop_after_first(*args, **kwargs):
            yield next(original(*args, **kwargs))
            raise MemoryPressureError("simulated pressure in round 2")
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "Test_3_carbon_flow.xlsx"
            with patch("carbon_flow_tools.config.available_memory", return_value=(12, 16)):
                with patch("carbon_flow_tools.pipeline.propagate_rounds", stop_after_first):
                    with self.assertRaises(MemoryPressureError):
                        trace_analysis("Test", {"carbon_sources": ["A[c]"]}, [rxn],
                                       metabolites, 3, Workbook(), output)
            with pd.ExcelFile(output, engine="openpyxl") as excel:
                self.assertEqual(set(excel.sheet_names), {"Summary", "Carbon_Flow_R1"})
                self.assertEqual(pd.read_excel(excel, "Summary").Round.tolist(), [1])
            previous = output.read_bytes()
            with patch("carbon_flow_tools.config.available_memory", return_value=(.5, 16)):
                with self.assertRaises(MemoryPressureError):
                    trace_analysis("Test", {"carbon_sources": ["A[c]"]}, [rxn],
                                   metabolites, 3, Workbook(), output)
            self.assertEqual(output.read_bytes(), previous)

    def test_excel_writing_and_step3_settings(self):
        from carbon_flow_tools.excel_io import write_sheet
        calls = []
        class StopWriting:
            def check(self, context="", force=False):
                calls.append(context)
                if len(calls) == 3:
                    raise MemoryPressureError("simulated pressure during Excel write")
        workbook = Workbook()
        with self.assertRaises(MemoryPressureError):
            write_sheet(workbook, "Carbon_Flow_R2",
                        pd.DataFrame({"value": range(3000)}), StopWriting())
        self.assertEqual(workbook["Carbon_Flow_R2"].max_row, 1025)
        workbook.close()
        with patch.object(step3, "MIN_AVAILABLE_MEMORY_GB", 3):
            self.assertEqual(step3.settings()["min_available_memory_gb"], 3)


if __name__ == "__main__":
    unittest.main()
