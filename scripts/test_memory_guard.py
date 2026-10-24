"""Check memory safety stops without exhausting real machine memory."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from openpyxl import Workbook

from carbon_flow_tools.config import MemoryGuard, MemoryPressureError, available_memory
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
            write_sheet(workbook, "Carbon_Flow",
                        pd.DataFrame({"value": range(3000)}), StopWriting())
        self.assertEqual(workbook["Carbon_Flow"].max_row, 1025)
        workbook.close()
        with patch.object(step3, "MIN_AVAILABLE_MEMORY_GB", 3):
            self.assertEqual(step3.settings()["min_available_memory_gb"], 3)


if __name__ == "__main__":
    unittest.main()
