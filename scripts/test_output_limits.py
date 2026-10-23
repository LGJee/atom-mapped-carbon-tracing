"""Check batch-level size stops and preservation of complete per-model results."""

from scripts.test_carbon_tracing import read_carbon_round

from pathlib import Path
import json
import os
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from openpyxl import Workbook

from carbon_flow_tools.excel_io import (
    save_output, write_sheet, workbook_limit, WorkbookSizeLimitError,
)
from carbon_flow_tools.pipeline import trace_analysis, run_carbon_flow
from carbon_flow_tools.flow_output import FlowOutput
from carbon_flow_tools.config import MemoryPressureError
from scripts.test_carbon_tracing import reaction
import step3_carbon_flow as step3


def summary(counts):
    return pd.DataFrame({"Round": range(1, len(counts) + 1),
                         "TotalCarbonMetabolites": 100,
                         "glucose[e]_LabeledCarbonMetabolites": counts})


def output_workbook(path, table):
    workbook = Workbook()
    del workbook[workbook.active.title]
    write_sheet(workbook, "Summary", table)
    write_sheet(workbook, "Carbon_Flow", pd.DataFrame({"Carbon_Flow_Round": [1, 1], "ProductionID": [1, 2]}))
    save_output(workbook, path)
    workbook.close()


def flow_fixture():
    from carbon_flow_tools.carbon_trace import propagate_rounds

    structures = {"A[c]": "C", "B[c]": "C"}
    rxn, _ = reaction("A[c]", "B[c]", structures, "[CH4:1]>>[CH4:1]")
    rxn.status = "mapped"
    metabolites = pd.DataFrame([
        dict(MetaboliteID=k, MetaboliteName=k, SMILES_Model=v)
        for k, v in structures.items()
    ])
    first = next(propagate_rounds("Test", [rxn], metabolites, 1, ["A[c]"]))
    return rxn, metabolites, first[1]


class OutputLimitTests(unittest.TestCase):
    def test_limit_and_shared_settings(self):
        self.assertEqual(workbook_limit(50), 50000000)
        self.assertEqual(step3.settings()["max_excel_size_mb"], 50)
        self.assertEqual(step3.settings()["rounds"], 10)
        self.assertEqual(workbook_limit(.001), 1000)
        for value in (None, True, 0, -.1, float("inf"), float("nan"), .0000001):
            with self.assertRaises(ValueError):
                workbook_limit(value)
        with patch.object(step3, "MAX_EXCEL_SIZE_MB", 75):
            self.assertEqual(step3.settings()["max_excel_size_mb"], 75)

    def test_rejected_save_preserves_exact_previous_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "flow.xlsx"
            workbook = Workbook()
            save_output(workbook, path)
            previous = path.read_bytes()
            with self.assertRaises(WorkbookSizeLimitError):
                save_output(workbook, path, max_bytes=1)
            self.assertEqual(path.read_bytes(), previous)
            self.assertFalse(path.with_suffix(".tmp.xlsx").exists())
            save_output(workbook, path, max_bytes=len(previous))
            self.assertEqual(path.stat().st_size, len(previous))
            workbook.close()

    def test_checkpoint_size_discards_exact_threshold_overshoot_and_first_round(self):
        rxn, metabolites, table = flow_fixture()
        original_checkpoint = FlowOutput.checkpoint
        original_stat = Path.stat
        for sizes in ((49999999, 50000000), (49999999, 80000000), (60000000,)):
            for suffix in ("", "_AIGenerated"):
                with self.subTest(sizes=sizes, suffix=suffix), tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / ("Test_3_carbon_flow" + suffix + ".xlsx")
                    output_workbook(path, summary([90]))
                    historical = path.read_bytes()
                    status, current_round = {}, [0]

                    def checkpoint(store):
                        current_round[0] = store.number
                        return original_checkpoint(store)

                    def measured_stat(candidate, *args, **kwargs):
                        result = original_stat(candidate, *args, **kwargs)
                        if candidate.name == "candidate.xlsx":
                            values = list(result)
                            values[6] = sizes[min(current_round[0] - 1, len(sizes) - 1)]
                            return os.stat_result(values)
                        return result

                    with patch.object(FlowOutput, "checkpoint", checkpoint), patch.object(Path, "stat", measured_stat):
                        result = trace_analysis("Test", {"carbon_sources": ["A[c]"],
                            "output_batch_rows": 1, "_run_status": status},
                            [rxn], metabolites, 10, Workbook(), path)
                    completed = len(sizes) - 1
                    self.assertEqual(status["stop_reason"], "file_size_limit")
                    self.assertEqual(status["completed_rounds"], completed)
                    self.assertEqual(result, path if completed else None)
                    self.assertEqual(status["output_file"], str(path) if completed else None)
                    if completed:
                        with pd.ExcelFile(path, engine="openpyxl") as excel:
                            self.assertEqual(excel.sheet_names, ["Summary", "Carbon_Flow"])
                            self.assertEqual(pd.read_excel(excel, "Summary").Round.tolist(), [1])
                            pd.testing.assert_frame_equal(read_carbon_round(excel, 1, dtype={"SourceIDs": str}), table)
                    else:
                        self.assertEqual(path.read_bytes(), historical)
                    self.assertFalse(list(Path(directory).glob(".carbon-flow-*")))

    def test_below_threshold_preserves_record_order_and_convergence(self):
        from carbon_flow_tools.carbon_trace import propagate_rounds
        rxn, metabolites, _ = flow_fixture()
        expected = list(propagate_rounds("Test", [rxn], metabolites, 3, ["A[c]"]))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "Test_3_carbon_flow.xlsx"
            status = {}
            trace_analysis("Test", {"carbon_sources": ["A[c]"], "output_batch_rows": 1,
                "_run_status": status}, [rxn], metabolites, 3, Workbook(), path)
            self.assertEqual(status["stop_reason"], "converged")
            self.assertEqual(status["completed_rounds"], len(expected))
            with pd.ExcelFile(path, engine="openpyxl") as excel:
                self.assertEqual(pd.read_excel(excel, "Summary").Round.tolist(), [r[0] for r in expected])
                for number, table, _, _ in expected:
                    pd.testing.assert_frame_equal(read_carbon_round(excel, number, dtype={"SourceIDs": str}), table)

    def test_failed_checkpoint_preserves_previous_complete_round(self):
        rxn, metabolites, _ = flow_fixture()
        original_checkpoint, original_complete = FlowOutput.checkpoint, FlowOutput.complete_round
        saved = []
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "Test_3_carbon_flow.xlsx"
            status = {}

            def checkpoint(store):
                original_checkpoint(store)
                if store.number == 2:
                    raise OSError("simulated disk failure")

            def complete(store):
                original_complete(store)
                saved.append(path.read_bytes())

            with patch.object(FlowOutput, "checkpoint", checkpoint), patch.object(FlowOutput, "complete_round", complete):
                with self.assertRaisesRegex(OSError, "simulated disk failure"):
                    trace_analysis("Test", {"carbon_sources": ["A[c]"], "output_batch_rows": 1,
                        "_run_status": status}, [rxn], metabolites, 3, Workbook(), path)
            self.assertEqual(path.read_bytes(), saved[0])
            self.assertEqual(status["completed_rounds"], 1)
            self.assertFalse(list(Path(directory).glob(".carbon-flow-*")))

    def test_row_limit_stops_in_batch_without_paginating_or_starting_next_round(self):
        rxn, metabolites, table = flow_fixture()
        calls = []
        def rounds(*args, **kwargs):
            store = kwargs["row_store"]
            row = table.to_dict("records")[0]
            for number in (1, 2):
                store.start_round(number)
                store.start_layer()
                store.append(row)
                if number == 2:
                    calls.append("combined sheet overflow should stop the first row")
                    store.append(dict(row, SourceIDs="9"))
                    calls.append("must not be reached")
                store.flush()
                yield number, None, 1, False
            raise AssertionError("No subsequent round after overflow")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "Test_3_carbon_flow.xlsx"
            status = {}
            with patch("carbon_flow_tools.pipeline.propagate_rounds", rounds), patch.object(FlowOutput, "MAX_DATA_ROWS", 1):
                trace_analysis("Test", {"carbon_sources": ["A[c]"], "output_batch_rows": 1,
                    "_run_status": status}, [rxn], metabolites, 5, Workbook(), path)
            self.assertEqual(calls, [])
            self.assertEqual(status["stop_reason"], "row_limit")
            self.assertEqual(status["completed_rounds"], 1)
            with pd.ExcelFile(path, engine="openpyxl") as excel:
                self.assertEqual(excel.sheet_names, ["Summary", "Carbon_Flow"])
                self.assertEqual(pd.read_excel(excel, "Summary").Round.tolist(), [1])

    def test_first_batch_overflow_produces_no_current_output(self):
        rxn, metabolites, _ = flow_fixture()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "Test_3_carbon_flow.xlsx"
            status = {}
            result = trace_analysis("Test", {"carbon_sources": ["A[c]"], "max_excel_size_mb": .000001,
                "output_batch_rows": 1, "_run_status": status}, [rxn], metabolites, 5, Workbook(), path)
            self.assertIsNone(result)
            self.assertFalse(path.exists())
            self.assertEqual(status["completed_rounds"], 0)
            self.assertEqual(status["stop_reason"], "file_size_limit")

    def test_error_status_distinguishes_versions_of_the_same_model(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            normal = root / "A_2_atom_mapping.xlsx"
            ai = root / "A_2_atom_mapping_AIGenerated.xlsx"
            normal.write_bytes(b"normal")
            ai.write_bytes(b"reviewed")

            def stopped(strain, job, rounds, mapping_file, output):
                if mapping_file == ai:
                    raise MemoryPressureError("AI worker stopped")
                output_workbook(output, summary([20]))

            jobs = {"A": {"mapping_files": [normal, ai], "carbon_sources": ["A[c]"]}}
            with patch("carbon_flow_tools.pipeline.trace_mapping_file", stopped):
                with self.assertRaises(MemoryPressureError):
                    run_carbon_flow(jobs, 10, root, max_processes=1)
            report = json.loads((root / "carbon_flow_run.json").read_text())
            statuses = {entry["mapping_file"]: entry for entry in report["models"]}
            self.assertEqual(statuses[str(normal)]["completed_rounds"], 1)
            self.assertEqual(statuses[str(normal)]["stop_reason"], "batch_interrupted")
            self.assertEqual(statuses[str(ai)]["completed_rounds"], 0)
            self.assertEqual(statuses[str(ai)]["stop_reason"], "memory_limit")
            self.assertIsNone(statuses[str(ai)]["output_file"])

    def test_memory_stop_preserves_saved_rounds_without_aggregate(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mapping = root / "A_2_atom_mapping.xlsx"
            mapping.touch()
            saved_hash = []
            def stopped(strain, job, rounds, mapping_file, output):
                output_workbook(output, summary([20, 30]))
                saved_hash.append(output.read_bytes())
                raise MemoryPressureError("simulated stop")
            jobs = {"A": {"mapping_file": mapping, "carbon_sources": ["A[c]"]}}
            with patch("carbon_flow_tools.pipeline.trace_mapping_file", stopped):
                with self.assertRaisesRegex(MemoryPressureError, "simulated stop"):
                    run_carbon_flow(jobs, 10, root, max_processes=1)
            result = root / "A_3_carbon_flow.xlsx"
            self.assertEqual(result.read_bytes(), saved_hash[0])
            self.assertEqual(pd.read_excel(result, sheet_name="Summary", engine="openpyxl")
                             .columns.tolist(), summary([20, 30]).columns.tolist())
            aggregate = root / "carbon_flow_summary.xlsx"
            report = json.loads((root / "carbon_flow_run.json").read_text())
            self.assertEqual(report["status"], "error")
            self.assertEqual(report["models"][0]["stop_reason"], "memory_limit")
            self.assertEqual(report["models"][0]["completed_rounds"], 2)
            self.assertFalse(aggregate.exists())
            with patch("carbon_flow_tools.pipeline.trace_mapping_file",
                       side_effect=MemoryPressureError("no new rounds")):
                with self.assertRaises(MemoryPressureError):
                    run_carbon_flow(jobs, 10, root, max_processes=1)
            self.assertFalse(aggregate.exists())
            report = json.loads((root / "carbon_flow_run.json").read_text())
            self.assertEqual(report["models"][0]["completed_rounds"], 0)
            self.assertIsNone(report["models"][0]["output_file"])
            self.assertEqual(result.read_bytes(), saved_hash[0])


if __name__ == "__main__":
    unittest.main()
