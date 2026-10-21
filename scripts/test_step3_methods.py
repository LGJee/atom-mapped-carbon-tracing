"""Verify combined method output and bounded Enumeration detail export."""
from hashlib import sha256
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from openpyxl import load_workbook

from carbon_flow_tools.carbon_trace import COLUMNS
from carbon_flow_tools.flow_output import FlowOutput
from carbon_flow_tools.pipeline import finish_analysis, run_carbon_flow
from carbon_flow_tools.reactions import prepare_reactions


class CombinedMethodTests(unittest.TestCase):
    def test_round_details_stop_at_three_but_summaries_and_emu_remain(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "flow.xlsx"
            summaries = []
            emu = {"EMU_Metabolites": pd.DataFrame([["P", 1.0]], columns=["MetaboliteID", "G"])}
            store = FlowOutput(path, 50000000, 1, summaries, None, detail_rounds=3, extra_sheets=emu)
            try:
                for number in range(1, 6):
                    store.start_round(number)
                    store.start_layer()
                    if number > 3:
                        store.MAX_DATA_ROWS = 0  # Non-exported rows have no worksheet row limit.
                    store.append(dict(zip(COLUMNS, [1, "r", "G", "1", "P", number, "C1", "[CH4:1]"])))
                    self.assertEqual(list(store.source_records()), [("P", "C1")])
                    summaries.append(dict(Round=number, All_sources_LabeledCarbonMetabolites=1))
                    store.complete_round()
                workbook = load_workbook(path, read_only=True)
                try:
                    self.assertEqual(set(workbook.sheetnames), {"Summary", "Carbon_Flow_R1", "Carbon_Flow_R2", "Carbon_Flow_R3", "EMU_Metabolites"})
                    self.assertEqual(workbook["Summary"].max_row, None)  # Read-only streaming dimensions.
                    self.assertEqual(len(list(workbook["Summary"].values)), 6)
                finally:
                    workbook.close()
                self.assertEqual(store.db.execute("SELECT DISTINCT round_no FROM records ORDER BY round_no").fetchall(), [(1,), (2,), (3,), (5,)])
            finally:
                store.close()

    def test_both_versions_read_their_own_mapping_and_emit_one_workbook_each(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = pd.DataFrame([["in", "Input", 1, " => A[e]"], ["r", "Transfer", 1, "A[e] => B[c]"], ["out", "Output", 1, "B[c] => "]],
                columns=["ReactionID", "ReactionName", "Flux", "Equation"])
            metabolites = pd.DataFrame([[n, n, "CCO"] for n in ("A[e]", "B[c]")],
                columns=["MetaboliteID", "MetaboliteName", "SMILES_Manual"])
            original = root / "input.xlsx"
            with pd.ExcelWriter(original, engine="openpyxl") as writer:
                source.to_excel(writer, sheet_name="Nonzero", index=False)
                metabolites.to_excel(writer, sheet_name="Metabolites", index=False)
            reactions = prepare_reactions(source, metabolites)
            reactions[1].status = "mapped"
            reactions[1].mapped_smiles = "[CH3:1][CH2:2][OH:3]>>[CH3:1][CH2:2][OH:3]"
            finish_analysis("Test", dict(input_file=original, metabolite_sheet="Metabolites"), metabolites,
                {"RXNMapper": reactions}, [], 1, ["RXNMapper"], root, stage="map")
            mapping = root / "Test_2_atom_mapping.xlsx"
            ai = root / "Test_2_atom_mapping_AIGenerated.xlsx"
            shutil.copy2(str(mapping), str(ai))
            before = {path: sha256(path.read_bytes()).hexdigest() for path in (original, mapping, ai)}
            jobs = {"Test": dict(mapping_files=[mapping, ai], carbon_sources=["A[e]"], analysis_methods=["Enumeration", "EMU"],
                export_detail_rounds=1, min_available_memory_gb=0.01)}
            with patch("carbon_flow_tools.pipeline.run_mapper", side_effect=AssertionError("Mapper called during step 3")):
                results = run_carbon_flow(jobs, 3, root / "combined", max_processes=1)
            self.assertTrue(all(row["emu_completed"] for row in results))
            self.assertEqual({row["emu"]["mapping_version"] for row in results}, {"Original", "AI"})
            for suffix in ("", "_AIGenerated"):
                output = root / "combined" / ("Test_3_carbon_flow" + suffix + ".xlsx")
                with pd.ExcelFile(output, engine="openpyxl") as excel:
                    self.assertIn("EMU_Metabolites", excel.sheet_names)
                    self.assertIn("EMU_BoundaryInputs", excel.sheet_names)
                    self.assertIn("Carbon_Flow_R1", excel.sheet_names)
                    self.assertNotIn("Carbon_Flow_R2", excel.sheet_names)
                    self.assertGreaterEqual(pd.read_excel(excel, sheet_name="Summary").Round.max(), 2)
                    table = pd.read_excel(excel, sheet_name="EMU_Metabolites").set_index("MetaboliteName")
                    self.assertAlmostEqual(table.loc["B[c]", "A[e]"], 1)
                    self.assertAlmostEqual(table.loc["B[c]", "Unknown"], 0)
            self.assertFalse(list((root / "combined").glob("*_EMU*.xlsx")))
            self.assertTrue((root / "combined/carbon_flow_summary.xlsx").is_file())
            from carbon_flow_tools.excel_io import WorkbookSizeLimitError
            previous = {p: p.read_bytes() for p in (root / "combined").glob("Test_3_carbon_flow*.xlsx")}
            with patch("carbon_flow_tools.flow_output.FlowOutput.checkpoint", side_effect=WorkbookSizeLimitError("test first-round rejection")):
                rejected = run_carbon_flow(jobs, 3, root / "combined", max_processes=1)
            self.assertTrue(all(row["emu_solved"] and not row["emu_completed"] for row in rejected))
            self.assertTrue(all(row["output_file"] is None and row["completed_rounds"] == 0 for row in rejected))
            self.assertEqual(previous, {p: p.read_bytes() for p in previous})
            jobs["Test"]["analysis_methods"] = ["EMU"]
            run_carbon_flow(jobs, 3, root / "emu_only", max_processes=1)
            with pd.ExcelFile(root / "emu_only/Test_3_carbon_flow.xlsx", engine="openpyxl") as excel:
                self.assertEqual(set(excel.sheet_names), {"EMU_Metabolites", "EMU_BoundaryInputs"})
            self.assertFalse((root / "emu_only/carbon_flow_summary.xlsx").exists())
            self.assertEqual(before, {p: sha256(p.read_bytes()).hexdigest() for p in before})


if __name__ == "__main__":
    unittest.main()
