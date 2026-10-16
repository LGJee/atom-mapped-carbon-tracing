"""Check independent Enumeration/size-1 EMU dispatch, provenance and exports."""
from hashlib import sha256
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from openpyxl import Workbook

from carbon_flow_tools.cli import main
from carbon_flow_tools.excel_io import write_sheet, save_output
import step3_carbon_flow as step3


def fixture(root):
    source = root / "Tiny_1_input.xlsx"
    mapping = root / "Tiny_2_atom_mapping.xlsx"
    rows = pd.DataFrame([
        dict(ReactionID="inF", ReactionName="Formate input", Flux=1, Equation=" -> F[c]"),
        dict(ReactionID="inA", ReactionName="Second carbon input", Flux=1, Equation=" -> A[c]"),
        dict(ReactionID="rF", ReactionName="CO2 release", Flux=1, Equation="F[c] -> co2[c]"),
        dict(ReactionID="rP", ReactionName="CO2 fixation", Flux=1, Equation="co2[c] + A[c] -> P[c]"),
        dict(ReactionID="out", ReactionName="Product export", Flux=1, Equation="P[c] -> ")])
    metabolites = pd.DataFrame([dict(MetaboliteID=n, MetaboliteName=n, SMILES_Manual=s)
        for n, s in (("F[c]", "O=CO"), ("A[c]", "C"), ("co2[c]", "O=C=O"), ("P[c]", "CC(=O)O"))])
    workbook = Workbook()
    del workbook[workbook.active.title]
    write_sheet(workbook, "Nonzero", rows)
    write_sheet(workbook, "Metabolites", metabolites)
    save_output(workbook, source)
    digest = sha256(source.read_bytes()).hexdigest()
    final = pd.DataFrame([
        dict(ReactionID="rF", Equation="F[c] -> co2[c]", MappingDirection="forward",
             MappedReactionSMILES="[CH:1](=[O:2])[OH:3]>>[O:2]=[C:1]=[O:3]", ReviewRequired=False, AIGenerated=False),
        dict(ReactionID="rP", Equation="co2[c] + A[c] -> P[c]", MappingDirection="forward",
             MappedReactionSMILES="[O:2]=[C:1]=[O:3].[CH4:4]>>[CH3:4][C:1](=[O:2])[OH:3]", ReviewRequired=False, AIGenerated=False)])
    write_sheet(workbook, "Mapping_Info", pd.DataFrame([dict(SourceInputSHA256=digest)]))
    write_sheet(workbook, "Final_Mapping", final)
    save_output(workbook, mapping)
    ai = root / "Tiny_2_atom_mapping_AIGenerated.xlsx"
    final["AIGenerated"] = True
    write_sheet(workbook, "Final_Mapping", final)
    save_output(workbook, ai)
    workbook.close()
    plan = root / "analysis_plan.json"
    plan.write_text(json.dumps(dict(version=3, output_dir=str(root), models=dict(Tiny=dict(
        input_file=str(source), input_sha256=digest, carbon_sources=["F[c]", "A[c]"])))), encoding="utf-8")
    options = step3.settings()
    options.update(analysis_plan=str(plan), tracing_python=sys.executable, rounds=2, max_carbon_processes=1,
                   output_batch_rows=2)
    return source, mapping, ai, plan, options


class Step3MethodTests(unittest.TestCase):
    def test_default_both_methods_preserve_inputs_and_different_co2_policies(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source, mapping, ai, plan, options = fixture(root)
            before = {p: p.read_bytes() for p in (source, mapping, ai, plan)}
            self.assertEqual(step3.ANALYSIS_METHODS, ["Enumeration", "EMU"])
            with patch("carbon_flow_tools.pipeline.run_mapper", side_effect=AssertionError("Mapper called")):
                self.assertEqual(main("trace", options, []), 0)
            expected = [root / ("Tiny_3_carbon_flow" + s + ".xlsx")
                        for s in ("_Enumeration", "_Enumeration_AIGenerated", "_EMU", "_EMU_AIGenerated")]
            self.assertTrue(all(p.is_file() for p in expected))
            for p, contents in before.items():
                self.assertEqual(p.read_bytes(), contents)
            enum = pd.read_excel(expected[0], sheet_name="Carbon_Flow_R1", engine="openpyxl")
            self.assertEqual(set(enum.loc[enum.Production == "P[c]", "CarbonLabels"]), {"C2"})
            summary = pd.read_excel(expected[0], sheet_name="Summary", engine="openpyxl")
            self.assertTrue(summary.Method.eq("Enumeration").all())
            self.assertTrue(summary.InorganicCarbonRecycling.eq(False).all())
            emu = pd.read_excel(expected[-1], sheet_name="MetaboliteContributions", engine="openpyxl").set_index("MetaboliteID")
            self.assertAlmostEqual(emu.loc["P[c]", "F[c]"], .5)
            self.assertAlmostEqual(emu.loc["P[c]", "A[c]"], .5)
            self.assertAlmostEqual(emu.loc["P[c]", "Unknown"], 0)
            with pd.ExcelFile(expected[-1], engine="openpyxl") as excel:
                self.assertEqual(excel.sheet_names, ["MetaboliteContributions", "BoundaryInputs"])
            self.assertNotIn("FluxResidual", emu.columns)
            self.assertNotIn("Status", emu.columns)
            report = json.loads((root/"source_contributions_run.json").read_text())
            self.assertIn("Size-1", report["models"][0]["emu_scope"])
            self.assertEqual(report["models"][0]["allow_inorganic_carbon_recycling"], True)
            self.assertEqual(report["method"], "EMU")
            self.assertEqual(len(report["models"]), 2)
            self.assertEqual({m["mapping_version"] for m in report["models"]}, {"Original", "AI"})
            self.assertEqual({m["mapping_file"] for m in report["models"]}, {str(mapping), str(ai)})
            self.assertEqual(report["models"][0]["mapping_file"], str(mapping))
            self.assertFalse((root/"carbon_flow_summary.xlsx").exists())

    def test_emu_only_reads_both_copies_ignores_rounds_and_honors_output_directory(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            _, _, ai, plan, options = fixture(root)
            destination = root / "separate"
            with patch("carbon_flow_tools.pipeline.run_carbon_flow", side_effect=AssertionError("Enumeration called")):
                self.assertEqual(main("trace", options, ["--method", "EMU", "--rounds", "99", "--output-dir", str(destination)]), 0)
            self.assertTrue((destination/"Tiny_3_carbon_flow_EMU.xlsx").is_file())
            self.assertTrue((destination/"Tiny_3_carbon_flow_EMU_AIGenerated.xlsx").is_file())
            self.assertFalse((destination/"carbon_flow_run.json").exists())
            self.assertFalse((destination/"carbon_flow_summary.xlsx").exists())
            self.assertFalse(any(destination.glob("*Enumeration*.xlsx")))

    def test_emu_accepts_ai_only_mapping(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            _, mapping, ai, _, options = fixture(root)
            mapping.unlink()
            self.assertEqual(main("trace", options, ["--method", "EMU"]), 0)
            self.assertFalse((root/"Tiny_3_carbon_flow_EMU.xlsx").exists())
            output = root/"Tiny_3_carbon_flow_EMU_AIGenerated.xlsx"
            self.assertTrue(output.is_file())
            report = json.loads((root/"source_contributions_run.json").read_text())
            self.assertEqual(len(report["models"]), 1)
            self.assertEqual(report["models"][0]["mapping_file"], str(ai))
            self.assertEqual(report["models"][0]["mapping_version"], "AI")

    def test_emu_versions_use_their_own_final_mapping_without_modifying_inputs(self):
        from openpyxl import load_workbook
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source, mapping, ai, plan, options = fixture(root)
            w = load_workbook(mapping)
            sheet = w["Final_Mapping"]
            sheet.delete_rows(2)
            save_output(w, mapping)
            w.close()
            before = {p: p.read_bytes() for p in (source, mapping, ai, plan)}
            self.assertEqual(main("trace", options, ["--method", "EMU"]), 0)
            for suffix, unknown in (("", .5), ("_AIGenerated", 0)):
                output = root/("Tiny_3_carbon_flow_EMU" + suffix + ".xlsx")
                d = pd.read_excel(output, sheet_name="MetaboliteContributions", engine="openpyxl").set_index("MetaboliteID")
                self.assertAlmostEqual(d.loc["P[c]", "Unknown"], unknown)
                self.assertAlmostEqual(d.loc["P[c]", "F[c]"], .5 - unknown)
            for p, contents in before.items():
                self.assertEqual(p.read_bytes(), contents)
            report = json.loads((root/"source_contributions_run.json").read_text())
            for m in report["models"]:
                self.assertEqual(m["mapping_sha256"], sha256(before[Path(m["mapping_file"])]).hexdigest())

    def test_ai_provenance_and_embedded_flux_are_still_validated(self):
        from openpyxl import load_workbook
        for invalid in ("hash", "flux"):
            with self.subTest(invalid=invalid), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                _, _, ai, _, options = fixture(root)
                w = load_workbook(ai)
                if invalid == "hash":
                    w["Mapping_Info"].cell(2, 1, "wrong hash")
                else:
                    sheet = w["Nonzero"]
                    headers = {c.value: c.column for c in sheet[1]}
                    sheet.cell(2, headers["Flux"], 2)
                save_output(w, ai)
                w.close()
                error = ValueError if invalid == "hash" else AssertionError
                with self.assertRaises(error):
                    main("trace", options, ["--method", "EMU"])
                self.assertFalse((root/"Tiny_3_carbon_flow_EMU_AIGenerated.xlsx").exists())

    def test_enumeration_only_needs_no_plan_or_original_input(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source, _, _, _, options = fixture(root)
            source.unlink()
            options.pop("analysis_plan")
            options.update(input_dir=str(root), output_dir=str(root), carbon_sources={"Tiny": ["F[c]", "A[c]"]})
            with patch("carbon_flow_tools.source_contributions.run_source_contributions", side_effect=AssertionError("EMU called")):
                self.assertEqual(main("trace", options, ["--method", "Enumeration"]), 0)
            self.assertTrue((root/"Tiny_3_carbon_flow_Enumeration.xlsx").is_file())
            self.assertFalse((root/"Tiny_3_carbon_flow_EMU.xlsx").exists())
            self.assertFalse((root/"carbon_flow_summary.xlsx").exists())

    def test_failed_emu_does_not_cancel_enumeration_or_report_success(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            _, _, _, _, options = fixture(root)
            with patch("carbon_flow_tools.source_contributions.solve_source_contributions", side_effect=ValueError("test solver failure")):
                with self.assertRaisesRegex(ValueError, "test solver failure"):
                    main("trace", options, [])
            self.assertTrue((root/"Tiny_3_carbon_flow_Enumeration.xlsx").is_file())
            self.assertFalse((root/"Tiny_3_carbon_flow_EMU.xlsx").exists())
            self.assertFalse((root/"carbon_flow_summary.xlsx").exists())
            self.assertEqual(json.loads((root/"source_contributions_run.json").read_text())["status"], "error")
            self.assertEqual(json.loads((root/"carbon_flow_run.json").read_text())["status"], "finished")


if __name__ == "__main__":
    unittest.main()
