"""Three-column manual mappings override selection and survive standalone handoff."""
import csv
from pathlib import Path
import tempfile
import unittest
import pandas as pd
from carbon_flow_tools.carbon_mapping_rules import apply_manual_mappings, apply_rules
from carbon_flow_tools.consensus import select_mappings
from carbon_flow_tools.reactions import prepare_reactions
from carbon_flow_tools.pipeline import finish_analysis, load_final_mappings
from carbon_flow_tools.source_contributions import solve_source_contributions
from scripts.test_reduced_mapper_inputs import fixture, MAPPED

UNIT = "[CH3:1][CH2:2]O>>[CH3:1][CH2:2]O"

class ManualTableTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.manual = self.root / "manual.csv"
    def tearDown(self):
        self.directory.cleanup()
    def write(self, rows):
        with self.manual.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["Species", "ReactionID", "MappedReactionSMILES"])
            writer.writerows(rows)
    def test_species_conditions_and_source_directory(self):
        self.write([["Escherichia_coli", "R", UNIT]])
        for model, source in (("Escherichia_coli_non_IUPD", None),
                              ("Escherichia_coli_IUPD", None),
                              ("GlucoseOnly", dict(path="data/Escherichia_coli/iML1515.mat"))):
            _, mets, _, selected, choices = fixture()
            apply_manual_mappings(model, selected, choices, mets, self.manual, source)
            self.assertEqual(choices.iloc[0].SelectedMapper, "Manual")
            self.assertEqual(choices.iloc[0].SelectedMappedReactionSMILES, UNIT)
            self.assertEqual(selected[0].partial_maps[0]["TransferCoefficient"], 2)
        _, mets, _, selected, choices = fixture()
        apply_manual_mappings("Saccharomyces_cerevisiae_IUPD", selected, choices, mets, self.manual)
        self.assertEqual(choices.iloc[0].SelectedMapper, "EnzymeMap")
    def test_complete_mapping_written_unchanged_to_both_sheets(self):
        source, mets, _, _, _ = fixture()
        source.loc[0, "Equation"] = "A[c] => B[c]"
        predictions = {"RXNMapper": prepare_reactions(source, mets)}
        mapped = "[OH:3][CH2:2][CH3:1]>>[CH3:1][CH2:2][OH:3]"
        self.write([["Escherichia_coli", "R", mapped]])
        path = self.root / "input.xlsx"
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            source.to_excel(writer, sheet_name="Nonzero", index=False)
            mets.to_excel(writer, sheet_name="Metabolites", index=False)
        model = "Escherichia_coli_non_IUPD"
        finish_analysis(model, dict(input_file=path, metabolite_sheet="Metabolites"), mets,
                        predictions, [], ["RXNMapper"], self.root, stage="map", manual_mapping_file=self.manual)
        output = self.root / (model + "_2_atom_mapping.xlsx")
        selection = pd.read_excel(str(output), sheet_name="Mapping_Selection", engine="openpyxl")
        final = pd.read_excel(str(output), sheet_name="Final_Mapping", engine="openpyxl")
        self.assertEqual(selection.iloc[0].SelectedMappedReactionSMILES, mapped)
        self.assertEqual(final.iloc[0].MappedReactionSMILES, mapped)
        self.assertEqual(final.iloc[0].MappingCoverage, "Complete")
        self.assertFalse(final.iloc[0].ReviewRequired)
        self.assertFalse(final.iloc[0].AIGenerated)
        loaded, _, _ = load_final_mappings(output, {})
        self.assertEqual(loaded[0].status, "mapped")
    def test_partial_copies_and_fractional_amounts_roundtrip(self):
        for amount, mapped in ((2, MAPPED), (55.3, UNIT)):
            source, mets, _, _, _ = fixture()
            source.loc[0, "Equation"] = "{} A[c] + ReducedProtein[c] => {} B[c] + OxidizedProtein[c]".format(amount, amount)
            predictions = {"RXNMapper": prepare_reactions(source, mets)}
            self.write([["Escherichia_coli", "R", mapped]])
            path = self.root / "input.xlsx"
            with pd.ExcelWriter(path, engine="openpyxl") as writer:
                source.to_excel(writer, sheet_name="Nonzero", index=False)
                mets.to_excel(writer, sheet_name="Metabolites", index=False)
            model = "Escherichia_coli_IUPD"
            finish_analysis(model, dict(input_file=path, metabolite_sheet="Metabolites"), mets,
                            predictions, [], ["RXNMapper"], self.root, stage="map", manual_mapping_file=self.manual)
            output = self.root / (model + "_2_atom_mapping.xlsx")
            self.manual.unlink()
            loaded, _, enriched = load_final_mappings(output, {})
            self.assertEqual(loaded[0].mapped_smiles, mapped)
            self.assertEqual(len(loaded[0].partial_maps), 1)
            self.assertAlmostEqual(float(loaded[0].partial_maps[0]["TransferCoefficient"]), amount)
            final = pd.read_excel(str(output), sheet_name="Final_Mapping", engine="openpyxl")
            self.assertTrue(final.iloc[0].ReviewRequired)
            self.assertEqual(final.iloc[0].MappingCoverage, "Partial")
            boundaries = pd.DataFrame([["in", "", amount, "=> A[c]"], ["out", "", amount, "B[c] =>"]], columns=source.columns)
            tables, _ = solve_source_contributions(loaded + prepare_reactions(boundaries, mets), enriched, ["A[c]"])
            row = tables["MetaboliteContributions"].set_index("MetaboliteName").loc["B[c]"]
            self.assertAlmostEqual(row["A[c]"], 1)
    def test_negative_flux_reverses_equation_direction_mapping(self):
        source, mets, _, _, _ = fixture(-1)
        mets.loc[mets.MetaboliteID.eq("b"), "SMILES_Model"] = "CC=O"
        selected, choices = select_mappings({"RXNMapper": prepare_reactions(source, mets)})
        mapped = "[CH3:1][CH2:2]O>>[CH3:1][CH:2]=O"
        self.write([["Escherichia_coli", "R", mapped]])
        apply_manual_mappings("Escherichia_coli_IUPD", selected, choices, mets, self.manual)
        self.assertEqual(selected[0].mapped_smiles, ">>".join(reversed(mapped.split(">>"))))
        self.assertEqual(selected[0].partial_maps[0]["SourceMetabolite"], "B[c]")
    def test_reject_invalid_duplicate_and_ambiguous_records(self):
        _, mets, _, selected, choices = fixture()
        for rows in ([ ["Escherichia_coli", "R", UNIT] ] * 2,
                     [["Escherichia_coli", "R", "not-smiles>>CC"]],
                     [["Escherichia_coli", "R", "[CH3:1][CH2:1]O>>CCO"]],
                     [["Escherichia_coli", "R", "CCO>>CCO"]]):
            self.write(rows)
            with self.assertRaises(ValueError):
                apply_manual_mappings("Escherichia_coli_IUPD", selected, choices, mets, self.manual)
        self.write([["Escherichia_coli", "R", UNIT]])
        source, mets, _, _, _ = fixture()
        source.loc[0, "Equation"] = "A[c] + D[c] => B[c]"
        mets.loc[len(mets)] = ["d", "D[c]", "CCO"]
        selected, choices = select_mappings({"RXNMapper": prepare_reactions(source, mets)})
        with self.assertRaisesRegex(ValueError, "one original metabolite"):
            apply_manual_mappings("Escherichia_coli_IUPD", selected, choices, mets, self.manual)
    def test_reject_inconsistent_stoichiometry_and_legacy_schema(self):
        source, mets, _, _, _ = fixture()
        source.loc[0, "Equation"] = "2 A[c] => 3 B[c]"
        selected, choices = select_mappings({"RXNMapper": prepare_reactions(source, mets)})
        self.write([["Escherichia_coli", "R", UNIT]])
        with self.assertRaisesRegex(ValueError, "stoichiometry"):
            apply_manual_mappings("Escherichia_coli_IUPD", selected, choices, mets, self.manual)
        self.manual.write_text("Model,ReactionID,MappedReactionSMILES\n")
        with self.assertRaisesRegex(ValueError, "exactly"):
            apply_manual_mappings("Escherichia_coli_IUPD", selected, choices, mets, self.manual)

if __name__ == "__main__":
    unittest.main()
