"""Formula spelling differences must not bypass real composition/charge checks."""
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import pandas as pd

from carbon_flow_tools.carbon_mapping_rules import formula_counts, formulas_match, fixed_formula_carbons, protonation_compatible
from carbon_flow_tools.model_structures import compatible
from carbon_flow_tools.smiles_completion import complete_smiles, select_candidate


class FormulaMatchingTests(unittest.TestCase):
    def test_equivalent_explicit_one_order_and_repeated_elements(self):
        for first, second in [("C35H67O8P1", "C35H67O8P"), ("C1H4", "CH4"),
                              ("H6OC2", "C2H6O"), ("CH3CH2OH", "C2H6O")]:
            self.assertTrue(formulas_match(first, second))
        self.assertEqual(formula_counts("C35H67O8P1"), {"C":35,"H":67,"O":8,"P":1})
        self.assertEqual(fixed_formula_carbons("CH3CH2OH"), 2)
        self.assertEqual(fixed_formula_carbons("Co"), 0)

    def test_real_composition_differences_remain_conflicts(self):
        for first, second in [("C35H67O8P", "C35H69O8P"), ("C35H67O8P", "C35H67O8P2"),
                              ("Co", "CO"), ("Cl", "C")]:
            self.assertFalse(formulas_match(first, second))

    def test_unsupported_formulas_do_not_compare_equal_from_failed_parses(self):
        for first, second in [("C20H37O4P(C5H8)n", "C25H45O4P"), ("C2H4R", "C2H4X"),
                              ("C2H4R", "C2H4R2"), ("bad formula", "other text"), ("NH4+", "NH4")]:
            self.assertFalse(formulas_match(first, second))
        self.assertIsNone(formula_counts("C20H37O4P(C5H8)n"))
        self.assertIsNone(fixed_formula_carbons("C20H37O4P(C5H8)n"))
        self.assertIsNone(formula_counts("C2H4R"))

    def test_manual_completion_accepts_one_but_checks_formula_and_charge(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"manual.csv"
            path.write_text("Species,MetaboliteName,SMILES\nEscherichia_coli,phosphate,OP(=O)(O)O\n")
            original = path.read_bytes()
            options = dict(manual_file=path,cache_file=Path(folder)/"cache.json",lookup=False,timeout=1,min_interval=0)
            table = pd.DataFrame([dict(MetaboliteID="p",MetaboliteName="phosphate[c]",MolecularFormula="H3O4P1",Charge=0)])
            unchanged = table.copy(deep=True)
            enriched = complete_smiles(table,"Escherichia_coli_Test",options)
            self.assertEqual(enriched.iloc[0].SMILES_Manual,"OP(=O)(O)O")
            pd.testing.assert_frame_equal(table,unchanged)
            self.assertEqual(path.read_bytes(),original)
            for formula,charge,error in [("H3O4P2",0,"structure/formula conflict"),
                                          ("H3O4P1",-1,"structure/charge conflict")]:
                wrong = table.assign(MolecularFormula=formula,Charge=charge)
                with self.assertRaisesRegex(ValueError,error):
                    complete_smiles(wrong,"Escherichia_coli_Test",options)

    def test_protonation_compatibility_requires_hydrogen_charge_balance(self):
        self.assertTrue(protonation_compatible("C19H37O7P", -2, "C19H39O7P", 0))
        self.assertTrue(protonation_compatible("C19H39O7P", 0, "C19H37O7P", -2))
        self.assertTrue(protonation_compatible("C2H8N", 1, "C2H7N", 0))
        for formula, charge in [("C19H37O7P", 0), ("C20H37O7P", -2),
                                ("C19H37O8P", -2), ("C19H39O7P", -2),
                                ("C19H37O7P", ""), ("C19H37O7P", float("nan")),
                                ("C19H37O7P", -1.5), ("C19H37O7P(C5H8)n", -2)]:
            self.assertFalse(protonation_compatible(formula, charge, "C19H39O7P", 0))

    def test_protonation_mismatch_warns_without_editing_inputs(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"manual.csv"
            path.write_text("Species,MetaboliteName,SMILES\nEscherichia_coli,lipid,CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(=O)(O)O\n")
            original = path.read_bytes()
            options = dict(manual_file=path,cache_file=Path(folder)/"cache.json",lookup=False,timeout=1,min_interval=0)
            table = pd.DataFrame([dict(MetaboliteID="p",MetaboliteName="lipid[c]",MolecularFormula="C19H37O7P1",Charge=-2)])
            unchanged = table.copy(deep=True)
            enriched = complete_smiles(table,"Escherichia_coli_Test",options)
            self.assertTrue(enriched.iloc[0].StructureReviewRequired)
            self.assertIn("Protonation-compatible mismatch", enriched.iloc[0].StructureStatus)
            self.assertEqual(enriched.iloc[0].Charge,-2)
            self.assertEqual(enriched.iloc[0].MolecularFormula,"C19H37O7P1")
            self.assertTrue(enriched.iloc[0].SMILES_Manual.endswith("P(=O)(O)O"))
            pd.testing.assert_frame_equal(table,unchanged)
            self.assertEqual(path.read_bytes(),original)
            # Hydrogen changes without the corresponding charge change may be redox chemistry.
            with self.assertRaisesRegex(ValueError,"structure/formula conflict"):
                complete_smiles(table.assign(Charge=0),"Escherichia_coli_Test",options)

    def test_protonation_review_flag_reaches_final_mapping_and_loader(self):
        from carbon_flow_tools.reactions import prepare_reactions
        from carbon_flow_tools.pipeline import finish_analysis, load_final_mappings
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            manual = root/"manual.csv"
            manual.write_text("Species,MetaboliteName,SMILES\nEscherichia_coli,acetate,CC(=O)O\n")
            metabolites = pd.DataFrame([dict(MetaboliteID=name,MetaboliteName=name,
                MolecularFormula="C2H3O2",Charge=-1) for name in ["acetate[e]","acetate[c]"]])
            source = pd.DataFrame([["T","Transport",1,"acetate[e] => acetate[c]"]],
                columns=["ReactionID","ReactionName","Flux","Equation"])
            options = dict(manual_file=manual,cache_file=root/"cache.json",lookup=False,timeout=1,min_interval=0)
            enriched = complete_smiles(metabolites,"Escherichia_coli_Test",options)
            path = root/"input.xlsx"
            with pd.ExcelWriter(path,engine="openpyxl") as writer:
                source.to_excel(writer,sheet_name="Nonzero",index=False)
                metabolites.to_excel(writer,sheet_name="Metabolites",index=False)
            before = path.read_bytes()
            reactions = prepare_reactions(source,enriched)
            finish_analysis("Escherichia_coli_Test",dict(input_file=path,metabolite_sheet="Metabolites"),
                enriched,{"LocalMapper":reactions},[],["LocalMapper"],root,stage="map")
            output = root/"Escherichia_coli_Test_2_atom_mapping.xlsx"
            final = pd.read_excel(output,sheet_name="Final_Mapping",engine="openpyxl")
            self.assertTrue(final.iloc[0].ReviewRequired)
            self.assertTrue(final.iloc[0].MappedReactionSMILES)
            loaded, _, saved = load_final_mappings(output,{})
            self.assertEqual(loaded[0].status,"mapped")
            self.assertTrue(saved.StructureReviewRequired.all())
            self.assertEqual(path.read_bytes(),before)

    def test_database_candidates_use_same_composition_rule(self):
        row = SimpleNamespace(MolecularFormula="H3O4P1",Charge=0)
        record = dict(CID=1,SMILES="OP(=O)(O)O",MolecularFormula="H3O4P",Charge=0)
        self.assertEqual(len(select_candidate([record],row)),1)
        self.assertEqual(select_candidate([dict(record,MolecularFormula="H3O4P2")],row),[])
        self.assertEqual(select_candidate([dict(record,Charge=-1)],row),[])

    def test_model_annotations_use_same_composition_rule(self):
        row = SimpleNamespace(MolecularFormula="C35H67O8P1",MetaboliteName="compound[c]")
        annotation = dict(name="compound",compartment="c",formula="C35H67O8P")
        self.assertTrue(compatible(annotation,row))
        self.assertFalse(compatible(dict(annotation,formula="C35H69O8P"),row))
        self.assertFalse(compatible(dict(annotation,name="other"),row))


if __name__ == "__main__":
    unittest.main()
