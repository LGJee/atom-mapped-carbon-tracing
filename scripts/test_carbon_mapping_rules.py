"""Check partial carbon rules, workbook handoff and actual label propagation."""

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd
from rdkit import Chem

from carbon_flow_tools.carbon_mapping_rules import (
    _REFERENCES, carbon_pairs, apply_rules, load_partials, rule_transfers)
from carbon_flow_tools.carbon_trace import prepare_mappings, map_one_reaction
from carbon_flow_tools.consensus import select_mappings
from carbon_flow_tools.pipeline import (finish_analysis, load_final_mappings,
                                        carbon_metabolites)
from carbon_flow_tools.reactions import prepare_reactions

NADPH = "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O"


def fixture():
    metabolites = pd.DataFrame([
        [name, name, smiles] for name, smiles in [
            ("ATP[c]", _REFERENCES["ATP"]), ("ADP[c]", _REFERENCES["ADP"]),
            ("NADPH[c]", NADPH), ("NADP[c]", _REFERENCES["NADP"]),
            ("protein[c]", ""), ("biomass[c]", "")]
    ], columns=["MetaboliteID", "MetaboliteName", "SMILES_Manual"])
    source = pd.DataFrame([
        ["B", "Biomass", 1, "55.3 ATP[c] + 3 NADPH[c] + protein[c] => 55.3 ADP[c] + 3 NADP[c] + biomass[c]"],
        ["U", "Unresolved", 1, "protein[c] => biomass[c]"],
        ["E", "Exchange", 1, "ATP[c] => "]
    ], columns=["ReactionID", "ReactionName", "Flux", "Equation"])
    reactions = prepare_reactions(source, metabolites)
    selected, selection = select_mappings({"LocalMapper": deepcopy(reactions)})
    partials = apply_rules(selected, selection, metabolites)
    return source, metabolites, reactions, selected, selection, partials


class CarbonRuleTests(unittest.TestCase):
    def test_adenylates_and_redox(self):
        for a, b in [("ATP", "ADP"), ("ATP", "AMP"), ("ADP", "AMP")]:
            rule, pairs = carbon_pairs(_REFERENCES[a], _REFERENCES[b])
            self.assertEqual(rule, "adenylate_skeleton")
            self.assertEqual(len(pairs), 10)
        self.assertEqual(len(carbon_pairs(NADPH, _REFERENCES["NADP"])[1]), 21)
        self.assertIsNone(carbon_pairs(NADPH, _REFERENCES["NAD"]))
        self.assertIsNone(carbon_pairs(_REFERENCES["ATP"], "CC"))
        self.assertIsNone(carbon_pairs("", _REFERENCES["ATP"]))
        self.assertIsNone(carbon_pairs("CC", "CC"))  # Symmetry not uniquely aligned.
        self.assertIsNone(carbon_pairs("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"))

    def test_carbon_presence_does_not_auto_classify_unmapped_reactions(self):
        from types import SimpleNamespace
        from carbon_flow_tools.carbon_mapping_rules import carbon_presence
        for formula, expected in [("H2O", False), ("Co", False), ("CO2", True), ("C2H4R", True), ("R", None), ("", None)]:
            self.assertIs(carbon_presence(SimpleNamespace(MolecularFormula=formula)), expected)
        metabolites = pd.DataFrame([["a", "Water", "H2O"], ["b", "Hydrogen", "H2"], ["c", "Oxygen", "O2"]],
            columns=["MetaboliteID", "MetaboliteName", "MolecularFormula"])
        source = pd.DataFrame([["R", "No carbon", 1, "2 Water => 2 Hydrogen + Oxygen"]],
            columns=["ReactionID", "ReactionName", "Flux", "Equation"])
        selected, choices = select_mappings({"LocalMapper": prepare_reactions(source, metabolites)})
        self.assertTrue(apply_rules(selected, choices, metabolites).empty)
        self.assertEqual(choices.iloc[0].MappingCoverage, "Unresolved")
        self.assertTrue(choices.iloc[0].ReviewRequired)
        self.assertEqual(selected[0].mapped_smiles, "")

    def test_withdrawn_flavin_and_ispg_fallbacks_are_not_generated(self):
        fmn = "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C[C@H](O)[C@H](O)[C@H](O)COP(=O)(O)O)c2cc1C"
        reduced = "Cc1cc2c(cc1C)N(C[C@H](O)[C@H](O)[C@H](O)COP(=O)(O)O)c1[nH]c(=O)[nH]c(=O)c1N2"
        self.assertIsNone(carbon_pairs(reduced, fmn))
        self.assertNotIn("FAD", _REFERENCES)
        self.assertNotIn("FMN", _REFERENCES)
        a = "C[C@@]1([C@@H](COP(=O)(OP(=O)(O1)[O-])[O-])O)CO"
        b = r"C/C(=C\COP(=O)([O-])OP(=O)([O-])[O-])/CO"
        self.assertIsNone(carbon_pairs(a, b))
        source = pd.DataFrame([["R", "IspG", 1,
            "2 ReducedProtein + MEcPP => 2 OxidizedProtein + HMBPP"]],
            columns=["ReactionID", "ReactionName", "Flux", "Equation"])
        metabolites = pd.DataFrame([["a", "MEcPP", a], ["b", "HMBPP", b],
            ["c", "ReducedProtein", ""], ["d", "OxidizedProtein", ""]],
            columns=["MetaboliteID", "MetaboliteName", "SMILES_Manual"])
        selected, choices = select_mappings({"LocalMapper": prepare_reactions(source, metabolites)})
        self.assertTrue(apply_rules(selected, choices, metabolites).empty)
        self.assertEqual(choices.iloc[0].MappingCoverage, "Unresolved")
        self.assertTrue(choices.iloc[0].ReviewRequired)

    def test_malformed_equations_remain_reviewable_rows(self):
        source, metabolites, original, _, _, _ = fixture()
        original[0].equation = "not a reaction"
        original[0].status = "Equation parsing failed"
        selected, choices = select_mappings({"LocalMapper": original})
        apply_rules(selected, choices, metabolites)
        self.assertEqual(choices.iloc[0].MappingCoverage, "Unresolved")
        self.assertTrue(choices.iloc[0].ReviewRequired)
        self.assertIn("cannot be parsed", choices.iloc[0].SelectionReason)

    def test_formula_mismatch_cannot_be_overruled_by_mapper_priority(self):
        from carbon_flow_tools.carbon_mapping_rules import carbon_structure_issues
        from carbon_flow_tools.reactions import Component, PreparedReaction
        metabolites = pd.DataFrame([["thf", "THF[c]", "C19H21N7O6", "C1CCOC1"]],
                     columns=["MetaboliteID", "MetaboliteName", "MolecularFormula", "SMILES_Manual"])
        issues = carbon_structure_issues(metabolites)
        self.assertIn("THF[c]", issues)
        mol = Chem.MolFromSmiles("C1CCOC1")
        for a in mol.GetAtoms():
            a.SetAtomMapNum(a.GetIdx() + 1)
        marked = Chem.MolToSmiles(mol)
        r = PreparedReaction("bad", "Wrong THF annotation", 1, "THF[c] => THF[m]", "C1CCOC1>>C1CCOC1",
                             [Component("THF[c]", 0, "C1CCOC1")], [Component("THF[m]", 0, "C1CCOC1")],
                             source_text="THF[c]", status="mapped", mapped_smiles=marked + ">>" + marked,
                             is_transport=True)
        selected, choices = select_mappings({"EnzymeMap": [r]}, issues)
        self.assertEqual(selected[0].mapped_smiles, r.mapped_smiles)
        self.assertEqual(selected[0].status, "candidate")
        self.assertTrue(choices.iloc[0].ReviewRequired)
        self.assertIn("formula C19 versus SMILES C4", choices.iloc[0].SelectionReason)
        partials = apply_rules(selected, choices, metabolites)
        self.assertTrue(partials.empty)
        self.assertEqual(choices.iloc[0].MappingCoverage, "Unresolved")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "inconsistent.xlsx"
            both = pd.concat([metabolites, metabolites.assign(MetaboliteID="thfm", MetaboliteName="THF[m]")])
            with pd.ExcelWriter(path, engine="openpyxl") as writer:
                pd.DataFrame([["bad", "Bad annotation", 1, r.equation]],
                    columns=["ReactionID", "ReactionName", "Flux", "Equation"]).to_excel(writer, sheet_name="Reactions", index=False)
                both.to_excel(writer, sheet_name="Metabolites", index=False)
                pd.DataFrame([["bad", "forward", r.equation, r.mapped_smiles]],
                    columns=["ReactionID", "MappingDirection", "Equation", "MappedReactionSMILES"]).to_excel(writer, sheet_name="Final_Mapping", index=False)
            with self.assertRaisesRegex(ValueError, "carbon-inconsistent structures"):
                load_final_mappings(path, {})

    def test_generic_pools_do_not_receive_name_lookup_structures(self):
        from unittest.mock import patch
        from carbon_flow_tools.smiles_completion import complete_smiles
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manual = root / "manual.csv"
            manual.write_text("Species,MetaboliteName,SMILES\n")
            table = pd.DataFrame([[name, name] for name in
                       ("biomass[c]", "protein[c]", "DNA[c]", "L-alanyl-tRNA[c]")],
                       columns=["MetaboliteID", "MetaboliteName"])
            with patch("carbon_flow_tools.smiles_completion.fetch_properties",
                       side_effect=AssertionError("Generic pools must not use name lookup")):
                enriched = complete_smiles(table, "Test", dict(manual_file=manual,
                          cache_file=root / "cache.json", lookup=True, timeout=1, min_interval=0))
            self.assertTrue(enriched.StructureStatus.str.startswith("Generic biomass/").all())
            self.assertTrue(enriched.StructureReviewRequired.all())

    def test_other_nucleotides_and_coa_core(self):
        for prefix, family in (("G", "guanylate"), ("C", "cytidylate"), ("U", "uridylate")):
            rule, pairs = carbon_pairs(_REFERENCES[prefix + "TP"], _REFERENCES[prefix + "DP"])
            self.assertEqual(rule, family + "_skeleton")
            self.assertEqual(len(pairs), 10 if prefix == "G" else 9)
        self.assertIsNone(carbon_pairs(_REFERENCES["ATP"], _REFERENCES["GTP"]))
        coa = _REFERENCES["CoA"]
        acetyl_coa = coa.replace("NCCS", "NCCSC(=O)C")
        rule, pairs = carbon_pairs(acetyl_coa, coa)
        self.assertEqual(rule, "CoA_core")
        self.assertEqual(len(pairs), 19)
        molecule = Chem.MolFromSmiles(acetyl_coa)
        carbonyl = [a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum() == 6
                    and any(n.GetAtomicNum() == 16 for n in a.GetNeighbors())
                    and any(n.GetAtomicNum() == 8 and molecule.GetBondBetweenAtoms(
                          a.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.DOUBLE for n in a.GetNeighbors())][0]
        self.assertNotIn(carbonyl, {i for i, _ in pairs})
        self.assertIsNone(carbon_pairs(acetyl_coa, _REFERENCES["ATP"]))

    def test_mapper_partial_precedes_rules_and_preserves_raw(self):
        from carbon_flow_tools.reactions import PreparedReaction, Component
        from carbon_flow_tools.consensus import carbon_candidate
        raw = "[CH3:1][CH2:2][OH:3]>>[CH3:1][CH2:8][OH:3]"
        r = PreparedReaction("R", "Test", 1, "A[c] => B[c]", "CCO>>CCO",
                             [Component("A[c]", 0, "CCO")], [Component("B[c]", 0, "CCO")],
                             source_text="A[c]", status="mapped", mapped_smiles=raw)
        predictions = {name: [deepcopy(r)] for name in ("LocalMapper", "RXNMapper", "EnzymeMap")}
        selected, choices = select_mappings(predictions)
        self.assertEqual(choices.iloc[0].SelectedMapper, "EnzymeMap")
        self.assertEqual(selected[0].status, "partial")
        self.assertTrue(choices.iloc[0].ReviewRequired)
        self.assertEqual(predictions["EnzymeMap"][0].mapped_smiles, raw)
        metabolites = pd.DataFrame([["a", "A[c]", "CCO"], ["b", "B[c]", "CCO"]],
                     columns=["MetaboliteID", "MetaboliteName", "SMILES_Manual"])
        self.assertTrue(apply_rules(selected, choices, metabolites).empty)
        self.assertEqual(choices.iloc[0].MappingCoverage, "Partial")
        self.assertTrue(selected[0].mapped_smiles)
        self.assertFalse(selected[0].partial_maps)
        prepared = prepare_mappings(selected)
        _, rows, _ = map_one_reaction(*prepared[0], {"A[c]": {1: ((0, 0, "C1"), (0, 1, "C2"))}},
                                     1, {}, 2, {"A[c]": "CCO", "B[c]": "CCO"})
        self.assertEqual(rows[0]["CarbonLabels"], "C1")
        full_carbon = deepcopy(r)
        full_carbon.mapped_smiles = "[CH3:1][CH2:2][OH:1]>>[CH3:1][CH2:2][OH:1]"
        recovered, record = select_mappings({"RXNMapper": [full_carbon]})
        self.assertEqual(recovered[0].status, "mapped")
        self.assertTrue(record.iloc[0].ReviewRequired)
        self.assertIn("carbon-only view", record.iloc[0].SelectionReason)
        for invalid in ["[CH4:1]>>[CH4:1]", "[CH3:1][CH2:1]O>>[CH3:1][CH2:1]O"]:
            bad = deepcopy(r)
            bad.mapped_smiles = invalid
            with self.assertRaises(ValueError):
                carbon_candidate(bad)

    def test_any_nonempty_mapper_candidate_survives_selection_rules_and_export(self):
        source, metabolites, original, _, _, _ = fixture()
        bad_maps = {
            "LocalMapper": "[CH4:11]>>[CH4:11]",
            "RXNMapper": "[CH4:22]>>[CH4:22]",
            "EnzymeMap": "[CH4:33]>>[CH4:33]",
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "input.xlsx"
            with pd.ExcelWriter(path, engine="openpyxl") as writer:
                source.iloc[:1].to_excel(writer, sheet_name="Reactions", index=False)
                metabolites.to_excel(writer, sheet_name="Metabolites", index=False)
            for supplied in ([name] for name in bad_maps):
                self.check_candidate_export(supplied, bad_maps, original, metabolites, root, path)
            self.check_candidate_export(list(bad_maps), bad_maps, original, metabolites, root, path)

    def check_candidate_export(self, supplied, bad_maps, original, metabolites, root, path):
        from unittest.mock import patch
        predictions = {name: deepcopy(original[:1]) for name in bad_maps}
        for name in supplied:
            predictions[name][0].mapped_smiles = bad_maps[name]
            predictions[name][0].status = "Invalid mapping: test diagnostic"
        winner = next(name for name in ("EnzymeMap", "RXNMapper", "LocalMapper") if name in supplied)
        expected = bad_maps[winner]
        # This biomass fixture would otherwise trigger ATP/NADPH fallback rules.
        with patch("carbon_flow_tools.carbon_mapping_rules.rule_transfers",
                   side_effect=AssertionError("A nonempty mapper candidate must precede rules")):
            selected, choices = select_mappings(predictions)
            self.assertTrue(apply_rules(selected, choices, metabolites).empty)
            self.assertEqual(choices.iloc[0].SelectedMappedReactionSMILES, expected)
            self.assertEqual(choices.iloc[0].SelectedMapper, winner)
            self.assertTrue(choices.iloc[0].ReviewRequired)
            self.assertEqual(choices.iloc[0].MappingCoverage, "Unresolved")
            self.assertIn("test diagnostic", choices.iloc[0].SelectionReason)
            self.assertEqual(prepare_mappings(selected), [])
            finish_analysis("Candidate", dict(input_file=path, metabolite_sheet="Metabolites"),
                            metabolites, predictions, [], 1, list(predictions), root, stage="map")
        output = root / "Candidate_2_atom_mapping.xlsx"
        with pd.ExcelFile(output, engine="openpyxl") as excel:
            chosen = pd.read_excel(excel, sheet_name="Mapping_Selection").iloc[0]
            final = pd.read_excel(excel, sheet_name="Final_Mapping").iloc[0]
        self.assertEqual(chosen.SelectedMappedReactionSMILES, expected)
        self.assertEqual(final.MappedReactionSMILES, expected)
        self.assertEqual(final.SelectedMapper, winner)
        self.assertTrue(final.ReviewRequired)
        self.assertEqual(final.MappingCoverage, "Unresolved")
        loaded, _, _ = load_final_mappings(output, {})
        self.assertFalse(loaded[0].mapped_smiles)
        for name in supplied:
            self.assertEqual(predictions[name][0].mapped_smiles, bad_maps[name])

    def test_partial_rules_preserve_fractional_amount_and_flags(self):
        _, metabolites, _, selected, selection, partials = fixture()
        self.assertEqual(selection.MappingCoverage.tolist(), ["Partial", "Unresolved", "Boundary"])
        self.assertTrue(selection.iloc[0].ReviewRequired)
        self.assertEqual(partials.TransferCoefficient.tolist(), [55.3, 3.0])
        self.assertEqual(len(selected[0].partial_maps), 2)
        self.assertIn(">>", selected[0].mapped_smiles)
        self.assertEqual(selection.iloc[0].SelectedMappedReactionSMILES, selected[0].mapped_smiles)
        self.assertNotIn("protein", selected[0].mapped_smiles)
        self.assertEqual(sum(map(len, carbon_metabolites(metabolites, selected)[0].values())), 4)
        self.assertEqual(partials.ReactionID.tolist(), ["B", "B"])
        # A supplied monomer is not a generic protein/RNA pool's atomic structure.
        selected[0].equation = "protein[c] => RNA[c]"
        self.assertEqual(rule_transfers(selected[0], {"protein[c]": "CCO", "RNA[c]": "CCO"}), [])

    def test_real_propagation_and_original_reaction_id(self):
        _, metabolites, _, selected, _, _ = fixture()
        selected[1].reaction_id = "B::carbon_rule:0"  # User IDs cannot collide with internal pairs.
        prepared = prepare_mappings(selected)
        self.assertNotEqual(prepared[0][0].reaction_id, selected[1].reaction_id)
        self.assertEqual(len(prepared), 2)
        self.assertEqual(len(set(r.reaction_id for r, _ in prepared)), 2)
        smiles = dict(zip(metabolites.MetaboliteName, metabolites.SMILES_Manual))
        for reaction, components in prepared:
            name = reaction.left[0].metabolite
            mol = Chem.MolFromSmiles(smiles[name])
            assignment = tuple((0, a.GetIdx(), "C" + str(i + 1))
                        for i, a in enumerate(a for a in mol.GetAtoms() if a.GetAtomicNum() == 6))
            _, rows, _ = map_one_reaction(reaction, components, {name: {1: assignment}},
                                         1, {}, 2, smiles)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["ReactionID"], "B")
            self.assertEqual(len(rows[0]["CarbonLabels"].split(",")), len(assignment))

    def test_reverse_and_ambiguity(self):
        source, metabolites, original, _, _, _ = fixture()
        original[0].flux = -1
        transfers = rule_transfers(original[0], dict(zip(metabolites.MetaboliteName, metabolites.SMILES_Manual)))
        self.assertEqual(transfers[0]["SourceMetabolite"], "ADP[c]")
        self.assertEqual(transfers[0]["MappingDirection"], "reverse")
        original[0].flux = 1
        original[0].equation = "ATP[c] => 2 ADP[c]"
        self.assertEqual(rule_transfers(original[0], dict(zip(metabolites.MetaboliteName, metabolites.SMILES_Manual))), [])
        smiles = dict(zip(metabolites.MetaboliteName, metabolites.SMILES_Manual))
        smiles["AMP[c]"] = _REFERENCES["AMP"]
        original[0].equation = "ATP[c] => ADP[c] + AMP[c]"
        self.assertEqual(rule_transfers(original[0], smiles), [])

    def test_no_rule_overwrites_complete_mapping(self):
        _, metabolites, original, _, _, _ = fixture()
        original[0].status, original[0].mapped_smiles = "mapped", "sentinel"
        _, selection = select_mappings({"LocalMapper": original})
        original[0].status = "mapped"
        partials = apply_rules(original, selection, metabolites)
        self.assertTrue(partials.empty)
        self.assertEqual(original[0].mapped_smiles, "sentinel")
        self.assertEqual(selection.iloc[0].MappingCoverage, "Complete")

    def test_workbook_roundtrip_rejection_and_tampering(self):
        source, metabolites, original, _, _, _ = fixture()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_file = root / "input.xlsx"
            with pd.ExcelWriter(input_file, engine="openpyxl") as writer:
                source.to_excel(writer, sheet_name="Reactions", index=False)
                metabolites.to_excel(writer, sheet_name="Metabolites", index=False)
            finish_analysis("Test", dict(input_file=input_file, metabolite_sheet="Metabolites"),
                            metabolites, {"LocalMapper": original}, [], 1,
                            ["LocalMapper"], root, stage="map")
            output = root / "Test_2_atom_mapping.xlsx"
            saved, _, _ = load_final_mappings(output, {})
            self.assertEqual(saved[0].status, "partial")
            with pd.ExcelFile(output, engine="openpyxl") as excel:
                final = pd.read_excel(excel, sheet_name="Final_Mapping").fillna("")
                from carbon_flow_tools.excel_io import read_partial_mapping_records
                self.assertNotIn("Partial_Mapping", excel.sheet_names)
                self.assertEqual(excel.book["Mapping_Info"].sheet_state, "hidden")
                partials = read_partial_mapping_records(excel.book)
            self.assertEqual([row["TransferCoefficient"] for row in saved[0].partial_maps], [55.3, 3.0])
            self.assertEqual(len(final), len(source))
            self.assertFalse(final.AIGenerated.any())
            partials.loc[0, "TransferCoefficient"] = 55
            with self.assertRaises(ValueError):
                load_partials(deepcopy(original), final, partials, metabolites)
            final.loc[0, "MappingCoverage"] = "Unresolved"
            disabled = deepcopy(original)
            load_partials(disabled, final, partials, metabolites)
            self.assertFalse(disabled[0].partial_maps)
            final.loc[0, "MappingCoverage"] = "Partial"
            with self.assertRaises(ValueError):
                load_partials(deepcopy(original), final, pd.DataFrame(), metabolites)


if __name__ == "__main__":
    unittest.main()
