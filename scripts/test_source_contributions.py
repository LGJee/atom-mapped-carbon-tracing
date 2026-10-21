"""Hand-solvable checks for the optional single-carbon origin balances."""

import unittest

import numpy as np
import pandas as pd

from carbon_flow_tools.reactions import prepare_reactions
from carbon_flow_tools.source_contributions import solve_source_contributions


def solve(structures, specifications, sources):
    metabolites = pd.DataFrame([dict(MetaboliteID=n, MetaboliteName=n, SMILES_Manual=s)
                                for n, s in structures.items()])
    table = pd.DataFrame([dict(ReactionID="r" + str(i), ReactionName="", Flux=v, Equation=e)
                          for i, (e, v, m) in enumerate(specifications)])
    reactions = prepare_reactions(table, metabolites)
    for r, (_, _, mapped) in zip(reactions, specifications):
        r.status = "mapped" if mapped else "missing"
        r.mapped_smiles = mapped
    return solve_source_contributions(reactions, metabolites, sources)


class SourceContributionTests(unittest.TestCase):
    def test_fractional_rule_coefficients_preserve_cofactor_origins(self):
        from scripts.test_carbon_mapping_rules import fixture
        from carbon_flow_tools.carbon_mapping_rules import apply_rules
        from carbon_flow_tools.consensus import select_mappings
        from carbon_flow_tools.reactions import prepare_reactions
        source, metabolites, _, _, _, _ = fixture()
        source = pd.DataFrame([["in", "Import", 1, " => ATP[c]"],
            ["transfer", "ATP turnover with unknown pool", 1, "0.5 ATP[c] + protein[c] => 0.5 ADP[c] + biomass[c]"],
            ["atp_out", "Export", 0.5, "ATP[c] => "],
            ["adp_out", "Export", 0.5, "ADP[c] => "]], columns=source.columns)
        reactions, choices = select_mappings({"LocalMapper": prepare_reactions(source, metabolites)})
        apply_rules(reactions, choices, metabolites)
        tables, _ = solve_source_contributions(reactions, metabolites, ["ATP[c]"])
        row = tables["MetaboliteContributions"].set_index("MetaboliteName").loc["ADP[c]"]
        self.assertAlmostEqual(row["ATP[c]"], 1)
        self.assertAlmostEqual(row.Unknown, 0)
        self.assertAlmostEqual(row.Inflow, .5)

    def test_partial_mapper_accounts_only_unmapped_product_positions_as_unknown(self):
        from carbon_flow_tools.consensus import carbon_candidate
        structures = pd.DataFrame([[n, n, "CCO"] for n in ("G", "P")],
            columns=["MetaboliteID", "MetaboliteName", "SMILES_Manual"])
        source = pd.DataFrame([["in", "", 1, " => G"], ["r", "", 1, "G => P"], ["out", "", 1, "P => "]],
            columns=["ReactionID", "ReactionName", "Flux", "Equation"])
        reactions = prepare_reactions(source, structures)
        reactions[1].mapped_smiles = "[CH3:1][CH2:2]O>>[CH3:1][CH2:8]O"
        reactions[1] = carbon_candidate(reactions[1])
        tables, _ = solve_source_contributions(reactions, structures, ["G"])
        row = tables["MetaboliteContributions"].set_index("MetaboliteName").loc["P"]
        self.assertAlmostEqual(row.G, .5)
        self.assertAlmostEqual(row.Unknown, .5)

    def test_emu_inorganic_carbon_policy_is_independent(self):
        structures = pd.DataFrame([[n, n, text] for n, text in [("G", "C"), ("CO2", "O=C=O"), ("P", "CO")]],
            columns=["MetaboliteID", "MetaboliteName", "SMILES_Manual"])
        source = pd.DataFrame([["in", "", 1, " => G"], ["r1", "", 1, "G => CO2"],
            ["r2", "", 1, "CO2 => P"], ["out", "", 1, "P => "]],
            columns=["ReactionID", "ReactionName", "Flux", "Equation"])
        reactions = prepare_reactions(source, structures)
        reactions[1].status, reactions[1].mapped_smiles = "mapped", "[CH4:1]>>O=[C:1]=O"
        reactions[2].status, reactions[2].mapped_smiles = "mapped", "O=[C:1]=O>>[CH3:1]O"
        for enabled in (True, False):
            tables, _ = solve_source_contributions(reactions, structures, ["G"], allow_inorganic_carbon_recycling=enabled)
            row = tables["MetaboliteContributions"].set_index("MetaboliteName").loc["P"]
            self.assertAlmostEqual(row.G, float(enabled))
            self.assertAlmostEqual(row.Unknown, float(not enabled))

    def test_flux_weighted_merge_and_partial_unknown(self):
        for mapping, expected, unknown in [("[CH4:1]>>[CH4:1]", .6, 0), ("", 0, .6)]:
            tables, info = solve({"G": "C", "X": "C", "P": "C"}, [
                (" -> G", 6, ""), (" -> X", 4, ""),
                ("G -> P", 6, mapping), ("X -> P", 4, "[CH4:1]>>[CH4:1]"),
                ("P -> ", 10, "")], ["G", "X"])
            p = tables["MetaboliteContributions"].set_index("MetaboliteID").loc["P"]
            self.assertAlmostEqual(p.G, expected)
            self.assertAlmostEqual(p.X, .4)
            self.assertAlmostEqual(p.Unknown, unknown)
            self.assertLess(info["max_scaled_residual"], 1e-10)

    def test_position_specific_condensation_and_cleavage(self):
        tables, _ = solve({"G": "C", "X": "CO", "P": "CCO"}, [
            (" -> G", 2, ""), (" -> X", 2, ""),
            ("G + X -> P", 3, "[CH4:1].[CH3:2][OH:3]>>[CH3:1][CH2:2][OH:3]"),
            ("P -> G + X", 1, "[CH3:1][CH2:2][OH:3]>>[CH4:1].[CH3:2][OH:3]"),
            ("P -> ", 2, "")], ["G", "X"])
        atoms = tables["AtomContributions"].query("MetaboliteID == 'P'")
        np.testing.assert_allclose(atoms.G, [1, 0])
        np.testing.assert_allclose(atoms.X, [0, 1])
        p = tables["MetaboliteContributions"].set_index("MetaboliteID").loc["P"]
        self.assertAlmostEqual(p.G, .5)
        self.assertAlmostEqual(p.X, .5)

    def test_circulation_does_not_create_extra_source(self):
        tables, _ = solve(dict.fromkeys(["G", "X", "A", "B"], "C"), [
            (" -> G", 6, ""), (" -> X", 4, ""),
            ("G -> A", 6, "[CH4:1]>>[CH4:1]"),
            ("X -> A", 4, "[CH4:1]>>[CH4:1]"),
            ("A -> B", 15, "[CH4:1]>>[CH4:1]"),
            ("B -> A", 5, "[CH4:1]>>[CH4:1]"), ("B -> ", 10, "")], ["G", "X"])
        p = tables["MetaboliteContributions"].set_index("MetaboliteID")
        np.testing.assert_allclose(p.loc[["A", "B"], "G"], [.6, .6])

    def test_unanchored_cycle_and_downstream(self):
        tables, info = solve(dict.fromkeys(["G", "A", "B"], "C"), [
            (" -> G", 1, ""), ("G -> ", 1, ""),
            ("A -> B", 3, "[CH4:1]>>[CH4:1]"),
            ("B -> A", 3, "[CH4:1]>>[CH4:1]")], ["G"])
        p = tables["MetaboliteContributions"].set_index("MetaboliteID")
        self.assertEqual(info["unanchored_positions"], 2)
        self.assertEqual(p.loc["A", "Status"], "unanchored_cycle")
        self.assertAlmostEqual(p.loc["A", "Unknown"], 1)

    def test_negative_flux_and_repeated_molecules(self):
        tables, _ = solve({"G": "C", "P": "C"}, [
            ("G -> ", -2, ""),
            ("2 P -> 2 G", -1, "[CH4:1].[CH4:2]>>[CH4:1].[CH4:2]"),
            ("P -> ", 2, "")], ["G"])
        p = tables["MetaboliteContributions"].set_index("MetaboliteID")
        self.assertAlmostEqual(p.loc["P", "G"], 1)
        self.assertEqual(p.loc["P", "Status"], "complete")

    def test_external_other_and_unbalanced_pool(self):
        for output_rate in (1, 2):
            tables, _ = solve(dict.fromkeys(["G", "Z", "P"], "C"), [
                (" -> G", 1, ""), ("G -> ", 1, ""),
                (" -> Z", 1, ""), ("Z -> P", 1, "[CH4:1]>>[CH4:1]"),
                ("P -> ", output_rate, "")], ["G"])
            p = tables["MetaboliteContributions"].set_index("MetaboliteID")
            if output_rate == 1:
                self.assertAlmostEqual(p.loc["P", "External:Z"], 1)
            else:
                self.assertEqual(p.loc["P", "Status"], "invalid_flux_balance")
                self.assertAlmostEqual(p.loc["P", "Unknown"], 1)

    def test_same_origin_symmetry_preserves_known_source(self):
        tables, info = solve({"G": "CC", "P": "CC"}, [
            (" -> G", 1, ""), ("G -> P", 1, "[CH3:1][CH3:2]>>[CH3:1][CH3:2]"),
            ("P -> ", 1, "")], ["G"])
        p = tables["MetaboliteContributions"].set_index("MetaboliteID")
        self.assertEqual(info["ambiguous_carbon_positions"], 4)
        self.assertAlmostEqual(p.loc["P", "G"], 1)
        self.assertAlmostEqual(p.loc["P", "Unknown"], 0)

    def test_same_origin_survives_symmetric_intermediate_and_cleavage(self):
        tables, info = solve({"G": "CCCCO", "P": "CC(C)CO", "Q": "C", "R": "CCCO"}, [
            (" -> G", 1, ""),
            ("G -> P", 1, "[CH3:1][CH2:2][CH2:3][CH2:4][OH:5]>>[CH3:1][CH:2]([CH3:3])[CH2:4][OH:5]"),
            ("P -> Q + R", 1, "[CH3:1][CH:2]([CH3:3])[CH2:4][OH:5]>>[CH4:1].[CH3:3][CH2:2][CH2:4][OH:5]"),
            ("Q -> ", 1, ""), ("R -> ", 1, "")], ["G"])
        p = tables["MetaboliteContributions"].set_index("MetaboliteID")
        np.testing.assert_allclose(p.loc[["P", "Q", "R"], "G"], 1)
        np.testing.assert_allclose(p.loc[["P", "Q", "R"], "Unknown"], 0)
        self.assertGreater(info["symmetry_mapped_transfers"], 0)

    def test_symmetric_positions_follow_the_selected_mapping_without_scrambling(self):
        for cleavage, expected_g in (("[CH3:1][CH3:2]>>[CH4:1].[CH4:2]", [1, 0]),
                                     ("[CH3:1][CH3:2]>>[CH4:2].[CH4:1]", [0, 1])):
            tables, _ = solve({"G": "C", "X": "C", "P": "CC", "Q": "C", "R": "C"}, [
                (" -> G", 1, ""), (" -> X", 1, ""),
                ("G + X -> P", 1, "[CH4:1].[CH4:2]>>[CH3:1][CH3:2]"),
                ("P -> Q + R", 1, cleavage),
                ("Q -> ", 1, ""), ("R -> ", 1, "")], ["G", "X"])
            p = tables["MetaboliteContributions"].set_index("MetaboliteID")
            self.assertAlmostEqual(p.loc["P", "G"], .5)
            self.assertAlmostEqual(p.loc["P", "X"], .5)
            np.testing.assert_allclose(p.loc[["Q", "R"], "G"], expected_g)
            np.testing.assert_allclose(p.loc[["P", "Q", "R"], "Unknown"], 0)

    def test_symmetric_transfer_preserves_partial_mixed_origins(self):
        tables, _ = solve(dict.fromkeys(["G", "X", "Z", "P", "Q"], "CC"), [
            (" -> G", 4, ""), (" -> X", 2, ""), (" -> Z", 4, ""),
            ("G -> P", 4, "[CH3:1][CH3:2]>>[CH3:1][CH3:2]"),
            ("X -> P", 2, "[CH3:1][CH3:2]>>[CH3:1][CH3:2]"),
            ("Z -> P", 4, ""),
            ("P -> Q", 10, "[CH3:1][CH3:2]>>[CH3:1][CH3:2]"),
            ("Q -> ", 10, "")], ["G", "X"])
        p = tables["MetaboliteContributions"].set_index("MetaboliteID")
        for name in ("P", "Q"):
            np.testing.assert_allclose(p.loc[name, ["G", "X", "Unknown"]].astype(float), [.4, .2, .4])

    def test_two_precursor_copies_preserve_mean_origins_in_a_symmetric_product(self):
        structures = dict.fromkeys(["G", "X", "Z", "P"], "CCO")
        structures["Q"] = "OCCCCO"
        tables, _ = solve(structures, [
            (" -> G", 4, ""), (" -> X", 2, ""), (" -> Z", 4, ""),
            ("G -> P", 4, "[CH3:1][CH2:2][OH:3]>>[CH3:1][CH2:2][OH:3]"),
            ("X -> P", 2, "[CH3:1][CH2:2][OH:3]>>[CH3:1][CH2:2][OH:3]"),
            ("Z -> P", 4, ""),
            ("2 P -> Q", 5, "[CH3:1][CH2:2][OH:3].[CH3:4][CH2:5][OH:6]>>[OH:3][CH2:2][CH2:1][CH2:4][CH2:5][OH:6]"),
            ("Q -> ", 5, "")], ["G", "X"])
        p = tables["MetaboliteContributions"].set_index("MetaboliteID")
        np.testing.assert_allclose(p.loc["Q", ["G", "X", "Unknown"]].astype(float), [.4, .2, .4])

    def test_unknown_precursor_is_not_resolved_just_because_one_source_is_configured(self):
        tables, _ = solve({"G": "CC", "P": "CC", "Q": "CC"}, [
            (" -> G", 1, ""), ("G -> P", 1, ""),
            ("P -> Q", 1, "[CH3:1][CH3:2]>>[CH3:1][CH3:2]"),
            ("Q -> ", 1, "")], ["G"])
        p = tables["MetaboliteContributions"].set_index("MetaboliteID")
        self.assertAlmostEqual(p.loc["Q", "Unknown"], 1)

    def test_selected_map_turns_over_both_symmetric_positions(self):
        tables, _ = solve({"G": "C", "P": "CC", "W": "C"}, [
            (" -> G", 1, ""),
            ("G + P -> P + W", 1, "[CH4:1].[CH3:2][CH3:3]>>[CH3:1][CH3:2].[CH4:3]"),
            ("W -> ", 1, "")], ["G"])
        p = tables["MetaboliteContributions"].set_index("MetaboliteID")
        # The explicit map replaces P.C1, shifts old C1 to C2, and exports old C2.
        self.assertAlmostEqual(p.loc["P", "G"], 1)
        self.assertAlmostEqual(p.loc["P", "Unknown"], 0)
        self.assertAlmostEqual(p.loc["W", "G"], 1)

    def test_identical_inflow_composition_at_configured_source_pool(self):
        tables, _ = solve(dict.fromkeys(["G", "X"], "C"), [
            (" -> G", 1, ""), (" -> X", 1, ""),
            ("G -> X", 1, "[CH4:1]>>[CH4:1]"),
            ("X -> ", 2, "")], ["G", "X"])
        p = tables["MetaboliteContributions"].set_index("MetaboliteID")
        self.assertAlmostEqual(p.loc["X", "G"], .5)
        self.assertAlmostEqual(p.loc["X", "X"], .5)

    def test_scientific_notation_stoichiometry(self):
        from carbon_flow_tools.reactions import parse_side
        self.assertEqual(parse_side("2.5e-05 Co2+[c] + .5 A[c] + 2E+1 B[c]"),
                         [(2.5e-5, "Co2+[c]"), (.5, "A[c]"), (20., "B[c]")])

    def test_each_position_can_have_multiple_origins(self):
        tables, _ = solve({"G": "C", "X": "CO", "P": "CCO"}, [
            (" -> G", 10, ""), (" -> X", 10, ""),
            ("G + X -> P", 6, "[CH4:1].[CH3:2][OH:3]>>[CH3:1][CH2:2][OH:3]"),
            ("G + X -> P", 4, "[CH4:1].[CH3:2][OH:3]>>[CH3:2][CH2:1][OH:3]"),
            ("P -> ", 10, "")], ["G", "X"])
        atoms = tables["AtomContributions"].query("MetaboliteID == 'P'")
        np.testing.assert_allclose(atoms.G, [.6, .4])
        np.testing.assert_allclose(atoms.X, [.4, .6])

    def test_structure_formula_mismatch_is_unavailable(self):
        structures = pd.DataFrame([dict(MetaboliteID="G", MetaboliteName="G", SMILES_Manual="C"),
            dict(MetaboliteID="P", MetaboliteName="P", SMILES_Manual="CO", MolecularFormula="C2H6O")])
        source = pd.DataFrame([dict(ReactionID="in", ReactionName="", Flux=1, Equation=" -> G"),
            dict(ReactionID="out", ReactionName="", Flux=1, Equation="G -> ")])
        tables, _ = solve_source_contributions(prepare_reactions(source, structures), structures, ["G"])
        p = tables["MetaboliteContributions"].set_index("MetaboliteID").loc["P"]
        self.assertEqual(p.Status, "Structure/formula carbon count mismatch")
        self.assertTrue(pd.isna(p.G))

    def test_unavailable_positions_do_not_break_mapped_transfers(self):
        rows = pd.DataFrame([
            dict(ReactionID="in", ReactionName="", Flux=1, Equation=" -> G"),
            dict(ReactionID="r", ReactionName="", Flux=1, Equation="G -> P"),
            dict(ReactionID="out", ReactionName="", Flux=1, Equation="P -> ")])
        for unavailable in ({"G"}, {"P"}, {"G", "P"}):
            with self.subTest(unavailable=unavailable):
                structures = pd.DataFrame([dict(MetaboliteID=n, MetaboliteName=n,
                    SMILES_Manual="C", MolecularFormula="C2" if n in unavailable else "CH4")
                    for n in ("G", "P")])
                reactions = prepare_reactions(rows, structures)
                reactions[1].mapped_smiles = "[CH4:1]>>[CH4:1]"
                reactions[1].status = "mapped"
                # A separate usable source keeps this case's other inputs unchanged.
                structures = pd.concat([structures, pd.DataFrame([dict(MetaboliteID="X", MetaboliteName="X",
                    SMILES_Manual="C", MolecularFormula="CH4")])], ignore_index=True)
                tables, _ = solve_source_contributions(reactions, structures, ["X"])
                result = tables["MetaboliteContributions"].set_index("MetaboliteID")
                if "P" in unavailable:
                    self.assertTrue(pd.isna(result.loc["P", "Unknown"]))
                    self.assertEqual(result.loc["P", "Status"], "Structure/formula carbon count mismatch")
                else:
                    self.assertAlmostEqual(result.loc["P", "Unknown"], 1)
                for name in unavailable:
                    self.assertFalse(tables["AtomContributions"].MetaboliteID.eq(name).any())

    def test_invalid_mapping_is_rejected(self):
        with self.assertRaises(ValueError):
            solve({"G": "C", "P": "C"}, [(" -> G", 1, ""),
                ("G -> P", 1, "[CH4:1]>>[CH4:2]"), ("P -> ", 1, "")], ["G"])


if __name__ == "__main__":
    unittest.main()
