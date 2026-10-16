"""Verify bounded output batches without changing the scientific propagation rules."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from carbon_flow_tools.carbon_trace import COLUMNS, propagate_rounds, map_one_reaction
from carbon_flow_tools.flow_output import FlowOutput
from carbon_flow_tools.pipeline import summarize_sources
from scripts.test_carbon_tracing import reaction


class BatchedOutputTests(unittest.TestCase):
    def test_staged_rounds_match_original_records_ids_and_convergence(self):
        structures = {"A[c]": "CC", "X[c]": "C", "P[c]": "CCC"}
        forward, parts = reaction("A[c] + X[c]", "P[c]", structures,
                                 "[CH3:1][CH3:2].[CH4:3]>>[CH3:1][CH2:2][CH3:3]")
        reverse, _ = reaction("P[c]", "A[c] + X[c]", structures,
                              "[CH3:1][CH2:2][CH3:3]>>[CH3:1][CH3:2].[CH4:3]")
        forward.reaction_id, reverse.reaction_id = "r1", "r2"
        for item in (forward, reverse):
            item.status = "mapped"
        metabolites = pd.DataFrame([
            dict(MetaboliteID=k, MetaboliteName=k, SMILES_Model=v)
            for k, v in structures.items()
        ])
        for sources in (["A[c]"], ["A[c]", "X[c]"], ["X[c]", "A[c]"]):
            for recycling in (False, True):
                with self.subTest(sources=sources, recycling=recycling), tempfile.TemporaryDirectory() as directory:
                    expected = list(propagate_rounds("Test", [forward, reverse], metabolites,
                                                     3, sources, allow_inorganic_carbon_recycling=recycling))
                    summaries = []
                    store = FlowOutput(Path(directory) / "flow.xlsx", 50000000, 2, summaries, None)
                    try:
                        actual = propagate_rounds("Test", [forward, reverse], metabolites, 3, sources,
                                                  allow_inorganic_carbon_recycling=recycling, row_store=store)
                        observed = []
                        for number, table, added, converged in actual:
                            self.assertIsNone(table)
                            data = pd.DataFrame(store.db.execute(
                                "SELECT " + store.fields + " FROM records WHERE round_no=? "
                                "ORDER BY phase, Step, ReactionID, SourceIDs, ProductionID, seq", (number,)),
                                columns=COLUMNS)
                            observed.append((number, data, added, converged))
                            self.assertLess(len(store.buffer), store.batch_rows)
                            store.complete_round()
                        self.assertEqual(len(expected), len(observed))
                        for old, new in zip(expected, observed):
                            self.assertEqual((old[0], old[2], old[3]), (new[0], new[2], new[3]))
                            pd.testing.assert_frame_equal(old[1], new[1])
                        self.assertFalse(any('_p' in name for name in pd.ExcelFile(store.output, engine="openpyxl").sheet_names))
                    finally:
                        store.close()
                    self.assertFalse(list(Path(directory).glob('.carbon-flow-*')))

        emitted = []
        produced, rows, _ = map_one_reaction(
            forward, parts, {"A[c]": {0: ((0, 0, "C1"),)}}, 1, {}, 1, structures,
            record_sink=emitted.append)
        self.assertTrue(emitted)
        self.assertEqual(produced, [])
        self.assertEqual(rows, [])

    def test_staged_inorganic_endpoints_and_other_substrates_match(self):
        from scripts.test_carbon_tracing import check_inorganic_carbon_termination

        def staged(*args, **kwargs):
            expected = list(propagate_rounds(*args, **kwargs))
            with tempfile.TemporaryDirectory() as directory:
                store = FlowOutput(Path(directory) / "flow.xlsx", 50000000, 2, [], None)
                try:
                    observed = []
                    for number, _, added, converged in propagate_rounds(*args, row_store=store, **kwargs):
                        table = pd.DataFrame(store.db.execute(
                            "SELECT " + store.fields + " FROM records WHERE round_no=? "
                            "ORDER BY phase, Step, ReactionID, SourceIDs, ProductionID, seq", (number,)),
                            columns=COLUMNS)
                        observed.append((number, table, added, converged))
                    self.assertEqual(len(expected), len(observed))
                    for old, new in zip(expected, observed):
                        self.assertEqual((old[0], old[2], old[3]), (new[0], new[2], new[3]))
                        pd.testing.assert_frame_equal(old[1], new[1])
                        yield new
                finally:
                    store.close()

        with patch("scripts.test_carbon_tracing.propagate_rounds", staged):
            check_inorganic_carbon_termination()

    def test_streaming_summary_keeps_mixed_state_definition(self):
        sources = [("A", (), {"C1"}), ("X", (), {"C7"})]
        records = [("A[c]", "C1"), ("A[c]", "C7"), ("P[c]", "C1,C7")]
        ids = {"A[c]": {"a"}, "P[c]": {"p"}}
        expected = summarize_sources(pd.DataFrame(records, columns=["Production", "CarbonLabels"]),
                                     sources, ids, 1, {})
        actual = summarize_sources(iter(records), sources, ids, 1, {})
        self.assertEqual(actual, expected)
        self.assertEqual(actual["MixedSources_LabeledCarbonMetabolites"], 1)

    def test_invalid_batch_size_fails_before_creating_output(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "flow.xlsx"
            for size in (0, -1, True, 2.5, None):
                with self.assertRaises(ValueError):
                    FlowOutput(path, 50000000, size, [], None)
            self.assertFalse(path.exists())
            self.assertFalse(list(Path(directory).iterdir()))


if __name__ == '__main__':
    unittest.main()
