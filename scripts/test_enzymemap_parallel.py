"""Check exact prediction reuse, shared lookup jobs and CPU/memory limits."""

import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from carbon_flow_tools.config import available_cpus
from carbon_flow_tools.mapper_runner import group_mapper_jobs, run_mapper
from carbon_flow_tools.mappers.enzymemap_adapter import rule_capacity, map_reactions
from carbon_flow_tools.reactions import PreparedReaction
from carbon_flow_tools.task_scheduler import run_tasks


ENZYME = "carbon_flow_tools.mappers.enzymemap_adapter"


def reaction(identifier, smiles):
    return PreparedReaction(identifier, "test", 1, "A => B", smiles, [], [])


def initialize_fixture(value):
    global fixture_value
    fixture_value = value


def initialized_fixture(index):
    return fixture_value + index


class EnzymeMapParallelTests(unittest.TestCase):
    def test_identical_settings_share_one_job_and_keep_original_objects(self):
        settings = dict(module=ENZYME, python=sys.executable,
                        options={"dataset": "test.csv"}, log_file="logs/A_EnzymeMap.log")
        first, second = reaction("a", "C>>C"), reaction("b", "CC>>CC")
        other = dict(settings, module="other.adapter")
        jobs = group_mapper_jobs([("A", "EnzymeMap", [first], settings),
            ("A", "Other", [first], other), ("B", "EnzymeMap", [second], settings)])
        self.assertEqual(len(jobs), 2)
        self.assertIs(jobs[0][1][0], first)
        self.assertIs(jobs[0][1][1], second)
        self.assertEqual(jobs[0][3]["log_file"], "logs/EnzymeMap.log")
        different = dict(settings, options={"dataset": "different.csv"})
        self.assertEqual(len(group_mapper_jobs([("A", "EnzymeMap", [], settings),
            ("B", "EnzymeMap", [], different)])), 2)

    def test_exact_duplicates_are_predicted_once_and_restored_to_each_row(self):
        inputs = [reaction("a", "C[C@H](O)F>>C[C@H](O)F"),
                  reaction("b", "C[C@H](O)F>>C[C@H](O)F"),
                  reaction("c", "C[C@@H](O)F>>C[C@@H](O)F")]
        settings = dict(module=ENZYME, python=sys.executable)
        calls = []
        def execute(command, **kwargs):
            payload = json.loads(Path(command[3]).read_text())
            calls.append(payload)
            Path(command[4]).write_text(json.dumps([
                dict(mapped_rxn="[CH4:1]>>[CH4:1]", template=str(i), confidence=.9)
                for i in range(len(payload["reactions"]))]))
            self.assertEqual(kwargs["env"]["OMP_NUM_THREADS"], "1")
            return SimpleNamespace(returncode=0)
        with patch("carbon_flow_tools.mapper_runner.subprocess.run", side_effect=execute):
            run_mapper("All models", inputs, "EnzymeMap", settings, 8, None)
        self.assertEqual(len(calls), 1)
        self.assertEqual(len(calls[0]["reactions"]), 2)
        self.assertEqual([r.template for r in inputs], ["0", "0", "1"])
        self.assertEqual([r.reaction_id for r in inputs], ["a", "b", "c"])
        self.assertTrue(all(r.status == "mapped" for r in inputs))
        self.assertEqual(calls[0]["cpu_threads"], 8)

    def test_cpu_capacity_scales_to_servers_and_respects_task_and_memory_limits(self):
        with patch("carbon_flow_tools.config.available_memory", return_value=(512,512)), \
             patch("carbon_flow_tools.config.available_cpus", return_value=128):
            self.assertEqual(rule_capacity(200, 128, None), 128)
            self.assertEqual(rule_capacity(200, 32, None), 32)
            self.assertEqual(rule_capacity(3, 128, None), 3)
            self.assertEqual(rule_capacity(200, 2, None), 2)
        with patch("carbon_flow_tools.config.available_memory", return_value=(10,10)), \
             patch("carbon_flow_tools.config.available_cpus", return_value=128):
            self.assertEqual(rule_capacity(200,128,None), 14)
            self.assertEqual(rule_capacity(200,128,2), 2)
        with patch("carbon_flow_tools.config.available_memory", side_effect=OSError):
            self.assertEqual(rule_capacity(200,128,None), 1)

    def test_affinity_and_cpu_quota_are_both_honored(self):
        def read(path, *args, **kwargs):
            if str(path) == "/proc/self/cgroup":
                return "0::/test\n"
            if str(path).endswith("/cpu.max"):
                return "200000 100000" if str(path) == "/sys/fs/cgroup/cpu.max" else "max 100000"
            raise FileNotFoundError(str(path))
        with patch("carbon_flow_tools.config.os.sched_getaffinity", return_value=set(range(32))), \
             patch.object(Path, "read_text", read):
            self.assertEqual(available_cpus(), 2)
        with patch("carbon_flow_tools.config.os.sched_getaffinity", return_value={1}), \
             patch.object(Path, "read_text", read):
            self.assertEqual(available_cpus(), 1)

    def test_scheduler_initializes_once_and_reports_each_completion(self):
        for workers in (1,2):
            completed, waiting = [], []
            with patch("carbon_flow_tools.task_scheduler.available_memory", return_value=(32,32)):
                values = run_tasks(initialized_fixture, [(i,) for i in range(4)], workers, .5,
                    initializer=initialize_fixture, initargs=(10,),
                    on_result=lambda: completed.append(True),
                    on_wait=lambda: waiting.append(len(completed)))
            self.assertEqual(values, [10,11,12,13])
            self.assertEqual(len(completed), 4)
            if workers == 2:
                self.assertTrue(waiting)
                self.assertEqual(waiting, sorted(waiting))
                self.assertEqual(waiting[0], 0)
            else:
                self.assertEqual(waiting, [])

    def test_database_and_identity_hits_bypass_rule_pool(self):
        progress = SimpleNamespace(label="EnzymeMap", begin=lambda *a, **kw: None,
                                   update=lambda *a, **kw: None)
        with tempfile.TemporaryDirectory() as folder:
            with patch("carbon_flow_tools.task_scheduler.run_tasks") as pool:
                values = map_reactions(dict(matches={}, predictor=object()), ["C>>C"],
                    32, None, str(Path(folder)/"mapper.log"), progress)
            self.assertTrue(values[0].get("mapped_rxn"))
            pool.assert_not_called()

    def test_lookup_failure_does_not_discard_other_reactions(self):
        progress = SimpleNamespace(label="EnzymeMap", begin=lambda *a, **kw: None,
                                   update=lambda *a, **kw: None)
        with tempfile.TemporaryDirectory() as folder:
            with patch("carbon_flow_tools.mappers.enzymemap_adapter.map_batch",
                       side_effect=[ValueError("fixture lookup failure"), [{"mapped_rxn": "C>>C"}]]):
                values = map_reactions(dict(matches={}, predictor=object()), ["bad", "C>>C"],
                    32, None, str(Path(folder)/"mapper.log"), progress)
        self.assertEqual(values, [{"error": "fixture lookup failure"}, {"mapped_rxn": "C>>C"}])

    def test_lookup_and_parallel_rules_have_independent_completed_counters(self):
        events, phases = [], []
        progress = SimpleNamespace(label="Shared / EnzymeMap")
        def begin(label, total, estimate=True):
            phases.append((label, total, estimate))
            progress.label = label
        progress.begin = begin
        progress.update = lambda count, **kw: events.append((progress.label, count))
        def parallel_rules(function, tasks, workers, memory, **kwargs):
            self.assertEqual(list(tasks), [("slow",), ("fast",)])
            kwargs["on_wait"]()
            kwargs["on_result"]()  # The second task may finish first.
            kwargs["on_wait"]()
            kwargs["on_result"]()
            return [{"mapped_rxn": "slow result"}, {"error": "fast failure"}]
        with tempfile.TemporaryDirectory() as folder, \
             patch(ENZYME + ".map_batch", side_effect=[[{"mapped_rxn": "database result"}], [{"error": "missing"}], [{"error": "missing"}]]), \
             patch(ENZYME + ".rule_capacity", return_value=2), \
             patch("carbon_flow_tools.task_scheduler.run_tasks", side_effect=parallel_rules):
            values = map_reactions(dict(matches={}, predictor=object(), options={}),
                ["hit", "slow", "fast"], 2, None, str(Path(folder)/"mapper.log"), progress)
        self.assertEqual([phase[1:] for phase in phases], [(3, False), (2, True)])
        self.assertTrue(phases[0][0].endswith("database matching"))
        self.assertTrue(phases[1][0].endswith("rule prediction"))
        self.assertEqual([count for label, count in events if label == phases[1][0]], [0, 1, 1, 2])
        self.assertEqual(values, [{"mapped_rxn": "database result"},
                                 {"mapped_rxn": "slow result"}, {"error": "fast failure"}])

    def test_single_scan_preserves_database_quality_and_tie_breaking(self):
        from carbon_flow_tools.mappers.enzymemap_adapter import load_model, key
        import csv
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder)/"database.csv"
            source.write_text("unmapped,mapped,quality\n"
                "C>>CO,[CH4:1]>>[CH3:1][OH:2],0.1\n"
                "C>>CO,[CH4:7]>>[CH3:7][OH:8],0.9\n"
                "C>>CO,[CH4:9]>>[CH3:9][OH:10],0.9\n")
            with patch("carbon_flow_tools.mappers.enzymemap_adapter.csv.DictReader", wraps=csv.DictReader) as reader:
                model = load_model(dict(dataset=str(source), reactions=["C>>CO", "CO>>C"]))
            self.assertEqual(reader.call_count, 1)
            self.assertEqual(model["matches"][key("C>>CO")]["mapped"], "[CH4:7]>>[CH3:7][OH:8]")
            self.assertIn("line 3; forward", model["matches"][key("C>>CO")]["template"])
            self.assertIn("line 3; reverse", model["matches"][key("CO>>C")]["template"])

    def test_failed_child_initialization_is_reported_without_respawn_loop(self):
        from carbon_flow_tools.mappers.enzymemap_adapter import _initialize_rules, _predict_rule
        with tempfile.TemporaryDirectory() as folder:
            _initialize_rules(dict(rules=str(Path(folder)/"missing.pkl")), str(Path(folder)/"mapper.log"))
            with self.assertRaisesRegex(RuntimeError, "initialization failed"):
                _predict_rule("CCO>>CC=O")


if __name__ == "__main__":
    unittest.main()
