"""Check bounded admission, deterministic results and isolated process limits."""
from pathlib import Path
import os
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

from carbon_flow_tools.config import MemoryPressureError, parallel_capacity
from carbon_flow_tools.task_scheduler import run_tasks


def process_probe(index, folder):
    Path(folder, str(index)).write_text(str(os.getpid())+":"+os.environ.get("OPENBLAS_NUM_THREADS", "missing"))
    time.sleep(0.03)
    return index * 2


def raise_error():
    raise ValueError("scheduler fixture failure")


class SchedulerTests(unittest.TestCase):
    def test_largest_first_and_results_remain_in_original_order(self):
        started = []
        def record(value):
            started.append(value)
            return value * 2
        with patch("carbon_flow_tools.task_scheduler.available_memory", return_value=(32, 32)):
            result = run_tasks(record, [(1,), (2,), (3,)], 1, 1, weights=[1, 3, 2])
        self.assertEqual(started, [2, 3, 1])
        self.assertEqual(result, [2, 4, 6])

    def test_thread_jobs_overlap_and_do_not_exceed_limit(self):
        lock = threading.Lock()
        active, peak = [0], [0]
        def work(value):
            with lock:
                active[0] += 1
                peak[0] = max(peak[0], active[0])
            time.sleep(0.05)
            with lock:
                active[0] -= 1
            return value
        with patch("carbon_flow_tools.task_scheduler.available_memory", return_value=(32,32)):
            result = run_tasks(work, [(n,) for n in range(5)], 2, 1, processes=False)
        self.assertEqual(result, list(range(5)))
        self.assertEqual(peak[0], 2)

    def test_active_budgets_are_reserved_before_usage_is_visible(self):
        lock = threading.Lock()
        active, peak = [0], [0]
        def work():
            with lock:
                active[0] += 1
                peak[0] = max(peak[0], active[0])
            time.sleep(0.05)
            with lock:
                active[0] -= 1
        with patch("carbon_flow_tools.task_scheduler.available_memory", return_value=(10,10)):
            run_tasks(work, [(), (), ()], 3, 3, processes=False)
        self.assertEqual(peak[0], 2)

    def test_failure_allows_running_jobs_to_finish_but_stops_new_jobs(self):
        finished = []
        def work(index):
            if index == 1:
                time.sleep(0.03)
                raise ValueError("failed task")
            time.sleep(0.10)
            finished.append(index)
        with patch("carbon_flow_tools.task_scheduler.available_memory", return_value=(32,32)):
            with self.assertRaisesRegex(ValueError, "failed task"):
                run_tasks(work, [(0,), (1,), (2,)], 2, 1, processes=False)
        self.assertEqual(finished, [0])

    def test_low_memory_holds_pending_jobs_until_memory_is_released(self):
        running = threading.Event()
        started = []
        def memory():
            return (2.5 if running.is_set() else 20, 20)
        def work(index):
            started.append(index)
            running.set()
            time.sleep(0.05)
            running.clear()
            return index
        with patch("carbon_flow_tools.task_scheduler.available_memory", side_effect=memory):
            self.assertEqual(run_tasks(work, [(0,), (1,)], 2, 3, processes=False), [0, 1])
        self.assertEqual(started, [0, 1])

    def test_unstartable_idle_queue_fails_without_hanging(self):
        with patch("carbon_flow_tools.task_scheduler.available_memory", return_value=(1,32)):
            with self.assertRaisesRegex(MemoryPressureError, "Not enough free memory"):
                run_tasks(lambda: None, [()], 1, 6)

    def test_missing_memory_readings_use_serial_admission(self):
        running, peak = [0], [0]
        def work():
            running[0] += 1
            peak[0] = max(peak[0], running[0])
            time.sleep(0.03)
            running[0] -= 1
        with patch("carbon_flow_tools.task_scheduler.available_memory", side_effect=OSError):
            run_tasks(work, [(), (), ()], 3, 1, processes=False)
        self.assertEqual(peak[0], 1)

    def test_spawned_jobs_receive_thread_limits_and_parent_is_restored(self):
        before = os.environ.get("OPENBLAS_NUM_THREADS")
        with tempfile.TemporaryDirectory() as folder:
            with patch("carbon_flow_tools.task_scheduler.available_memory", return_value=(32,32)):
                self.assertEqual(run_tasks(process_probe, [(n,folder) for n in range(3)],
                    2, 1, cpu_threads=2), [0, 2, 4])
            for path in Path(folder).iterdir():
                pid, threads = path.read_text().split(":")
                self.assertNotEqual(int(pid), os.getpid())
                self.assertEqual(threads, "2")
        self.assertEqual(os.environ.get("OPENBLAS_NUM_THREADS"), before)

    def test_worker_exceptions_keep_their_original_type_and_message(self):
        for processes in (False, True):
            with patch("carbon_flow_tools.task_scheduler.available_memory", return_value=(32,32)):
                with self.assertRaisesRegex(ValueError, "scheduler fixture failure"):
                    run_tasks(raise_error, [(), ()], 2, 1, processes=processes)

    def test_capacity_uses_cgroup_aware_memory(self):
        with patch("carbon_flow_tools.config.os.sched_getaffinity", return_value=set(range(32))), \
                patch("carbon_flow_tools.config.available_memory", return_value=(8,8)):
            self.assertEqual(parallel_capacity(8, memory_per_process_gb=2), (3,10))

    def test_empty_queue_and_invalid_parameters(self):
        self.assertEqual(run_tasks(None, [], 1, 1), [])
        for kwargs in [dict(workers=0), dict(memory_gb=-1), dict(weights=[1,2]), dict(cpu_threads=0)]:
            params=dict(workers=1,memory_gb=1); params.update(kwargs)
            with self.assertRaises(ValueError):
                run_tasks(lambda: None, [()], **params)


if __name__ == "__main__":
    unittest.main()
