"""Bounded scheduling for independent mapper/model jobs, not shared atom states.

Use existing threads to launch isolated mappers and spawned processes for analyses.
Admit one job at a time after checking live cgroup-aware memory. A task's memory
budget is an estimate, not a peak guarantee; each analysis retains its MemoryGuard
and owns its mapping input, source-balance solver and output workbook.
"""

from concurrent.futures import ThreadPoolExecutor
from multiprocessing import get_context
import math
import os
import signal
import time

from .config import available_memory, task_memory_budget, MemoryPressureError
from .console import detail

_THREAD_VARIABLES = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS")


def _cancel_analysis(signum, frame):
    raise KeyboardInterrupt("Analysis process cancelled")


def _analysis_process(function, args, connection, initializer, initargs):
    # Detect abrupt worker exits instead of waiting forever for a missing result.
    signal.signal(signal.SIGTERM, _cancel_analysis)
    try:
        if initializer is not None:
            initializer(*initargs)
        connection.send((True, function(*args)))
    except BaseException as error:
        try:
            connection.send((False, error))
        except Exception:
            connection.send((False, RuntimeError("Analysis failed: " + str(error))))
    finally:
        connection.close()


class _AnalysisResult:
    def __init__(self, function, args, initializer=None, initargs=()):
        context = get_context("spawn")
        self.connection, sender = context.Pipe(duplex=False)
        self.process = context.Process(target=_analysis_process, args=(function, args, sender, initializer, initargs))
        try:
            self.process.start()
        except BaseException:
            self.connection.close()
            raise
        finally:
            sender.close()

    def ready(self):
        return self.connection.poll() or not self.process.is_alive()

    def get(self):
        try:
            if not self.connection.poll() and not self.process.is_alive():
                raise RuntimeError("Analysis process exited without a result (exit {}).".format(self.process.exitcode))
            try:
                success, value = self.connection.recv()
            except EOFError:
                raise RuntimeError("Analysis process closed without a result (exit {}).".format(self.process.exitcode))
            if not success:
                raise value
            return value
        finally:
            self.process.join()
            self.connection.close()

    def terminate(self):
        if self.process.is_alive():
            self.process.terminate()
        self.process.join()
        self.connection.close()


def run_tasks(function, tasks, workers, memory_gb, processes=True, weights=None,
              cpu_threads=1, initializer=None, initargs=(), on_result=None,
              on_wait=None):
    """Run expensive jobs first, but return results in the original input order.

    Never enqueue the whole workload into a pool: only admitted jobs are submitted.
    Process workers inherit native-library thread limits before importing NumPy;
    this process mode is called from the main thread, after mapper threads finish.
    Mapper orchestration uses thread mode and its existing subprocess environments.
    """
    tasks = list(tasks)
    if not tasks:
        return []
    if type(workers) is not int or workers < 1 or type(cpu_threads) is not int or cpu_threads < 1:
        raise ValueError("Scheduler workers and CPU threads must be positive integers")
    if memory_gb is not None and (type(memory_gb) not in (int, float)
            or not math.isfinite(memory_gb) or memory_gb <= 0):
        raise ValueError("Scheduler memory budget must be positive or None")
    weights = list(weights) if weights is not None else [0] * len(tasks)
    if len(weights) != len(tasks) or any(not math.isfinite(value) or value < 0 for value in weights):
        raise ValueError("Scheduler requires one finite nonnegative weight per task")
    pending = sorted(range(len(tasks)), key=lambda index: (-weights[index], index))
    workers = min(workers, len(tasks))
    results, active, reservations = [None] * len(tasks), {}, {}
    failure = None
    previous_env = {}
    pool = None
    try:
        if workers > 1:
            if processes:
                # Keep limits for replacement workers too; restore the parent on exit.
                previous_env = {key: os.environ.get(key) for key in _THREAD_VARIABLES}
                for key in _THREAD_VARIABLES:
                    os.environ[key] = str(cpu_threads)
            else:
                pool = ThreadPoolExecutor(max_workers=workers)
        waiting_notice = False
        while pending or active:
            for index, result in list(active.items()):
                ready = result.ready() if processes else result.done()
                if ready:
                    try:
                        results[index] = result.get() if processes else result.result()
                        if on_result is not None:
                            on_result()
                    except Exception as error:
                        if failure is None:
                            failure = error
                        pending.clear()  # Let running analyses publish/clean up normally.
                    del active[index]
                    reservations.pop(index, None)
            if failure is not None and not active:
                raise failure
            while pending and len(active) < workers:
                try:
                    available, total = available_memory()
                    budget = task_memory_budget(memory_gb, total, available)
                    reserve = max(2, available * 0.2, total * 0.1)
                    # Reserve full active estimates even before child allocations appear.
                    # This is intentionally conservative on memory-constrained hosts.
                    enough = available - sum(reservations.values()) >= reserve + budget
                except (OSError, ValueError, KeyError):
                    # Unknown memory permits only one job, matching serial fallback.
                    enough = not active
                    available = total = None
                if not enough:
                    if not active:
                        raise MemoryPressureError(
                            "Not enough free memory to start the next calculation. "
                            "Free memory, or lower the per-task memory estimate only if measured usage permits it. "
                            "No new output was written by this pending task.")
                    if not waiting_notice:
                        detail("Scheduler: waiting for memory before starting another task.", flush=True)
                        waiting_notice = True
                    break
                index = pending.pop(0)
                waiting_notice = False
                detail("Scheduler: starting task {}/{}; {} already running.".format(
                    index + 1, len(tasks), len(active)), flush=True)
                if workers == 1:
                    if initializer is not None:
                        initializer(*initargs)
                        initializer = None
                    results[index] = function(*tasks[index])
                    if on_result is not None:
                        on_result()
                elif processes:
                    active[index] = _AnalysisResult(function, tasks[index], initializer, initargs)
                    reservations[index] = budget if available is not None else 0
                else:
                    active[index] = pool.submit(function, *tasks[index])
                    reservations[index] = budget if available is not None else 0
                if available is None:
                    break
            if active:
                if on_wait is not None:
                    on_wait()
                time.sleep(0.05)
        return results
    except BaseException:
        if processes:
            for result in active.values():
                result.terminate()
        raise
    finally:
        if pool is not None:
            pool.shutdown(wait=True)
        for key, value in previous_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
