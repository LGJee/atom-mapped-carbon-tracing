"""Pure candidate computation; IDs, state merges and output stay in the caller."""
from collections import defaultdict, deque
from itertools import islice
from multiprocessing import current_process, get_context
import os
import time

from rdkit import RDLogger

from .atom_helpers import labeled_smiles
from .config import available_cpus, available_memory, MemoryGuard
from .console import detail

CHUNK_SIZE = 256
MIN_PARALLEL_COMBINATIONS = 10000
MAX_CHUNK_BYTES = 2 * 1024 ** 2
WORKER_MEMORY_GB = 0.5  # Estimate, not a hard bound; live memory checks remain active.
_THREAD_VARIABLES = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS")


def candidate_rows(left, right, selections, smiles, step, use_first_steps, required, unordered, compact=False):
    """Read immutable projected states in their original Cartesian-product order."""
    for selected in selections:
        if not any(choice is not None for choice in selected):
            continue
        if required and not any(choice and choice[3] for choice in selected):
            continue
        labels, terms, steps = {}, defaultdict(list), []
        for (term, _, _), choice in zip(left, selected):
            if not compact:
                terms[term].append(0 if choice is None else choice[0])
            if choice is not None:
                labels.update(choice[1])
                if use_first_steps:
                    steps.append(choice[2])
        source_ids = "" if compact else " + ".join("&".join(map(str, sorted(terms[i]) if unordered else terms[i]))
                                for i in sorted(terms))
        actual_step = max(steps, default=0) + 1 if use_first_steps else step
        for _, name, transfers in right:
            assignment = tuple(sorted((part, atom, labels[number])
                for part, atom, number in transfers if number in labels))
            if assignment:
                marked = labeled_smiles(smiles[name], assignment)
                carbon_labels = ",".join(sorted((a[2] for a in assignment), key=lambda value: int(value[1:])))
                yield actual_step, source_ids, name, carbon_labels, marked


def _initialize_candidates():
    RDLogger.DisableLog("rdApp.warning")


def _candidate_batch(left, right, selections, smiles, step, use_first_steps, required, unordered, compact=False):
    guard = MemoryGuard()
    guard.check("starting candidate batch", force=True)
    def checked():
        for index, value in enumerate(selections):
            if index % 16 == 0:
                guard.check("calculating candidate batch")
            yield value
    values = list(candidate_rows(left, right, checked(), smiles, step, use_first_steps, required, unordered, compact))
    guard.check("returning candidate batch", force=True)
    return values


class EnumerationPool:
    """Reuse a bounded process pool; consume contiguous chunks strictly in order."""
    def __init__(self, cpu_budget=1, memory_budget_gb=None, max_processes=None):
        if max_processes is not None and (type(max_processes) is not int or max_processes < 1):
            raise ValueError("Maximum enumeration processes must be a positive integer or None")
        self.workers = min(cpu_budget, available_cpus(), max_processes or cpu_budget)
        try:
            available, total = available_memory()
            spare = available - max(2, available * .2, total * .1) - 1
            if memory_budget_gb is not None:
                spare = min(spare, memory_budget_gb - 1)
            self.workers = min(self.workers, max(1, int(spare / WORKER_MEMORY_GB)))
        except (OSError, ValueError, KeyError):
            self.workers = 1
        if current_process().daemon:
            self.workers = 1  # Legacy daemonic callers retain the serial implementation.
        self.pool, self.previous_env = None, {}

    def eligible(self, combinations):
        return self.workers > 1 and combinations >= MIN_PARALLEL_COMBINATIONS

    def rows(self, left, right, selections, smiles, step, use_first_steps, required,
             unordered, memory_guard=None, progress=None, compact=False):
        guard = memory_guard or MemoryGuard()
        if self.pool is None:
            guard.check("starting enumeration workers", force=True,
                        extra_bytes=int(self.workers * WORKER_MEMORY_GB * 1024 ** 3))
            self.previous_env = {key: os.environ.get(key) for key in _THREAD_VARIABLES}
            for key in _THREAD_VARIABLES:
                os.environ[key] = "1"
            try:
                self.pool = get_context("spawn").Pool(self.workers, initializer=_initialize_candidates)
            except BaseException:
                self.close()
                raise
            detail("Enumeration: {} candidate workers; ordered merge stays in the parent."
                   .format(self.workers), flush=True)
        worker_ids = {worker.pid for worker in self.pool._pool}
        pending, completed = deque(), 0
        product_smiles = {name: smiles[name] for _, name, _ in right}

        # Account for molecular size as well as combination count when prefetching.
        bytes_per_choice = (sum(256 + len(smiles[name]) * 12 for _, name, _ in right)
                            + sum(128 + len(transfers) * 96 for _, _, transfers in left))
        chunk_size = max(1, min(CHUNK_SIZE, MAX_CHUNK_BYTES // max(1, bytes_per_choice)))
        def submit():
            guard.check("reserving enumeration chunk", force=True,
                        extra_bytes=(self.workers + 1) * chunk_size * bytes_per_choice * 3)
            chunk = tuple(islice(selections, chunk_size))
            if not chunk:
                return False
            guard.check("submitting enumeration chunk")
            pending.append((self.pool.apply_async(_candidate_batch,
                (left, right, chunk, product_smiles, step, use_first_steps, required, unordered, compact)), len(chunk)))
            return True

        try:
            for _ in range(self.workers):
                if not submit():
                    break
            while pending:
                result, count = pending.popleft()
                while not result.ready():
                    # Pool can replace a crashed worker without resolving its lost task.
                    if ({worker.pid for worker in self.pool._pool} != worker_ids
                            or any(worker.exitcode is not None for worker in self.pool._pool)):
                        raise RuntimeError("Enumeration worker exited unexpectedly; incomplete round was not published")
                    guard.check("waiting for enumeration chunk")
                    if progress is not None:
                        progress.update(completed)
                    time.sleep(.02)
                values = result.get()
                submit()  # At most one pending chunk per worker, plus this ordered result.
                for value in values:
                    yield value
                completed += count
                if progress is not None:
                    progress.update(completed)
        except BaseException:
            self.close(terminate=True)
            raise

    def close(self, terminate=False):
        if self.pool is not None:
            if terminate:
                self.pool.terminate()
            else:
                self.pool.close()
            self.pool.join()
            self.pool = None
        for key, value in self.previous_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
        self.previous_env = {}
