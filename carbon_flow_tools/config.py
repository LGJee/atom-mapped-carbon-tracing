import math
import os
from pathlib import Path


BATCH_SIZE = 1
MAX_REACTION_ATOMS = 500


def configure(batch_size, max_reaction_atoms):
    """Apply validated mapping limits supplied by the entry point."""
    if any(
        type(value) is not int or value < 1
        for value in (batch_size, max_reaction_atoms)
    ):
        raise ValueError(
            "Batch size and maximum reaction atoms must be positive integers."
        )
    global BATCH_SIZE, MAX_REACTION_ATOMS
    BATCH_SIZE, MAX_REACTION_ATOMS = batch_size, max_reaction_atoms
    print(
        f"Configuration: batch size {BATCH_SIZE}; maximum {MAX_REACTION_ATOMS} atoms per reaction (both sides combined).",
        flush=True,
    )


def parallel_capacity(
    task_count, max_processes=None, memory_per_process_gb=6, task_kind="prediction"
):
    """Estimate concurrency from CPU and available memory; budgets are not peak guarantees."""
    if max_processes is not None and (
        type(max_processes) is not int or max_processes < 1
    ):
        raise ValueError("Maximum processes must be a positive integer or None.")
    if memory_per_process_gb is not None and (
        type(memory_per_process_gb) not in (int, float)
        or not math.isfinite(memory_per_process_gb)
        or memory_per_process_gb <= 0
    ):
        raise ValueError("Memory budget must be a finite positive number or None.")
    if task_count == 0:
        print("No pending predictions; no mapper processes needed.", flush=True)
        return 1, 1
    cpus = (
        len(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else os.cpu_count() or 1
    )
    available, total = None, None
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemTotal:"):
                total = int(line.split()[1]) / (1024**2)
            elif line.startswith("MemAvailable:"):
                available = int(line.split()[1]) / (1024**2)
                break
    except (OSError, ValueError):
        pass
    automatic = memory_per_process_gb is None
    if automatic:
        # Budget one quarter of total memory per strain; this is not a peak bound.
        memory_per_process_gb = max(2, (total or available or 0) / 4)

    # Reserve 2 GiB or 20%; use serial execution when memory cannot be detected.
    memory_limit = (
        max(1, int((available - max(2, available * 0.2)) / memory_per_process_gb))
        if available
        else 1
    )
    workers = max(1, min(task_count, cpus, memory_limit, max_processes or task_count))
    threads = max(1, cpus // workers)
    memory_text = f"{available:.1f} GiB" if available is not None else "unavailable"
    budget_kind = "estimated" if automatic else "configured"
    print(
        f"Resources: {cpus} CPU threads, {memory_text} available memory; "
        f"{task_count} {task_kind} tasks, up to {workers} concurrent processes; "
        f"{budget_kind} memory budget {memory_per_process_gb:.1f} GiB/process.",
        flush=True,
    )
    return workers, threads
