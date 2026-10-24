import math
import os
import shutil
from zipfile import ZipFile, BadZipFile
import time
from pathlib import Path

from .console import detail


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
    detail(
        f"Configuration: batch size {BATCH_SIZE}; maximum {MAX_REACTION_ATOMS} atoms per reaction (both sides combined).",
        flush=True,
    )


def task_memory_budget(memory_gb, total_gb, available_gb=None):
    """Use the same per-task estimate in capacity planning and live admission."""
    return max(2, (total_gb or available_gb or 0) / 4) if memory_gb is None else memory_gb


def available_cpus():
    """Respect CPU affinity and tighter cgroup-v2 quotas on shared servers."""
    count = (len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity")
             else os.cpu_count() or 1)
    root = Path("/sys/fs/cgroup")
    try:
        for line in Path("/proc/self/cgroup").read_text().splitlines():
            if not line.startswith("0::"):
                continue
            group = root / line.split("::", 1)[1].lstrip("/")
            for directory in (group,) + tuple(group.parents):
                if directory != root and root not in directory.parents:
                    break
                try:
                    quota, period = (directory / "cpu.max").read_text().split()
                    if quota != "max":
                        count = min(count, max(1, int(quota) // int(period)))
                except (OSError, ValueError, ZeroDivisionError):
                    continue
            break
    except OSError:
        pass
    return max(1, count)


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
        detail("No pending predictions; no mapper processes needed.", flush=True)
        return 1, 1
    cpus = available_cpus()
    available, total = None, None
    try:
        available, total = available_memory()
    except (OSError, ValueError, KeyError):
        pass
    automatic = memory_per_process_gb is None
    memory_per_process_gb = task_memory_budget(memory_per_process_gb, total, available)

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
    detail("{}: {} task(s), {} parallel worker(s).".format(task_kind.capitalize(), task_count, workers), flush=True)
    detail(
        f"Resources: {cpus} CPU threads, {memory_text} available memory; "
        f"{task_count} {task_kind} tasks, up to {workers} concurrent processes; "
        f"{budget_kind} memory budget {memory_per_process_gb:.1f} GiB/process.",
        flush=True,
    )
    return workers, threads


class MemoryPressureError(RuntimeError):
    """Stop carbon-source analysis before the system runs out of memory."""


def available_memory():
    """Read Linux/WSL memory and any tighter cgroup-v2 ancestor limits, in GiB."""
    values = {}
    for line in Path("/proc/meminfo").read_text().splitlines():
        name, _, value = line.partition(":")
        if name in ("MemTotal", "MemAvailable"):
            values[name] = int(value.split()[0]) / (1024 ** 2)
    available, total = values["MemAvailable"], values["MemTotal"]
    cgroup_root = Path("/sys/fs/cgroup")
    try:
        for line in Path("/proc/self/cgroup").read_text().splitlines():
            if line.startswith("0::"):
                group = cgroup_root / line.split("::", 1)[1].lstrip("/")
                for directory in (group,) + tuple(group.parents):
                    if directory != cgroup_root and cgroup_root not in directory.parents:
                        break
                    try:
                        limit = (directory / "memory.max").read_text().strip()
                        if limit != "max":
                            limit = int(limit) / (1024 ** 3)
                            used = int(
                                (directory / "memory.current").read_text()
                            ) / (1024 ** 3)
                            total = min(total, limit)
                            available = min(available, max(0, limit - used))
                    except (OSError, ValueError):
                        continue
                break
    except OSError:
        pass
    return available, total


class MemoryGuard:
    """Check a shared memory reserve before large allocations and exports."""

    def __init__(self, minimum_gb=2):
        if (
            type(minimum_gb) not in (int, float)
            or not math.isfinite(minimum_gb)
            or minimum_gb <= 0
        ):
            raise ValueError("Minimum available memory must be a finite positive GiB value.")
        self.minimum_gb = minimum_gb
        self.last_check = None

    def check_file(self, path, context="loading input"):
        """Reserve for decoded tables, not only a compressed workbook's file size."""
        path = Path(path)
        size = path.stat().st_size
        if path.suffix.lower() == ".xlsx":
            try:
                with ZipFile(str(path)) as archive:
                    size = sum(info.file_size for info in archive.infolist())
            except BadZipFile:
                pass  # The actual reader supplies the format error.
        self.check(context, force=True, extra_bytes=size * 12)

    def check(self, context="", force=False, extra_bytes=0):
        now = time.monotonic()
        if not force and self.last_check is not None and now - self.last_check < 0.5:
            return
        self.last_check = now
        try:
            available, total = available_memory()
        except (OSError, ValueError, KeyError) as error:
            raise MemoryPressureError(
                "Cannot read available memory; stopping propagation safely: " + str(error)
            )
        reserve = max(self.minimum_gb, total * 0.10)
        required = reserve + extra_bytes / (1024 ** 3)
        if available < required:
            raise MemoryPressureError(
                "Low memory during {}: {:.2f} GiB available; {:.2f} GiB required "
                "(safety reserve plus estimated temporary allocation). "
                "Propagation stopped, not converged. No states were truncated. "
                "Free memory or reduce MAX_CARBON_PROCESSES before rerunning."
                .format(context or "carbon propagation", available, required)
            )


class DiskSpaceLimitError(OSError):
    """Keep the last complete workbook when working storage is nearly full."""


class DiskGuard:
    def __init__(self, directory):
        self.directory, self.last_check = str(directory), None

    def check(self, force=False):
        now = time.monotonic()
        if not force and self.last_check is not None and now - self.last_check < 1:
            return
        self.last_check = now
        total, _, free = shutil.disk_usage(self.directory)
        reserve = max(1024 ** 3, total * .01)
        if free < reserve:
            raise DiskSpaceLimitError(
                "Low working-disk space: {:.2f} GiB free; {:.2f} GiB reserve required. "
                "The incomplete export is discarded; existing outputs remain unchanged."
                .format(free / 1024 ** 3, reserve / 1024 ** 3))
