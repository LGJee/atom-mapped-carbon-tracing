"""Concise progress shared by mapping workers and carbon propagation."""

import os
import math
from collections import deque
import sys
import time
import traceback
from contextlib import contextmanager
from pathlib import Path


def verbose():
    return os.environ.get("CARBON_TRACING_VERBOSE") == "1"


def detail(*args, **kwargs):
    if verbose():
        print(*args, **kwargs)


@contextmanager
def library_output(path):
    """Capture library Python/C output inside one isolated mapper process only."""
    if verbose():
        yield
        return
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    sys.stdout.flush()
    sys.stderr.flush()
    with Path(path).open("a", encoding="utf-8") as log:
        saved = [os.dup(1), os.dup(2)]
        try:
            # RDKit may write directly to stderr, bypassing Python redirection.
            os.dup2(log.fileno(), 1)
            os.dup2(log.fileno(), 2)
            # Keep stream objects alive: library logging handlers may retain them.
            try:
                yield
            except Exception:
                traceback.print_exc(file=log)
                raise
            finally:
                sys.stdout.flush()
                sys.stderr.flush()
                log.flush()
        finally:
            for descriptor, original in zip((1, 2), saved):
                os.dup2(original, descriptor)
                os.close(original)


def duration(seconds):
    seconds = max(0, int(seconds + 0.5))
    if seconds >= 3600:
        return "{}h {:02d}m".format(seconds // 3600, seconds % 3600 // 60)
    if seconds >= 60:
        return "{}m {:02d}s".format(seconds // 60, seconds % 60)
    return "{}s".format(seconds)


class Progress:
    """Estimate the current fixed phase from completed work in a recent window."""

    def __init__(self, label, total=None, unit="items", estimate=True):
        self.unit, self.estimate = unit, estimate
        self.last_output = time.monotonic()
        self.begin(label, total)

    def begin(self, label, total=None, estimate=None):
        # Keep the output throttle, but never reuse speed across phases/reactions.
        self.label, self.total = label, total
        if estimate is not None:
            self.estimate = estimate
        self.started = self.last_advance = time.monotonic()
        self.completed, self.rate = 0, None
        self.samples = deque([(self.started, 0)])

    def update(self, completed=0, force=False):
        now = time.monotonic()
        if completed < self.completed:
            self.begin(self.label, self.total)
        if completed > self.completed:
            self.last_advance = now
        self.completed = completed
        # Sample independently of printing; batched callbacks must not imply zero-time work.
        if now - self.samples[-1][0] >= 1:
            self.samples.append((now, completed))
            while len(self.samples) > 2 and self.samples[1][0] <= now - 10:
                self.samples.popleft()
        self.rate = None
        elapsed = now - self.started
        if self.estimate and elapsed >= 5 and completed >= 2 and now - self.last_advance < 15:
            sample_time, sample_done = self.samples[0]
            sample_elapsed = now - sample_time
            if sample_elapsed > 0:
                try:
                    rate = (completed - sample_done) / sample_elapsed
                    if rate > 0 and math.isfinite(rate):
                        self.rate = rate
                except OverflowError:
                    pass
        if not force and now - self.last_output < 15:
            return
        self.last_output = now
        completed_text = ">10^300" if completed.bit_length() > 1000 else str(completed)
        if self.total is None:
            count, eta = "{} {}".format(completed_text, self.unit), "Unavailable"
        else:
            total_text = ">10^300" if self.total.bit_length() > 1000 else str(self.total)
            count = "{}/{} {}".format(completed_text, total_text, self.unit)
            eta = "Estimating..." if self.estimate else "Unavailable"
            if completed >= self.total:
                eta = "0s"
            elif self.rate:
                try:
                    remaining = (self.total - completed) / self.rate
                    # Coarse estimates communicate uncertainty; elapsed time stays exact.
                    quantum = 5 if remaining < 60 else 15 if remaining < 600 else 60
                    eta = "~" + duration(math.ceil(remaining / quantum) * quantum)
                except (OverflowError, ValueError):
                    eta = "Unavailable"
        if verbose():
            print("{} | {} | elapsed {} | ETA {}".format(
                self.label, count, duration(elapsed), eta), flush=True)
        else:
            if self.total is None:
                message = "still working ({} elapsed)".format(duration(elapsed))
            elif completed >= self.total:
                message = "finished ({} {})".format(total_text, self.unit)
            else:
                message = "{}% ({})".format(int(100 * completed // self.total), count)
                if eta.startswith("~"):
                    message += "; about {} left".format(eta[1:])
            print("{}: {}.".format(self.label.replace(" | ", ": "), message), flush=True)
