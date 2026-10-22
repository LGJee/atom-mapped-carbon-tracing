"""Concise progress shared by mapping workers and carbon propagation."""

import os
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
    """Re-estimate one fixed work unit from recent samples, at most every 15 s."""

    def __init__(self, label, total=None, unit="items"):
        self.unit = unit
        self.last_output = time.monotonic()
        self.begin(label, total)

    def begin(self, label, total=None):
        # Reuse the output throttle while changing the current reaction.
        self.label, self.total = label, total
        self.started = self.sample_time = time.monotonic()
        self.sample_done, self.rate = 0, None

    def update(self, completed=0, force=False):
        now = time.monotonic()
        if not force and now - self.last_output < 15:
            return
        elapsed = now - self.started
        sample_elapsed = now - self.sample_time
        if sample_elapsed > 0 and completed > self.sample_done:
            recent = (completed - self.sample_done) / sample_elapsed
            self.rate = recent if self.rate is None else 0.5 * recent + 0.5 * self.rate
        self.sample_time, self.sample_done = now, completed
        self.last_output = now
        if self.total is None:
            count, eta = "{} {}".format(completed, self.unit), "Unavailable"
        else:
            total_text = ">10^300" if self.total.bit_length() > 1000 else str(self.total)
            count = "{}/{} {}".format(completed, total_text, self.unit)
            if completed >= self.total:
                eta = "0s"
            elif self.rate and elapsed >= 3:
                try:
                    eta = "~" + duration((self.total - completed) / self.rate)
                except (OverflowError, ValueError):
                    eta = "Unavailable"
            else:
                eta = "Estimating..."
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
