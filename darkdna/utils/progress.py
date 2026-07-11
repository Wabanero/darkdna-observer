"""Small console progress helpers for long CLI runs."""

from __future__ import annotations

import sys
import time
from typing import TextIO


def _elapsed(seconds: float) -> str:
    seconds = max(0, int(seconds))
    minutes, sec = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return f"{hours:02d}:{minutes:02d}:{sec:02d}"
    return f"{minutes:02d}:{sec:02d}"


def progress_message(label: str, message: str, stream: TextIO | None = None) -> None:
    """Print a flushed status line suitable for console and redirected logs."""

    out = stream or sys.stdout
    print(f"[{label}] {message}", file=out, flush=True)


class ProgressReporter:
    """Line-based ASCII progress reporter.

    It intentionally prints complete lines rather than carriage-return updates,
    so redirected log files remain readable.
    """

    def __init__(
        self,
        label: str,
        total: int | None = None,
        *,
        width: int = 24,
        step_percent: float = 5.0,
        min_interval: float = 5.0,
        stream: TextIO | None = None,
    ) -> None:
        self.label = label
        self.total = int(total) if total is not None else None
        self.width = width
        self.step_fraction = max(0.001, step_percent / 100.0)
        self.min_interval = min_interval
        self.stream = stream or sys.stdout
        self.current = 0
        self.started = time.monotonic()
        self.last_print = 0.0
        self.next_fraction = 0.0
        self._finished = False

    def start(self, message: str | None = None) -> None:
        suffix = f": {message}" if message else ""
        total = f" total={self.total}" if self.total is not None else ""
        progress_message(self.label, f"start{total}{suffix}", self.stream)
        self.last_print = time.monotonic()

    def update(self, current: int, *, message: str | None = None, force: bool = False) -> None:
        if self._finished:
            return
        self.current = int(current)
        now = time.monotonic()
        fraction = None
        if self.total and self.total > 0:
            fraction = min(1.0, max(0.0, self.current / self.total))
        should_print = force or now - self.last_print >= self.min_interval
        if fraction is not None and fraction >= self.next_fraction:
            should_print = True
            while self.next_fraction <= fraction:
                self.next_fraction += self.step_fraction
        if fraction == 1.0:
            should_print = True
        if not should_print:
            return

        elapsed = _elapsed(now - self.started)
        suffix = f" {message}" if message else ""
        if fraction is None:
            print(f"[{self.label}] {self.current} elapsed={elapsed}{suffix}", file=self.stream, flush=True)
        else:
            filled = int(round(fraction * self.width))
            bar = "#" * filled + "-" * (self.width - filled)
            percent = fraction * 100.0
            print(
                f"[{self.label}] {percent:5.1f}% [{bar}] {self.current}/{self.total} elapsed={elapsed}{suffix}",
                file=self.stream,
                flush=True,
            )
        self.last_print = now

    def step(self, amount: int = 1, *, message: str | None = None, force: bool = False) -> None:
        self.update(self.current + amount, message=message, force=force)

    def finish(self, message: str | None = None) -> None:
        if self._finished:
            return
        if self.total is not None:
            self.update(self.total, message=message, force=True)
        else:
            self.update(self.current, message=message, force=True)
        progress_message(self.label, f"done elapsed={_elapsed(time.monotonic() - self.started)}", self.stream)
        self._finished = True
