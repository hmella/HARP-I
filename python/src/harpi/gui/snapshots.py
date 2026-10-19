"""Keeping earlier runs so two can be compared.

Once every parameter is exposed, it becomes easy to lose track of which settings
produced a result. A snapshot records a result together with the settings behind
it, so two runs can be put side by side and the difference read against the
settings that differ.

Snapshots are held outside ``Session._results`` on purpose. ``invalidate_from``
sweeps results forward when an upstream setting changes, and a snapshot must
survive exactly that, since its whole point is to outlive the settings that made
it. Nothing here imports Qt.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

__all__ = ["Snapshot", "SnapshotStore", "diff_settings"]


@dataclass
class Snapshot:
    """One retained result and the settings that produced it."""

    label: str
    settings: dict = field(default_factory=dict)
    strain: object = None
    segments: object = None
    taken: float = field(default_factory=time.time)
    notes: str = ""

    @property
    def when(self) -> str:
        return time.strftime("%H:%M:%S", time.localtime(self.taken))

    def curve(self, component="CC"):
        """Global curve, or None when that component was not computed."""
        if self.strain is None:
            return None
        try:
            return self.strain.global_curve(component)
        except (AttributeError, ValueError):
            return None

    def field(self, component, frame):
        """One frame of a pixelwise component, or None."""
        v = getattr(self.strain, component, None) if self.strain else None
        if v is None:
            return None
        frame = min(max(int(frame), 0), v.shape[-1] - 1)
        return v[..., frame]

    def summary(self) -> dict:
        """The few numbers worth showing next to the settings."""
        out = {"pixels": None, "peak CC": None, "peak RR": None}
        if self.strain is not None:
            mask = getattr(self.strain, "mask", None)
            if mask is not None:
                out["pixels"] = int(np.asarray(mask).sum())
            for comp in ("CC", "RR"):
                c = self.curve(comp)
                if c is not None and np.isfinite(c).any():
                    out[f"peak {comp}"] = float(np.nanmin(c))
        return out


class SnapshotStore:
    """A capped list of snapshots, newest last."""

    def __init__(self, limit=6):
        self.limit = int(limit)
        self._items: list[Snapshot] = []

    def __len__(self):
        return len(self._items)

    def __iter__(self):
        return iter(self._items)

    def __getitem__(self, i):
        return self._items[i]

    @property
    def labels(self) -> list[str]:
        return [f"{s.label} ({s.when})" for s in self._items]

    def add(self, snapshot: Snapshot) -> Snapshot:
        """Store one, dropping the oldest when the cap is reached."""
        self._items.append(snapshot)
        while len(self._items) > self.limit:
            self._items.pop(0)
        return snapshot

    def remove(self, index: int) -> None:
        if 0 <= index < len(self._items):
            self._items.pop(index)

    def clear(self) -> None:
        self._items.clear()


def diff_settings(a: dict, b: dict) -> list[tuple[str, object, object]]:
    """Settings that differ, as ``(name, a value, b value)``.

    Keys missing from one side are reported with None, so a setting that only
    became available later still shows up rather than being silently skipped.
    """
    out = []
    for key in sorted(set(a) | set(b)):
        va, vb = a.get(key), b.get(key)
        if not _same(va, vb):
            out.append((key, va, vb))
    return out


def _same(x, y) -> bool:
    if isinstance(x, float) and isinstance(y, float):
        return bool(np.isclose(x, y, equal_nan=True))
    if isinstance(x, np.ndarray) or isinstance(y, np.ndarray):
        try:
            return bool(np.array_equal(np.asarray(x), np.asarray(y)))
        except (TypeError, ValueError):
            return False
    return x == y
