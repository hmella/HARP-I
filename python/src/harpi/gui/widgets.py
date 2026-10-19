"""Small shared widgets."""

from __future__ import annotations

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSlider,
    QSpinBox,
    QWidget,
)

__all__ = ["Navigator", "MplCanvas", "StageBanner", "mask_bbox", "nan_range"]


class Navigator(QWidget):
    """Frame slider plus optional slice spinner, shared by every tab."""

    changed = Signal(int, int)       # slice, frame

    def __init__(self, parent=None):
        super().__init__(parent)
        self.n_frames = 1
        self.n_slices = 1

        self.frame = QSlider(Qt.Horizontal)
        self.frame.setMinimum(0)
        self.frame.valueChanged.connect(self._emit)
        self.frame_label = QLabel("frame 1/1")
        self.frame_label.setMinimumWidth(90)

        self.slice = QSpinBox()
        self.slice.setMinimum(0)
        self.slice.valueChanged.connect(self._emit)
        self.slice_label = QLabel("slice")

        self.play = QPushButton("▶")
        self.play.setCheckable(True)
        self.play.setFixedWidth(32)
        self.play.setToolTip("Play / pause (Space)")

        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.play)
        lay.addWidget(self.frame_label)
        lay.addWidget(self.frame, 1)
        lay.addWidget(self.slice_label)
        lay.addWidget(self.slice)
        self.set_extent(1, 1)

    def set_extent(self, n_frames, n_slices=1):
        self.n_frames = max(int(n_frames), 1)
        self.n_slices = max(int(n_slices), 1)
        self.frame.setMaximum(self.n_frames - 1)
        self.slice.setMaximum(self.n_slices - 1)
        has_slices = self.n_slices > 1
        self.slice.setVisible(has_slices)
        self.slice_label.setVisible(has_slices)
        self._update_label()

    def _update_label(self):
        self.frame_label.setText(f"frame {self.frame.value() + 1}/{self.n_frames}")

    def _emit(self):
        self._update_label()
        self.changed.emit(self.slice.value(), self.frame.value())

    def step_frame(self, delta):
        self.frame.setValue((self.frame.value() + delta) % self.n_frames)

    def step_slice(self, delta):
        if self.n_slices > 1:
            self.slice.setValue((self.slice.value() + delta) % self.n_slices)


class MplCanvas(FigureCanvasQTAgg):
    """A matplotlib canvas for the plots."""

    def __init__(self, width=5, height=3, tight=True):
        self.fig = Figure(figsize=(width, height), layout="constrained" if tight else None)
        super().__init__(self.fig)

    def clear(self):
        self.fig.clear()

    def axes(self, *args, **kwargs):
        return self.fig.add_subplot(*args, **kwargs)


class StageBanner(QLabel):
    """Status text shown alongside a tab's controls."""

    def __init__(self, text="", parent=None):
        super().__init__(text, parent)
        self.setAlignment(Qt.AlignCenter)
        self.setWordWrap(True)
        self.setStyleSheet(
            "QLabel { background: #2c2f36; color: #ddd; padding: 10px; "
            "border-radius: 4px; }"
        )


def mask_bbox(mask, margin=6):
    """Bounding box of a mask as ``(r0, r1, c0, c1)``, padded by ``margin``."""
    m = np.asarray(mask, dtype=bool)
    if not m.any():
        return 0, m.shape[0], 0, m.shape[1]
    rows = np.flatnonzero(m.any(axis=1))
    cols = np.flatnonzero(m.any(axis=0))
    r0 = max(int(rows[0]) - margin, 0)
    r1 = min(int(rows[-1]) + margin + 1, m.shape[0])
    c0 = max(int(cols[0]) - margin, 0)
    c1 = min(int(cols[-1]) + margin + 1, m.shape[1])
    return r0, r1, c0, c1


def nan_range(a, lo=2, hi=98):
    """Percentile display range, ignoring NaN."""
    a = np.asarray(a, dtype=float)
    f = a[np.isfinite(a)]
    if f.size == 0:
        return 0.0, 1.0
    return float(np.percentile(f, lo)), float(np.percentile(f, hi))
