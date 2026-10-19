# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

"""Small shared widgets."""

from __future__ import annotations

from collections import namedtuple

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QValidator
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSlider,
    QSpinBox,
    QWidget,
)

__all__ = ["Navigator", "FrameSpin", "MplCanvas", "StageBanner",
           "mask_bbox", "nan_range", "BACKGROUND_CHOICES", "Background",
           "background_image", "BackgroundControls"]


class FrameSpin(QSpinBox):
    """A frame number: counted from zero inside, shown from one.

    Every array axis here is zero based, and so is ``Navigator.frame.value()``,
    but a frame called 0 on screen reads as a mistake.
    """

    def textFromValue(self, value):
        return str(int(value) + 1)

    def valueFromText(self, text):
        try:
            return int(text) - 1
        except ValueError:
            return self.value()

    def validate(self, text, pos):
        if not text:
            return QValidator.Intermediate, text, pos
        if not text.isdigit():
            return QValidator.Invalid, text, pos
        n = int(text)
        if 1 <= n <= self.maximum() + 1:
            return QValidator.Acceptable, text, pos
        return QValidator.Intermediate, text, pos


class Navigator(QWidget):
    """Frame buttons plus optional slice spinner, shared by every tab."""

    changed = Signal(int, int)       # slice, frame

    def __init__(self, parent=None):
        super().__init__(parent)
        self.n_frames = 1
        self.n_slices = 1

        self.frame = FrameSpin()
        self.frame.setMinimum(0)
        self.frame.setToolTip("Frame number. Type one to jump to it.")
        self.frame.valueChanged.connect(self._emit)
        self.frame_label = QLabel("frame 1/1")
        self.frame_label.setMinimumWidth(90)

        # Stepping one frame at a time is what this is mostly used for, so it
        # gets buttons rather than a slider, which was easy to overshoot.
        # The ends get a double triangle against a bar, the steps a single
        # triangle. An earlier set used |◀ for first and ▶| for next, which are
        # mirror images of each other and read as the same control reversed.
        self.first = QPushButton("|◀◀")
        self.prev = QPushButton("◀")
        self.next = QPushButton("▶")
        self.last = QPushButton("▶▶|")
        for b, tip, step in ((self.first, "First frame (Home)", None),
                             (self.prev, "Previous frame ([)", -1),
                             (self.next, "Next frame (])", +1),
                             (self.last, "Last frame (End)", None)):
            b.setFixedWidth(42)
            b.setToolTip(tip)
            b.setAutoRepeat(True)
            b.setAutoRepeatDelay(300)
            b.setAutoRepeatInterval(80)
            if step is not None:
                b.clicked.connect(lambda _=False, d=step: self.step_frame(d))
        self.first.clicked.connect(lambda: self.frame.setValue(0))
        self.last.clicked.connect(lambda: self.frame.setValue(self.n_frames - 1))

        self.slice = QSpinBox()
        self.slice.setMinimum(0)
        self.slice.valueChanged.connect(self._emit)
        self.slice_label = QLabel("slice")

        # Spelled out, because a triangle here would be the next-frame glyph.
        self.play = QPushButton("Play")
        self.play.setCheckable(True)
        self.play.setFixedWidth(58)
        self.play.setToolTip("Play / pause (Space)")

        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.play)
        lay.addWidget(self.first)
        lay.addWidget(self.prev)
        lay.addWidget(self.frame)
        lay.addWidget(self.next)
        lay.addWidget(self.last)
        lay.addWidget(self.frame_label)
        lay.addStretch(1)
        lay.addWidget(self.slice_label)
        lay.addWidget(self.slice)
        self.set_extent(1, 1)

    def set_extent(self, n_frames, n_slices=1):
        self.n_frames = max(int(n_frames), 1)
        self.n_slices = max(int(n_slices), 1)
        self.frame.setMaximum(self.n_frames - 1)
        self.slice.setMaximum(self.n_slices - 1)
        single = self.n_frames <= 1
        for b in (self.first, self.prev, self.next, self.last, self.play):
            b.setEnabled(not single)
        self.frame.setEnabled(not single)
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


# What a tab can draw underneath its own data, in the order the boxes offer it.
BACKGROUND_CHOICES = (
    ("Magnitude (raw)", "magnitude"),
    ("Magnitude (raw, encoding 2)", "magnitude2"),
    ("Harmonic magnitude (filtered)", "harmonic"),
    ("Harmonic phase 1", "phase1"),
    ("Harmonic phase 2", "phase2"),
)

Background = namedtuple("Background", "array vmin vmax cmap title choice")
"""One frame ready to draw. ``cmap`` is None for the grey images, which is
what ImageCanvas expects; a matplotlib caller passes ``cmap or "gray"``."""


def background_image(dataset, filt, choice, frame, slice_index=0):
    """The tagged image a tab draws its own layer over.

    ``filt`` is the filter stage result, or None when it has not run, in which
    case a harmonic choice falls back to the raw magnitude. The returned
    ``choice`` is the one actually used, so a caller can say that it
    substituted. Returns None when there is no dataset.
    """
    if dataset is None or getattr(dataset, "image", None) is None:
        return None
    image = dataset.image
    frame = int(np.clip(frame, 0, image.shape[-1] - 1))

    def plane(array, encoding):
        """One frame of one encoding, whether or not there is a slice axis."""
        encoding = int(np.clip(encoding, 0, array.shape[-2] - 1))
        if array.ndim == 5:
            sl = int(np.clip(slice_index, 0, array.shape[2] - 1))
            return array[:, :, sl, encoding, frame]
        return array[:, :, encoding, frame]

    if choice not in ("magnitude", "magnitude2") and filt is None:
        choice = "magnitude"

    if choice == "magnitude2":
        a = np.real(plane(image, 1))
        return Background(a, *nan_range(a), None, "magnitude, encoding 2",
                          choice)
    if choice == "harmonic":
        a = np.abs(plane(filt["image"], 0))
        return Background(a, *nan_range(a), None, "harmonic magnitude", choice)
    if choice in ("phase1", "phase2"):
        d = 0 if choice == "phase1" else 1
        a = np.angle(plane(filt["image"], d))
        return Background(a, -np.pi, np.pi, "twilight",
                          f"harmonic phase {d + 1}", choice)
    a = np.real(plane(image, 0))
    return Background(a, *nan_range(a), None, "magnitude", "magnitude")


class BackgroundControls(QGroupBox):
    """Whether to draw the tagged image under a plot, which one, and how
    strongly the data layer covers it."""

    changed = Signal()

    def __init__(self, title="Background", parent=None):
        super().__init__(title, parent)
        # not called "show": QWidget already defines that
        self.active = QCheckBox("over the tagged image")
        self.active.setChecked(True)
        self.image = QComboBox()
        for label, key in BACKGROUND_CHOICES:
            self.image.addItem(label, key)
        self.alpha = QSlider(Qt.Horizontal)
        self.alpha.setRange(10, 100)
        self.alpha.setValue(75)
        self.alpha_lbl = QLabel("75%")

        self.active.toggled.connect(self._emit)
        self.image.currentIndexChanged.connect(self._emit)
        self.alpha.valueChanged.connect(self._alpha_changed)

        form = QFormLayout()
        form.addRow("", self.active)
        form.addRow("image", self.image)
        row = QHBoxLayout()
        row.addWidget(self.alpha); row.addWidget(self.alpha_lbl)
        holder = QWidget(); holder.setLayout(row)
        form.addRow("data opacity", holder)
        self.setLayout(form)

    def _alpha_changed(self, value):
        self.alpha_lbl.setText(f"{value}%")
        self._emit()

    def _emit(self, *_):
        self.changed.emit()

    def opacity(self):
        """Alpha for the data layer, 0 to 1."""
        return self.alpha.value() / 100.0

    def choice(self):
        return self.image.currentData()

    def wanted(self):
        return self.active.isChecked()
