# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

"""Comparing two retained runs."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .snapshots import diff_settings
from .widgets import MplCanvas, StageBanner, mask_bbox, nan_range

__all__ = ["CompareTab"]


class CompareTab(QWidget):
    """Two snapshots side by side, with their difference and their settings.

    The panels share one colour scale, since two maps on independent scales
    cannot be compared by eye, and the difference panel gets its own symmetric
    scale centred on zero.
    """

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session

        self.a = QComboBox()
        self.b = QComboBox()
        for w in (self.a, self.b):
            w.currentIndexChanged.connect(self.refresh)
        self.component = QComboBox()
        self.component.addItems(["CC", "RR", "p1", "p2"])
        self.component.currentTextChanged.connect(self.refresh)

        self.drop_btn = QPushButton("Forget A")
        self.drop_btn.clicked.connect(self._drop)
        self.clear_btn = QPushButton("Forget all")
        self.clear_btn.clicked.connect(self._clear)

        form = QFormLayout()
        form.addRow("A", self.a)
        form.addRow("B", self.b)
        form.addRow("component", self.component)
        box = QGroupBox("Snapshots"); box.setLayout(form)

        self.banner = StageBanner(
            "No snapshots yet.\n\n"
            "Compute strain, then press Take snapshot on the Strain tab. "
            "Snapshots survive a settings change, so a result is not lost by "
            "trying something else."
        )

        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["setting", "A", "B"])
        self.table.setAlternatingRowColors(True)
        self.table.setMaximumHeight(220)

        side = QVBoxLayout()
        side.addWidget(box)
        side.addWidget(self.drop_btn)
        side.addWidget(self.clear_btn)
        side.addWidget(self.banner)
        side.addStretch(1)
        panel = QWidget(); panel.setLayout(side); panel.setMaximumWidth(320)

        self.plot = MplCanvas(9, 5)
        right = QVBoxLayout()
        right.addWidget(self.plot, 3)
        right.addWidget(self.table, 1)
        right_w = QWidget(); right_w.setLayout(right)

        lay = QHBoxLayout(self)
        lay.addWidget(panel)
        lay.addWidget(right_w, 1)

    # ----------------------------------------------------------- selection
    def _drop(self):
        self.session.snapshots.remove(self.a.currentIndex())
        self.sync()

    def _clear(self):
        self.session.snapshots.clear()
        self.sync()

    def sync(self):
        """Repopulate the selectors, keeping the current choices where possible."""
        labels = self.session.snapshots.labels
        for w, default in ((self.a, 0), (self.b, len(labels) - 1)):
            want = w.currentText()
            blocked = w.blockSignals(True)
            w.clear()
            w.addItems(labels)
            if want in labels:
                w.setCurrentText(want)
            elif labels:
                w.setCurrentIndex(max(default, 0))
            w.blockSignals(blocked)
        enabled = len(labels) > 0
        for w in (self.a, self.b, self.component, self.drop_btn, self.clear_btn):
            w.setEnabled(enabled)
        self.refresh()

    def _pair(self):
        store = self.session.snapshots
        if len(store) == 0:
            return None, None
        ia = min(max(self.a.currentIndex(), 0), len(store) - 1)
        ib = min(max(self.b.currentIndex(), 0), len(store) - 1)
        return store[ia], store[ib]

    # ------------------------------------------------------------- display
    def refresh(self):
        self.plot.clear()
        sa, sb = self._pair()
        if sa is None:
            self.table.setRowCount(0)
            self.plot.draw()
            return
        if len(self.session.snapshots) == 1:
            self.banner.setText("One snapshot. Take another to compare.")
        else:
            self.banner.setText(f"A: {sa.label}\nB: {sb.label}")

        comp = self.component.currentText()
        frame = self.session.current_frame
        fa, fb = sa.field(comp, frame), sb.field(comp, frame)

        self._maps(sa, sb, fa, fb, comp, frame)
        self._curves(sa, sb, comp)
        self._segments(sa, sb)
        self._fill_table(sa, sb)
        self.plot.draw()

    def _maps(self, sa, sb, fa, fb, comp, frame):
        if fa is None or fb is None:
            return
        both = np.concatenate([fa[np.isfinite(fa)].ravel(),
                               fb[np.isfinite(fb)].ravel()])
        lo, hi = nan_range(both) if both.size else (-0.5, 0.25)
        mask = None
        for s in (sa, sb):
            m = getattr(s.strain, "mask", None)
            mask = m if mask is None else (mask | m)
        r0, r1, c0, c1 = mask_bbox(mask) if mask is not None else (
            0, fa.shape[0], 0, fa.shape[1])

        for i, (fld, title) in enumerate(((fa, f"A {comp}"), (fb, f"B {comp}"))):
            ax = self.plot.axes(2, 3, i + 1)
            im = ax.imshow(fld[r0:r1, c0:c1], cmap="RdBu_r", vmin=lo, vmax=hi)
            ax.set_title(f"{title} frame {frame + 1}", fontsize=9)
            ax.axis("off")
            self.plot.fig.colorbar(im, ax=ax, shrink=0.75)

        ax = self.plot.axes(2, 3, 3)
        if fa.shape == fb.shape:
            d = fa - fb
            finite = d[np.isfinite(d)]
            lim = float(np.nanpercentile(np.abs(finite), 98)) if finite.size else 0.05
            lim = lim or 0.05
            im = ax.imshow(d[r0:r1, c0:c1], cmap="PuOr_r", vmin=-lim, vmax=lim)
            ax.set_title("A minus B", fontsize=9)
            self.plot.fig.colorbar(im, ax=ax, shrink=0.75)
        else:
            ax.text(0.5, 0.5, "shapes differ", ha="center", va="center",
                    fontsize=8, transform=ax.transAxes)
            ax.set_title("A minus B", fontsize=9)
        ax.axis("off")

    def _curves(self, sa, sb, comp):
        ax = self.plot.axes(2, 3, 4)
        for s, style, name in ((sa, "-o", "A"), (sb, "--s", "B")):
            c = s.curve(comp)
            if c is not None:
                ax.plot(c, style, ms=3, label=name)
        ax.axhline(0, color="#888", lw=0.5)
        ax.set_xlabel("frame"); ax.set_ylabel(f"global {comp}")
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
        ax.set_title("global curves", fontsize=9)

    def _segments(self, sa, sb):
        ax = self.plot.axes(2, 3, 5)
        pa = _segment_peaks(sa)
        pb = _segment_peaks(sb)
        if pa is None or pb is None or pa.size != pb.size:
            ax.text(0.5, 0.5, "segment counts differ", ha="center", va="center",
                    fontsize=8, transform=ax.transAxes)
            ax.set_xticks([]); ax.set_yticks([])
        else:
            pos = np.arange(pa.size)
            ax.barh(pos - 0.2, pa, height=0.4, label="A", color="#4878a8")
            ax.barh(pos + 0.2, pb, height=0.4, label="B", color="#d08b26")
            ax.set_yticks(pos)
            ax.set_yticklabels([str(i + 1) for i in pos], fontsize=6)
            ax.invert_yaxis()
            ax.axvline(0, color="#888", lw=0.5)
            ax.legend(fontsize=7)
            ax.grid(alpha=0.3, axis="x")
        ax.set_title("peak per segment", fontsize=9)

        ax2 = self.plot.axes(2, 3, 6)
        rows = []
        for name, s in (("A", sa), ("B", sb)):
            d = s.summary()
            rows.append(f"{name}: {d['pixels']} px, peak CC "
                        f"{_fmt(d['peak CC'])}, peak RR {_fmt(d['peak RR'])}")
        ax2.text(0.02, 0.7, "\n".join(rows), fontsize=7, va="top",
                 transform=ax2.transAxes, family="monospace")
        ax2.axis("off")
        ax2.set_title("summary", fontsize=9)

    def _fill_table(self, sa, sb):
        rows = diff_settings(sa.settings, sb.settings)
        self.table.setRowCount(len(rows))
        if not rows:
            self.table.setRowCount(1)
            item = QTableWidgetItem("the settings are identical")
            item.setFlags(Qt.ItemIsEnabled)
            self.table.setItem(0, 0, item)
            for c in (1, 2):
                self.table.setItem(0, c, QTableWidgetItem(""))
            return
        for r, (name, va, vb) in enumerate(rows):
            for c, v in enumerate((name, _fmt(va), _fmt(vb))):
                item = QTableWidgetItem(v)
                item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                self.table.setItem(r, c, item)
        self.table.resizeColumnsToContents()


def _segment_peaks(snap):
    seg = snap.segments
    curves = getattr(seg, "CC", None) if seg is not None else None
    if curves is None:
        curves = getattr(seg, "RR", None) if seg is not None else None
    if curves is None:
        return None
    return np.nanmin(np.asarray(curves, dtype=float), axis=1)


def _fmt(v):
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.4g}"
    if isinstance(v, (tuple, list, np.ndarray)):
        return ", ".join(_fmt(x) for x in np.ravel(v))
    return str(v)
