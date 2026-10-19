"""Choosing a series from a DICOM folder."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QLabel,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)

__all__ = ["SeriesPicker", "choose_series"]

HEADERS = ("series", "slice", "frames", "size", "tag spacing", "tag angles",
           "tagged")


class SeriesPicker(QDialog):
    """A table of the series in a folder, the tagged ones first.

    A clinical export can hold a couple of hundred series of which a handful are
    tagged, so the likely ones are listed first and the rest are behind a
    checkbox. They are never removed: the detection is a guess, and a series it
    scores badly may still be the one wanted.
    """

    def __init__(self, catalogue, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Choose a series")
        self.resize(820, 440)
        self._all = list(catalogue)
        self._rows: list = []

        n_tagged = sum(1 for s in self._all if s.looks_tagged)
        self.summary = QLabel(
            f"{len(self._all)} series, {n_tagged} of them look tagged. "
            "The tag spacing and angles are measured from the harmonic peaks."
        )
        self.summary.setWordWrap(True)

        self.show_all = QCheckBox("show the series that do not look tagged")
        self.show_all.setChecked(n_tagged == 0)
        self.show_all.toggled.connect(self._fill)

        self.table = QTableWidget(0, len(HEADERS))
        self.table.setHorizontalHeaderLabels(list(HEADERS))
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.setAlternatingRowColors(True)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.doubleClicked.connect(lambda _: self.accept())

        self.buttons = QDialogButtonBox(
            QDialogButtonBox.Open | QDialogButtonBox.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

        lay = QVBoxLayout(self)
        lay.addWidget(self.summary)
        lay.addWidget(self.show_all)
        lay.addWidget(self.table, 1)
        lay.addWidget(self.buttons)

        self._fill()

    # ------------------------------------------------------------- contents
    def _visible(self):
        if self.show_all.isChecked():
            return list(self._all)
        tagged = [s for s in self._all if s.looks_tagged]
        return tagged or list(self._all)

    def _fill(self):
        self._rows = self._visible()
        self.table.setRowCount(len(self._rows))
        for r, s in enumerate(self._rows):
            spacing = ("-" if s.tag is None
                       else f"{s.tag.spacing_mm(s.pixel_spacing):.2f} mm")
            if s.tag is None:
                angles = "-"
            else:
                a = np.degrees(s.tag.angles) % 180
                angles = f"{a[0]:.0f}, {a[1]:.0f} deg"
            cells = (
                s.description or s.protocol or s.uid[-12:],
                "-" if s.slice_location is None else f"{s.slice_location:+.1f}",
                str(s.n_frames) + (" multiframe" if s.multiframe else ""),
                f"{s.shape[0]}x{s.shape[1]}",
                spacing,
                angles,
                "yes" if s.looks_tagged else "",
            )
            for c, text in enumerate(cells):
                item = QTableWidgetItem(text)
                item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                self.table.setItem(r, c, item)
        self.table.resizeColumnsToContents()
        if self._rows:
            self.table.selectRow(0)

    def chosen(self):
        """The selected series, or None."""
        rows = {i.row() for i in self.table.selectedIndexes()}
        if not rows:
            return None
        r = sorted(rows)[0]
        return self._rows[r] if 0 <= r < len(self._rows) else None


def choose_series(catalogue, parent=None):
    """Pick a series, skipping the dialog when there is nothing to choose.

    Returns the chosen :class:`~harpi.io.dicomio.SeriesInfo`, or None if the
    dialog was cancelled.
    """
    catalogue = list(catalogue)
    if not catalogue:
        return None
    tagged = [s for s in catalogue if s.looks_tagged]
    if len(catalogue) == 1:
        return catalogue[0]
    if len(tagged) == 1:
        return tagged[0]

    dlg = SeriesPicker(catalogue, parent)
    if dlg.exec() != QDialog.Accepted:
        return None
    return dlg.chosen()
