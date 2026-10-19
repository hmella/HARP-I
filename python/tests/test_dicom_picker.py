"""The series picker, over a synthetic catalogue."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication, QDialog  # noqa: E402

from harpi.harmonic import TagGeometry  # noqa: E402
from harpi.io.dicomio import SeriesInfo  # noqa: E402
from harpi.gui.dicom_picker import SeriesPicker, choose_series  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _geometry(spacing_px=6.0, angles=(0.0, np.pi / 2), confidence=50.0):
    return TagGeometry(
        centres=np.array([[128.0, 86.0], [86.0, 128.0]]),
        angles=np.asarray(angles, dtype=float),
        spacing_px=np.array([spacing_px, spacing_px]),
        confidence=confidence,
    )


def _series(description, tagged=True, n_frames=24, slice_location=None, uid=None):
    return SeriesInfo(
        uid=uid or f"uid-{description}",
        description=description,
        series_number=1,
        slice_location=slice_location,
        shape=(256, 256),
        pixel_spacing=(1.0, 1.0),
        n_files=n_frames,
        n_frames=n_frames,
        tag=_geometry() if tagged else _geometry(spacing_px=40.0, confidence=1.0),
    )


def test_a_plausible_geometry_marks_a_series_tagged():
    assert _series("tagged").looks_tagged
    assert not _series("cine", tagged=False).looks_tagged


def test_only_the_tagged_series_are_listed_at_first(qapp):
    cat = [_series("tag a"), _series("cine", tagged=False), _series("tag b")]
    dlg = SeriesPicker(cat)
    assert dlg.table.rowCount() == 2


def test_showing_all_reveals_the_rest(qapp):
    cat = [_series("tag a"), _series("cine", tagged=False)]
    dlg = SeriesPicker(cat)
    dlg.show_all.setChecked(True)
    assert dlg.table.rowCount() == 2


def test_nothing_tagged_falls_back_to_showing_everything(qapp):
    """A detection that scores nothing must not leave an empty table."""
    cat = [_series("a", tagged=False), _series("b", tagged=False)]
    dlg = SeriesPicker(cat)
    assert dlg.table.rowCount() == 2
    assert dlg.show_all.isChecked()


def test_the_summary_counts_the_series(qapp):
    cat = [_series("a"), _series("b"), _series("c", tagged=False)]
    dlg = SeriesPicker(cat)
    assert "3 series" in dlg.summary.text()
    assert "2 of them" in dlg.summary.text()


def test_the_table_shows_the_measured_geometry(qapp):
    dlg = SeriesPicker([_series("tag a", slice_location=-12.5)])
    row = [dlg.table.item(0, c).text() for c in range(dlg.table.columnCount())]
    assert "tag a" in row[0]
    assert "-12.5" in row[1]
    assert "24" in row[2]
    assert "256x256" in row[3]
    assert "mm" in row[4]
    assert "deg" in row[5]
    assert row[6] == "yes"


def test_a_multiframe_series_is_marked(qapp):
    s = _series("cine file")
    s.multiframe = True
    dlg = SeriesPicker([s])
    assert "multiframe" in dlg.table.item(0, 2).text()


def test_the_first_row_is_selected_so_open_always_works(qapp):
    dlg = SeriesPicker([_series("a"), _series("b")])
    assert dlg.chosen() is not None
    assert dlg.chosen().description == "a"


def test_selecting_a_row_returns_that_series(qapp):
    dlg = SeriesPicker([_series("a"), _series("b")])
    dlg.table.selectRow(1)
    assert dlg.chosen().description == "b"


def test_a_missing_slice_location_is_shown_as_a_dash(qapp):
    dlg = SeriesPicker([_series("a", slice_location=None)])
    assert dlg.table.item(0, 1).text() == "-"


def test_a_series_without_detection_is_shown_without_geometry(qapp):
    s = _series("unknown")
    s.tag = None
    dlg = SeriesPicker([s])
    assert dlg.table.item(0, 4).text() == "-"
    assert dlg.table.item(0, 5).text() == "-"


# ------------------------------------------------------------------ choosing

def test_one_series_needs_no_dialog(qapp):
    only = _series("the only one")
    assert choose_series([only]) is only


def test_one_tagged_series_among_many_needs_no_dialog(qapp):
    """The common case: a folder of cines with a single tagged acquisition."""
    tagged = _series("tagged")
    cat = [tagged, _series("cine a", tagged=False), _series("cine b", tagged=False)]
    assert choose_series(cat) is tagged


def test_an_empty_catalogue_returns_nothing(qapp):
    assert choose_series([]) is None


def test_several_tagged_series_open_the_dialog(qapp, monkeypatch):
    """Twelve short-axis slices is a clinical choice, not a technical one."""
    cat = [_series("slice a", slice_location=-10.0),
           _series("slice b", slice_location=-20.0)]
    seen = {}

    def fake_exec(self):
        seen["opened"] = True
        self.table.selectRow(1)
        return QDialog.Accepted

    monkeypatch.setattr(SeriesPicker, "exec", fake_exec)
    got = choose_series(cat)
    assert seen.get("opened")
    assert got.description == "slice b"


def test_cancelling_the_dialog_chooses_nothing(qapp, monkeypatch):
    cat = [_series("a"), _series("b")]
    monkeypatch.setattr(SeriesPicker, "exec", lambda self: QDialog.Rejected)
    assert choose_series(cat) is None
