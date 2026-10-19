"""Snapshots and the Compare tab."""

from __future__ import annotations

import numpy as np
import pytest

from harpi.gui.snapshots import Snapshot, SnapshotStore, diff_settings


# ------------------------------------------------------------------ the store

def test_the_store_keeps_order():
    st = SnapshotStore()
    st.add(Snapshot(label="one"))
    st.add(Snapshot(label="two"))
    assert [s.label for s in st] == ["one", "two"]


def test_the_cap_drops_the_oldest():
    st = SnapshotStore(limit=2)
    for n in ("a", "b", "c"):
        st.add(Snapshot(label=n))
    assert [s.label for s in st] == ["b", "c"]


def test_removing_and_clearing():
    st = SnapshotStore()
    st.add(Snapshot(label="a"))
    st.add(Snapshot(label="b"))
    st.remove(0)
    assert [s.label for s in st] == ["b"]
    st.clear()
    assert len(st) == 0


def test_removing_out_of_range_is_harmless():
    st = SnapshotStore()
    st.remove(3)
    assert len(st) == 0


def test_labels_carry_the_time():
    st = SnapshotStore()
    st.add(Snapshot(label="run 1"))
    assert st.labels[0].startswith("run 1 (")


# ------------------------------------------------------------- settings diff

def test_identical_settings_give_no_rows():
    assert diff_settings({"a": 1}, {"a": 1}) == []


def test_a_changed_value_is_reported():
    assert diff_settings({"a": 1}, {"a": 2}) == [("a", 1, 2)]


def test_a_missing_key_is_reported_rather_than_skipped():
    assert diff_settings({"a": 1}, {}) == [("a", 1, None)]
    assert diff_settings({}, {"b": 2}) == [("b", None, 2)]


def test_close_floats_count_as_equal():
    assert diff_settings({"eta": 15.0}, {"eta": 15.0}) == []
    assert diff_settings({"eta": 15.0}, {"eta": 15.5}) != []


def test_arrays_are_compared_by_value():
    a = {"a": np.array([1.0, 2.0])}
    b = {"a": np.array([1.0, 2.0])}
    assert diff_settings(a, b) == []
    assert diff_settings(a, {"a": np.array([1.0, 3.0])}) != []


def test_rows_are_sorted_for_a_stable_table():
    rows = diff_settings({"z": 1, "a": 1}, {"z": 2, "a": 2})
    assert [r[0] for r in rows] == ["a", "z"]


# ------------------------------------------------------------------- the tab

pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication  # noqa: E402

from harpi import MotionResult  # noqa: E402
from harpi.gui.session import Session, Stage  # noqa: E402
from harpi.gui.tabs_compare import CompareTab  # noqa: E402
from harpi.gui.tabs_motion import StrainTab  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def pair(qapp):
    """A Strain tab and a Compare tab over one contracting ring."""
    rows = cols = 40
    yy, xx = np.mgrid[0:rows, 0:cols]
    r = np.hypot(xx - 20.0, yy - 20.0)
    ring = (r > 7) & (r < 15)
    nfr = 5
    u = np.zeros((rows, cols, 2, nfr))
    for k in range(nfr):
        s = -0.10 * np.sin(np.pi * k / (nfr - 1))
        u[:, :, 0, k] = s * (xx - 20.0)
        u[:, :, 1, k] = s * (yy - 20.0)

    sess = Session()
    sess.mask = np.repeat(ring[:, :, None], nfr, 2)
    sess.set_result(Stage.MOTION, MotionResult(u=u, meta={"kernel": "multiquadric3",
                                                          "a": 37.5, "eta": 15.0}))
    return StrainTab(sess), CompareTab(sess)


def test_a_snapshot_needs_a_result_first(pair):
    strain, _ = pair
    strain.take_snapshot()
    assert len(strain.session.snapshots) == 0
    assert "before taking" in strain.banner.text()


def test_taking_a_snapshot_keeps_the_settings(pair):
    strain, _ = pair
    strain.nseg.setValue(6)
    strain.run()
    strain.take_snapshot()
    assert len(strain.session.snapshots) == 1
    snap = strain.session.snapshots[0]
    assert snap.settings["segments"] == 6
    assert snap.settings["kernel"] == "multiquadric3"
    assert snap.strain is not None and snap.segments is not None


def test_a_snapshot_survives_invalidation(pair):
    """The whole point: it must outlive the settings that produced it."""
    strain, _ = pair
    strain.run()
    strain.take_snapshot()
    strain.session.invalidate_from(Stage.DATA)
    assert not strain.session.has(Stage.STRAIN)
    assert len(strain.session.snapshots) == 1
    assert strain.session.snapshots[0].strain is not None


def test_the_summary_reports_shortening(pair):
    strain, _ = pair
    strain.run()
    strain.take_snapshot()
    d = strain.session.snapshots[0].summary()
    assert d["pixels"] > 0
    assert d["peak CC"] < 0


def test_the_compare_tab_is_empty_until_a_snapshot_exists(pair):
    _, compare = pair
    compare.sync()
    assert compare.a.count() == 0
    assert not compare.a.isEnabled()
    assert "No snapshots" in compare.banner.text()


def test_one_snapshot_asks_for_another(pair):
    strain, compare = pair
    strain.run()
    strain.take_snapshot()
    compare.sync()
    assert compare.a.count() == 1
    assert "another" in compare.banner.text()


def test_two_snapshots_are_compared(pair):
    strain, compare = pair
    strain.nseg.setValue(6)
    strain.run(); strain.take_snapshot()
    strain.nseg.setValue(4)
    strain.run(); strain.take_snapshot()

    compare.sync()
    assert compare.a.count() == 2
    compare.a.setCurrentIndex(0)
    compare.b.setCurrentIndex(1)
    compare.refresh()
    names = [compare.table.item(r, 0).text()
             for r in range(compare.table.rowCount())]
    assert "segments" in names


def test_identical_settings_are_stated_rather_than_left_blank(pair):
    strain, compare = pair
    strain.run(); strain.take_snapshot()
    strain.run(); strain.take_snapshot()
    compare.sync()
    compare.a.setCurrentIndex(0); compare.b.setCurrentIndex(1)
    compare.refresh()
    assert compare.table.rowCount() == 1
    assert "identical" in compare.table.item(0, 0).text()


def test_differing_segment_counts_do_not_crash_the_bar_chart(pair):
    strain, compare = pair
    strain.nseg.setValue(6); strain.run(); strain.take_snapshot()
    strain.nseg.setValue(4); strain.run(); strain.take_snapshot()
    compare.sync()
    compare.a.setCurrentIndex(0); compare.b.setCurrentIndex(1)
    compare.refresh()          # must not raise


def test_forgetting_a_snapshot(pair):
    strain, compare = pair
    strain.run(); strain.take_snapshot()
    strain.run(); strain.take_snapshot()
    compare.sync()
    compare.a.setCurrentIndex(0)
    compare._drop()
    assert len(strain.session.snapshots) == 1


def test_forgetting_all_snapshots(pair):
    strain, compare = pair
    strain.run(); strain.take_snapshot()
    compare.sync()
    compare._clear()
    assert len(strain.session.snapshots) == 0
    assert compare.a.count() == 0


def test_every_component_can_be_compared(pair):
    strain, compare = pair
    strain.run(); strain.take_snapshot()
    strain.run(); strain.take_snapshot()
    compare.sync()
    for i in range(compare.component.count()):
        compare.component.setCurrentIndex(i)
        compare.refresh()      # must not raise for any component
