"""The Strain tab: parameter exposure and regional analysis."""

from __future__ import annotations

import csv
from unittest import mock

import numpy as np
import pytest

pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication  # noqa: E402

from harpi import MotionResult  # noqa: E402
from harpi.gui.session import Session, Stage  # noqa: E402
from harpi.gui import tabs_motion  # noqa: E402
from harpi.gui.tabs_motion import StrainTab  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def tab(qapp):
    """A Strain tab holding a contracting ring, ready to compute."""
    rows = cols = 48
    yy, xx = np.mgrid[0:rows, 0:cols]
    r = np.hypot(xx - 24.0, yy - 24.0)
    ring = (r > 8) & (r < 17)
    nfr = 6
    u = np.zeros((rows, cols, 2, nfr))
    for k in range(nfr):
        s = -0.10 * np.sin(np.pi * k / (nfr - 1))
        u[:, :, 0, k] = s * (xx - 24.0)
        u[:, :, 1, k] = s * (yy - 24.0)

    sess = Session()
    sess.mask = np.repeat(ring[:, :, None], nfr, 2)
    sess.set_result(Stage.MOTION, MotionResult(u=u))
    return StrainTab(sess)


def test_running_with_defaults_produces_segments(tab):
    tab.run()
    res = tab.session.result(Stage.STRAIN)
    assert res is not None
    assert res["segments"].n_segments == 6


def test_only_computed_components_are_offered(tab):
    tab.run()
    offered = [tab.component.itemText(i) for i in range(tab.component.count())]
    assert "CC (circumferential)" in offered
    # LL is a 3-D component and must not appear for a 2-D result
    assert not any(o.startswith("LL") for o in offered)


def test_the_table_is_filled_with_aha_names(tab):
    tab.level.setCurrentText("mid")
    tab.run()
    assert tab.table.rowCount() == 6
    names = [tab.table.item(r, 1).text() for r in range(6)]
    assert names[0] == "mid anterior"


def test_choosing_apical_offers_four_segments(tab):
    tab.level.setCurrentText("apical")
    assert tab.nseg.value() == 4
    tab.run()
    assert tab.table.rowCount() == 4
    assert tab.table.item(0, 1).text() == "apical anterior"


def test_a_user_centre_reaches_the_result(tab):
    tab.origin_auto.setChecked(False)
    tab.origin_x.setValue(25.0)
    tab.origin_y.setValue(26.0)
    tab.run()
    st = tab.session.result(Stage.STRAIN)["strain"]
    assert st.origin == pytest.approx([25.0, 26.0])


def test_an_insertion_point_rotates_the_segments(tab):
    tab.run()
    before = tab.session.result(Stage.STRAIN)["segments"].segments.copy()
    tab.ins_auto.setChecked(False)
    tab.ins_x.setValue(24.0)
    tab.ins_y.setValue(8.0)
    tab.run()
    after = tab.session.result(Stage.STRAIN)["segments"].segments
    both = np.isfinite(before) & np.isfinite(after)
    assert not np.array_equal(before[both], after[both])


def test_clicking_the_map_sets_the_picked_point(tab):
    tab.run()
    tab.ins_auto.setChecked(False)
    tab._set_pick("insertion", True)

    class Event:
        inaxes = tab.plot.fig.axes[0]
        xdata, ydata = 5.0, 6.0

    tab._clicked(Event())
    r0, c0 = tab._crop
    # the library grid is 1-based and the panel shows a crop
    assert tab.ins_x.value() == pytest.approx(c0 + 5.0 + 1.0)
    assert tab.ins_y.value() == pytest.approx(r0 + 6.0 + 1.0)
    assert not tab.pick_ins.isChecked()


def test_picking_one_point_cancels_the_other(tab):
    tab.origin_auto.setChecked(False)
    tab.ins_auto.setChecked(False)
    tab.pick_origin.setChecked(True)
    tab.pick_ins.setChecked(True)
    assert not tab.pick_origin.isChecked()


def test_a_click_is_ignored_when_no_pick_is_armed(tab):
    tab.run()
    tab.ins_auto.setChecked(False)
    before = tab.ins_x.value()

    class Event:
        inaxes = tab.plot.fig.axes[0]
        xdata, ydata = 5.0, 6.0

    tab._clicked(Event())
    assert tab.ins_x.value() == before


def test_the_outlier_threshold_is_forwarded(tab):
    """The cut is one sided: it drops strain above the threshold.

    A contracting ring never exceeds it, so the tab's job is only to pass the
    value through. The rejection itself belongs to strain_by_segments.
    """
    tab.outlier.setValue(0.25)
    seen = {}
    real = tabs_motion.strain_by_segments

    def spy(st, **kw):
        seen.update(kw)
        return real(st, **kw)

    with mock.patch.object(tabs_motion, "strain_by_segments", spy):
        tab.run()
    assert seen["outlier_threshold"] == pytest.approx(0.25)


def test_the_time_axis_is_frames_without_an_acquisition(tab):
    """No dataset means no frame times, so the table counts frames."""
    tab.run()
    ttp = float(tab.table.item(0, 3).text())
    assert ttp == pytest.approx(round(ttp))        # a frame index, not a duration


def test_acquisition_times_reach_the_table(tab):
    """Frame times are an acquisition property, set on the Data tab."""
    from harpi import Dataset2D

    tab.run()
    ttp_frames = float(tab.table.item(0, 3).text())

    nfr = tab.session.result(Stage.MOTION).u.shape[-1]
    tab.session.dataset = Dataset2D(
        image=np.zeros((48, 48, 2, nfr)), mask=np.ones((48, 48, nfr)),
        pixel_size=[0.002, 0.002], tag_spacing=0.008,
        trigger_times=np.arange(nfr) * 40.0)
    tab.run()
    ttp_ms = float(tab.table.item(0, 3).text())
    assert ttp_ms == pytest.approx(ttp_frames * 40.0)


def test_a_manual_display_range_is_used(tab):
    tab.run()
    tab.auto_range.setChecked(False)
    tab.vmin.setValue(-0.3)
    tab.vmax.setValue(0.1)
    assert tab._range(np.zeros((4, 4))) == (-0.3, 0.1)


def test_the_percentile_range_follows_the_data(tab):
    tab.run()
    tab.auto_range.setChecked(True)
    field = np.concatenate([np.full(50, -0.2), np.full(50, 0.1)])
    lo, hi = tab._range(field)
    assert lo < hi
    assert -0.3 < lo < 0.0


def test_every_component_can_be_displayed(tab):
    tab.run()
    for i in range(tab.component.count()):
        tab.component.setCurrentIndex(i)
        tab.refresh()          # must not raise for any component


def test_a_contracting_ring_shortens_in_every_segment(tab):
    tab.run()
    peaks = [float(tab.table.item(r, 2).text()) for r in range(tab.table.rowCount())]
    assert all(p < 0 for p in peaks), peaks


def test_the_segment_table_exports_as_csv(tab, tmp_path):
    tab.level.setCurrentText("mid")
    tab.run()
    path = tmp_path / "segments.csv"
    with mock.patch("harpi.gui.tabs_motion.QFileDialog.getSaveFileName",
                    return_value=(str(path), "")):
        tab.export()
    rows = list(csv.reader(path.open()))
    assert rows[0][:2] == ["segment", "name"]
    assert len(rows) == 7                      # header plus six segments
    assert rows[1][1] == "mid anterior"


def test_changing_a_setting_warns_that_a_rerun_is_needed(tab):
    tab.run()
    tab.nseg.setValue(8)
    assert "again" in tab.banner.text().lower()


def test_refresh_without_a_result_is_harmless(tab):
    tab.refresh()
    assert tab.table.rowCount() == 0
