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
    tab.ins_ant_x.setValue(24.0)
    tab.ins_ant_y.setValue(8.0)
    tab.run()
    after = tab.session.result(Stage.STRAIN)["segments"].segments
    both = np.isfinite(before) & np.isfinite(after)
    assert not np.array_equal(before[both], after[both])


def test_clicking_the_map_sets_the_picked_point(tab):
    tab.run()
    tab.ins_auto.setChecked(False)
    tab._set_pick("anterior", True)

    class Event:
        inaxes = tab.plot.fig.axes[0]
        xdata, ydata = 5.0, 6.0

    tab._clicked(Event())
    r0, c0 = tab._crop
    # the library grid is 1-based and the panel shows a crop
    assert tab.ins_ant_x.value() == pytest.approx(c0 + 5.0 + 1.0)
    assert tab.ins_ant_y.value() == pytest.approx(r0 + 6.0 + 1.0)
    assert not tab.pick_ant.isChecked()


def test_picking_one_point_cancels_the_other(tab):
    tab.origin_auto.setChecked(False)
    tab.ins_auto.setChecked(False)
    tab.pick_origin.setChecked(True)
    tab.pick_ant.setChecked(True)
    assert not tab.pick_origin.isChecked()


def test_a_click_is_ignored_when_no_pick_is_armed(tab):
    tab.run()
    tab.ins_auto.setChecked(False)
    before = tab.ins_ant_x.value()

    class Event:
        inaxes = tab.plot.fig.axes[0]
        xdata, ydata = 5.0, 6.0

    tab._clicked(Event())
    assert tab.ins_ant_x.value() == before


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


# ------------------------------------------------- arming the pick buttons

def test_arming_the_second_picker_keeps_the_mode(tab):
    """Un-checking the first button lands back in _set_pick, which used to
    clear the mode just set: the button stayed lit and the click went
    nowhere."""
    tab.origin_auto.setChecked(False)
    tab.ins_auto.setChecked(False)
    tab.pick_origin.setChecked(True)
    assert tab._pick == "origin"
    tab.pick_ant.setChecked(True)
    assert tab._pick == "anterior", "the mode was cleared by the un-check"
    assert tab.pick_ant.isChecked()
    assert not tab.pick_origin.isChecked()


def test_a_click_after_switching_pickers_lands(tab):
    """The symptom of the bug above, from the user's side."""
    tab.run()
    tab.origin_auto.setChecked(False)
    tab.ins_auto.setChecked(False)
    tab.pick_origin.setChecked(True)
    tab.pick_ant.setChecked(True)

    class Event:
        inaxes = tab.plot.fig.axes[0]
        xdata, ydata = 5.0, 6.0

    tab._clicked(Event())
    r0, c0 = tab._crop
    assert tab.ins_ant_x.value() == pytest.approx(c0 + 6.0)
    assert tab.ins_ant_y.value() == pytest.approx(r0 + 7.0)


def test_un_toggling_a_picker_disarms_it(tab):
    tab.origin_auto.setChecked(False)
    tab.pick_origin.setChecked(True)
    tab.pick_origin.setChecked(False)
    assert tab._pick is None


def test_the_banner_names_the_point_that_was_set(tab):
    """It read 'None set to (6.0, 7.0)' because the mode was cleared first."""
    tab.run()
    tab.ins_auto.setChecked(False)
    tab.pick_ant.setChecked(True)

    class Event:
        inaxes = tab.plot.fig.axes[0]
        xdata, ydata = 5.0, 6.0

    tab._clicked(Event())
    assert "None" not in tab.banner.text(), tab.banner.text()
    assert "anterior" in tab.banner.text()


def test_a_click_outside_the_map_is_ignored(tab):
    """The curves would hand over a time, the bullseye an angle."""
    tab.run()
    tab.ins_auto.setChecked(False)
    tab.pick_ant.setChecked(True)
    before = (tab.ins_ant_x.value(), tab.ins_ant_y.value())

    class Event:
        inaxes = tab.plot.fig.axes[1]
        xdata, ydata = 5.0, 6.0

    tab._clicked(Event())
    assert (tab.ins_ant_x.value(), tab.ins_ant_y.value()) == before
    assert tab._pick == "anterior", "still armed, waiting for the map"


# ------------------------------------------------- the two insertion points

def _place_both(tab, anterior=(24.0, 8.0), inferior=(8.0, 32.0)):
    tab.ins_auto.setChecked(False)
    tab.ins_ant_x.setValue(anterior[0]); tab.ins_ant_y.setValue(anterior[1])
    tab.ins_inf_x.setValue(inferior[0]); tab.ins_inf_y.setValue(inferior[1])


def test_both_points_reach_the_library(tab):
    seen = {}
    real = tabs_motion.strain_by_segments

    def spy(*a, **kw):
        seen.update(kw)
        return real(*a, **kw)

    _place_both(tab)
    with mock.patch.object(tabs_motion, "strain_by_segments", spy):
        tab.run()
    assert seen["anterior"] == (24.0, 8.0)
    assert seen["inferior"] == (8.0, 32.0)


def test_the_inferior_point_is_optional(tab):
    """Leaving it unplaced keeps the single-point behaviour."""
    seen = {}
    real = tabs_motion.strain_by_segments

    def spy(*a, **kw):
        seen.update(kw)
        return real(*a, **kw)

    tab.ins_auto.setChecked(False)
    tab.ins_ant_x.setValue(24.0); tab.ins_ant_y.setValue(8.0)
    with mock.patch.object(tabs_motion, "strain_by_segments", spy):
        tab.run()
    assert seen["anterior"] == (24.0, 8.0)
    assert seen["inferior"] is None


def test_the_two_points_decide_the_direction(tab):
    _place_both(tab)
    tab.run()
    geom = tab.session.result(Stage.STRAIN)["segments"].geometry
    assert geom["direction_from"] == "insertions"


def test_the_clockwise_box_is_disabled_once_both_points_are_set(tab):
    """It cannot do anything then, so it should not look as if it can."""
    assert tab.clockwise.isEnabled()
    _place_both(tab)
    assert not tab.clockwise.isEnabled()
    tab.ins_inf_x.setValue(0.0); tab.ins_inf_y.setValue(0.0)
    assert tab.clockwise.isEnabled(), "one point again, so it applies"


def test_the_septum_readout_reports_the_measured_arc(tab):
    _place_both(tab)
    tab.run()
    text = tab.septum_lbl.text()
    assert "deg between the insertions" in text
    assert "120" in text, f"the model's arc should be named: {text}"


def test_the_septum_readout_warns_when_the_points_look_swapped(tab):
    """A departure that large usually means they are the wrong way round."""
    _place_both(tab, anterior=(24.0, 8.0), inferior=(25.0, 9.0))
    tab.run()
    assert "Check which point is which" in tab.septum_lbl.text()


def test_the_septum_readout_says_when_only_one_point_is_given(tab):
    tab.ins_auto.setChecked(False)
    tab.ins_ant_x.setValue(24.0); tab.ins_ant_y.setValue(8.0)
    tab.run()
    assert "clockwise" in tab.septum_lbl.text()


def test_the_run_banner_carries_the_septum(tab):
    _place_both(tab)
    tab.run()
    assert "Septum:" in tab.banner.text()


def test_both_points_are_drawn_on_the_map(tab):
    """With their rays, so the numbering can be checked by eye."""
    _place_both(tab)
    tab.run()
    ax = tab.plot.fig.axes[0]
    markers = [ln.get_marker() for ln in ax.lines]
    assert "x" in markers, "the anterior point"
    assert "*" in markers, "the inferior point"
    styles = [ln.get_linestyle() for ln in ax.lines]
    assert "--" in styles and ":" in styles, "the two rays"


def test_the_settings_record_both_points(tab):
    _place_both(tab)
    tab.run()
    s = tab.settings()
    assert s["anterior insertion"] == (24.0, 8.0)
    assert s["inferior insertion"] == (8.0, 32.0)
    assert "deg" in s["septal arc"]


def test_picking_the_inferior_point_cancels_the_others(tab):
    tab.run()
    tab.origin_auto.setChecked(False)
    tab.ins_auto.setChecked(False)
    tab.pick_origin.setChecked(True)
    tab.pick_inf.setChecked(True)
    assert tab._pick == "inferior"
    assert not tab.pick_origin.isChecked() and not tab.pick_ant.isChecked()


def test_clicking_sets_the_inferior_point(tab):
    tab.run()
    tab.ins_auto.setChecked(False)
    tab.pick_inf.setChecked(True)

    class Event:
        inaxes = tab.plot.fig.axes[0]
        xdata, ydata = 7.0, 9.0

    tab._clicked(Event())
    r0, c0 = tab._crop
    assert tab.ins_inf_x.value() == pytest.approx(c0 + 8.0)
    assert tab.ins_inf_y.value() == pytest.approx(r0 + 10.0)


# -------------------------------------------- the tagged image underneath

def _backgrounds(tab):
    """Background layers only, told from the data by their z order."""
    return [im for ax in tab.plot.fig.axes for im in ax.images
            if im.get_zorder() == 0]


def test_the_overlay_is_on_by_default(tab):
    assert tab.bg.active.isChecked()


def test_the_tab_works_without_a_dataset(tab):
    """It is built from a mask and a motion result alone, which is how the
    whole of this file drives it."""
    assert tab.session.dataset is None
    tab.run()
    assert not _backgrounds(tab)
    assert tab.session.result(Stage.STRAIN) is not None


def _with_dataset(tab):
    from harpi import Dataset2D

    rows, cols, nfr = 48, 48, 6
    rng = np.random.default_rng(0)
    tab.session.dataset = Dataset2D(
        image=rng.normal(size=(rows, cols, 2, nfr)),
        mask=np.ones((rows, cols, nfr), bool),
        pixel_size=[2.0, 2.0], tag_spacing=8.0)
    return tab


def test_the_image_goes_under_the_map_and_the_segments(tab):
    _with_dataset(tab)
    tab.run()
    assert len(_backgrounds(tab)) >= 2, "the component map and the wedges"


def test_turning_the_overlay_off_removes_it(tab):
    _with_dataset(tab)
    tab.bg.active.setChecked(False)
    tab.run()
    assert not _backgrounds(tab)


def test_the_opacity_reaches_the_data_layer(tab):
    _with_dataset(tab)
    tab.bg.alpha.setValue(40)
    tab.run()
    data = [im for ax in tab.plot.fig.axes for im in ax.images
            if im.get_zorder() == 1]
    assert data
    assert all(im.get_alpha() == pytest.approx(0.4) for im in data)


def test_the_colour_bar_stays_solid(tab):
    """It is built from a mappable, not from the faded image, or the scale
    is as washed out as the data."""
    _with_dataset(tab)
    tab.bg.alpha.setValue(30)
    tab.run()
    bars = [ax for ax in tab.plot.fig.axes
            if getattr(ax, "_colorbar", None) is not None]
    assert bars
    for ax in bars:
        assert ax._colorbar.mappable.get_alpha() in (None, 1.0)
