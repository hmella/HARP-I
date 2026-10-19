"""Controls driven the way a user drives them.

Every other GUI test calls a handler directly, so a control wired to the
wrong slot, or to nothing, passes them all. That is how Clear this frame
reached release never having worked: it raised the moment it was pressed,
and no test had pressed it.

These go through the widget: `button.click()`, `box.setChecked(...)`,
`combo.setCurrentIndex(...)`, so the `connect` is part of what is tested.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication  # noqa: E402

from harpi.gui.session import Stage  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def window(qapp, dataset):
    from harpi.gui.app import MainWindow

    w = MainWindow()
    w.data_tab._adopt(dataset, "test")
    return w


def _ring(window, frame=0, n=16, shrink=1.0):
    ds = window.session.dataset
    yy, xx = np.nonzero(ds.mask[:, :, 0])
    cx, cy = xx.mean(), yy.mean()
    r = float(np.median(np.hypot(yy - cy, xx - cx))) * shrink
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ring = np.column_stack([cx + r * np.cos(t), cy + r * np.sin(t)])
    window.session.contours[(0, frame)] = [ring]
    window.session.labels[(0, frame)] = ["lv_epi"]
    window.session.current_frame = frame
    window.segment_tab.refresh()
    return ring


# --------------------------------------------------- buttons that are pressed

def test_clearing_a_frame_from_the_button(window):
    """The one that was broken: pressed, not called."""
    _ring(window)
    window.segment_tab.clear_btn.click()
    assert not window.session.contours.get((0, 0))


def test_deleting_a_contour_from_the_button(window):
    """_update_list clears the selection on every refresh, so pressing this
    after a frame step used to do nothing at all."""
    _ring(window)
    st = window.segment_tab
    assert st.contour_list.count() == 1
    st.contour_list.setCurrentRow(-1)
    st.del_btn.click()
    assert not window.session.contours.get((0, 0))


def test_deleting_with_nothing_drawn_says_so(window):
    st = window.segment_tab
    st.del_btn.click()
    assert "delete" in st.banner.text().lower()


def test_propagating_from_the_button(window, monkeypatch):
    from harpi.gui import tabs_segment

    ring = _ring(window)
    st = window.segment_tab
    st.prop_from.setValue(1)
    st.prop_to.setValue(window.session.dataset.n_frames)
    monkeypatch.setattr(tabs_segment.QMessageBox, "information",
                        lambda *a, **k: None)
    st.prop_btn.click()
    assert window.session.contours.get((0, 2))
    assert np.allclose(window.session.contours[(0, 2)][0], ring)


def test_tracking_from_the_button(window):
    window.filter_tab.run()
    ring = _ring(window)
    window.segment_tab.track_btn.click()
    moved = max(
        np.hypot(*(window.session.contours[(0, fr)][0] - ring).T).max()
        for fr in range(1, window.session.dataset.n_frames))
    assert moved > 0.1, "the contours did not follow the tissue"


def test_the_mask_from_the_dataset_button(window):
    window.segment_tab.load_mask_btn.click()
    assert window.session.mask is not None
    assert window.session.mask_source == "file"


def test_exporting_strain_before_computing_says_so(window):
    """It returned without a word, while every sibling wrote a banner. The
    default banner already mentions strain, so the text has to be the one
    this press produces."""
    window.strain_tab.banner.setText("")
    window.strain_tab.export_btn.click()
    assert "export" in window.strain_tab.banner.text().lower()


def test_running_strain_before_motion_says_so(window):
    window.strain_tab.banner.setText("")
    window.strain_tab.run_btn.click()
    assert "motion" in window.strain_tab.banner.text().lower()


def test_taking_a_snapshot_before_computing_says_so(window):
    window.strain_tab.banner.setText("")
    window.strain_tab.snap_btn.click()
    assert window.strain_tab.banner.text()


# ------------------------------------------------------------- shortcuts

def test_the_mask_toggle_hides_the_overlay(window):
    """toggle_mask flipped a flag that refresh never read, so M did nothing
    until some unrelated event redrew."""
    window.tabs.setCurrentWidget(window.segment_tab)
    window.segment_tab.load_mask_btn.click()
    canvas = window.segment_tab.canvas
    assert not canvas._overlay.pixmap().isNull(), "no overlay to hide"
    window._toggle_mask()
    assert canvas._overlay.pixmap().isNull(), "M did not hide it"
    window._toggle_mask()
    assert not canvas._overlay.pixmap().isNull(), "M did not bring it back"


def test_the_mask_toggle_reaches_the_tab_in_front(window):
    """It only ever acted on the Segmentation tab's canvases."""
    window.tabs.setCurrentWidget(window.filter_tab)
    before = window.filter_tab.k_before._show_mask
    window._toggle_mask()
    assert window.filter_tab.k_before._show_mask is not before


def test_fit_does_not_raise_on_a_plot_tab(window):
    """F looked for canvases that Motion, Strain and Compare do not have."""
    for tab in (window.motion_tab, window.strain_tab):
        window.tabs.setCurrentWidget(tab)
        window._fit()


# ----------------------------------------------------------- the wheel zoom

def test_a_refresh_keeps_the_zoom(window):
    """Zooming in to place a node accurately was undone the moment the node
    was released, because every refresh refitted the view."""
    _ring(window)
    canvas = window.segment_tab.canvas
    canvas._user_zoomed = True
    before = canvas.transform().m11()
    window.segment_tab.refresh()
    assert canvas._user_zoomed, "the zoom was forgotten"
    assert canvas.transform().m11() == pytest.approx(before)


def test_fit_still_overrides_the_zoom(window):
    """Pressing F is the way back to the whole image."""
    _ring(window)
    canvas = window.segment_tab.canvas
    canvas._user_zoomed = True
    canvas.fit()
    assert not canvas._user_zoomed


# ------------------------------------------- settings reaching another tab

def test_changing_the_connectivity_reaches_the_motion_tab(window):
    """invalidate_from stays quiet when there is nothing to drop, which is
    the usual case before the first solve, so the Motion tab went on
    reporting the old setting while solving with the new one."""
    window.unwrap_tab.connectivity.setCurrentText("8")
    assert "8" in window.motion_tab.unwrap_lbl.text()


def test_choosing_a_seed_mode_reaches_the_motion_tab(window):
    window.unwrap_tab.seed.setCurrentText("chosen point")
    assert window.unwrap_tab.pick_seed.isEnabled()


# ------------------------------------------------------- the overlay boxes

def test_the_motion_overlay_box_redraws(window):
    import harpi

    window.filter_tab.run()
    ds = window.session.dataset
    res = harpi.estimate_motion_2d(
        window.session.result(Stage.FILTER)["image"], ds.mask,
        a=ds.suggested_a, eta=15.0, njumps=2)
    window.session.set_result(Stage.MOTION, res)
    mt = window.motion_tab
    mt.refresh()

    def backgrounds():
        return [im for ax in mt.plot.fig.axes for im in ax.images
                if im.get_zorder() == 0]

    assert backgrounds(), "on by default"
    mt.bg.active.setChecked(False)
    assert not backgrounds(), "the box did not reach refresh"
    mt.bg.active.setChecked(True)
    assert backgrounds()


def test_the_slice_thickness_row_is_not_shown(window):
    """It could never be enabled, the tab taking 2-D studies only."""
    assert not window.data_tab.pixel_size_slice.isVisible()


# ------------------------------------- a result whose controls have moved

def test_note_preview_says_so_without_dropping_anything():
    """The third verb: invalidate_from drops and tells, note_settings always
    tells, this one tells and keeps."""
    from harpi.gui.session import Session, Stage

    s = Session()
    s.set_result(Stage.FILTER, "f")
    s.set_result(Stage.MOTION, "m")
    s.dirty = False
    seen = []
    s.previewChanged.connect(seen.append)

    s.note_preview(Stage.FILTER)
    assert seen == [int(Stage.FILTER)]
    assert s.has(Stage.FILTER) and s.has(Stage.MOTION), "nothing was dropped"
    assert not s.dirty, "nothing persisted changed"


def test_the_header_marks_a_result_its_controls_have_moved_past(window):
    """Only a tab that offers settings_match_result can reach this state, so
    the rule is exercised through a stub rather than waiting for one."""
    from harpi.gui.app import _STATE_ICON
    from harpi.gui.session import Stage

    def header():
        for i, (stage, _) in enumerate(window._tabs):
            if stage is Stage.FILTER:
                return window.tabs.tabText(i)
        raise AssertionError("no filtering tab")

    window.filter_tab.run()
    assert _STATE_ICON["done"] in header()

    matches = [True]
    window.filter_tab.settings_match_result = lambda: matches[0]
    window._update_headers()
    assert _STATE_ICON["done"] in header()

    matches[0] = False
    window._update_headers()
    assert _STATE_ICON["edited"] in header(), header()

    del window.filter_tab.settings_match_result
    window._update_headers()
    assert _STATE_ICON["done"] in header(), "a tab without the hook is unchanged"


def test_a_tab_without_the_hook_is_never_marked_edited(window):
    from harpi.gui.app import _STATE_ICON
    from harpi.gui.session import Stage

    window.filter_tab.run()
    window.session.set_result(Stage.MOTION, "m")
    window._update_headers()
    for i, (stage, _) in enumerate(window._tabs):
        if stage is Stage.MOTION:
            assert _STATE_ICON["edited"] not in window.tabs.tabText(i)


def test_the_motion_controls_are_reachable_in_the_default_window(window):
    """They are taller than the window, so without a scroll area the banner
    at the bottom is simply cut off."""
    from PySide6.QtWidgets import QScrollArea

    assert window.motion_tab.findChildren(QScrollArea), \
        "the motion controls cannot be scrolled to"


def test_the_arrow_scale_is_reachable_and_redraws(window):
    """Pressed through the widget, so the connection is covered."""
    import harpi
    from harpi.gui.session import Stage

    window.filter_tab.run()
    ds = window.session.dataset
    res = harpi.estimate_motion_2d(
        window.session.result(Stage.FILTER)["image"], ds.mask,
        a=ds.suggested_a, eta=15.0, njumps=2)
    window.session.set_result(Stage.MOTION, res)
    mt = window.motion_tab
    mt.refresh()

    def quiver_scale():
        for ax in mt.plot.fig.axes:
            for coll in ax.collections:
                if hasattr(coll, "U"):
                    return coll.scale
        raise AssertionError("no arrows")

    assert quiver_scale() == pytest.approx(1.0), "not true length by default"
    mt.arrow_scale.setValue(3.0)
    assert quiver_scale() == pytest.approx(1 / 3.0)
