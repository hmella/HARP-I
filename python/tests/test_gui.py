"""GUI tests, run offscreen.

These cover the state machine and the wiring rather than the rendering:
that stale results are cleared, that undo works, and that the whole chain
runs.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PySide6", reason="GUI extra not installed")

from PySide6.QtCore import QEvent, Qt  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from harpi.gui.canvas import AddContourCommand, DeleteContourCommand  # noqa: E402
from harpi import frequency_to_centre  # noqa: E402
from harpi.gui.session import Session, Stage  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


# ------------------------------------------------------------------ session

def test_invalidation_sweeps_downstream():
    s = Session()
    for st in Stage:
        s.set_result(st, f"result-{st}")
    assert s.has(Stage.STRAIN)

    s.invalidate_from(Stage.FILTER)
    assert s.has(Stage.DATA)
    for st in (Stage.FILTER, Stage.SEGMENT, Stage.UNWRAP, Stage.MOTION, Stage.STRAIN):
        assert not s.has(st), st


def test_setting_a_config_invalidates_that_stage():
    s = Session()
    s.set_result(Stage.MOTION, "m")
    s.set_config(Stage.FILTER, {"cutoff": 17})
    assert not s.has(Stage.MOTION)


def test_setting_an_identical_config_does_not_invalidate():
    s = Session()
    cfg = {"cutoff": 17}
    s.set_config(Stage.FILTER, cfg)
    s.set_result(Stage.MOTION, "m")
    s.set_config(Stage.FILTER, {"cutoff": 17})
    assert s.has(Stage.MOTION)


def test_new_result_clears_later_stages():
    s = Session()
    s.set_result(Stage.MOTION, "m")
    s.set_result(Stage.STRAIN, "s")
    s.set_result(Stage.MOTION, "m2")
    assert not s.has(Stage.STRAIN)


# ------------------------------------------------------------------- window

@pytest.fixture
def window(qapp, dataset):
    from harpi.gui.app import MainWindow

    w = MainWindow()
    w.data_tab._adopt(dataset, "test")
    return w


def test_full_pipeline(window, dataset):
    from skimage import measure

    window.filter_tab.run()
    assert window.session.has(Stage.FILTER), window.filter_tab.banner.text()

    for fr in range(dataset.n_frames):
        cs = measure.find_contours(dataset.mask[:, :, fr].astype(float), 0.5)
        cs = sorted(cs, key=lambda c: -len(c))[:2]
        window.session.contours[(0, fr)] = [
            np.column_stack([c[:, 1], c[:, 0]])[::4] for c in cs
        ]
    window.segment_tab._contours_changed()
    assert window.session.mask is not None

    window.unwrap_tab.run()
    assert window.session.has(Stage.UNWRAP)
    assert window.session.result(Stage.UNWRAP)["complete"].all()

    from harpi import estimate_motion_2d

    res = estimate_motion_2d(
        window.session.result(Stage.FILTER)["image"], window.session.mask,
        a=dataset.suggested_a, eta=15.0, temporal_fit=True,
    )
    window.session.set_result(Stage.MOTION, res)

    window.strain_tab.run()
    st = window.session.result(Stage.STRAIN)["strain"]
    gcs = st.global_curve("CC")
    assert gcs[0] == pytest.approx(0, abs=1e-12)
    assert gcs.min() < -0.05, "a contracting ventricle shortens circumferentially"


def test_editing_a_filter_setting_keeps_motion(window):
    """Filtering previews its controls and commits on Apply, so twiddling a
    dial must not throw away a solved pipeline. It used to.

    This moves the kind, not the cutoff: cutoff does not apply to the
    transmission window and is left out of the recorded settings, so nudging
    it is correctly a no change, and setValue works on a disabled box, which
    would make this pass for the wrong reason.
    """
    ft = window.filter_tab
    ft.run()
    window.session.set_result(Stage.MOTION, "m")
    ft.kind.setCurrentText("butterworth")
    assert window.session.has(Stage.MOTION), "a solved field was dropped"
    assert window.session.has(Stage.FILTER)
    assert not ft.settings_match_result(), "the result is behind the controls"


def test_applying_a_filter_clears_motion(window):
    """What Apply still does, because the harmonic image really changed."""
    ft = window.filter_tab
    ft.run()
    window.session.set_result(Stage.MOTION, "m")
    ft.kind.setCurrentText("butterworth")
    ft.run()
    assert not window.session.has(Stage.MOTION)
    assert ft.settings_match_result()


def test_a_dead_setting_does_not_mark_the_result_behind(window):
    """Cutoff and order do nothing to the transmission window."""
    ft = window.filter_tab
    ft.run()
    ft.cutoff.setValue(ft.cutoff.value() + 3)
    assert ft.settings_match_result()


def test_tab_headers_show_state(window):
    titles = [window.tabs.tabText(i) for i in range(window.tabs.count())]
    assert titles[0].startswith("●"), "Data should be marked done"
    assert any("Strain" in t for t in titles)


def test_compare_is_not_a_pipeline_stage(window):
    """It carries no state marker, since it has no result of its own."""
    titles = [window.tabs.tabText(i) for i in range(window.tabs.count())]
    compare = [t for t in titles if "Compare" in t]
    assert len(compare) == 1
    assert not compare[0].startswith("●")
    # the stage tabs keep their index-to-stage pairing
    assert len(window._tabs) == window.tabs.count() - 1


def test_undo_add_and_delete_contour(window):
    c = window.segment_tab.canvases[0]
    c.key = (0, 0)
    tri = np.array([[10.0, 10.0], [20.0, 10.0], [15.0, 20.0]])

    window.undo_stack.push(AddContourCommand(c, tri))
    assert len(window.session.contours[(0, 0)]) == 1
    window.undo_stack.undo()
    assert len(window.session.contours[(0, 0)]) == 0
    window.undo_stack.redo()
    assert len(window.session.contours[(0, 0)]) == 1

    window.undo_stack.push(DeleteContourCommand(c, 0))
    assert len(window.session.contours[(0, 0)]) == 0
    window.undo_stack.undo()
    assert len(window.session.contours[(0, 0)]) == 1
    assert np.allclose(window.session.contours[(0, 0)][0], tri)


def test_propagation_prefers_explicit_range(window, dataset):
    c = window.segment_tab.canvases[0]
    c.key = (0, 0)
    tri = np.array([[10.0, 10.0], [30.0, 10.0], [20.0, 30.0]])
    window.session.contours[(0, 0)] = [tri]

    window.segment_tab.prop_from.setValue(1)
    window.segment_tab.prop_to.setValue(4)
    window.segment_tab._propagate()

    for fr in range(1, 4):
        assert np.allclose(window.session.contours[(0, fr)][0], tri)
    window.undo_stack.undo()
    assert not window.session.contours.get((0, 3))


def test_project_round_trip(window, tmp_path):
    window.session.contours[(0, 0)] = [np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 1.0]])]
    p = tmp_path / "s.harpi"
    window.session.project_path = str(p)
    window.save_project()

    from harpi.gui.app import MainWindow

    w2 = MainWindow()
    w2.load_project(p)
    assert (0, 0) in w2.session.contours
    assert np.allclose(w2.session.contours[(0, 0)][0],
                       window.session.contours[(0, 0)][0])


# ------------------------------------------------------- mask and wall UI

def test_mask_can_be_loaded_from_the_dataset(window):
    st = window.segment_tab
    st._load_mask_from_dataset()
    assert window.session.mask is not None
    assert window.session.mask_source == "file"
    assert "dataset" in st.source_lbl.text()


def test_a_loaded_mask_survives_a_refresh(window):
    st = window.segment_tab
    st._load_mask_from_dataset()
    before = window.session.mask.sum()
    st._contours_changed()
    assert window.session.mask is not None, "the rebuild discarded the loaded mask"
    assert window.session.mask.sum() == before


def test_drawing_replaces_a_loaded_mask(window):
    st = window.segment_tab
    st._load_mask_from_dataset()
    c = st.canvases[0]
    c.key = (0, 0)
    th = np.linspace(0, 2 * np.pi, 14, endpoint=False)
    outer = np.column_stack([80 + 20 * np.cos(th), 80 + 20 * np.sin(th)])
    window.undo_stack.push(AddContourCommand(c, outer, "lv_epi"))
    st._contours_changed()
    assert window.session.mask_source == "contours"


def test_mask_load_updates_the_tab_header(window):
    window.segment_tab._load_mask_from_dataset()
    titles = [window.tabs.tabText(i) for i in range(window.tabs.count())]
    assert titles[2].startswith("●"), "Segmentation should be marked done"


def test_wall_label_is_kept_through_undo(window):
    c = window.segment_tab.canvases[0]
    c.key = (0, 0)
    th = np.linspace(0, 2 * np.pi, 12, endpoint=False)
    pts = np.column_stack([80 + 15 * np.cos(th), 80 + 15 * np.sin(th)])

    window.undo_stack.push(AddContourCommand(c, pts, "rv_epi"))
    assert c.labels_for((0, 0)) == ["rv_epi"]

    window.undo_stack.push(DeleteContourCommand(c, 0))
    assert c.labels_for((0, 0)) == []
    window.undo_stack.undo()
    assert c.labels_for((0, 0)) == ["rv_epi"]


def test_walls_get_distinct_colours(window):
    c = window.segment_tab.canvases[0]
    c.key = (0, 0)
    th = np.linspace(0, 2 * np.pi, 12, endpoint=False)
    for i, wall in enumerate(("lv_epi", "lv_endo", "rv_epi", "rv_endo")):
        pts = np.column_stack([80 + (20 - 3 * i) * np.cos(th),
                               80 + (20 - 3 * i) * np.sin(th)])
        window.undo_stack.push(AddContourCommand(c, pts, wall))
    colours = {c.colour_for((0, 0), i).name() for i in range(4)}
    assert len(colours) == 4


def test_background_can_be_the_filtered_image(window):
    st = window.segment_tab
    window.filter_tab.run()
    for i in range(st.background.count()):
        st.background.setCurrentIndex(i)
        st.refresh()
        assert st.title.text(), st.background.itemText(i)
    st.background.setCurrentIndex(st.background.findData("harmonic"))
    st.refresh()
    assert "harmonic" in st.title.text()


def test_filtered_background_falls_back_before_filtering(window):
    """Choosing the harmonic image before filtering shows the magnitude."""
    st = window.segment_tab
    st.background.setCurrentIndex(st.background.findData("harmonic"))
    st.refresh()
    assert st.title.text() == "magnitude"


def test_project_round_trip_keeps_wall_labels(window, tmp_path):
    c = window.segment_tab.canvases[0]
    c.key = (0, 0)
    th = np.linspace(0, 2 * np.pi, 12, endpoint=False)
    for wall, r in (("lv_epi", 20), ("lv_endo", 11)):
        pts = np.column_stack([80 + r * np.cos(th), 80 + r * np.sin(th)])
        window.undo_stack.push(AddContourCommand(c, pts, wall))

    p = tmp_path / "walls.harpi"
    window.session.project_path = str(p)
    window.save_project()

    from harpi.gui.app import MainWindow

    w2 = MainWindow()
    w2.load_project(p)
    assert w2.session.labels[(0, 0)] == ["lv_epi", "lv_endo"]
    assert np.allclose(w2.session.contours[(0, 0)][0],
                       window.session.contours[(0, 0)][0])


def test_flip_controls_reuse_the_existing_solution(window, dataset):
    """Toggling a flip must not re-solve, and must update the result."""
    from harpi import estimate_motion_2d

    window.filter_tab.run()
    window.segment_tab._load_mask_from_dataset()
    res = estimate_motion_2d(window.session.result(Stage.FILTER)["image"],
                             window.session.mask, a=dataset.suggested_a, eta=15.0)
    mt = window.motion_tab
    mt._done(res)

    base = window.session.result(Stage.MOTION).best.copy()

    mt.flip_y.setChecked(True)
    flipped = window.session.result(Stage.MOTION).best
    assert np.allclose(flipped[..., 1, :], -base[..., 1, :], equal_nan=True)
    assert np.allclose(flipped[..., 0, :], base[..., 0, :], equal_nan=True)

    mt.flip_y.setChecked(False)
    assert np.allclose(window.session.result(Stage.MOTION).best, base, equal_nan=True)


def test_flipping_invalidates_strain(window, dataset):
    from harpi import estimate_motion_2d

    window.filter_tab.run()
    window.segment_tab._load_mask_from_dataset()
    res = estimate_motion_2d(window.session.result(Stage.FILTER)["image"],
                             window.session.mask, a=dataset.suggested_a, eta=15.0)
    window.motion_tab._done(res)
    window.strain_tab.run()
    assert window.session.has(Stage.STRAIN)

    window.motion_tab.swap_xy.setChecked(True)
    assert not window.session.has(Stage.STRAIN), "strain depends on the components"


def _solve(window, dataset):
    from harpi import estimate_motion_2d

    mt = window.motion_tab
    return estimate_motion_2d(
        window.session.result(Stage.FILTER)["image"], window.session.mask,
        kernel=mt.kernel.currentText(), a=mt.a.value(), eta=mt.eta.value(),
        temporal_fit=mt.tfit.isChecked(), temporal_order=mt.torder.value(),
        njumps=mt.njumps.value(), phase_jump_rule=mt.jump_rule.currentData(),
    )


@pytest.fixture
def solved(window, dataset):
    window.filter_tab.run()
    window.segment_tab._load_mask_from_dataset()
    window.motion_tab._done(_solve(window, dataset))
    return window


@pytest.mark.parametrize("control,value", [
    ("torder", 6),
    ("njumps", 3),
    ("downsample", 2),
    ("eta", 1.0),
])
def test_changing_a_solver_setting_clears_the_result(solved, control, value):
    """A stale result must not stay on screen looking current."""
    assert solved.session.has(Stage.MOTION)
    getattr(solved.motion_tab, control).setValue(value)
    assert not solved.session.has(Stage.MOTION), control
    assert "re-run" in solved.motion_tab.banner.text()


def test_changing_the_jump_rule_clears_the_result(solved):
    mt = solved.motion_tab
    mt.jump_rule.setCurrentIndex(mt.jump_rule.findData("published"))
    assert not solved.session.has(Stage.MOTION)
    assert "re-run" in mt.banner.text()


def test_the_jump_search_box_follows_the_rule(window):
    """The search budget belongs to the published rule only, so it starts off."""
    mt = window.motion_tab
    assert mt.jump_rule.currentData() == "nearest"
    assert not mt.njumps.isEnabled()
    mt.jump_rule.setCurrentIndex(mt.jump_rule.findData("published"))
    assert mt.njumps.isEnabled()
    mt.jump_rule.setCurrentIndex(mt.jump_rule.findData("nearest"))
    assert not mt.njumps.isEnabled()
    mt.jump_rule.setCurrentIndex(mt.jump_rule.findData("published"))
    assert mt.njumps.isEnabled()


def test_toggling_the_fit_clears_the_result(solved):
    solved.motion_tab.tfit.setChecked(not solved.motion_tab.tfit.isChecked())
    assert not solved.session.has(Stage.MOTION)


def test_flips_do_not_clear_the_result(solved):
    """Flips re-use the solution, so they must not throw it away."""
    solved.motion_tab.flip_x.setChecked(True)
    assert solved.session.has(Stage.MOTION)


def test_fit_order_box_follows_the_checkbox(window):
    mt = window.motion_tab
    mt.tfit.setChecked(True)
    assert mt.torder.isEnabled()
    mt.tfit.setChecked(False)
    assert not mt.torder.isEnabled()


def test_fit_order_changes_the_fitted_motion(window, dataset):
    window.filter_tab.run()
    window.segment_tab._load_mask_from_dataset()
    mt = window.motion_tab

    mt.tfit.setChecked(True)
    mt.torder.setValue(10)
    a = _solve(window, dataset).u_fitted
    mt.torder.setValue(6)
    b = _solve(window, dataset).u_fitted
    assert np.nanmax(np.abs(a - b)) > 1e-3, "the fit order had no effect"


def test_turning_the_fit_off_changes_what_is_used(window, dataset):
    window.filter_tab.run()
    window.segment_tab._load_mask_from_dataset()
    mt = window.motion_tab

    mt.tfit.setChecked(True)
    on = _solve(window, dataset)
    mt.tfit.setChecked(False)
    off = _solve(window, dataset)
    assert on.u_fitted is not None and off.u_fitted is None
    assert np.nanmax(np.abs(on.best - off.best)) > 1e-3


def test_jump_search_reports_whether_the_budget_binds(solved, dataset):
    mt = solved.motion_tab
    # the readout describes the published rule's search, so solve with that rule
    mt.jump_rule.setCurrentIndex(mt.jump_rule.findData("published"))
    mt._done(_solve(solved, dataset))
    mt.refresh()
    text = mt.jump_lbl.text()
    assert "search used" in text
    info = solved.session.result(Stage.MOTION).meta["jump_search"]
    used = int(np.abs(info["offset"]).max())
    assert f"used {used} of" in text
    if used < info["njumps"]:
        assert "will not change" in text


def test_the_budget_advice_is_for_the_published_rule_only(solved, dataset):
    """The nearest-coordinate rule has no budget to raise."""
    mt = solved.motion_tab
    mt.jump_rule.setCurrentIndex(mt.jump_rule.findData("nearest"))
    mt._done(_solve(solved, dataset))
    res = solved.session.result(Stage.MOTION)
    assert res.meta["phase_jump_rule"] == "nearest"
    text = mt.jump_lbl.text()
    assert "search used" in text
    assert "will not change" not in text


def test_settings_survive_running(solved, dataset):
    """Estimating must not reset the boxes to the suggested values."""
    mt = solved.motion_tab
    mt.a.setValue(55.0)
    mt.eta.setValue(3.5)
    mt.torder.setValue(6)
    mt.njumps.setValue(3)
    before = (mt.a.value(), mt.eta.value(), mt.torder.value(), mt.njumps.value())

    for _ in range(3):
        mt._done(_solve(solved, dataset))
        after = (mt.a.value(), mt.eta.value(), mt.torder.value(), mt.njumps.value())
        assert after == before, "a run reset the settings"


def test_settings_survive_any_stage_change(window, dataset):
    """Any stage change refreshes the headers, which must not touch the boxes."""
    window.filter_tab.run()
    window.segment_tab._load_mask_from_dataset()
    mt = window.motion_tab
    mt.a.setValue(42.0)
    mt.eta.setValue(7.0)

    window.segment_tab._load_mask_from_dataset()   # emits stageChanged
    window.filter_tab.run()
    assert mt.a.value() == 42.0
    assert mt.eta.value() == 7.0


def test_the_solver_uses_the_values_on_screen(solved, dataset):
    mt = solved.motion_tab
    mt.a.setValue(48.0)
    mt.eta.setValue(2.0)
    mt.njumps.setValue(4)
    mt.torder.setValue(8)
    mt._done(_solve(solved, dataset))

    meta = solved.session.result(Stage.MOTION).meta
    assert float(np.atleast_1d(meta["a"])[0]) == 48.0
    assert meta["eta"] == 2.0
    assert meta["njumps"] == 4
    assert meta["temporal_order"] == 8


def test_choosing_a_kernel_reloads_its_suggestions(window, dataset):
    """Picking a kernel is deliberate, so it may replace an edited value."""
    mt = window.motion_tab
    mt.a.setValue(99.0)
    mt.kernel.setCurrentText("gaussian")
    from harpi import get_kernel

    expected = get_kernel("gaussian").suggested_a_scale / dataset.pixels_per_wavelength
    assert mt.a.value() == pytest.approx(expected, rel=1e-6)
    assert mt.eta.value() == pytest.approx(get_kernel("gaussian").suggested_eta)


def test_an_edited_value_is_kept_when_the_dataset_reloads(window, dataset):
    mt = window.motion_tab
    mt.a.setValue(31.0)
    window.data_tab._adopt(dataset, "again.mat")
    assert mt.a.value() == 31.0


# --------------------------------------------------------- filter centre UI

def test_centres_are_detected_when_the_tab_is_shown(window):
    ft = window.filter_tab
    ft.refresh()
    assert ft.centres.any()
    assert ft.c_row.value() == pytest.approx(ft.centres[0][0])
    assert ft.c_col.value() == pytest.approx(ft.centres[0][1])


def test_dragging_a_marker_updates_the_boxes(window):
    ft = window.filter_tab
    ft.refresh()
    ft._centre_dragged(0, 123.5, 77.25)          # x is the column, y the row
    assert ft.c_row.value() == pytest.approx(77.25)
    assert ft.c_col.value() == pytest.approx(123.5)
    assert np.allclose(ft.centres[0], [77.25, 123.5])


def test_moving_a_centre_keeps_the_applied_filter(window):
    """The panels follow the new centre at once; what was applied stands
    until Apply is pressed again."""
    ft = window.filter_tab
    ft.refresh()
    ft.run()
    assert window.session.has(Stage.FILTER)
    ft.c_col.setValue(ft.c_col.value() + 2)
    assert window.session.has(Stage.FILTER)
    assert not ft.settings_match_result()
    assert "Apply" in ft.banner.text()


def test_dragging_keeps_the_applied_filter(window):
    ft = window.filter_tab
    ft.refresh()
    ft.run()
    ft._centre_dragged(0, ft.centres[0][1] + 4, ft.centres[0][0])
    assert window.session.has(Stage.FILTER), "the applied filter was dropped"
    assert not ft.settings_match_result()


def test_the_encoding_selector_shows_the_right_centre(window):
    ft = window.filter_tab
    ft.refresh()
    ft.encoding.setCurrentIndex(1)
    assert ft.c_row.value() == pytest.approx(ft.centres[1][0])
    assert ft.c_col.value() == pytest.approx(ft.centres[1][1])
    ft.encoding.setCurrentIndex(0)
    assert ft.c_row.value() == pytest.approx(ft.centres[0][0])


def test_a_centre_on_dc_is_refused(window):
    ft = window.filter_tab
    ft.refresh()
    rows, cols = window.session.dataset.image.shape[:2]
    ft._centre_dragged(0, cols / 2, rows / 2)
    ft.run()
    assert not window.session.has(Stage.FILTER)
    assert "DC" in ft.banner.text()


def test_reset_to_nominal_uses_the_tag_frequency(window, dataset):
    ft = window.filter_tab
    ft.refresh()
    ft._centre_dragged(0, 100, 100)
    ft.reset_centres()
    shape = dataset.image.shape[:2]
    expected = frequency_to_centre([dataset.ke[0], 0.0], shape)
    assert np.allclose(ft.centres[0], expected)


def test_both_markers_are_shown(window):
    ft = window.filter_tab
    ft.refresh()
    assert len(ft.k_before._markers) == 2
    assert ft.k_before.markers_draggable


def test_filtering_uses_the_moved_centre(window):
    ft = window.filter_tab
    ft.refresh()
    ft.run()
    before = np.abs(window.session.result(Stage.FILTER)["image"][:, :, 0, 5]).sum()
    ft._centre_dragged(0, ft.centres[0][1] + 6, ft.centres[0][0])
    ft.run()
    after = np.abs(window.session.result(Stage.FILTER)["image"][:, :, 0, 5]).sum()
    assert not np.isclose(before, after), "the moved centre was ignored"


# ------------------------------------------------- one segmentation canvas

def test_segmentation_shows_a_single_canvas(window):
    """Contours are placed by hand, so the panel gets the whole width."""
    st = window.segment_tab
    assert len(st.canvases) == 1
    assert st.canvases[0] is st.canvas


def test_the_single_canvas_is_the_editable_one(window):
    assert window.segment_tab.canvas.editable


def test_every_image_choice_can_be_drawn_on(window):
    """The selector replaces the read-only panels, so it must cover them."""
    st = window.segment_tab
    window.filter_tab.run()
    seen = set()
    for i in range(st.background.count()):
        st.background.setCurrentIndex(i)
        st.refresh()
        seen.add(st.title.text())
    # the two harmonic phases that used to have panels of their own
    assert any("phase 1" in t for t in seen)
    assert any("phase 2" in t for t in seen)
    assert len(seen) == st.background.count()


def test_the_second_raw_encoding_needs_no_filter(window):
    """It comes straight from the dataset, unlike the harmonic choices."""
    st = window.segment_tab
    st.background.setCurrentIndex(st.background.findData("magnitude2"))
    st.refresh()
    assert "encoding 2" in st.title.text()


def test_drawing_still_works_on_the_single_canvas(window):
    st = window.segment_tab
    c = st.canvas
    c.key = (0, 0)
    tri = np.array([[10.0, 10.0], [20.0, 10.0], [15.0, 20.0]])
    window.undo_stack.push(AddContourCommand(c, tri, "lv_epi"))
    assert len(window.session.contours[(0, 0)]) == 1
    window.undo_stack.undo()
    assert not window.session.contours.get((0, 0))


# ----------------------------------------------- tracking contours in the tab

def _draw_ring(window, frame=0, n=16):
    """A ring of nodes on one frame, as drawing produces."""
    ds = window.session.dataset
    yy, xx = np.nonzero(ds.mask[:, :, 0])
    cx, cy = xx.mean(), yy.mean()
    r = float(np.median(np.hypot(yy - cy, xx - cx)))
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ring = np.column_stack([cx + r * np.cos(t), cy + r * np.sin(t)])
    window.session.contours[(0, frame)] = [ring]
    window.session.labels[(0, frame)] = ["lv_epi"]
    window.session.current_frame = frame
    return ring


def test_tracking_needs_the_filter(window, dataset, monkeypatch):
    st = window.segment_tab
    _draw_ring(window)
    seen = {}
    monkeypatch.setattr("harpi.gui.tabs_segment.QMessageBox.information",
                        lambda *a, **k: seen.setdefault("msg", a[-1]))
    st._track_contours()
    assert "filter" in seen["msg"].lower()
    assert not window.session.contours.get((0, 1))


def test_tracking_fills_every_other_frame(window, dataset):
    window.filter_tab.run()
    st = window.segment_tab
    ring = _draw_ring(window)
    st._track_contours()
    nfr = dataset.n_frames
    for fr in range(nfr):
        assert window.session.contours.get((0, fr)), f"frame {fr} empty"
    assert np.array_equal(window.session.contours[(0, 0)][0], ring)


def test_tracked_contours_actually_move(window, dataset):
    """The point of it: copying would leave them identical."""
    window.filter_tab.run()
    st = window.segment_tab
    ring = _draw_ring(window)
    st._track_contours()
    moved = max(
        np.hypot(*(window.session.contours[(0, fr)][0] - ring).T).max()
        for fr in range(1, dataset.n_frames))
    assert moved > 1.0, "the contours did not follow the tissue"


def test_tracking_keeps_the_wall_label(window, dataset):
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    assert window.session.labels[(0, 3)] == ["lv_epi"]


def test_tracking_is_one_undo(window, dataset):
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    assert window.session.contours.get((0, 5))
    window.undo_stack.undo()
    assert not window.session.contours.get((0, 5))
    assert window.session.contours.get((0, 0)), "the drawn frame survives"


def test_tracking_respects_the_overwrite_checkbox(window, dataset, monkeypatch):
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    keep = np.array([[5.0, 5.0], [9.0, 5.0], [7.0, 9.0]])
    window.session.contours[(0, 4)] = [keep]
    window.session.labels[(0, 4)] = ["other"]

    st.prop_overwrite.setChecked(False)
    st._track_contours()
    assert np.array_equal(window.session.contours[(0, 4)][0], keep)

    st.prop_overwrite.setChecked(True)
    st._track_contours()
    assert not np.array_equal(window.session.contours[(0, 4)][0], keep)


def test_tracking_reports_nothing_to_track(window, dataset, monkeypatch):
    window.filter_tab.run()
    st = window.segment_tab
    window.session.current_frame = 2
    seen = {}
    monkeypatch.setattr("harpi.gui.tabs_segment.QMessageBox.information",
                        lambda *a, **k: seen.setdefault("msg", a[-1]))
    st._track_contours()
    assert "no contours" in seen["msg"].lower()


# ------------------------------------------- the window shown while tracking

class _FakeProgress:
    """Stands in for the dialog, and can report Cancel on a chosen call."""

    def __init__(self, total, cancel_after=None):
        self.total = total
        self.cancel_after = cancel_after
        self.values = []
        self.labels = []
        self.closed = False
        self._calls = 0

    def wasCanceled(self):
        self._calls += 1
        return (self.cancel_after is not None
                and self._calls > self.cancel_after)

    def setValue(self, v):
        self.values.append(v)

    def setLabelText(self, text):
        self.labels.append(text)

    def close(self):
        self.closed = True

    def deleteLater(self):
        self.deleted = True


def _fake_dialogs(st, monkeypatch, cancel_after=None):
    made = []

    def make(total):
        made.append(_FakeProgress(total, cancel_after))
        return made[-1]

    monkeypatch.setattr(st, "_make_progress", make)
    return made


def test_tracking_shows_a_progress_window(window, dataset, monkeypatch):
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    made = _fake_dialogs(st, monkeypatch)
    st._track_contours()
    assert len(made) == 1, "one window for the run"
    assert made[0].total == dataset.n_frames - 1, "one step per hop"


def test_the_progress_window_advances_to_the_end(window, dataset, monkeypatch):
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    made = _fake_dialogs(st, monkeypatch)
    st._track_contours()
    dialog = made[0]
    assert dialog.values == sorted(dialog.values), "progress must not go back"
    assert max(dialog.values) == dialog.total
    assert any("frame" in t for t in dialog.labels)


def test_the_progress_window_closes_when_the_run_ends(window, dataset,
                                                      monkeypatch):
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    made = _fake_dialogs(st, monkeypatch)
    st._track_contours()
    assert made[0].closed


def test_the_window_counts_every_contour(window, dataset, monkeypatch):
    window.filter_tab.run()
    st = window.segment_tab
    ring = _draw_ring(window)
    window.session.contours[(0, 0)] = [ring, ring * 0.9]
    window.session.labels[(0, 0)] = ["lv_epi", "lv_endo"]
    made = _fake_dialogs(st, monkeypatch)
    st._track_contours()
    assert made[0].total == 2 * (dataset.n_frames - 1)
    assert max(made[0].values) == made[0].total


def test_cancelling_changes_nothing(window, dataset, monkeypatch):
    """A half-tracked sequence is worse than none, so a cancel writes nothing."""
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    made = _fake_dialogs(st, monkeypatch, cancel_after=2)
    st._track_contours()
    assert made[0].closed
    assert not window.session.contours.get((0, 1))
    assert "cancel" in st.banner.text().lower()


def test_cancelling_leaves_the_undo_stack_alone(window, dataset, monkeypatch):
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    before = st.undo_stack.count()
    _fake_dialogs(st, monkeypatch, cancel_after=1)
    st._track_contours()
    assert st.undo_stack.count() == before


def test_the_window_closes_when_tracking_fails(window, dataset, monkeypatch):
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    made = _fake_dialogs(st, monkeypatch)
    monkeypatch.setattr("harpi.tracking.track_points_2d",
                        lambda *a, **k: (_ for _ in ()).throw(ValueError("boom")))
    monkeypatch.setattr("harpi.gui.tabs_segment.QMessageBox.warning",
                        lambda *a, **k: None)
    st._track_contours()
    assert made[0].closed
    assert not window.session.contours.get((0, 1))


def test_the_real_progress_window_is_modal_and_cancellable(window, dataset):
    """The fakes above bypass the widget, so build the real one once."""
    from PySide6.QtWidgets import QProgressDialog

    dialog = window.segment_tab._make_progress(10)
    try:
        assert isinstance(dialog, QProgressDialog)
        assert dialog.maximum() == 10
        assert dialog.windowModality() == Qt.WindowModal
        assert not dialog.wasCanceled()
        dialog.cancel()
        assert dialog.wasCanceled(), "the Cancel button must take effect"
    finally:
        dialog.close()


def test_the_window_stays_up_while_the_frames_are_written(window, dataset,
                                                          monkeypatch):
    """Writing is the slower half on a long sequence. Closing the window first
    leaves the interface looking frozen with nothing on screen to explain it.
    """
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    made = _fake_dialogs(st, monkeypatch)
    seen = {}
    real_push = st.undo_stack.push
    monkeypatch.setattr(st.undo_stack, "push",
                        lambda cmd: (seen.update(open=not made[0].closed),
                                     real_push(cmd))[1])
    st._track_contours()
    assert seen["open"], "the window had already closed"
    assert made[0].closed, "and it must close afterwards"
    assert any("Applying" in t for t in made[0].labels)


def test_writing_the_frames_refreshes_once(window, dataset, monkeypatch):
    """Per frame, the mask was rebuilt over the whole sequence: 17 frames cost
    170 rasterisations and five seconds."""
    from harpi.gui import tabs_segment

    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    counts = {"mask": 0, "refresh": 0}
    real_mask = tabs_segment.contours_to_mask
    monkeypatch.setattr(
        tabs_segment, "contours_to_mask",
        lambda *a, **k: (counts.update(mask=counts["mask"] + 1),
                         real_mask(*a, **k))[1])
    real_refresh = type(st.canvas).refresh
    monkeypatch.setattr(
        type(st.canvas), "refresh",
        lambda self, *a, **k: (counts.update(refresh=counts["refresh"] + 1),
                               real_refresh(self, *a, **k))[1])
    st._track_contours()
    assert counts["refresh"] <= 2, "one refresh for the whole run"
    assert counts["mask"] <= dataset.n_frames, "the mask is rebuilt once"


def test_undoing_a_tracked_run_also_refreshes_once(window, dataset, monkeypatch):
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    counts = {"refresh": 0}
    real_refresh = type(st.canvas).refresh
    monkeypatch.setattr(
        type(st.canvas), "refresh",
        lambda self, *a, **k: (counts.update(refresh=counts["refresh"] + 1),
                               real_refresh(self, *a, **k))[1])
    st.undo_stack.undo()
    assert counts["refresh"] <= 2
    assert not window.session.contours.get((0, 1))


def test_tracking_from_a_middle_frame_fills_both_sides(window, dataset):
    """The contour is drawn on whichever frame is clearest, which is rarely the
    first. The frames before it must be tracked too, not left empty."""
    window.filter_tab.run()
    st = window.segment_tab
    mid = dataset.n_frames // 2
    ring = _draw_ring(window, frame=mid)
    st._track_contours()
    for fr in range(dataset.n_frames):
        assert window.session.contours.get((0, fr)), f"frame {fr} empty"
    assert np.array_equal(window.session.contours[(0, mid)][0], ring)


def test_the_frames_before_the_drawn_one_actually_move(window, dataset):
    """Copying would leave them identical to the drawn contour."""
    window.filter_tab.run()
    st = window.segment_tab
    mid = dataset.n_frames // 2
    ring = _draw_ring(window, frame=mid)
    st._track_contours()
    for fr in range(mid):
        moved = np.hypot(*(window.session.contours[(0, fr)][0] - ring).T).max()
        assert moved > 0.1, f"frame {fr} was left where it was drawn"


def test_the_window_counts_both_directions(window, dataset, monkeypatch):
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window, frame=dataset.n_frames // 2)
    made = _fake_dialogs(st, monkeypatch)
    st._track_contours()
    assert made[0].total == dataset.n_frames - 1, "every hop, both ways"
    assert max(made[0].values) == made[0].total
    assert made[0].values == sorted(made[0].values)


# ------------------------------------------- editing after the whole sequence

def _count_rasters(st, monkeypatch):
    """How many frames get rasterised by the next rebuild."""
    from harpi.gui import tabs_segment

    calls = {"n": 0}
    real = tabs_segment.contours_to_mask
    monkeypatch.setattr(
        tabs_segment, "contours_to_mask",
        lambda *a, **k: (calls.__setitem__("n", calls["n"] + 1), real(*a, **k))[1])
    return calls


def _move_a_node(window, st, frame=0, by=(0.7, 0.3)):
    from harpi.gui.canvas import MovePointCommand

    key = (0, frame)
    window.session.current_frame = frame
    st.canvas.key = key
    old = np.array(window.session.contours[key][0][0])
    st.undo_stack.push(MovePointCommand(st.canvas, 0, 0, old, old + np.array(by)))
    return old


def test_editing_one_frame_rasterises_only_that_frame(window, dataset,
                                                      monkeypatch):
    """Tracking leaves every frame annotated, and rasterising all of them on
    each edit put half a second between dropping a node and the interface
    responding."""
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    calls = _count_rasters(st, monkeypatch)
    _move_a_node(window, st)
    assert calls["n"] == 1, f"{calls['n']} frames rasterised for a one frame edit"


def test_the_cached_mask_matches_one_built_from_scratch(window, dataset,
                                                        monkeypatch):
    """The check that matters: reusing frames must not change the answer."""
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    _move_a_node(window, st)
    cached = window.session.mask.copy()
    cached_structures = {k: v.copy()
                         for k, v in window.session.structure_masks.items()}

    st._mask_cache.clear()
    st._rebuild_mask()
    assert np.array_equal(cached, window.session.mask)
    assert set(cached_structures) == set(window.session.structure_masks)
    for name, m in cached_structures.items():
        assert np.array_equal(m, window.session.structure_masks[name]), name


def test_moving_a_node_changes_that_frames_mask(window, dataset):
    """A cache that never invalidated would pass every timing test and leave the
    mask showing the contour where it used to be."""
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    before = window.session.mask[:, :, 0].copy()
    _move_a_node(window, st, by=(4.0, 4.0))
    assert not np.array_equal(before, window.session.mask[:, :, 0])


def test_the_other_frames_keep_their_masks(window, dataset):
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    before = window.session.mask.copy()
    _move_a_node(window, st, frame=0, by=(4.0, 4.0))
    for fr in range(1, dataset.n_frames):
        assert np.array_equal(before[:, :, fr], window.session.mask[:, :, fr]), fr


def test_changing_the_resolution_redraws_every_frame(window, dataset,
                                                     monkeypatch):
    """It changes every rasterisation, so none of the cached frames survive."""
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    calls = _count_rasters(st, monkeypatch)
    st.resolution.setValue(st.resolution.value() - 0.1)
    assert calls["n"] == dataset.n_frames


def test_changing_the_origin_redraws_every_frame(window, dataset, monkeypatch):
    """The half pixel offset moves every mask, cached or not."""
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    calls = _count_rasters(st, monkeypatch)
    window.session.mask_origin = 0.0
    st._rebuild_mask()
    assert calls["n"] == dataset.n_frames


def test_undoing_an_edit_restores_the_mask(window, dataset):
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    before = window.session.mask.copy()
    _move_a_node(window, st, by=(4.0, 4.0))
    st.undo_stack.undo()
    assert np.array_equal(before, window.session.mask)


def _click_out(canvas, points, label="lv_epi"):
    """Place a contour node by node, as clicking does."""
    canvas.start_contour(label)
    for x, y in points:
        canvas._drawing.append([float(x), float(y)])
        canvas._undone_points.clear()
        canvas.session.contours[canvas.key][-1] = np.array(canvas._drawing)
        canvas.refresh()


def test_drawing_rasterises_nothing_until_the_contour_is_finished(
        window, dataset, monkeypatch):
    """Placing nodes only redraws the scene. Rebuilding the mask on each click
    would be the same trap as rebuilding it on each tracked frame."""
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    calls = _count_rasters(st, monkeypatch)
    t = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    _click_out(st.canvas, np.column_stack([20 + 5 * np.cos(t), 20 + 5 * np.sin(t)]))
    assert calls["n"] == 0, "the mask was rebuilt while drawing"


def test_finishing_a_contour_rasterises_only_its_frame(window, dataset,
                                                       monkeypatch):
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    t = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    _click_out(st.canvas, np.column_stack([20 + 5 * np.cos(t), 20 + 5 * np.sin(t)]))
    calls = _count_rasters(st, monkeypatch)
    assert st.canvas.finish_contour()
    assert calls["n"] == 1, f"{calls['n']} frames rasterised for one new contour"


def test_stepping_through_frames_rasterises_nothing(window, dataset,
                                                    monkeypatch):
    """Looking over what tracking produced is the common thing to do next."""
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    calls = _count_rasters(st, monkeypatch)
    for fr in range(dataset.n_frames):
        window._navigate(0, fr)
    assert calls["n"] == 0


def test_clearing_a_frame_is_not_an_error(window, dataset):
    """The button passed its text where the labels go, so it raised instead."""
    window.filter_tab.run()
    st = window.segment_tab
    _draw_ring(window)
    st._track_contours()
    window._navigate(0, 0)
    st._clear_frame()
    assert not window.session.contours.get((0, 0))
    assert window.session.contours.get((0, 1)), "only this frame"
    st.undo_stack.undo()
    assert window.session.contours.get((0, 0))


# ---------------------------------------------- what a project remembers

def _saved(window, tmp_path, name="g.harpi"):
    """Save, then load into a fresh window, as reopening does."""
    from harpi.gui.app import MainWindow

    p = tmp_path / name
    window.session.project_path = str(p)
    window.save_project()
    w2 = MainWindow()
    w2.load_project(p)
    return w2


def test_a_project_keeps_the_strain_geometry(window, tmp_path):
    """None of it was saved, so every reopen started from the defaults."""
    st = window.strain_tab
    st.origin_auto.setChecked(False)
    st.origin_x.setValue(23.5); st.origin_y.setValue(24.5)
    st.ins_auto.setChecked(False)
    st.ins_ant_x.setValue(16.0); st.ins_ant_y.setValue(10.0)
    st.ins_inf_x.setValue(10.0); st.ins_inf_y.setValue(34.0)
    st.level.setCurrentText("basal")

    back = _saved(window, tmp_path).strain_tab
    assert not back.origin_auto.isChecked()
    assert (back.origin_x.value(), back.origin_y.value()) == (23.5, 24.5)
    assert (back.ins_ant_x.value(), back.ins_ant_y.value()) == (16.0, 10.0)
    assert (back.ins_inf_x.value(), back.ins_inf_y.value()) == (10.0, 34.0)
    assert back.level.currentText() == "basal"


def test_a_project_without_a_saved_geometry_keeps_the_defaults(window,
                                                               tmp_path):
    """Every project written before this existed, which must not be read as
    a deliberate choice of anything."""
    from harpi.io import HarpiProject
    from harpi.gui.app import MainWindow

    p = HarpiProject(config={"segment": {"mask_origin": 0.5}})
    f = p.save(tmp_path / "old.harpi")
    w2 = MainWindow()
    w2.load_project(f)
    assert w2.strain_tab.ins_auto.isChecked()
    assert w2.strain_tab.origin_auto.isChecked()


def test_a_loaded_mask_survives_a_resolution_change(window, tmp_path):
    """The mask was saved but its source was not, so it read as drawn from
    contours and the next rebuild threw it away. Nudging contour resolution
    after opening a project destroyed it."""
    window.segment_tab._load_mask_from_dataset()
    assert window.session.mask is not None

    w2 = _saved(window, tmp_path, "m.harpi")
    assert w2.session.mask is not None, "the mask did not survive the load"
    assert w2.session.mask_source == "file"

    st = w2.segment_tab
    st.resolution.setValue(st.resolution.value() - 0.1)
    assert w2.session.mask is not None, "wiped by the resolution change"


# ------------------------------------------- the filter previews its controls

def test_the_panels_preview_before_anything_is_applied(window):
    """The shading used to appear only after Apply and vanish on the next
    edit, so it was absent exactly while the filter was being chosen."""
    ft = window.filter_tab
    ft.refresh()
    assert not window.session.has(Stage.FILTER)
    assert ft.k_before._mask is not None, "nothing shaded"
    assert not ft.k_after._pix.pixmap().isNull(), "the result panel is blank"
    assert not ft.mag._pix.pixmap().isNull()
    assert "preview" in ft.source_lbl.text()


def test_the_shading_is_the_response_not_a_threshold(window):
    ft = window.filter_tab
    ft.refresh()
    shaded = ft.k_before._mask
    assert shaded.dtype != bool, "a threshold again"
    assert np.array_equal(shaded, ft.preview()["response"])
    between = (shaded > 0.01) & (shaded < 0.99)
    assert between.any(), "no rolloff is visible, only a support"


def test_what_is_previewed_is_what_apply_stores(window):
    """Bitwise. _build_filter serves both, and filter_plane is the same
    arithmetic as filtering the series."""
    ft = window.filter_tab
    ft.refresh()
    before = ft.preview()["magnitude"]
    ft.run()
    res = window.session.result(Stage.FILTER)
    enc, fr = ft._current_encoding(), window.session.current_frame
    assert np.array_equal(before, np.abs(res["image"][:, :, enc, fr]))


def test_the_preview_follows_the_cutoff(window):
    ft = window.filter_tab
    ft.kind.setCurrentText("butterworth")
    ft.refresh()
    ft.cutoff.setValue(8)
    narrow = ft.preview()["response"].copy()
    ft.cutoff.setValue(40)
    wide = ft.preview()["response"]
    assert wide.sum() > narrow.sum(), "a wider cutoff passes no more"


def test_the_preview_follows_the_centre(window):
    ft = window.filter_tab
    ft.refresh()
    before = ft.preview()["response"].copy()
    ft.c_col.setValue(ft.c_col.value() + 6)
    assert not np.array_equal(before, ft.preview()["response"])


def test_the_preview_follows_the_frame(window, dataset):
    ft = window.filter_tab
    ft.refresh()
    first = ft.preview()["k_before"].copy()
    window._navigate(0, dataset.n_frames - 1)
    assert not np.array_equal(first, ft.preview()["k_before"])


def test_the_two_k_space_panels_share_a_range(window):
    """Giving the filtered panel its own percentiles renormalises it, so the
    suppression the panel exists to show is scaled away."""
    ft = window.filter_tab
    ft.refresh()
    pv = ft.preview()
    assert pv["k_after"].max() < pv["k_before"].max(), "nothing was suppressed"


def test_a_centre_on_dc_leaves_the_panels_alone(window, dataset):
    ft = window.filter_tab
    ft.refresh()
    before = ft.k_before._mask.copy()
    rows, cols = dataset.image.shape[:2]
    ft._centre_dragged(0, cols / 2, rows / 2)
    assert "DC" in ft.banner.text()
    assert np.array_equal(before, ft.k_before._mask), "the panels were redrawn"


def test_the_shading_slider_only_changes_the_strength(window):
    ft = window.filter_tab
    ft.run()
    ft.shading.setValue(90)
    assert ft.k_before._overlay.opacity() == pytest.approx(0.9)
    assert "90%" in ft.shading_lbl.text()
    assert ft.settings_match_result(), "a strength change is not a setting"


def test_cutoff_and_order_are_disabled_for_transmission(window):
    """Two of the three bandpass controls do nothing on this kind."""
    ft = window.filter_tab
    assert ft.kind.currentText() == "transmission"
    assert not ft.cutoff.isEnabled(), "enabled on a freshly built tab"
    assert not ft.order.isEnabled()
    assert "wave vector" in ft.kind_note.text()
    ft.kind.setCurrentText("butterworth")
    assert ft.cutoff.isEnabled() and ft.order.isEnabled()
    assert ft.kind_note.text() == ""
    ft.kind.setCurrentText("transmission")
    assert not ft.cutoff.isEnabled()


# ---------------------------------------------- following a marker as it moves

def _send(canvas, kind, x, y):
    """A real mouse event at a scene position, not a call to the handler."""
    from PySide6.QtCore import QPointF
    from PySide6.QtGui import QMouseEvent
    from PySide6.QtWidgets import QApplication as QA

    pos = canvas.mapFromScene(QPointF(float(x), float(y)))
    ev = QMouseEvent(kind, QPointF(pos), QPointF(pos), Qt.LeftButton,
                     Qt.LeftButton if kind != QEvent.MouseButtonRelease
                     else Qt.NoButton, Qt.NoModifier)
    QA.sendEvent(canvas.viewport(), ev)


def _ready_to_drag(window):
    ft = window.filter_tab
    ft.refresh()
    ft.k_before.resize(300, 300)
    ft.k_before.fit()
    return ft


def test_the_panels_follow_a_marker_before_it_is_dropped(window):
    ft = _ready_to_drag(window)
    row, col = ft.centres[0]
    before = ft.preview()["response"].copy()

    _send(ft.k_before, QEvent.MouseButtonPress, col, row)
    _send(ft.k_before, QEvent.MouseMove, col + 8, row)
    assert not np.array_equal(before, ft.preview()["response"]), \
        "the filter did not follow the marker"
    assert ft.c_col.value() == pytest.approx(col + 8, abs=1.0)
    _send(ft.k_before, QEvent.MouseButtonRelease, col + 8, row)


def test_dragging_does_not_rebuild_the_markers(window):
    """set_markers destroys every marker item, and the canvas holds an index
    into the list it is dragging. Rebuilding mid-drag is a dangling item."""
    ft = _ready_to_drag(window)
    row, col = ft.centres[0]
    _send(ft.k_before, QEvent.MouseButtonPress, col, row)
    held = ft.k_before._markers[0]
    for step in (3, 6, 9):
        _send(ft.k_before, QEvent.MouseMove, col + step, row)
        assert ft.k_before._markers[0] is held, "the markers were rebuilt"
    _send(ft.k_before, QEvent.MouseButtonRelease, col + 9, row)


def test_dragging_does_not_touch_the_stage_until_release(window):
    ft = _ready_to_drag(window)
    ft.run()
    row, col = ft.centres[0]
    seen = []
    window.session.previewChanged.connect(seen.append)

    _send(ft.k_before, QEvent.MouseButtonPress, col, row)
    _send(ft.k_before, QEvent.MouseMove, col + 5, row)
    _send(ft.k_before, QEvent.MouseMove, col + 9, row)
    assert not seen, "the session heard about an unfinished drag"

    _send(ft.k_before, QEvent.MouseButtonRelease, col + 9, row)
    assert seen, "the release was not reported"
    assert not ft.settings_match_result()


def test_the_magnitude_greys_hold_still_while_dragging(window):
    """Rescaling on every move makes the panel pump as the marker crosses
    the peak, which is the one thing the eye follows."""
    ft = _ready_to_drag(window)
    row, col = ft.centres[0]
    _send(ft.k_before, QEvent.MouseButtonPress, col, row)
    _send(ft.k_before, QEvent.MouseMove, col + 4, row)
    held = ft._mag_range
    _send(ft.k_before, QEvent.MouseMove, col + 12, row)
    assert ft._mag_range == held
    _send(ft.k_before, QEvent.MouseButtonRelease, col + 12, row)
    assert not ft._dragging


def test_a_new_dataset_is_not_shown_the_old_k_space(window, dataset):
    """The cache was keyed on id(), and a freed dataset can have its address
    reused, which would show the series before last."""
    ft = window.filter_tab
    ft.refresh()
    first = ft.preview()["k_before"].copy()

    from harpi import Dataset2D

    rng = np.random.default_rng(7)
    other = Dataset2D(image=rng.normal(size=dataset.image.shape),
                      mask=dataset.mask, pixel_size=[2.0, 2.0],
                      tag_spacing=dataset.tag_spacing)
    window.session.dataset = other
    assert not np.array_equal(first, ft.preview()["k_before"])


def test_the_cache_holds_the_dataset_it_belongs_to(window):
    ft = window.filter_tab
    ft.refresh()
    ft.preview()
    assert ft._k_cache[2] is window.session.dataset


def test_a_series_with_no_peak_is_refused_rather_than_looping(window):
    """_build_filter asks detect_centres for a centre, and detect_centres
    draws, which comes back here. It has to stop on the second miss."""
    ft = window.filter_tab
    ft.refresh()
    ft.centres[:] = 0.0
    ft.detect_centres = lambda: None          # a search that finds nothing
    with pytest.raises(ValueError, match="by hand"):
        ft._build_filter()


# ------------------------------------------- where the filter effectively ends

def test_an_edge_line_is_drawn_on_the_response(window):
    """The shading says how much is passed; the line says where that has
    fallen to almost nothing."""
    ft = window.filter_tab
    ft.refresh()
    assert ft.k_before._outlines, "no edge was drawn"


def test_the_edge_follows_its_level(window):
    """A lower level is further out, so it encloses more."""
    from harpi.gui.tabs_data import _response_edge

    ft = window.filter_tab
    ft.refresh()
    h = ft.preview()["response"]

    def extent(level):
        lines = _response_edge(h, level)
        assert lines, f"nothing at {level}"
        pts = np.vstack(lines)
        return np.ptp(pts[:, 0]) * np.ptp(pts[:, 1])

    assert extent(0.01) > extent(0.20)


def test_changing_the_level_redraws_the_edge(window):
    ft = window.filter_tab
    ft.refresh()
    before = [i.path().elementCount() for i in ft.k_before._outlines]
    ft.edge.setValue(25.0)
    after = [i.path().elementCount() for i in ft.k_before._outlines]
    assert before != after


def test_the_edge_level_is_not_a_filter_setting(window):
    """It is a thing drawn, so it cannot mark a result out of date."""
    ft = window.filter_tab
    ft.run()
    ft.edge.setValue(12.0)
    assert ft.settings_match_result()


def test_the_edge_sits_where_the_response_crosses(window):
    """Sampled back out of the response, the drawn line is at the level
    asked for.

    Except where it runs off the edge of k-space: a contour that leaves the
    domain is closed along its border, and those points are on the boundary
    rather than on the level.
    """
    from scipy.ndimage import map_coordinates

    ft = window.filter_tab
    ft.refresh()
    h = ft.preview()["response"]
    peak = float(np.nanmax(h))
    level = 0.05 * peak

    pts = []
    for item in ft.k_before._outlines:
        path = item.path()
        pts += [[path.elementAt(k).x, path.elementAt(k).y]
                for k in range(path.elementCount())]
    pts = np.array(pts)
    assert len(pts) > 50

    rows, cols = h.shape
    on_border = ((pts[:, 0] <= 0) | (pts[:, 0] >= cols - 1)
                 | (pts[:, 1] <= 0) | (pts[:, 1] >= rows - 1))
    inside = pts[~on_border]
    assert len(inside) > 0.9 * len(pts), "the contour is mostly on the border"

    got = map_coordinates(h, np.vstack([inside[:, 1], inside[:, 0]]), order=1)
    assert np.allclose(got, level, atol=0.02 * peak)


def test_a_flat_response_draws_no_edge():
    from harpi.gui.tabs_data import _response_edge

    assert _response_edge(np.zeros((8, 8)), 0.05) == []
