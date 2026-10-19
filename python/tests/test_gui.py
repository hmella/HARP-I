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


def test_editing_a_filter_setting_clears_motion(window):
    window.filter_tab.run()
    window.session.set_result(Stage.MOTION, "m")
    window.filter_tab.cutoff.setValue(window.filter_tab.cutoff.value() + 1)
    assert not window.session.has(Stage.MOTION)


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


def test_moving_a_centre_invalidates_the_filter(window):
    ft = window.filter_tab
    ft.refresh()
    ft.run()
    assert window.session.has(Stage.FILTER)
    ft.c_col.setValue(ft.c_col.value() + 2)
    assert not window.session.has(Stage.FILTER)


def test_dragging_invalidates_the_filter(window):
    ft = window.filter_tab
    ft.refresh()
    ft.run()
    ft._centre_dragged(0, ft.centres[0][1] + 4, ft.centres[0][0])
    assert not window.session.has(Stage.FILTER)


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
