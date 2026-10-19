"""The half-pixel offset between drawn contours and the rasterised mask.

The canvas draws the image with pixel ``i`` covering scene ``[i, i+1)``, so the
centre of that pixel is at ``i + 0.5``. Rasterising at integer coordinates
therefore tests pixel corners, and the resulting mask sits half a pixel off the
shape that was drawn.
"""

from __future__ import annotations

import numpy as np
import pytest

from harpi.io.project import SCHEMA_VERSION, migrate
from harpi.segmentation import Contour, contours_to_mask


def _circle(cx, cy, r, n=512):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([cx + r * np.cos(t), cy + r * np.sin(t)])


SHAPE = (40, 40)


def _centroid(mask):
    return np.argwhere(mask).mean(0)


def test_sampling_pixel_centres_is_symmetric_about_the_drawn_centre():
    """A circle centred on a pixel corner must straddle it evenly."""
    c = Contour(_circle(20.0, 20.0, 6.0))
    mask = contours_to_mask([c], SHAPE, 0.2, origin=0.5)
    row, col = _centroid(mask)
    # scene 20.0 is the corner between pixels 19 and 20
    assert row == pytest.approx(19.5, abs=1e-6)
    assert col == pytest.approx(19.5, abs=1e-6)


def test_sampling_corners_shifts_the_mask_by_half_a_pixel():
    c = Contour(_circle(20.0, 20.0, 6.0))
    corners = contours_to_mask([c], SHAPE, 0.2, origin=0)
    centres = contours_to_mask([c], SHAPE, 0.2, origin=0.5)
    assert _centroid(corners)[0] - _centroid(centres)[0] == pytest.approx(0.5, abs=1e-6)


def test_the_two_conventions_disagree_on_boundary_pixels():
    c = Contour(_circle(20.0, 20.0, 6.0))
    corners = contours_to_mask([c], SHAPE, 0.2, origin=0)
    centres = contours_to_mask([c], SHAPE, 0.2, origin=0.5)
    differ = int((corners ^ centres).sum())
    assert differ > 0
    # the disagreement is confined to the rim, not the interior
    assert differ < int((corners | centres).sum()) // 2


def test_pixel_centres_give_the_better_area():
    r = 6.0
    c = Contour(_circle(20.0, 20.0, r))
    truth = np.pi * r ** 2
    err = {o: abs(contours_to_mask([c], SHAPE, 0.2, origin=o).sum() - truth)
           for o in (0, 0.5)}
    assert err[0.5] <= err[0]


def test_a_mask_in_scene_coordinates_recovers_the_drawn_centre():
    """Index centroid plus the origin is where the contour actually was."""
    c = Contour(_circle(20.0, 20.0, 6.0))
    for origin in (0.0, 0.5):
        mask = contours_to_mask([c], SHAPE, 0.2, origin=origin)
        row, col = _centroid(mask)
        assert col + origin == pytest.approx(20.0, abs=0.1)


# ---------------------------------------------------------------- migration

def test_an_old_project_is_pinned_to_the_old_convention():
    meta = migrate({"schema_version": 1})
    assert meta["config"]["segment"]["mask_origin"] == 0.0
    assert meta["schema_version"] == SCHEMA_VERSION


def test_migration_keeps_the_rest_of_the_config():
    meta = migrate({"schema_version": 1,
                    "config": {"motion": {"eta": 15.0}}})
    assert meta["config"]["motion"]["eta"] == 15.0
    assert meta["config"]["segment"]["mask_origin"] == 0.0


def test_migration_does_not_override_an_explicit_value():
    meta = migrate({"schema_version": 1,
                    "config": {"segment": {"mask_origin": 0.25}}})
    assert meta["config"]["segment"]["mask_origin"] == 0.25


def test_a_current_project_is_left_alone():
    meta = {"schema_version": SCHEMA_VERSION, "config": {}}
    assert migrate(dict(meta)) == meta


# ------------------------------------------------------------------- the gui

pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication  # noqa: E402

from harpi import Dataset2D  # noqa: E402
from harpi.gui.session import Session  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def test_the_session_defaults_to_pixel_centres(qapp):
    assert Session().mask_origin == 0.5


@pytest.fixture
def window(qapp, tmp_path):
    from harpi.gui.app import MainWindow

    w = MainWindow()
    nfr = 3
    w.data_tab._adopt(Dataset2D(
        image=np.zeros((40, 40, 2, nfr)), mask=np.ones((40, 40, nfr)),
        pixel_size=[0.002, 0.002], tag_spacing=0.008), "synthetic")
    return w


def test_the_segment_tab_uses_the_session_convention(window):
    sess = window.session
    sess.contours[(0, 0)] = [_circle(20.0, 20.0, 8.0, 64)]
    sess.labels[(0, 0)] = ["lv_epi"]

    sess.mask_origin = 0.5
    window.segment_tab._rebuild_mask()
    centres = sess.mask[:, :, 0].copy()

    sess.mask_origin = 0.0
    window.segment_tab._rebuild_mask()
    corners = sess.mask[:, :, 0]

    assert not np.array_equal(centres, corners)
    assert _centroid(corners)[0] - _centroid(centres)[0] == pytest.approx(0.5, abs=0.05)


def test_an_older_convention_is_announced(window):
    sess = window.session
    sess.contours[(0, 0)] = [_circle(20.0, 20.0, 8.0, 64)]
    sess.labels[(0, 0)] = ["lv_epi"]
    sess.mask_origin = 0.0
    window.segment_tab._rebuild_mask()
    assert "older project" in window.segment_tab.source_lbl.text()
    sess.mask_origin = 0.5
    window.segment_tab._rebuild_mask()
    assert "older project" not in window.segment_tab.source_lbl.text()


def test_a_project_round_trip_keeps_the_convention(window, tmp_path):
    sess = window.session
    sess.contours[(0, 0)] = [_circle(20.0, 20.0, 8.0, 64)]
    sess.labels[(0, 0)] = ["lv_epi"]
    window.segment_tab._rebuild_mask()

    path = tmp_path / "p.harpi"
    # save_project takes an "ask" flag, not a path; passing one opens a dialog
    sess.project_path = str(path)
    window.save_project()
    sess.mask_origin = 0.0            # pretend something changed it
    window.load_project(path)
    assert window.session.mask_origin == 0.5
