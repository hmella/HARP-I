"""Restricting analysis to one outlined structure."""

from __future__ import annotations

import numpy as np
import pytest

from harpi.segmentation import Contour, contours_to_mask, structure_masks


def _ring(cx, cy, r, n=72):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([cx + r * np.cos(t), cy + r * np.sin(t)])


SHAPE = (48, 64)


def _two_ventricles():
    return [Contour(_ring(20, 24, 10), label="lv_epi"),
            Contour(_ring(20, 24, 5), label="lv_endo"),
            Contour(_ring(44, 24, 8), label="rv_epi")]


def test_structures_are_returned_separately():
    per = structure_masks(_two_ventricles(), SHAPE, 0.5, origin=0)
    assert sorted(per) == ["lv", "rv"]


def test_the_parts_reconstruct_the_union():
    cs = _two_ventricles()
    union = contours_to_mask(cs, SHAPE, 0.5, origin=0)
    per = structure_masks(cs, SHAPE, 0.5, origin=0)
    assert np.array_equal(per["lv"] | per["rv"], union)


def test_the_parts_do_not_overlap():
    per = structure_masks(_two_ventricles(), SHAPE, 0.5, origin=0)
    assert not (per["lv"] & per["rv"]).any()


def test_the_endocardial_hole_is_respected():
    per = structure_masks(_two_ventricles(), SHAPE, 0.5, origin=0)
    solid = structure_masks([Contour(_ring(20, 24, 10), label="lv_epi")],
                            SHAPE, 0.5, origin=0)
    assert per["lv"].sum() < solid["lv"].sum()


def test_unlabelled_contours_give_one_region_matching_the_union():
    cs = [_ring(20, 24, 10), _ring(20, 24, 5)]
    per = structure_masks(cs, SHAPE, 0.5, origin=0)
    assert list(per) == ["other"]
    assert np.array_equal(per["other"], contours_to_mask(cs, SHAPE, 0.5, origin=0))


def test_an_empty_list_is_refused():
    with pytest.raises(ValueError, match="at least one contour"):
        structure_masks([], SHAPE)


# --------------------------------------------------------------------- the tab

pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication  # noqa: E402

from harpi import MotionResult  # noqa: E402
from harpi.gui.session import Session, Stage  # noqa: E402
from harpi.gui.tabs_motion import StrainTab  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def tab(qapp):
    rows, cols = SHAPE
    nfr = 5
    per = structure_masks(_two_ventricles(), SHAPE, 0.5, origin=0)
    union = per["lv"] | per["rv"]

    yy, xx = np.mgrid[0:rows, 0:cols]
    u = np.zeros((rows, cols, 2, nfr))
    for k in range(nfr):
        s = -0.10 * np.sin(np.pi * k / (nfr - 1))
        u[:, :, 0, k] = s * (xx - 20.0)
        u[:, :, 1, k] = s * (yy - 24.0)

    sess = Session()
    sess.mask = np.repeat(union[:, :, None], nfr, 2)
    sess.structure_masks = {k: np.repeat(m[:, :, None], nfr, 2)
                            for k, m in per.items()}
    sess.set_result(Stage.MOTION, MotionResult(u=u))
    return StrainTab(sess)


def test_the_selector_lists_the_outlined_structures(tab):
    names = [tab.region.itemText(i) for i in range(tab.region.count())]
    assert names[0] == "all outlined tissue"
    assert "left ventricle" in names and "right ventricle" in names


def test_the_selector_is_disabled_when_nothing_is_labelled(qapp):
    sess = Session()
    sess.mask = np.ones((8, 8, 3), dtype=bool)
    t = StrainTab(sess)
    assert t.region.count() == 1
    assert not t.region.isEnabled()


def test_restricting_to_one_ventricle_shrinks_the_analysed_area(tab):
    tab.region.setCurrentText("all outlined tissue")
    tab.run()
    everything = int(tab.session.result(Stage.STRAIN)["strain"].mask.sum())
    tab.region.setCurrentText("left ventricle")
    tab.run()
    lv_only = int(tab.session.result(Stage.STRAIN)["strain"].mask.sum())
    assert 0 < lv_only < everything


def test_the_two_ventricles_analyse_different_pixels(tab):
    tab.region.setCurrentText("left ventricle")
    tab.run()
    lv = tab.session.result(Stage.STRAIN)["strain"].mask.copy()
    tab.region.setCurrentText("right ventricle")
    tab.run()
    rv = tab.session.result(Stage.STRAIN)["strain"].mask
    assert not (lv & rv).any()


def test_the_banner_names_the_region(tab):
    tab.region.setCurrentText("right ventricle")
    tab.run()
    assert "right ventricle" in tab.banner.text()


def test_changing_the_region_warns_that_a_rerun_is_needed(tab):
    tab.run()
    tab.region.setCurrentText("left ventricle")
    assert "again" in tab.banner.text().lower()


def test_an_empty_region_is_reported_rather_than_crashing(tab):
    nfr = tab.session.mask.shape[-1]
    tab.session.structure_masks["rv"] = np.zeros(SHAPE + (nfr,), dtype=bool)
    tab.region.setCurrentText("right ventricle")
    tab.run()
    assert "empty" in tab.banner.text().lower()
