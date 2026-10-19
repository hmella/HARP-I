"""Parameters that the tabs offered but never passed to the solver."""

from __future__ import annotations

from unittest import mock

import numpy as np
import pytest

pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication  # noqa: E402

from harpi import Dataset2D  # noqa: E402
from harpi.gui import tabs_motion  # noqa: E402
from harpi.gui.session import Session, Stage  # noqa: E402
from harpi.gui.tabs_motion import MotionTab  # noqa: E402
from harpi.gui.tabs_segment import UnwrapTab  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def session(qapp):
    nfr = 5
    sess = Session()
    sess.dataset = Dataset2D(
        image=np.zeros((16, 16, 2, nfr)), mask=np.ones((16, 16, nfr)),
        pixel_size=[0.002, 0.002], tag_spacing=0.008)
    sess.mask = np.ones((16, 16, nfr), dtype=bool)
    sess.set_result(Stage.FILTER, {"image": np.zeros((16, 16, 2, nfr), complex)})
    return sess


@pytest.fixture
def motion(session):
    return MotionTab(session)


@pytest.fixture
def unwrap(session):
    return UnwrapTab(session)


def _captured(tab):
    """Run the tab's solve with the solver stubbed, returning its kwargs."""
    seen = {}

    def fake(Ih, mask, **kw):
        seen.update(kw)
        raise RuntimeError("stop here, the arguments are what matter")

    with mock.patch.object(tabs_motion, "estimate_motion_2d", fake):
        tab.run()
        # the work runs on a thread; wait for it to finish
        if tab._thread is not None:
            tab._thread.quit()
            tab._thread.wait(5000)
    return seen


def test_the_reference_frame_reaches_the_solver(motion):
    motion.reference.setValue(3)
    seen = _captured(motion)
    # the control is one based, the library is zero based
    assert seen["reference"] == 2


def test_the_reference_frame_defaults_to_the_first(motion):
    seen = _captured(motion)
    assert seen["reference"] == 0


def test_the_reference_frame_is_capped_at_the_frame_count(motion):
    motion.sync_defaults()
    assert motion.reference.maximum() == motion.session.dataset.n_frames


def test_the_published_jump_rule_is_the_default(motion):
    seen = _captured(motion)
    assert seen["phase_jump_rule"] == "published"


def test_the_jump_rule_reaches_the_solver(motion):
    motion.jump_rule.setCurrentIndex(motion.jump_rule.findData("nearest"))
    seen = _captured(motion)
    assert seen["phase_jump_rule"] == "nearest"


def test_connectivity_reaches_the_solver(motion, unwrap):
    unwrap.connectivity.setCurrentText("8")
    seen = _captured(motion)
    assert seen["connectivity"] == 8


def test_the_seed_reaches_the_solver(motion, unwrap):
    unwrap.seed.setCurrentText("chosen point")
    unwrap.pick_seed.setChecked(True)
    unwrap._seed_picked(6.4, 9.6)          # scene coordinates
    seen = _captured(motion)
    assert seen["seed"] == (10, 6)          # rounded to a (row, col) index


def test_the_default_seed_is_automatic(motion):
    seen = _captured(motion)
    assert seen["seed"] == "auto"


def test_the_preview_and_the_solver_share_one_setting(unwrap):
    unwrap.connectivity.setCurrentText("8")
    assert unwrap.session.unwrap["connectivity"] == 8
    unwrap.connectivity.setCurrentText("4")
    assert unwrap.session.unwrap["connectivity"] == 4


def test_choosing_automatic_clears_a_picked_seed(unwrap):
    unwrap.seed.setCurrentText("chosen point")
    unwrap.pick_seed.setChecked(True)
    unwrap._seed_picked(3.0, 4.0)
    assert unwrap.session.unwrap["seed"] == (4, 3)
    unwrap.seed.setCurrentText("auto")
    assert unwrap.session.unwrap["seed"] == "auto"


def test_the_pick_button_only_works_in_chosen_mode(unwrap):
    assert not unwrap.pick_seed.isEnabled()
    unwrap.seed.setCurrentText("chosen point")
    assert unwrap.pick_seed.isEnabled()


def test_a_click_without_arming_is_ignored(unwrap):
    unwrap.seed.setCurrentText("chosen point")
    unwrap.pick_seed.setChecked(False)
    unwrap._seed_picked(3.0, 4.0)
    assert unwrap.session.unwrap["seed"] == "auto"


def test_changing_connectivity_invalidates_motion(unwrap):
    unwrap.session.set_result(Stage.MOTION, object())
    unwrap.connectivity.setCurrentText("8")
    assert not unwrap.session.has(Stage.MOTION)


def test_the_motion_tab_reports_the_shared_settings(motion, unwrap):
    unwrap.connectivity.setCurrentText("8")
    motion.sync_defaults()
    assert "8-connected" in motion.unwrap_lbl.text()
    unwrap.seed.setCurrentText("chosen point")
    unwrap.pick_seed.setChecked(True)
    unwrap._seed_picked(2.0, 5.0)
    motion.sync_defaults()
    assert "row 5" in motion.unwrap_lbl.text()
