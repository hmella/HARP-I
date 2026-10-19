"""The Unwrapping tab is a preview, not a prerequisite.

The solver unwraps internally, so requiring this tab to run first forced the
same work twice. Worse, the two could disagree: the solver works on a
downsampled grid and seeds later frames from the reference, which the preview
does not, so the per-frame flags shown here could differ from the ones the solve
acted on.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication  # noqa: E402

from harpi import Dataset2D, MotionResult  # noqa: E402
from harpi.gui.session import Session, Stage  # noqa: E402
from harpi.gui.tabs_segment import UnwrapTab  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


ROWS = COLS = 16
NFR = 4


@pytest.fixture
def session(qapp):
    sess = Session()
    sess.dataset = Dataset2D(
        image=np.zeros((ROWS, COLS, 2, NFR)), mask=np.ones((ROWS, COLS, NFR)),
        pixel_size=[0.002, 0.002], tag_spacing=0.008)
    sess.mask = np.ones((ROWS, COLS, NFR), dtype=bool)
    sess.set_result(Stage.DATA, sess.dataset)
    sess.set_result(Stage.FILTER, {"image": np.zeros((ROWS, COLS, 2, NFR), complex)})
    return sess


@pytest.fixture
def tab(session):
    return UnwrapTab(session)


def _motion(nfr=NFR, complete=None):
    phase = np.zeros((ROWS, COLS, 2, nfr))
    for i in range(nfr):
        phase[:, :, 0, i] = i
        phase[:, :, 1, i] = -i
    return MotionResult(
        u=np.zeros((ROWS, COLS, 2, nfr)),
        unwrapped_phase=phase,
        quality=np.full((ROWS, COLS, nfr), 0.5),
        unwrap_complete=(np.ones(nfr, bool) if complete is None
                         else np.asarray(complete, bool)),
    )


# ------------------------------------------------------------- not a gate

def test_motion_can_run_without_the_preview(session):
    assert not session.has(Stage.UNWRAP)
    assert session.ready_for(Stage.MOTION)


def test_motion_still_needs_the_filter(session):
    session.invalidate_from(Stage.FILTER)
    assert not session.ready_for(Stage.MOTION)


def test_motion_still_needs_a_mask(session):
    session.mask = None
    assert not session.ready_for(Stage.MOTION)


def test_strain_still_needs_motion(session):
    assert not session.ready_for(Stage.STRAIN)
    session.set_result(Stage.MOTION, _motion())
    assert session.ready_for(Stage.STRAIN)


# --------------------------------------------------- which source is shown

def test_nothing_unwrapped_is_stated(tab):
    tab.refresh()
    assert "nothing unwrapped" in tab.source_lbl.text()


def test_the_preview_is_shown_when_there_is_no_solve(tab):
    tab.session.set_result(Stage.UNWRAP, {
        "phase": [[np.zeros((ROWS, COLS))] * 2] * NFR,
        "quality": [np.zeros((ROWS, COLS))] * NFR,
        "complete": np.ones(NFR, bool),
    })
    tab.refresh()
    assert "preview" in tab.source_lbl.text()


def test_the_solve_takes_precedence_over_the_preview(tab):
    tab.session.set_result(Stage.UNWRAP, {
        "phase": [[np.zeros((ROWS, COLS))] * 2] * NFR,
        "quality": [np.zeros((ROWS, COLS))] * NFR,
        "complete": np.ones(NFR, bool),
    })
    tab.session.set_result(Stage.MOTION, _motion())
    tab.refresh()
    assert "motion solve" in tab.source_lbl.text()


def test_the_flags_come_from_the_solve(tab):
    """The point of the change: what is listed is what the solve acted on."""
    tab.session.set_result(Stage.UNWRAP, {
        "phase": [[np.zeros((ROWS, COLS))] * 2] * NFR,
        "quality": [np.zeros((ROWS, COLS))] * NFR,
        "complete": np.ones(NFR, bool),          # preview says all complete
    })
    tab.session.set_result(Stage.MOTION, _motion(complete=[True, False, True, True]))
    tab.refresh()
    listed = [tab.status.item(i).text() for i in range(tab.status.count())]
    assert "INCOMPLETE" in listed[1]
    assert "complete" in listed[0]


def test_a_solve_without_stored_phase_falls_back_to_the_preview(tab):
    tab.session.set_result(Stage.UNWRAP, {
        "phase": [[np.zeros((ROWS, COLS))] * 2] * NFR,
        "quality": [np.zeros((ROWS, COLS))] * NFR,
        "complete": np.ones(NFR, bool),
    })
    tab.session.set_result(Stage.MOTION, MotionResult(u=np.zeros((ROWS, COLS, 2, NFR))))
    tab.refresh()
    assert "preview" in tab.source_lbl.text()


def test_a_frame_count_from_the_solve_is_respected(tab):
    """The solver may hold fewer frames than the preview listed."""
    tab.session.set_result(Stage.MOTION, _motion(nfr=2))
    tab.refresh()
    assert tab.status.count() == 2


def test_refresh_is_harmless_with_a_frame_past_the_end(tab):
    tab.session.set_result(Stage.MOTION, _motion(nfr=2))
    tab.session.current_frame = 7
    tab.refresh()          # must not raise
    assert tab.status.count() == 2
