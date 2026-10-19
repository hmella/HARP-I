"""Frame navigation by buttons rather than a slider.

The frame number is counted from zero everywhere in the arrays and in
``Navigator.frame.value()``, but shown from one, so the two conventions are
pinned here.
"""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication  # noqa: E402

from harpi.gui.widgets import FrameSpin, Navigator  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def nav(qapp):
    n = Navigator()
    n.set_extent(28)
    return n


# ------------------------------------------------------------- the frame box

def test_the_box_counts_from_zero_but_shows_from_one(qapp):
    box = FrameSpin()
    box.setMaximum(27)
    box.setValue(0)
    assert box.value() == 0
    assert box.text() == "1"
    box.setValue(9)
    assert box.value() == 9
    assert box.text() == "10"


def test_typing_a_frame_number_is_one_based(qapp):
    box = FrameSpin()
    box.setMaximum(27)
    assert box.valueFromText("10") == 9


def test_a_number_past_the_end_is_not_accepted(qapp):
    from PySide6.QtGui import QValidator

    box = FrameSpin()
    box.setMaximum(27)                      # 28 frames, so 1 to 28 on screen
    assert box.validate("28", 2)[0] == QValidator.Acceptable
    assert box.validate("29", 2)[0] == QValidator.Intermediate
    assert box.validate("x", 1)[0] == QValidator.Invalid


# --------------------------------------------------------------- the buttons

def test_there_is_no_slider_any_more(nav):
    from PySide6.QtWidgets import QSlider

    assert not nav.findChildren(QSlider)


def test_next_advances_one_frame(nav):
    nav.next.click()
    assert nav.frame.value() == 1
    nav.next.click()
    assert nav.frame.value() == 2


def test_previous_goes_back_one_frame(nav):
    nav.frame.setValue(5)
    nav.prev.click()
    assert nav.frame.value() == 4


def test_first_and_last_jump_to_the_ends(nav):
    nav.last.click()
    assert nav.frame.value() == 27
    nav.first.click()
    assert nav.frame.value() == 0


def test_stepping_wraps_round_the_cycle(nav):
    """A cine loops, so the ends join up."""
    nav.last.click()
    nav.next.click()
    assert nav.frame.value() == 0
    nav.prev.click()
    assert nav.frame.value() == 27


def test_the_label_follows_the_buttons(nav):
    assert nav.frame_label.text() == "frame 1/28"
    nav.next.click()
    assert nav.frame_label.text() == "frame 2/28"
    nav.last.click()
    assert nav.frame_label.text() == "frame 28/28"


def test_every_step_announces_the_change(nav):
    seen = []
    nav.changed.connect(lambda s, f: seen.append(f))
    nav.next.click()
    nav.next.click()
    nav.prev.click()
    assert seen == [1, 2, 1]


def test_the_buttons_repeat_when_held(nav):
    """Holding one should scrub, which is what the slider was used for."""
    for b in (nav.first, nav.prev, nav.next, nav.last):
        assert b.autoRepeat()


def test_a_single_frame_dataset_disables_the_controls(nav):
    nav.set_extent(1)
    for b in (nav.first, nav.prev, nav.next, nav.last, nav.play):
        assert not b.isEnabled()
    assert not nav.frame.isEnabled()


def test_the_controls_come_back_for_a_cine(nav):
    nav.set_extent(1)
    nav.set_extent(10)
    assert nav.next.isEnabled()
    assert nav.frame.isEnabled()
    assert nav.frame.maximum() == 9


def test_the_tooltips_name_the_real_shortcuts(nav):
    """The menu binds [ and ], not the arrow keys."""
    assert "[" in nav.prev.toolTip()
    assert "]" in nav.next.toolTip()
    assert "Home" in nav.first.toolTip()
    assert "End" in nav.last.toolTip()


def test_setting_the_value_directly_still_works(nav):
    """app.py connects session.frameRequested straight to frame.setValue."""
    nav.frame.setValue(13)
    assert nav.frame.value() == 13
    assert nav.frame_label.text() == "frame 14/28"


def test_the_slice_spinner_appears_only_when_there_are_slices(nav):
    assert not nav.slice.isVisible()
    nav.set_extent(10, 4)
    assert nav.slice.maximum() == 3


# ------------------------------------------------------------- the labels

def test_no_two_buttons_carry_the_same_label(nav):
    labels = [b.text() for b in (nav.first, nav.prev, nav.next, nav.last, nav.play)]
    assert len(set(labels)) == len(labels), labels


def test_first_is_not_the_mirror_of_next(nav):
    """|◀ and ▶| read as one control reversed, which is what this replaced."""
    assert nav.first.text() != nav.next.text()[::-1]
    assert nav.last.text() != nav.prev.text()[::-1]


def test_the_ends_are_double_triangles_and_the_steps_single(nav):
    assert nav.first.text().count("◀") == 2
    assert nav.last.text().count("▶") == 2
    assert nav.prev.text().count("◀") == 1
    assert nav.next.text().count("▶") == 1


def test_the_ends_carry_a_bar_and_the_steps_do_not(nav):
    assert "|" in nav.first.text() and "|" in nav.last.text()
    assert "|" not in nav.prev.text() and "|" not in nav.next.text()


def test_play_is_spelled_out_so_it_is_not_a_triangle(nav):
    """A triangle here would be indistinguishable from next frame."""
    assert nav.play.text() == "Play"
    assert "▶" not in nav.play.text()


def test_playing_relabels_the_button(qapp):
    """The relabelling lives in app.py, so call that rather than restate it."""
    from harpi.gui.app import MainWindow

    from PySide6.QtCore import QTimer

    class Holder:
        def __init__(self):
            self.nav = Navigator()
            self._timer = QTimer()

    h = Holder()
    MainWindow._toggle_play(h, True)
    assert h.nav.play.text() == "Pause"
    MainWindow._toggle_play(h, False)
    assert h.nav.play.text() == "Play"
    assert "▶" not in h.nav.play.text()
