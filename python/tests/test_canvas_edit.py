"""Editing contours: undo while drawing, and hover highlighting."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("PySide6")

from PySide6.QtCore import QPointF, Qt  # noqa: E402
from PySide6.QtGui import QUndoStack  # noqa: E402
from PySide6.QtWidgets import (  # noqa: E402
    QApplication,
    QGraphicsEllipseItem,
    QGraphicsPathItem,
)

from harpi import Dataset2D  # noqa: E402
from harpi.gui.canvas import AddContourCommand, ImageCanvas  # noqa: E402
from harpi.gui.session import Session, Stage  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def canvas(qapp):
    sess = Session()
    sess.dataset = Dataset2D(image=np.zeros((40, 40, 2, 3)),
                             mask=np.ones((40, 40, 3)),
                             pixel_size=[2.0, 2.0], tag_spacing=8.0)
    sess.set_result(Stage.DATA, sess.dataset)
    c = ImageCanvas(sess, "magnitude")
    c.undo_stack = QUndoStack()
    c.editable = True
    c.key = (0, 0)
    c.set_image(np.zeros((40, 40)), 0, 1)
    return c


def _place(c, points, label="lv_epi"):
    """Draw a contour node by node, as clicking does."""
    c.start_contour(label)
    for x, y in points:
        c._drawing.append([float(x), float(y)])
        c._undone_points.clear()
        c.session.contours[c.key][-1] = np.array(c._drawing)
    c.refresh()


TRI = [(10, 10), (20, 10), (15, 20)]


# ------------------------------------------------------ undo while drawing

def test_undo_removes_the_last_node_placed(canvas):
    _place(canvas, TRI)
    assert len(canvas._drawing) == 3
    assert canvas.undo_last_point()
    assert len(canvas._drawing) == 2
    assert canvas.session.contours[canvas.key][-1].shape[0] == 2


def test_undo_steps_back_through_every_node(canvas):
    _place(canvas, TRI)
    for expected in (2, 1, 0):
        assert canvas.undo_last_point()
        assert len(canvas._drawing) == expected
    assert not canvas.undo_last_point(), "nothing left to take off"


def test_redo_puts_a_node_back(canvas):
    _place(canvas, TRI)
    canvas.undo_last_point()
    assert canvas.redo_last_point()
    assert len(canvas._drawing) == 3
    assert canvas._drawing[-1] == [15.0, 20.0]


def test_redo_does_nothing_when_no_node_was_removed(canvas):
    _place(canvas, TRI)
    assert not canvas.redo_last_point()


def test_placing_a_node_discards_the_redo_list(canvas):
    """Branching from an undone state drops the old branch, as editors do."""
    _place(canvas, TRI)
    canvas.undo_last_point()
    canvas._drawing.append([30.0, 30.0])
    canvas._undone_points.clear()
    assert not canvas.redo_last_point()


def test_undo_reports_nothing_to_do_when_not_drawing(canvas):
    """So the caller falls through to the committed history."""
    assert canvas._drawing is None
    assert not canvas.undo_last_point()
    assert not canvas.redo_last_point()


def test_starting_a_contour_clears_the_redo_list(canvas):
    _place(canvas, TRI)
    canvas.undo_last_point()
    canvas.start_contour("lv_endo")
    assert not canvas.redo_last_point()


def test_a_committed_contour_goes_on_the_stack_instead(canvas):
    """Once finished, undo is the ordinary one and removes the whole contour."""
    _place(canvas, TRI)
    assert canvas.finish_contour()
    assert canvas._drawing is None
    assert not canvas.undo_last_point(), "nothing in progress any more"
    assert canvas.undo_stack.canUndo()
    canvas.undo_stack.undo()
    assert not canvas.session.contours.get(canvas.key)


def test_cancelling_clears_the_nodes(canvas):
    _place(canvas, TRI)
    canvas.cancel_contour()
    assert canvas._drawing is None
    assert not canvas.redo_last_point()


# ----------------------------------------------------------------- hover

def _items(canvas):
    nodes = [i for i in canvas._items
             if isinstance(i, QGraphicsEllipseItem) and i.data(0) is not None]
    paths = [i for i in canvas._items if isinstance(i, QGraphicsPathItem)]
    return nodes, paths


def _committed(canvas):
    canvas.undo_stack.push(
        AddContourCommand(canvas, np.array(TRI, dtype=float), "lv_epi"))
    return _items(canvas)


def test_the_canvas_tracks_the_mouse(canvas):
    """Without this, hover only fires while a button is held."""
    assert canvas.hasMouseTracking()


def test_hovering_a_node_enlarges_it(canvas):
    nodes, _ = _committed(canvas)
    before = nodes[0].rect()
    canvas._set_hover(nodes[0])
    assert nodes[0].rect().width() > before.width()


def test_leaving_a_node_restores_it(canvas):
    nodes, _ = _committed(canvas)
    before = nodes[0].rect()
    canvas._set_hover(nodes[0])
    canvas._set_hover(None)
    assert nodes[0].rect() == before


def test_moving_between_nodes_restores_the_first(canvas):
    nodes, _ = _committed(canvas)
    first = nodes[0].rect()
    canvas._set_hover(nodes[0])
    canvas._set_hover(nodes[1])
    assert nodes[0].rect() == first
    assert nodes[1].rect().width() > first.width()


def test_hovering_a_contour_thickens_it(canvas):
    _, paths = _committed(canvas)
    before = paths[0].pen().widthF()
    canvas._set_hover(paths[0])
    assert paths[0].pen().widthF() > before


def test_leaving_a_contour_restores_its_width(canvas):
    _, paths = _committed(canvas)
    before = paths[0].pen().widthF()
    canvas._set_hover(paths[0])
    canvas._set_hover(None)
    assert paths[0].pen().widthF() == before


def test_hovering_the_same_item_twice_does_not_compound(canvas):
    """A second call must not thicken on top of the thickened pen."""
    _, paths = _committed(canvas)
    canvas._set_hover(paths[0])
    once = paths[0].pen().widthF()
    canvas._set_hover(paths[0])
    assert paths[0].pen().widthF() == once


def test_a_refresh_forgets_the_highlight(canvas):
    """The items are destroyed, so holding a reference would be a dangling one."""
    nodes, _ = _committed(canvas)
    canvas._set_hover(nodes[0])
    canvas.refresh()
    assert canvas._hover is None
    assert canvas._hover_style is None


def test_leaving_the_canvas_clears_the_highlight(canvas):
    nodes, _ = _committed(canvas)
    canvas._set_hover(nodes[0])
    canvas.leaveEvent(None)
    assert canvas._hover is None


# ------------------------------------------- undo routed from the window

def _holder(canvas):
    """Stand in for MainWindow: _undo touches only these two attributes."""
    from harpi.gui.app import MainWindow

    class Holder:
        pass

    h = Holder()
    h.segment_tab = type("T", (), {"canvas": canvas})()
    h.undo_stack = canvas.undo_stack
    return h, MainWindow


def test_undo_takes_a_node_before_touching_the_stack(canvas):
    """A drawn-but-uncommitted contour must not be skipped past."""
    canvas.undo_stack.push(
        AddContourCommand(canvas, np.array(TRI, dtype=float), "lv_epi"))
    _place(canvas, [(30, 30), (34, 30), (32, 34)])
    h, MainWindow = _holder(canvas)

    MainWindow._undo(h)
    assert len(canvas._drawing) == 2
    assert canvas.undo_stack.count() == 1, "the committed contour is untouched"


def test_undo_falls_through_to_the_stack_when_nothing_is_in_progress(canvas):
    canvas.undo_stack.push(
        AddContourCommand(canvas, np.array(TRI, dtype=float), "lv_epi"))
    h, MainWindow = _holder(canvas)

    assert canvas.session.contours.get(canvas.key)
    MainWindow._undo(h)
    assert not canvas.session.contours.get(canvas.key)


def test_redo_puts_the_node_back_before_touching_the_stack(canvas):
    _place(canvas, TRI)
    h, MainWindow = _holder(canvas)
    MainWindow._undo(h)
    assert len(canvas._drawing) == 2
    MainWindow._redo(h)
    assert len(canvas._drawing) == 3


def test_redo_falls_through_to_the_stack(canvas):
    canvas.undo_stack.push(
        AddContourCommand(canvas, np.array(TRI, dtype=float), "lv_epi"))
    h, MainWindow = _holder(canvas)
    MainWindow._undo(h)
    assert not canvas.session.contours.get(canvas.key)
    MainWindow._redo(h)
    assert len(canvas.session.contours[canvas.key]) == 1
