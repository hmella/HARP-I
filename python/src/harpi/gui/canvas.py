"""Image canvas with interactive contour editing.

Built on QGraphicsView, which gives native hit testing and cursor handling
for dragging control points. Matplotlib is used for the plots elsewhere.

Supports undo and redo, deleting points, scroll wheel zoom, drag to pan,
keyboard frame and slice navigation, and a live mask overlay.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import (
    QBrush,
    QColor,
    QImage,
    QPainter,
    QPainterPath,
    QPen,
    QPixmap,
    QUndoCommand,
)
from PySide6.QtWidgets import (
    QGraphicsEllipseItem,
    QGraphicsPathItem,
    QGraphicsPixmapItem,
    QGraphicsScene,
    QGraphicsView,
)

from ..segmentation import WALLS, resample_contour

__all__ = ["ImageCanvas", "MovePointCommand", "AddPointCommand", "DeletePointCommand",
           "AddContourCommand", "DeleteContourCommand", "SetContoursCommand",
           "to_qimage"]

_COLOURS = ["#e8503a", "#2fb463", "#3a7fe8", "#e8c53a", "#a03ae8"]

# One colour per wall, so the structure is readable at a glance.
_WALL_COLOURS = {
    "lv_epi": "#e8503a",
    "lv_endo": "#2fb463",
    "rv_epi": "#3a7fe8",
    "rv_endo": "#e8c53a",
    "other": "#a03ae8",
}
_WALL_KEYS = [k for k, _ in WALLS]


def to_qimage(arr, vmin=None, vmax=None, cmap=None) -> QImage:
    """Grayscale (or colour-mapped) QImage from a 2-D array."""
    a = np.asarray(arr, dtype=float)
    finite = np.isfinite(a)
    if vmin is None:
        vmin = float(a[finite].min()) if finite.any() else 0.0
    if vmax is None:
        vmax = float(a[finite].max()) if finite.any() else 1.0
    if vmax <= vmin:
        vmax = vmin + 1.0
    norm = np.clip((a - vmin) / (vmax - vmin), 0, 1)
    norm[~finite] = 0.0

    if cmap is None:
        rgb = np.repeat((norm * 255).astype(np.uint8)[:, :, None], 3, axis=2)
    else:
        import matplotlib as mpl

        rgb = (mpl.colormaps[cmap](norm)[:, :, :3] * 255).astype(np.uint8)

    rgb = np.ascontiguousarray(rgb)
    h, w, _ = rgb.shape
    return QImage(rgb.data, w, h, 3 * w, QImage.Format_RGB888).copy()


# --------------------------------------------------------------- undo stack

class _ContourCommand(QUndoCommand):
    def __init__(self, canvas, text):
        super().__init__(text)
        self.canvas = canvas
        self.key = canvas.key


class MovePointCommand(_ContourCommand):
    def __init__(self, canvas, ci, pi, old, new):
        super().__init__(canvas, "move point")
        self.ci, self.pi, self.old, self.new = ci, pi, np.array(old), np.array(new)

    def undo(self):
        self.canvas.set_point(self.key, self.ci, self.pi, self.old)

    def redo(self):
        self.canvas.set_point(self.key, self.ci, self.pi, self.new)


class AddPointCommand(_ContourCommand):
    def __init__(self, canvas, ci, pi, point):
        super().__init__(canvas, "add point")
        self.ci, self.pi, self.point = ci, pi, np.array(point)

    def undo(self):
        self.canvas.remove_point(self.key, self.ci, self.pi)

    def redo(self):
        self.canvas.insert_point(self.key, self.ci, self.pi, self.point)


class DeletePointCommand(_ContourCommand):
    def __init__(self, canvas, ci, pi, point):
        super().__init__(canvas, "delete point")
        self.ci, self.pi, self.point = ci, pi, np.array(point)

    def undo(self):
        self.canvas.insert_point(self.key, self.ci, self.pi, self.point)

    def redo(self):
        self.canvas.remove_point(self.key, self.ci, self.pi)


class AddContourCommand(_ContourCommand):
    def __init__(self, canvas, points, label="other"):
        super().__init__(canvas, "add contour")
        self.points = np.asarray(points, dtype=float)
        self.label = label or "other"
        self.ci = None

    def undo(self):
        self.canvas.remove_contour(self.key, self.ci)

    def redo(self):
        self.ci = self.canvas.append_contour(self.key, self.points, self.label)


class DeleteContourCommand(_ContourCommand):
    def __init__(self, canvas, ci):
        super().__init__(canvas, "delete contour")
        self.ci = ci
        self.saved = None
        self.label = "other"

    def undo(self):
        if self.saved is not None:
            self.canvas.insert_contour(self.key, self.ci, self.saved, self.label)

    def redo(self):
        removed = self.canvas.remove_contour(self.key, self.ci)
        if removed is not None:
            self.saved, self.label = removed


class SetContoursCommand(QUndoCommand):
    """Replace several frames at once. Used by propagation."""

    def __init__(self, canvas, changes, labels=None, text="propagate contours"):
        super().__init__(text)
        self.canvas = canvas
        self.changes = changes                       # {key: [points, ...]}
        self.labels = labels or {}                   # {key: [label, ...]}
        self.previous = {k: canvas.contours_for(k) for k in changes}
        self.previous_labels = {k: list(canvas.labels_for(k)) for k in changes}

    def undo(self):
        for k, v in self.previous.items():
            self.canvas.replace_contours(k, v, self.previous_labels.get(k))

    def redo(self):
        for k, v in self.changes.items():
            self.canvas.replace_contours(
                k, [np.array(p, dtype=float) for p in v], self.labels.get(k)
            )


# -------------------------------------------------------------------- view

class ImageCanvas(QGraphicsView):
    """One linked image panel with editable contours."""

    contoursChanged = Signal()
    pointClicked = Signal(float, float)
    markerMoved = Signal(int, float, float)
    """A draggable marker was released at a new scene position."""

    HANDLE = 3.0

    def __init__(self, session, title="", parent=None):
        super().__init__(parent)
        self.session = session
        self.title = title
        self.key = (0, 0)
        self.editable = False
        self._drawing = None
        self._drag = None
        self._show_mask = True
        self._user_zoomed = False
        self.markers_draggable = False
        self._markers = []
        self._marker_drag = None

        self.setScene(QGraphicsScene(self))
        self.setRenderHint(QPainter.Antialiasing)
        self.setDragMode(QGraphicsView.NoDrag)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setBackgroundBrush(QBrush(QColor("#111")))
        self.setMinimumSize(220, 220)

        self._pix = QGraphicsPixmapItem()
        self.scene().addItem(self._pix)
        self._overlay = QGraphicsPixmapItem()
        self._overlay.setOpacity(0.35)
        self._overlay.setZValue(1)
        self.scene().addItem(self._overlay)
        self._items: list = []

    # ------------------------------------------------------------- display
    def set_image(self, arr, vmin=None, vmax=None, cmap=None):
        if arr is None:
            self._pix.setPixmap(QPixmap())
            return
        img = to_qimage(arr, vmin, vmax, cmap)
        self._pix.setPixmap(QPixmap.fromImage(img))
        rect = QRectF(0, 0, img.width(), img.height())
        if rect != self.scene().sceneRect():
            self.scene().setSceneRect(rect)
            self._user_zoomed = False

    def set_mask_overlay(self, mask):
        if mask is None or not self._show_mask:
            self._overlay.setPixmap(QPixmap())
            return
        m = np.asarray(mask, dtype=bool)
        rgba = np.zeros(m.shape + (4,), dtype=np.uint8)
        rgba[m] = (60, 180, 255, 255)
        rgba = np.ascontiguousarray(rgba)
        h, w, _ = rgba.shape
        img = QImage(rgba.data, w, h, 4 * w, QImage.Format_RGBA8888).copy()
        self._overlay.setPixmap(QPixmap.fromImage(img))

    def toggle_mask(self):
        self._show_mask = not self._show_mask
        self.refresh()

    def fit(self):
        if not self.scene().sceneRect().isEmpty():
            self.fitInView(self.scene().sceneRect(), Qt.KeepAspectRatio)
            self._user_zoomed = False

    def resizeEvent(self, event):
        super().resizeEvent(event)
        # Keep the image framed as the window changes, until the user zooms.
        if not self._user_zoomed:
            self.fit()

    def showEvent(self, event):
        super().showEvent(event)
        if not self._user_zoomed:
            self.fit()

    # ------------------------------------------------------------- markers
    def set_markers(self, markers, radius=4.0):
        """Show draggable points as ``(x, y, colour, label)`` tuples."""
        for item in self._markers:
            self.scene().removeItem(item)
        self._markers.clear()

        for i, m in enumerate(markers):
            x, y = float(m[0]), float(m[1])
            colour = QColor(m[2] if len(m) > 2 else "#e8503a")
            item = QGraphicsEllipseItem(
                QRectF(x - radius, y - radius, 2 * radius, 2 * radius)
            )
            item.setPen(QPen(colour, max(radius / 3.0, 0.6)))
            item.setBrush(QBrush(Qt.NoBrush))
            item.setZValue(5)
            item.setData(2, i)
            if len(m) > 3 and m[3]:
                item.setToolTip(str(m[3]))
            self.scene().addItem(item)
            self._markers.append(item)

    def move_marker(self, index, x, y, radius=4.0):
        if 0 <= index < len(self._markers):
            self._markers[index].setRect(
                QRectF(x - radius, y - radius, 2 * radius, 2 * radius)
            )

    def _marker_at(self, pos):
        for item in self.scene().items(pos):
            if item in self._markers:
                return item
        return None

    # ------------------------------------------------------------ contours
    def contours_for(self, key):
        return [np.array(p, dtype=float) for p in self.session.contours.get(key, [])]

    def labels_for(self, key):
        labels = self.session.labels.setdefault(key, [])
        n = len(self.session.contours.get(key, []))
        while len(labels) < n:
            labels.append("other")
        del labels[n:]
        return labels

    def colour_for(self, key, ci):
        labels = self.labels_for(key)
        lb = labels[ci] if ci < len(labels) else ""
        return QColor(_WALL_COLOURS.get(lb, _COLOURS[ci % len(_COLOURS)]))

    def replace_contours(self, key, value, labels=None):
        self.session.contours[key] = [np.array(p, dtype=float) for p in value]
        if labels is not None:
            self.session.labels[key] = list(labels)
        else:
            self.session.labels[key] = ["other"] * len(value)
        self.refresh()
        self.contoursChanged.emit()

    def append_contour(self, key, points, label=""):
        lst = self.session.contours.setdefault(key, [])
        lst.append(np.asarray(points, dtype=float))
        self.session.labels.setdefault(key, []).append(label or "other")
        self.refresh()
        self.contoursChanged.emit()
        return len(lst) - 1

    def insert_contour(self, key, ci, points, label="other"):
        self.session.contours.setdefault(key, []).insert(ci, np.asarray(points, float))
        self.session.labels.setdefault(key, []).insert(ci, label)
        self.refresh()
        self.contoursChanged.emit()

    def remove_contour(self, key, ci):
        lst = self.session.contours.get(key, [])
        if 0 <= ci < len(lst):
            out = lst.pop(ci)
            labels = self.session.labels.get(key, [])
            label = labels.pop(ci) if ci < len(labels) else "other"
            self.refresh()
            self.contoursChanged.emit()
            return out, label
        return None

    def set_point(self, key, ci, pi, value):
        self.session.contours[key][ci][pi] = value
        self.refresh()
        self.contoursChanged.emit()

    def insert_point(self, key, ci, pi, value):
        c = self.session.contours[key][ci]
        self.session.contours[key][ci] = np.insert(c, pi, value, axis=0)
        self.refresh()
        self.contoursChanged.emit()

    def remove_point(self, key, ci, pi):
        c = self.session.contours[key][ci]
        self.session.contours[key][ci] = np.delete(c, pi, axis=0)
        self.refresh()
        self.contoursChanged.emit()

    # -------------------------------------------------------------- render
    def refresh(self):
        for it in self._items:
            self.scene().removeItem(it)
        self._items.clear()

        contours = self.session.contours.get(self.key, [])
        self.labels_for(self.key)
        for ci, pts in enumerate(contours):
            colour = self.colour_for(self.key, ci)
            if len(pts) >= 3:
                dense = resample_contour(pts, closed=True, resolution=0.4)
                path = QPainterPath(QPointF(*dense[0]))
                for p in dense[1:]:
                    path.lineTo(QPointF(*p))
                item = QGraphicsPathItem(path)
                item.setPen(QPen(colour, 0.6))
                item.setZValue(2)
                self.scene().addItem(item)
                self._items.append(item)
            elif len(pts) == 2:
                path = QPainterPath(QPointF(*pts[0]))
                path.lineTo(QPointF(*pts[1]))
                item = QGraphicsPathItem(path)
                item.setPen(QPen(colour, 0.6, Qt.DashLine))
                item.setZValue(2)
                self.scene().addItem(item)
                self._items.append(item)

            if self.editable:
                for pi, p in enumerate(pts):
                    h = QGraphicsEllipseItem(
                        QRectF(p[0] - self.HANDLE / 2, p[1] - self.HANDLE / 2,
                               self.HANDLE, self.HANDLE)
                    )
                    h.setBrush(QBrush(colour))
                    h.setPen(QPen(Qt.black, 0.2))
                    h.setZValue(3)
                    h.setData(0, ci)
                    h.setData(1, pi)
                    self.scene().addItem(h)
                    self._items.append(h)

    # -------------------------------------------------------------- events
    def wheelEvent(self, event):
        factor = 1.15 if event.angleDelta().y() > 0 else 1 / 1.15
        self.scale(factor, factor)
        self._user_zoomed = True

    def mousePressEvent(self, event):
        pos = self.mapToScene(event.position().toPoint())
        if event.button() == Qt.MiddleButton:
            self.setDragMode(QGraphicsView.ScrollHandDrag)
            return super().mousePressEvent(event)

        if self.markers_draggable and event.button() == Qt.LeftButton:
            item = self._marker_at(pos)
            if item is not None:
                self._marker_drag = int(item.data(2))
                return

        if not self.editable:
            self.pointClicked.emit(pos.x(), pos.y())
            return super().mousePressEvent(event)

        item = self.scene().itemAt(pos, self.transform())
        if event.button() == Qt.LeftButton:
            if isinstance(item, QGraphicsEllipseItem) and item.data(0) is not None:
                ci, pi = item.data(0), item.data(1)
                self._drag = (ci, pi, np.array(self.session.contours[self.key][ci][pi]))
                return
            if self._drawing is not None:
                self._drawing.append([pos.x(), pos.y()])
                self.session.contours[self.key][-1] = np.array(self._drawing)
                self.refresh()
                return
        elif event.button() == Qt.RightButton:
            if isinstance(item, QGraphicsEllipseItem) and item.data(0) is not None:
                ci, pi = item.data(0), item.data(1)
                pts = self.session.contours[self.key][ci]
                if len(pts) > 3:
                    self.undo_stack.push(
                        DeletePointCommand(self, ci, pi, pts[pi].copy())
                    )
                return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._marker_drag is not None:
            pos = self.mapToScene(event.position().toPoint())
            self.move_marker(self._marker_drag, pos.x(), pos.y())
            return
        if self._drag is not None:
            pos = self.mapToScene(event.position().toPoint())
            ci, pi, _ = self._drag
            self.session.contours[self.key][ci][pi] = [pos.x(), pos.y()]
            self.refresh()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if self._marker_drag is not None:
            pos = self.mapToScene(event.position().toPoint())
            index, self._marker_drag = self._marker_drag, None
            self.markerMoved.emit(index, pos.x(), pos.y())
            return
        if self._drag is not None:
            ci, pi, old = self._drag
            new = np.array(self.session.contours[self.key][ci][pi])
            self._drag = None
            if not np.allclose(old, new):
                # Restore the old position so the command can reapply it.
                self.session.contours[self.key][ci][pi] = old
                self.undo_stack.push(MovePointCommand(self, ci, pi, old, new))
            return
        self.setDragMode(QGraphicsView.NoDrag)
        super().mouseReleaseEvent(event)

    # -------------------------------------------------------------- drawing
    def start_contour(self, label="other"):
        self._drawing = []
        self._drawing_label = label or "other"
        self.session.contours.setdefault(self.key, []).append(np.zeros((0, 2)))
        self.session.labels.setdefault(self.key, []).append(self._drawing_label)
        self.refresh()

    def finish_contour(self):
        """Commit the contour being drawn. Returns True if one was kept."""
        if self._drawing is None:
            return False
        pts = np.array(self._drawing, dtype=float)
        label = getattr(self, "_drawing_label", "other")
        lst = self.session.contours.setdefault(self.key, [])
        if lst:
            lst.pop()                     # drop the provisional entry
        labels = self.session.labels.setdefault(self.key, [])
        if labels:
            labels.pop()
        self._drawing = None
        if len(pts) >= 3:
            self.undo_stack.push(AddContourCommand(self, pts, label))
            return True
        self.refresh()
        return False

    def cancel_contour(self):
        if self._drawing is not None:
            lst = self.session.contours.get(self.key, [])
            if lst:
                lst.pop()
            labels = self.session.labels.get(self.key, [])
            if labels:
                labels.pop()
            self._drawing = None
            self.refresh()
