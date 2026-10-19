"""The HARP-I desktop application.

One window with six tabs, in the order the data flows:

    Data -> Filtering -> Segmentation -> Unwrapping -> Motion -> Strain

Segmentation comes after filtering because the borders are drawn on the
harmonic phase, which only exists once the filter has run.

Tab headers show the state of each stage, and changing a setting clears
every stage after it.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction, QKeySequence, QUndoStack
from PySide6.QtWidgets import (
    QApplication,
    QFileDialog,
    QLabel,
    QMainWindow,
    QMessageBox,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from .session import Session, Stage
from .tabs_compare import CompareTab
from .tabs_data import DataTab, FilterTab
from .tabs_motion import MotionTab, StrainTab
from .tabs_segment import SegmentTab, UnwrapTab
from .widgets import Navigator

__all__ = ["MainWindow", "main"]

_STATE_ICON = {"empty": "○", "done": "●", "stale": "◌"}


class MainWindow(QMainWindow):
    def __init__(self, project=None):
        super().__init__()
        self.setWindowTitle("HARP-I")
        self.resize(1400, 860)

        self.session = Session(self)
        self.undo_stack = QUndoStack(self)

        self.tabs = QTabWidget()
        self.data_tab = DataTab(self.session)
        self.filter_tab = FilterTab(self.session)
        self.segment_tab = SegmentTab(self.session, self.undo_stack)
        self.unwrap_tab = UnwrapTab(self.session)
        self.motion_tab = MotionTab(self.session)
        self.strain_tab = StrainTab(self.session)

        self._tabs = [
            (Stage.DATA, self.data_tab),
            (Stage.FILTER, self.filter_tab),
            (Stage.SEGMENT, self.segment_tab),
            (Stage.UNWRAP, self.unwrap_tab),
            (Stage.MOTION, self.motion_tab),
            (Stage.STRAIN, self.strain_tab),
        ]
        for stage, tab in self._tabs:
            self.tabs.addTab(tab, stage.title)

        # Compare is not a pipeline stage. It stays out of _tabs, which the
        # header logic walks assuming tab index equals stage, and it carries no
        # state icon.
        self.compare_tab = CompareTab(self.session)
        self.tabs.addTab(self.compare_tab, "Compare")
        self.session.snapshotsChanged.connect(self.compare_tab.sync)
        self.compare_tab.sync()

        self.nav = Navigator()
        self.nav.changed.connect(self._navigate)
        self.session.frameRequested.connect(self.nav.frame.setValue)
        self.nav.play.toggled.connect(self._toggle_play)
        self._timer = QTimer(self)
        self._timer.setInterval(120)
        self._timer.timeout.connect(lambda: self.nav.step_frame(1))

        central = QWidget()
        lay = QVBoxLayout(central)
        lay.addWidget(self.tabs, 1)
        lay.addWidget(self.nav)
        self.setCentralWidget(central)

        self.status = QLabel("Open a dataset to begin.")
        self.statusBar().addWidget(self.status)

        self._build_menus()
        self.session.stageChanged.connect(self._stage_changed)
        self.session.busyChanged.connect(self._busy_changed)
        self.tabs.currentChanged.connect(lambda _: self._refresh_current())

        if project:
            self.load_project(Path(project))
        self._update_headers()

    # --------------------------------------------------------------- menus
    def _build_menus(self):
        m = self.menuBar().addMenu("&File")
        for text, slot, key in (
            ("&Open dataset…", self.data_tab._open_mat, QKeySequence.Open),
            ("&Save project", self.save_project, QKeySequence.Save),
            ("Save project &as…", lambda: self.save_project(True), None),
            ("&Load project…", self.load_project_dialog, None),
        ):
            a = QAction(text, self)
            a.triggered.connect(lambda _=False, s=slot: s())
            if key:
                a.setShortcut(key)
            m.addAction(a)
        m.addSeparator()
        quit_action = QAction("&Quit", self)
        quit_action.setShortcut(QKeySequence.Quit)
        quit_action.triggered.connect(self.close)
        m.addAction(quit_action)

        e = self.menuBar().addMenu("&Edit")
        # Not createUndoAction: a contour being drawn is not on the stack yet,
        # so undo has to take its nodes off one at a time before falling
        # through to the committed history.
        undo = QAction("&Undo", self)
        undo.setShortcut(QKeySequence.Undo)
        undo.triggered.connect(self._undo)
        redo = QAction("&Redo", self)
        redo.setShortcut(QKeySequence.Redo)
        redo.triggered.connect(self._redo)
        e.addAction(undo)
        e.addAction(redo)

        v = self.menuBar().addMenu("&View")
        for text, slot, key in (
            ("Next frame", lambda: self.nav.step_frame(1), "]"),
            ("Previous frame", lambda: self.nav.step_frame(-1), "["),
            ("First frame", lambda: self.nav.frame.setValue(0), "Home"),
            ("Last frame",
             lambda: self.nav.frame.setValue(self.nav.n_frames - 1), "End"),
            ("Next slice", lambda: self.nav.step_slice(1), "."),
            ("Previous slice", lambda: self.nav.step_slice(-1), ","),
            ("Play / pause", lambda: self.nav.play.toggle(), "Space"),
            ("Toggle mask overlay", self._toggle_mask, "M"),
            ("Fit view", self._fit, "F"),
        ):
            a = QAction(text, self)
            a.setShortcut(QKeySequence(key))
            a.triggered.connect(slot)
            v.addAction(a)

        h = self.menuBar().addMenu("&Help")
        a = QAction("Keyboard shortcuts", self)
        a.triggered.connect(self._show_help)
        h.addAction(a)

    def _show_help(self):
        QMessageBox.information(self, "Shortcuts", (
            "[ ]      previous / next frame\n"
            "Home End first / last frame\n"
            ", .      previous / next slice\n"
            "Space    play / pause\n"
            "M        toggle mask overlay\n"
            "F        fit view\n"
            "Ctrl+Z   undo      Ctrl+Shift+Z  redo\n"
            "Ctrl+S   save project\n\n"
            "Segmentation tab:\n"
            "  click        place a contour point\n"
            "  Enter        close the contour\n"
            "  Esc          cancel drawing\n"
            "  right-click  delete a point\n"
            "  drag         move a point\n"
            "  Ctrl+Z       while drawing, take off the last point;\n"
            "               otherwise undo the last change\n"
            "  wheel        zoom      middle-drag  pan\n"
            "The point or contour under the cursor is highlighted."
        ))

    # ----------------------------------------------------------- navigation
    def _navigate(self, sl, fr):
        self.session.current_slice = sl
        self.session.current_frame = fr
        self._refresh_current()

    def _undo(self):
        """A node of the contour being drawn, else the last committed change."""
        if self.segment_tab.canvas.undo_last_point():
            return
        self.undo_stack.undo()

    def _redo(self):
        if self.segment_tab.canvas.redo_last_point():
            return
        self.undo_stack.redo()

    def _toggle_play(self, on):
        self.nav.play.setText("Pause" if on else "Play")
        self._timer.start() if on else self._timer.stop()

    def _toggle_mask(self):
        for c in self.segment_tab.canvases:
            c.toggle_mask()

    def _fit(self):
        w = self.tabs.currentWidget()
        for c in getattr(w, "canvases", []):
            c.fit()
        for name in ("canvas", "k_before", "k_after", "mag", "qual", "u1", "u2"):
            c = getattr(w, name, None)
            if c is not None:
                c.fit()

    def keyPressEvent(self, event):
        if self.tabs.currentWidget() is self.segment_tab:
            if event.key() in (Qt.Key_Return, Qt.Key_Enter):
                self.segment_tab.finish_drawing()
                return
            if event.key() == Qt.Key_Escape:
                self.segment_tab.cancel_drawing()
                return
        super().keyPressEvent(event)

    # --------------------------------------------------------------- state
    def _stage_changed(self, _stage):
        self._update_headers()
        self._refresh_current()
        if not self.session.busy:
            self.status.setText(self._summary())

    def _busy_changed(self, busy):
        self.tabs.setEnabled(not busy)
        self.status.setText("Working." if busy else self._summary())

    def _summary(self):
        ds = self.session.dataset
        if ds is None:
            return "Open a dataset to begin."
        done = [s.title for s in Stage if s != Stage.SEGMENT and self.session.has(s)]
        if self.session.mask is not None:
            done.append(Stage.SEGMENT.title)
        name = Path(ds.meta.get("source", "dataset")).name
        return f"{name}: {' | '.join(done) if done else 'loaded'}"

    def _update_headers(self):
        for i, (stage, _) in enumerate(self._tabs):
            if stage == Stage.SEGMENT:
                state = "done" if self.session.mask is not None else "empty"
            elif self.session.has(stage):
                state = "done"
            elif self.session.ready_for(stage):
                state = "empty"
            else:
                state = "stale"
            self.tabs.setTabText(i, f"{_STATE_ICON[state]}  {stage.title}")

        ds = self.session.dataset
        if ds is not None:
            n_slices = ds.image.shape[2] if ds.image.ndim == 5 else 1
            self.nav.set_extent(ds.n_frames, n_slices)
            self.motion_tab.sync_defaults()

    def _refresh_current(self):
        w = self.tabs.currentWidget()
        if hasattr(w, "refresh"):
            try:
                w.refresh()
            except Exception as exc:  # noqa: BLE001 - never kill the UI
                self.status.setText(f"refresh failed: {exc}")

    # ------------------------------------------------------------ projects
    def save_project(self, ask=False):
        from ..io import HarpiProject

        path = self.session.project_path
        if ask or path is None:
            path, _ = QFileDialog.getSaveFileName(self, "Save project", "session.harpi",
                                                  "HARP-I project (*.harpi)")
            if not path:
                return
        p = HarpiProject()
        ds = self.session.dataset
        if ds is not None and ds.meta.get("source"):
            src = Path(ds.meta["source"])
            if src.exists():
                p.set_source(src)
        p.contours = {
            f"{k[0]},{k[1]}": [
                {"points": np.asarray(c).tolist(),
                 "wall": (self.session.labels.get(k, []) + ["other"] * len(v))[i]}
                for i, c in enumerate(v)
            ]
            for k, v in self.session.contours.items() if v
        }
        if self.session.mask is not None:
            p.arrays["mask"] = self.session.mask
        p.config = {"segment": {"mask_origin": self.session.mask_origin}}
        p.save(path)
        self.session.project_path = path
        self.session.dirty = False
        self.status.setText(f"Saved {path}")

    def load_project_dialog(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load project", "",
                                              "HARP-I project (*.harpi)")
        if path:
            self.load_project(Path(path))

    def load_project(self, path):
        from ..io import HarpiProject

        p = HarpiProject.load(path)
        unchanged = p.source_unchanged()
        if p.source and unchanged is False:
            QMessageBox.warning(self, "Project", (
                f"The source data has changed since this project was saved:\n"
                f"{p.source}\nResults may not correspond."
            ))
        if p.source and Path(p.source).exists():
            from .. import load_mat_2d

            try:
                self.data_tab._adopt(load_mat_2d(p.source), p.source)
            except Exception:
                pass
        contours, labels = {}, {}
        for k, v in p.contours.items():
            key = tuple(int(x) for x in k.split(","))
            pts, walls = [], []
            for item in v:
                if isinstance(item, dict):          # current format
                    pts.append(np.asarray(item["points"], float))
                    walls.append(item.get("wall", "other"))
                else:                               # older projects
                    pts.append(np.asarray(item, float))
                    walls.append("other")
            contours[key] = pts
            labels[key] = walls
        self.session.contours = contours
        self.session.labels = labels
        # A project written before the sampling was corrected carries 0.0, put
        # there by the schema migration, so its mask is reproduced rather than
        # silently shifted.
        origin = ((p.config or {}).get("segment") or {}).get("mask_origin")
        self.session.mask_origin = 0.5 if origin is None else float(origin)
        if "mask" in p.arrays:
            self.session.mask = p.arrays["mask"].astype(bool)
        self.session.project_path = str(path)
        self.session.dirty = False
        self._update_headers()
        self._refresh_current()
        self.status.setText(f"Loaded {path}")

    def closeEvent(self, event):
        if not self.session.dirty:
            return super().closeEvent(event)
        r = QMessageBox.question(
            self, "Unsaved work", "Save this project before closing?",
            QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel,
        )
        if r == QMessageBox.Cancel:
            event.ignore()
            return
        if r == QMessageBox.Save:
            self.save_project()
        super().closeEvent(event)


def main(project=None) -> int:
    app = QApplication.instance() or QApplication(sys.argv)
    win = MainWindow(project)
    win.show()
    return app.exec()


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else None))
