"""Segmentation and Unwrapping tabs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from ..segmentation import (WALLS, Contour, contours_to_mask,
                            structure_masks, wall_label)
from ..unwrap import phase_quality_2d, unwrap2
from .canvas import DeleteContourCommand, ImageCanvas, SetContoursCommand
from .session import Stage
from .widgets import StageBanner, nan_range

__all__ = ["SegmentTab", "UnwrapTab"]


def _read_mask(path):
    """Pull a boolean mask out of a .mat, .npy or .npz file.

    Takes the first array that looks like a mask: 2-D or 3-D, and either
    boolean or holding only zeros and ones.
    """
    path = Path(path)
    if path.suffix.lower() in (".npy", ".npz"):
        data = np.load(path, allow_pickle=False)
        if path.suffix.lower() == ".npy":
            return data
        candidates = {k: data[k] for k in data.files}
    else:
        import scipy.io

        raw = scipy.io.loadmat(path, struct_as_record=False, squeeze_me=True)
        candidates = {}
        for key, value in raw.items():
            if key.startswith("__"):
                continue
            if hasattr(value, "_fieldnames"):        # a struct, e.g. `data`
                for field in value._fieldnames:
                    candidates[f"{key}.{field}"] = getattr(value, field)
            else:
                candidates[key] = value

    def looks_like_mask(a):
        a = np.asarray(a)
        if a.ndim not in (2, 3) or a.size == 0 or not np.issubdtype(a.dtype, np.number) \
                and a.dtype != bool:
            return False
        if a.dtype == bool:
            return True
        vals = np.unique(a)
        return vals.size <= 2 and np.all(np.isin(vals, (0, 1)))

    # Prefer anything actually called "mask".
    for key in sorted(candidates, key=lambda k: (0 if "mask" in k.lower() else 1, k)):
        value = candidates[key]
        try:
            if looks_like_mask(value):
                return np.asarray(value).astype(bool)
        except Exception:
            continue
    raise ValueError("no 2-D or 3-D binary array found in this file")


class SegmentTab(QWidget):
    """Draw endocardial and epicardial contours.

    Three linked panels showing the magnitude and the two phase images.
    """

    def __init__(self, session, undo_stack, parent=None):
        super().__init__(parent)
        self.session = session
        self.undo_stack = undo_stack

        self.canvases = []
        self.titles = []
        grid = QHBoxLayout()
        for title in ("magnitude", "phase 1", "phase 2"):
            c = ImageCanvas(session, title)
            c.undo_stack = undo_stack
            c.contoursChanged.connect(self._contours_changed)
            self.canvases.append(c)
            col = QVBoxLayout()
            lbl = QLabel(title); lbl.setAlignment(Qt.AlignCenter)
            self.titles.append(lbl)
            col.addWidget(lbl); col.addWidget(c)
            w = QWidget(); w.setLayout(col)
            grid.addWidget(w, 1)
        self.canvases[0].editable = True

        # Which wall the next contour outlines.
        self.wall = QComboBox()
        for key, label in WALLS:
            self.wall.addItem(label, key)
        self.wall.setCurrentIndex(0)

        # What the editable panel shows. Drawing on the harmonic magnitude
        # is often easier than on the raw image, since the tag pattern is
        # suppressed and the wall stands out.
        self.background = QComboBox()
        self.background.addItem("Magnitude (raw)", "magnitude")
        self.background.addItem("Harmonic magnitude (filtered)", "harmonic")
        self.background.addItem("Harmonic phase 1", "phase1")
        self.background.addItem("Harmonic phase 2", "phase2")
        self.background.currentIndexChanged.connect(lambda _: self.refresh())

        self.draw_btn = QPushButton("Draw contour (Enter to finish)")
        self.draw_btn.setCheckable(True)
        self.draw_btn.clicked.connect(self._toggle_draw)
        self.del_btn = QPushButton("Delete selected contour")
        self.del_btn.clicked.connect(self._delete_contour)
        self.clear_btn = QPushButton("Clear this frame")
        self.clear_btn.clicked.connect(self._clear_frame)

        self.load_mask_btn = QPushButton("Use mask from dataset")
        self.load_mask_btn.clicked.connect(self._load_mask_from_dataset)
        self.load_mask_file_btn = QPushButton("Load mask from file\u2026")
        self.load_mask_file_btn.clicked.connect(self._load_mask_from_file)

        self.contour_list = QListWidget()
        self.contour_list.setMaximumHeight(110)

        self.prop_from = QSpinBox(); self.prop_from.setMinimum(1)
        self.prop_to = QSpinBox(); self.prop_to.setMinimum(1)
        self.prop_overwrite = QCheckBox("overwrite existing")
        self.prop_btn = QPushButton("Propagate…")
        self.prop_btn.clicked.connect(self._propagate)

        pform = QFormLayout()
        pform.addRow("from frame", self.prop_from)
        pform.addRow("through frame", self.prop_to)
        pform.addRow("", self.prop_overwrite)
        pform.addRow("", self.prop_btn)
        pbox = QGroupBox("Propagate"); pbox.setLayout(pform)

        self.resolution = QDoubleSpinBox()
        self.resolution.setRange(0.05, 2.0); self.resolution.setValue(0.5)
        self.resolution.setSingleStep(0.05)
        self.resolution.valueChanged.connect(self._contours_changed)

        self.coverage = QLabel("0 / 0 frames annotated")
        self.source_lbl = QLabel("mask: none")
        self.source_lbl.setWordWrap(True)
        self.banner = StageBanner(
            "Draw the outer (epicardial) contour first, then the inner ones.\n"
            "Click to place points, Enter to close, right click a point to "
            "delete it, Ctrl+Z to undo."
        )

        form = QFormLayout()
        form.addRow("wall", self.wall)
        form.addRow("draw on", self.background)
        form.addRow("contour resolution", self.resolution)

        mbox = QGroupBox("Mask")
        mlay = QVBoxLayout()
        mlay.addWidget(self.load_mask_btn)
        mlay.addWidget(self.load_mask_file_btn)
        mlay.addWidget(self.source_lbl)
        mbox.setLayout(mlay)

        side = QVBoxLayout()
        side.addLayout(form)
        side.addWidget(self.draw_btn)
        side.addWidget(self.del_btn)
        side.addWidget(self.clear_btn)
        side.addWidget(QLabel("contours on this frame:"))
        side.addWidget(self.contour_list)
        side.addWidget(mbox)
        side.addWidget(pbox)
        side.addWidget(self.coverage)
        side.addWidget(self.banner)
        side.addStretch(1)

        panel = QWidget(); panel.setLayout(side); panel.setMaximumWidth(300)
        lay = QHBoxLayout(self)
        lay.addWidget(panel)
        lay.addLayout(grid, 1)

    # -------------------------------------------------------------- actions
    def _toggle_draw(self, checked):
        c = self.canvases[0]
        if checked:
            c.start_contour(self.wall.currentData())
            self.banner.setText("Click to place points. Enter to close, Esc to cancel.")
        else:
            c.finish_contour()

    def finish_drawing(self):
        if self.draw_btn.isChecked():
            self.draw_btn.setChecked(False)
            self.canvases[0].finish_contour()

    def cancel_drawing(self):
        if self.draw_btn.isChecked():
            self.draw_btn.setChecked(False)
            self.canvases[0].cancel_contour()

    def _delete_contour(self):
        ci = self.contour_list.currentRow()
        if ci < 0:
            return
        self.undo_stack.push(DeleteContourCommand(self.canvases[0], ci))

    def _clear_frame(self):
        key = self.canvases[0].key
        if not self.session.contours.get(key):
            return
        self.undo_stack.push(
            SetContoursCommand(self.canvases[0], {key: []}, "clear frame")
        )

    def _adopt_mask(self, mask, origin):
        """Take a mask that was not drawn here."""
        ds = self.session.dataset
        mask = np.asarray(mask).astype(bool)
        if mask.ndim == 2:
            mask = np.repeat(mask[:, :, None], ds.n_frames, axis=2)
        if mask.shape[:2] != ds.image.shape[:2]:
            QMessageBox.warning(self, "Mask", (
                f"The mask is {mask.shape[:2]} but the images are "
                f"{ds.image.shape[:2]}."
            ))
            return
        if mask.shape[2] < ds.n_frames:
            pad = np.repeat(mask[:, :, -1:], ds.n_frames - mask.shape[2], axis=2)
            mask = np.concatenate([mask, pad], axis=2)
        elif mask.shape[2] > ds.n_frames:
            mask = mask[:, :, :ds.n_frames]

        self.session.mask = mask
        self.session.mask_source = "file"
        self.session.invalidate_from(Stage.UNWRAP)
        self.session.stageChanged.emit(int(Stage.SEGMENT))
        nframes = int((mask.reshape(-1, mask.shape[2]).any(axis=0)).sum())
        self.source_lbl.setText(f"mask: {origin} ({nframes} frames)")
        self.banner.setText(
            "Using a mask loaded from file. Drawing a contour replaces it "
            "with one built from the contours."
        )
        self.refresh()

    def _load_mask_from_dataset(self):
        ds = self.session.dataset
        if ds is None:
            self.banner.setText("Load a dataset first.")
            return
        if getattr(ds, "mask", None) is None or not np.any(ds.mask):
            self.banner.setText("This dataset has no mask stored in it.")
            return
        self._adopt_mask(ds.mask, "from the dataset")

    def _load_mask_from_file(self):
        ds = self.session.dataset
        if ds is None:
            self.banner.setText("Load a dataset first.")
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Load mask", "", "MATLAB (*.mat);;NumPy (*.npy *.npz);;All files (*)"
        )
        if not path:
            return
        try:
            mask = _read_mask(path)
        except Exception as exc:
            QMessageBox.warning(self, "Mask", f"Could not read a mask from\n{path}\n\n{exc}")
            return
        self._adopt_mask(mask, Path(path).name)

    def _propagate(self):
        ds = self.session.dataset
        if ds is None:
            return
        src = self.prop_from.value() - 1
        dst = self.prop_to.value() - 1
        sl = self.session.current_slice
        source = self.session.contours.get((sl, src))
        source_labels = list(self.session.labels.get((sl, src), []))
        if not source:
            QMessageBox.information(self, "Propagate",
                                    f"Frame {src + 1} has no contours to copy.")
            return
        lo, hi = (src, dst) if src <= dst else (dst, src)
        changes, labels = {}, {}
        for fr in range(lo, hi + 1):
            if fr == src:
                continue
            key = (sl, fr)
            if self.session.contours.get(key) and not self.prop_overwrite.isChecked():
                continue
            changes[key] = [np.array(p, dtype=float) for p in source]
            labels[key] = list(source_labels)
        if not changes:
            QMessageBox.information(self, "Propagate",
                                    "Nothing to do (tick 'overwrite existing'?).")
            return
        self.undo_stack.push(SetContoursCommand(self.canvases[0], changes, labels))

    # ------------------------------------------------------------- refresh
    def _contours_changed(self):
        self._rebuild_mask()
        self._update_list()
        self.refresh()

    def _rebuild_mask(self):
        ds = self.session.dataset
        if ds is None:
            return
        shape = ds.image.shape[:2]
        nfr = ds.n_frames
        mask = np.zeros(shape + (nfr,), dtype=bool)
        per_structure: dict[str, np.ndarray] = {}
        annotated = 0
        for fr in range(nfr):
            key = (self.session.current_slice, fr)
            pts = self.session.contours.get(key)
            if not pts:
                continue
            labels = self.session.labels.get(key, [])
            usable = [(p, labels[i] if i < len(labels) else "other")
                      for i, p in enumerate(pts) if len(p) >= 3]
            if not usable:
                continue
            cs = [Contour(p, label=lb) for p, lb in usable]
            origin = self.session.mask_origin
            try:
                mask[:, :, fr] = contours_to_mask(
                    cs, shape, self.resolution.value(), origin=origin,
                )
                annotated += 1
            except Exception:
                continue
            try:
                for key, m in structure_masks(
                        cs, shape, self.resolution.value(), origin=origin).items():
                    per_structure.setdefault(
                        key, np.zeros(shape + (nfr,), dtype=bool))[:, :, fr] = m
            except Exception:
                # the union is what the pipeline needs, so a failure here only
                # costs the regional breakdown
                pass

        self.coverage.setText(f"{annotated} / {nfr} frames annotated")

        if annotated:
            self.session.mask = mask
            self.session.structure_masks = per_structure
            self.session.mask_source = "contours"
            note = ("" if self.session.mask_origin == 0.5
                    else f", sampled at {self.session.mask_origin} for an older project")
            self.source_lbl.setText(f"mask: from the contours drawn here{note}")
        elif self.session.mask_source == "file":
            # A mask was loaded and nothing has been drawn over it yet.
            return
        else:
            self.session.mask = None
            self.session.structure_masks = {}
            self.session.mask_source = "none"
            self.source_lbl.setText("mask: none")

        self.session.invalidate_from(Stage.UNWRAP)
        # The mask is not a stage result, so announce it explicitly.
        self.session.stageChanged.emit(int(Stage.SEGMENT))

    def _update_list(self):
        key = self.canvases[0].key
        self.contour_list.clear()
        labels = self.canvases[0].labels_for(key)
        for i, p in enumerate(self.session.contours.get(key, [])):
            lb = labels[i] if i < len(labels) else "other"
            self.contour_list.addItem(f"{wall_label(lb)}, {len(p)} points")

    def refresh(self):
        ds = self.session.dataset
        if ds is None:
            return
        fr = min(self.session.current_frame, ds.n_frames - 1)
        key = (self.session.current_slice, fr)
        for c in self.canvases:
            c.key = key

        self.prop_from.setMaximum(ds.n_frames)
        self.prop_to.setMaximum(ds.n_frames)

        filt = self.session.result(Stage.FILTER)
        mag = np.real(ds.image[:, :, 0, fr])

        # What the editable panel shows.
        choice = self.background.currentData()
        if choice != "magnitude" and filt is None:
            choice = "magnitude"
            if self.background.currentIndex() != 0:
                self.banner.setText("Run the filter first to draw on the harmonic image.")
        if choice == "magnitude":
            self.canvases[0].set_image(mag, *nan_range(mag))
            self.titles[0].setText("magnitude")
        elif choice == "harmonic":
            h = np.abs(filt["image"][:, :, 0, fr])
            self.canvases[0].set_image(h, *nan_range(h))
            self.titles[0].setText("harmonic magnitude")
        else:
            d = 0 if choice == "phase1" else 1
            pha = np.angle(filt["image"][:, :, d, fr])
            self.canvases[0].set_image(pha, -np.pi, np.pi, cmap="twilight")
            self.titles[0].setText(f"harmonic phase {d + 1}")

        for i in (1, 2):
            if filt is None:
                self.canvases[i].set_image(np.zeros_like(mag), 0, 1)
            else:
                pha = np.angle(filt["image"][:, :, i - 1, fr])
                self.canvases[i].set_image(pha, -np.pi, np.pi, cmap="twilight")

        m = self.session.mask
        overlay = m[:, :, fr] if m is not None else None
        for c in self.canvases:
            c.set_mask_overlay(overlay)
            c.refresh()
            c.fit()
        self._update_list()


class UnwrapTab(QWidget):
    """Show the phase quality and the unwrapped phase, frame by frame."""

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session

        self.connectivity = QComboBox(); self.connectivity.addItems(["4", "8"])
        self.connectivity.currentTextChanged.connect(self._settings_changed)
        self.seed = QComboBox(); self.seed.addItems(["auto", "chosen point"])
        self.seed.currentTextChanged.connect(self._seed_mode)
        self.pick_seed = QPushButton("pick on the quality map")
        self.pick_seed.setCheckable(True)
        self.pick_seed.setEnabled(False)
        self.seed_lbl = QLabel("highest quality pixel")
        self.run_btn = QPushButton("Unwrap all frames")
        self.run_btn.clicked.connect(self.run)

        form = QFormLayout()
        form.addRow("connectivity", self.connectivity)
        form.addRow("seed", self.seed)
        form.addRow("", self.pick_seed)
        form.addRow("seed point", self.seed_lbl)
        box = QGroupBox("Unwrapping"); box.setLayout(form)
        box.setToolTip("These settings are used by the motion solver too.")

        self.status = QListWidget()
        self.banner = StageBanner("Segment the myocardium first.")

        side = QVBoxLayout()
        side.addWidget(box)
        side.addWidget(self.run_btn)
        side.addWidget(QLabel("per-frame result:"))
        side.addWidget(self.status)
        side.addWidget(self.banner)
        side.addStretch(1)

        self.qual = ImageCanvas(session, "quality")
        self.qual.pointClicked.connect(self._seed_picked)
        self.u1 = ImageCanvas(session, "unwrapped 1")
        self.u2 = ImageCanvas(session, "unwrapped 2")
        grid = QHBoxLayout()
        for c, lbl in ((self.qual, "phase quality"), (self.u1, "unwrapped 1"),
                       (self.u2, "unwrapped 2")):
            col = QVBoxLayout()
            t = QLabel(lbl); t.setAlignment(Qt.AlignCenter)
            col.addWidget(t); col.addWidget(c)
            w = QWidget(); w.setLayout(col)
            grid.addWidget(w, 1)

        panel = QWidget(); panel.setLayout(side); panel.setMaximumWidth(300)
        lay = QHBoxLayout(self)
        lay.addWidget(panel)
        lay.addLayout(grid, 1)

    # ----------------------------------------------------------- settings
    def _settings_changed(self, *_):
        self.session.unwrap["connectivity"] = int(self.connectivity.currentText())
        self.session.invalidate_from(Stage.MOTION)

    def _seed_mode(self, text):
        chosen = text == "chosen point"
        self.pick_seed.setEnabled(chosen)
        if not chosen:
            self.session.unwrap["seed"] = "auto"
            self.seed_lbl.setText("highest quality pixel")
            self.pick_seed.setChecked(False)
        elif self.session.unwrap["seed"] == "auto":
            self.seed_lbl.setText("click the quality map")
        self.session.invalidate_from(Stage.MOTION)

    def _seed_picked(self, x, y):
        """A seed is a (row, col) index, so the scene point is rounded to one."""
        if not self.pick_seed.isChecked():
            return
        row, col = int(round(y)), int(round(x))
        self.session.unwrap["seed"] = (row, col)
        self.seed_lbl.setText(f"row {row}, column {col}")
        self.pick_seed.setChecked(False)
        self.session.invalidate_from(Stage.MOTION)

    def run(self):
        filt = self.session.result(Stage.FILTER)
        mask = self.session.mask
        if filt is None:
            self.banner.setText("Run the filter first.")
            return
        if mask is None:
            self.banner.setText("Segment at least one frame first.")
            return

        Ih = filt["image"]
        nfr = Ih.shape[-1]
        conn = self.session.unwrap["connectivity"]
        seed = self.session.unwrap["seed"]
        out = {"phase": [], "quality": [], "complete": np.zeros(nfr, bool)}

        for fr in range(nfr):
            m = mask[:, :, fr]
            if not m.any():
                out["phase"].append(None); out["quality"].append(None)
                continue
            frames, ok = [], True
            q = phase_quality_2d(np.angle(Ih[:, :, 0, fr]), m, (1.0, 1.0), conn)
            for d in range(2):
                u, c = unwrap2(np.angle(Ih[:, :, d, fr]), m, (1.0, 1.0), seed,
                               conn, warn=False)
                frames.append(u); ok &= c
            out["phase"].append(frames)
            out["quality"].append(q)
            out["complete"][fr] = ok

        self.session.set_result(Stage.UNWRAP, out)
        bad = np.flatnonzero(~out["complete"] &
                             np.array([p is not None for p in out["phase"]]))
        self.banner.setText(
            "All frames unwrapped completely." if bad.size == 0 else
            f"Incomplete on frame(s) {[int(b) + 1 for b in bad]}. The mask is "
            "probably disconnected there."
        )
        self.refresh()

    def refresh(self):
        res = self.session.result(Stage.UNWRAP)
        self.status.clear()
        if res is None:
            return
        for fr, ph in enumerate(res["phase"]):
            if ph is None:
                self.status.addItem(f"frame {fr + 1}: not segmented")
            elif res["complete"][fr]:
                self.status.addItem(f"frame {fr + 1}: complete")
            else:
                self.status.addItem(f"frame {fr + 1}: INCOMPLETE")

        fr = min(self.session.current_frame, len(res["phase"]) - 1)
        q = res["quality"][fr]
        if q is not None:
            self.qual.set_image(q, 0, 1)
            self.qual.fit()
        ph = res["phase"][fr]
        if ph is not None:
            for c, p in ((self.u1, ph[0]), (self.u2, ph[1])):
                c.set_image(p, *nan_range(p))
                c.fit()
