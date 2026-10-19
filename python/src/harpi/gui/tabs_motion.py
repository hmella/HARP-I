"""Motion and Strain tabs."""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from .. import (
    KERNELS,
    aha_labels,
    bullseye_rings,
    estimate_motion_2d,
    flip_components,
    get_kernel,
    pixel_strain_2d,
    segment_metrics,
    segment_table,
    strain_by_segments,
)
from ..analysis import NATIVE_SEGMENTS
from .session import Stage
from .snapshots import Snapshot
from .widgets import MplCanvas, StageBanner, mask_bbox, nan_range
from .worker import run_in_background

__all__ = ["MotionTab", "StrainTab"]

_REGION_NAMES = {"lv": "left ventricle", "rv": "right ventricle",
                 "other": "other"}


class MotionTab(QWidget):
    """Set up and run the RBF motion estimation."""

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session
        self._thread = None
        self._worker = None
        self._raw = None          # result as solved, before any flip
        self._syncing = False     # applying defaults, not a user edit
        self._edited = set()      # settings the user has changed by hand

        self.kernel = QComboBox(); self.kernel.addItems(sorted(KERNELS))
        self.kernel.setCurrentText("multiquadric3")
        self.kernel.currentTextChanged.connect(self._kernel_changed)

        self.a = QDoubleSpinBox(); self.a.setRange(0.001, 1e5); self.a.setDecimals(4)
        self.eta = QDoubleSpinBox(); self.eta.setRange(0.0, 1e6); self.eta.setDecimals(8)
        self.a_dist = QComboBox()
        self.a_dist.addItems(["constant", "linear", "decreasinglinear", "quadratic"])
        self.downsample = QSpinBox(); self.downsample.setRange(1, 8)
        self.smooth = QCheckBox("smooth reference phase")
        self.tfit = QCheckBox("temporal polynomial fit"); self.tfit.setChecked(True)
        self.torder = QSpinBox(); self.torder.setRange(2, 20); self.torder.setValue(10)
        self.njumps = QSpinBox(); self.njumps.setRange(1, 16); self.njumps.setValue(8)
        self.jump_rule = QComboBox()
        self.jump_rule.addItem("published", "published")
        self.jump_rule.addItem("nearest coordinate", "nearest")
        self.jump_rule.setToolTip(
            "How the 2π offset of each frame is chosen. The published rule "
            "compares the largest distance between virtual coordinates and "
            "searches the steps set in 'jump search'. The nearest-coordinate "
            "rule compares every coordinate with the nearest reference one.")
        self.reference = QSpinBox(); self.reference.setMinimum(1); self.reference.setValue(1)
        self.reference.setToolTip(
            "Frame every displacement is measured from, usually end diastole.")
        self.unwrap_lbl = QLabel("-")
        self.unwrap_lbl.setToolTip("Set on the Unwrapping tab.")

        # Changing any of these means the current solution no longer matches
        # the settings, so drop it rather than leaving a stale result on
        # screen looking as though the control did nothing.
        # Mark the value as edited before anything else reacts: invalidating
        # refreshes the headers, which reapplies the defaults, and that would
        # overwrite the new value if the flag were not set yet.
        self.a.valueChanged.connect(lambda _: self._mark_edited("a"))
        self.eta.valueChanged.connect(lambda _: self._mark_edited("eta"))
        for w in (self.a, self.eta, self.downsample, self.njumps, self.torder,
                  self.reference):
            w.valueChanged.connect(self._settings_changed)
        for w in (self.kernel, self.a_dist, self.jump_rule):
            w.currentIndexChanged.connect(self._settings_changed)
        for w in (self.smooth, self.tfit):
            w.toggled.connect(self._settings_changed)
        self.tfit.toggled.connect(self.torder.setEnabled)
        self.torder.setEnabled(self.tfit.isChecked())
        # Only the published rule has a search budget.
        self.jump_rule.currentIndexChanged.connect(
            lambda _: self.njumps.setEnabled(self.jump_rule.currentData() == "published"))

        form = QFormLayout()
        form.addRow("kernel", self.kernel)
        form.addRow("shape a", self.a)
        form.addRow("a schedule", self.a_dist)
        form.addRow("eta", self.eta)
        form.addRow("downsample", self.downsample)
        form.addRow("jump rule", self.jump_rule)
        form.addRow("jump search", self.njumps)
        form.addRow("reference frame", self.reference)
        form.addRow("unwrapping", self.unwrap_lbl)
        form.addRow("", self.smooth)
        form.addRow("", self.tfit)
        form.addRow("fit order", self.torder)
        box = QGroupBox("RBF interpolation"); box.setLayout(form)

        # The sign and order of the components follow the tag directions and
        # the row and column convention of the images. These correct a
        # mismatch without solving again.
        self.swap_xy = QCheckBox("swap x and y")
        self.flip_x = QCheckBox("flip x")
        self.flip_y = QCheckBox("flip y")
        for cb in (self.swap_xy, self.flip_x, self.flip_y):
            cb.toggled.connect(self._apply_flips)
        fbox = QGroupBox("Displacement components")
        flay = QVBoxLayout()
        for cb in (self.swap_xy, self.flip_x, self.flip_y):
            flay.addWidget(cb)
        fbox.setLayout(flay)
        self.flip_box = fbox

        self.run_btn = QPushButton("Estimate motion")
        self.run_btn.clicked.connect(self.run)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self._cancel)
        self.cancel_btn.setEnabled(False)
        self.progress = QProgressBar()

        self.jump_lbl = QLabel("-")
        self.jump_lbl.setWordWrap(True)
        self.banner = StageBanner("Unwrap the phase first.")

        side = QVBoxLayout()
        side.addWidget(box)
        side.addWidget(self.run_btn)
        side.addWidget(self.cancel_btn)
        side.addWidget(self.progress)
        side.addWidget(self.flip_box)
        side.addWidget(QLabel("phase jumps (x2π), per frame:"))
        side.addWidget(self.jump_lbl)
        side.addWidget(self.banner)
        side.addStretch(1)

        self.plot = MplCanvas(7, 5)
        panel = QWidget(); panel.setLayout(side); panel.setMaximumWidth(320)
        lay = QHBoxLayout(self)
        lay.addWidget(panel)
        lay.addWidget(self.plot, 1)

    def _mark_edited(self, name):
        """Remember that a value was set by hand, so defaults leave it alone."""
        if not self._syncing:
            self._edited.add(name)

    def _settings_changed(self, *_):
        """A solver setting changed, so the current solution is out of date."""
        if self._syncing:
            return
        if self.session.has(Stage.MOTION):
            self._raw = None
            self.session.invalidate_from(Stage.MOTION)
            self.banner.setText("Settings changed. Press Estimate motion to re-run.")

    def _flipped(self, res):
        return flip_components(
            res,
            swap_xy=self.swap_xy.isChecked(),
            flip_x=self.flip_x.isChecked(),
            flip_y=self.flip_y.isChecked(),
        )

    def _apply_flips(self):
        """Re-apply the component options to the existing solution."""
        if self._raw is None:
            return
        self.session.set_result(Stage.MOTION, self._flipped(self._raw))
        parts = [t for t, c in (("x and y swapped", self.swap_xy),
                                ("x flipped", self.flip_x),
                                ("y flipped", self.flip_y)) if c.isChecked()]
        self.banner.setText(
            "Components: " + (", ".join(parts) if parts else "unchanged")
        )
        self.refresh()

    def _apply_kernel_defaults(self, force):
        """Fill in the suggested shape and eta for the current kernel.

        With ``force`` false, anything the user has typed is left alone.
        """
        k = get_kernel(self.kernel.currentText())
        ds = self.session.dataset
        self._syncing = True
        try:
            if (force or "a" not in self._edited) \
                    and k.suggested_a_scale is not None and ds is not None:
                self.a.setValue(k.suggested_a_scale / ds.pixels_per_wavelength)
            if force or "eta" not in self._edited:
                self.eta.setValue(k.suggested_eta)
            self.a.setEnabled(k.uses_shape)
        finally:
            self._syncing = False

    def _kernel_changed(self, name):
        # Picking a kernel is a deliberate choice, so its suggested values
        # replace whatever was there.
        self._edited.discard("a")
        self._edited.discard("eta")
        self._apply_kernel_defaults(force=True)

    def sync_defaults(self):
        # Called whenever a stage changes, so it must never overwrite a value
        # the user typed.
        self._apply_kernel_defaults(force=False)
        ds = self.session.dataset
        if ds is not None:
            self.reference.setMaximum(max(ds.n_frames, 1))
        conn = self.session.unwrap["connectivity"]
        seed = self.session.unwrap["seed"]
        where = "highest quality" if seed == "auto" else f"row {seed[0]}, col {seed[1]}"
        self.unwrap_lbl.setText(f"{conn}-connected, seed {where}")

    # ------------------------------------------------------------------ run
    def run(self):
        filt = self.session.result(Stage.FILTER)
        mask = self.session.mask
        if filt is None or mask is None:
            self.banner.setText("Need a filtered image and a mask.")
            return

        Ih = filt["image"]
        kwargs = dict(
            kernel=self.kernel.currentText(),
            a=self.a.value(),
            a_dist=self.a_dist.currentText(),
            eta=self.eta.value(),
            smooth_reference=self.smooth.isChecked(),
            temporal_fit=self.tfit.isChecked(),
            temporal_order=self.torder.value(),
            downsample=self.downsample.value(),
            njumps=self.njumps.value(),
            phase_jump_rule=self.jump_rule.currentData(),
            reference=self.reference.value() - 1,
            connectivity=self.session.unwrap["connectivity"],
            seed=self.session.unwrap["seed"],
        )

        def job(progress, cancelled):
            def cb(i, n, stage):
                if cancelled():
                    raise RuntimeError("cancelled at a frame boundary")
                progress(i, n, stage)
            return estimate_motion_2d(Ih, mask, progress=cb, **kwargs)

        self.session.set_busy(True)
        self.run_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.banner.setText("Estimating. Cancel takes effect at the next frame.")

        self._thread, self._worker = run_in_background(
            self, job, self._done, self._failed, self._progress
        )

    def _progress(self, i, n, stage):
        self.progress.setMaximum(max(n, 1))
        self.progress.setValue(i + 1)
        self.progress.setFormat(f"{stage} %v/%m")

    def _cancel(self):
        if self._worker:
            self._worker.cancel()
            self.banner.setText("Cancelling at the next frame boundary.")

    def _done(self, res):
        self.session.set_busy(False)
        self.run_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.progress.reset()
        self._raw = res
        self.session.set_result(Stage.MOTION, self._flipped(res))

        # Frame 1 is the reference, where the displacement is zero by
        # definition, so show the frame that actually moved.
        mag = np.nanmean(np.abs(res.best), axis=tuple(range(res.best.ndim - 1)))
        if np.isfinite(mag).any():
            self.session.frameRequested.emit(int(np.nanargmax(mag)))

        bad = np.flatnonzero(~res.unwrap_complete)
        msg = "Motion estimated."
        if bad.size:
            msg += f" Unwrapping was incomplete on frame(s) {[int(b) + 1 for b in bad]}."
        self.banner.setText(msg)
        self.refresh()

    def _failed(self, msg, tb):
        self.session.set_busy(False)
        self.run_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.progress.reset()
        self.banner.setText(f"Failed: {msg}")

    # -------------------------------------------------------------- display
    def refresh(self):
        res = self.session.result(Stage.MOTION)
        self.plot.clear()
        if res is None:
            self.plot.draw()
            return

        jumps = res.phase_jumps / (2 * np.pi)
        text = ("x: " + np.array2string(jumps[0].astype(int), separator=" ") +
                "\ny: " + np.array2string(jumps[1].astype(int), separator=" "))

        # Say whether the search budget is doing anything. The candidates are
        # centred on a coarse estimate, so if the chosen jump never sits more
        # than a step or two away, raising the budget cannot change the
        # result.
        info = res.meta.get("jump_search")
        if info is not None:
            used = int(np.abs(info["offset"]).max())
            budget = int(info["njumps"])
            text += f"\nsearch used {used} of {budget} steps"
            if used < budget and res.meta.get("phase_jump_rule", "published") == "published":
                text += " (raising it will not change the result)"
        self.jump_lbl.setText(text)

        fr = min(self.session.current_frame, res.u.shape[-1] - 1)
        u = res.best
        mag = np.hypot(u[:, :, 0, fr], u[:, :, 1, fr])

        valid = np.isfinite(mag)
        r0, r1, c0, c1 = mask_bbox(valid)
        crop = mag[r0:r1, c0:c1]

        ax = self.plot.axes(1, 2, 1)
        im = ax.imshow(crop, cmap="magma")
        ax.set_title(f"|u| frame {fr + 1}")
        ax.axis("off")
        self.plot.fig.colorbar(im, ax=ax, shrink=0.8, label="px")

        ax2 = self.plot.axes(1, 2, 2)
        ux = u[r0:r1, c0:c1, 0, fr]
        uy = u[r0:r1, c0:c1, 1, fr]
        step = max(crop.shape[0] // 30, 1)
        ys, xs = np.mgrid[0:crop.shape[0]:step, 0:crop.shape[1]:step]
        ax2.quiver(xs, ys, ux[::step, ::step], -uy[::step, ::step],
                   color="#3a7fe8", scale_units="xy", angles="xy", scale=0.4)
        ax2.set_aspect("equal"); ax2.invert_yaxis(); ax2.axis("off")
        ax2.set_title("displacement")
        self.plot.draw()

class StrainTab(QWidget):
    """Strain maps, global curves and regional analysis.

    Every parameter the strain computation accepts is exposed here, because the
    ones that matter are properties of the subject rather than of the method:
    where the cavity centre is, where the right-ventricular insertion sits, how
    many segments to report and at which short-axis level.
    """

    COMPONENTS = ("CC", "RR", "XX", "YY", "XY", "RC", "p1", "p2", "p1or")
    LABELS = {
        "CC": "CC (circumferential)", "RR": "RR (radial)",
        "XX": "XX", "YY": "YY", "XY": "XY", "RC": "RC (shear)",
        "p1": "p1 (1st principal)", "p2": "p2 (2nd principal)",
        "p1or": "p1 orientation",
    }
    SHORTENING = {"CC", "RR", "LL"}

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session
        self._pick = None          # None, "origin" or "insertion"
        self._crop = (0, 0)        # row, col offset of the displayed crop

        # -------------------------------------------------- geometry
        self.origin_x = QDoubleSpinBox(); self.origin_x.setRange(0, 4096)
        self.origin_y = QDoubleSpinBox(); self.origin_y.setRange(0, 4096)
        for w in (self.origin_x, self.origin_y):
            w.setDecimals(1); w.setSingleStep(0.5)
        self.origin_auto = QCheckBox("from the mask centroid")
        self.origin_auto.setChecked(True)
        self.origin_auto.toggled.connect(self._origin_mode)
        self.pick_origin = QPushButton("pick on the map")
        self.pick_origin.setCheckable(True)
        self.pick_origin.toggled.connect(lambda on: self._set_pick("origin", on))

        self.ins_x = QDoubleSpinBox(); self.ins_x.setRange(0, 4096)
        self.ins_y = QDoubleSpinBox(); self.ins_y.setRange(0, 4096)
        for w in (self.ins_x, self.ins_y):
            w.setDecimals(1); w.setSingleStep(0.5)
        self.ins_auto = QCheckBox("default reference angle")
        self.ins_auto.setChecked(True)
        self.ins_auto.toggled.connect(self._ins_mode)
        self.pick_ins = QPushButton("pick on the map")
        self.pick_ins.setCheckable(True)
        self.pick_ins.toggled.connect(lambda on: self._set_pick("insertion", on))

        # -------------------------------------------------- segments
        self.region = QComboBox()
        self.region.setToolTip(
            "Which outlined structure to analyse. Segment wedges spanning both "
            "ventricles do not mean anything, so pick one when both are drawn.")
        self.region.currentTextChanged.connect(self._mark_stale)
        self.nseg = QSpinBox(); self.nseg.setRange(1, 32); self.nseg.setValue(6)
        self.clockwise = QCheckBox("clockwise numbering")
        self.level = QComboBox()
        self.level.addItems(["none", "basal", "mid", "apical"])
        self.level.setCurrentText("mid")
        self.level.currentTextChanged.connect(self._level_changed)
        self.outlier = QDoubleSpinBox()
        self.outlier.setRange(0.01, 100.0); self.outlier.setDecimals(2)
        self.outlier.setValue(0.5)
        self.outlier.setToolTip(
            "Pixels whose strain exceeds this are left out of the segment mean.")

        # -------------------------------------------------- display
        self.component = QComboBox()
        self.component.currentTextChanged.connect(self.refresh)
        self.auto_range = QCheckBox("percentile range")
        self.auto_range.setChecked(True)
        self.auto_range.toggled.connect(self._range_mode)
        self.pct_lo = QSpinBox(); self.pct_lo.setRange(0, 49); self.pct_lo.setValue(2)
        self.pct_hi = QSpinBox(); self.pct_hi.setRange(51, 100); self.pct_hi.setValue(98)
        self.vmin = QDoubleSpinBox(); self.vmin.setRange(-10, 10)
        self.vmax = QDoubleSpinBox(); self.vmax.setRange(-10, 10)
        for w in (self.vmin, self.vmax):
            w.setDecimals(3); w.setSingleStep(0.05)
        self.vmin.setValue(-0.5); self.vmax.setValue(0.25)
        for w in (self.pct_lo, self.pct_hi, self.vmin, self.vmax):
            w.valueChanged.connect(self.refresh)

        self.run_btn = QPushButton("Compute strain")
        self.run_btn.clicked.connect(self.run)
        self.export_btn = QPushButton("Export…")
        self.export_btn.clicked.connect(self.export)
        self.snap_btn = QPushButton("Take snapshot")
        self.snap_btn.setToolTip(
            "Keep this result and its settings so it can be compared later. "
            "Snapshots survive a settings change.")
        self.snap_btn.clicked.connect(self.take_snapshot)

        # a setting change only invalidates the analysis, not the motion fit,
        # so recomputing is cheap and the button stays enabled throughout
        for w in (self.nseg, self.outlier, self.origin_x, self.origin_y,
                  self.ins_x, self.ins_y):
            w.valueChanged.connect(self._mark_stale)
        for w in (self.clockwise, self.origin_auto, self.ins_auto):
            w.toggled.connect(self._mark_stale)

        geo = QFormLayout()
        geo.addRow("", self.origin_auto)
        geo.addRow("centre x", self.origin_x)
        geo.addRow("centre y", self.origin_y)
        geo.addRow("", self.pick_origin)
        geo.addRow("", self.ins_auto)
        geo.addRow("insertion x", self.ins_x)
        geo.addRow("insertion y", self.ins_y)
        geo.addRow("", self.pick_ins)
        geo_box = QGroupBox("Geometry"); geo_box.setLayout(geo)

        seg = QFormLayout()
        seg.addRow("region", self.region)
        seg.addRow("segments", self.nseg)
        seg.addRow("level", self.level)
        seg.addRow("", self.clockwise)
        seg.addRow("outlier above", self.outlier)
        seg_box = QGroupBox("Segments"); seg_box.setLayout(seg)

        disp = QFormLayout()
        disp.addRow("component", self.component)
        disp.addRow("", self.auto_range)
        disp.addRow("low percentile", self.pct_lo)
        disp.addRow("high percentile", self.pct_hi)
        disp.addRow("vmin", self.vmin)
        disp.addRow("vmax", self.vmax)
        disp_box = QGroupBox("Display"); disp_box.setLayout(disp)

        self.banner = StageBanner(
            "Estimate motion first.\n\n"
            "The circumferential and radial components are the ones reported "
            "in the paper. The principal strains are a different quantity."
        )

        side = QVBoxLayout()
        for w in (geo_box, seg_box, disp_box, self.run_btn, self.snap_btn,
                  self.export_btn, self.banner):
            side.addWidget(w)
        side.addStretch(1)
        panel = QWidget(); panel.setLayout(side)
        panel.setMaximumWidth(340)
        scroll = QScrollArea()
        scroll.setWidget(panel); scroll.setWidgetResizable(True)
        scroll.setMaximumWidth(360)

        self.plot = MplCanvas(9, 6)
        self.toolbar = NavigationToolbar2QT(self.plot, self)
        self.plot.mpl_connect("button_press_event", self._clicked)

        self.table = QTableWidget(0, 0)
        self.table.setMaximumHeight(190)
        self.table.setAlternatingRowColors(True)

        right = QVBoxLayout()
        right.addWidget(self.toolbar)
        right.addWidget(self.plot, 3)
        right.addWidget(self.table, 1)
        right_w = QWidget(); right_w.setLayout(right)

        lay = QHBoxLayout(self)
        lay.addWidget(scroll)
        lay.addWidget(right_w, 1)

        self._origin_mode(True)
        self._ins_mode(True)
        self._range_mode(True)
        self._sync_components()
        self._sync_regions()

    # ------------------------------------------------------------- controls

    def _origin_mode(self, auto):
        for w in (self.origin_x, self.origin_y, self.pick_origin):
            w.setEnabled(not auto)

    def _ins_mode(self, auto):
        for w in (self.ins_x, self.ins_y, self.pick_ins):
            w.setEnabled(not auto)

    def _range_mode(self, auto):
        for w in (self.pct_lo, self.pct_hi):
            w.setEnabled(auto)
        for w in (self.vmin, self.vmax):
            w.setEnabled(not auto)
        self.refresh()

    def _level_changed(self, _text):
        """Offer the segment count the chosen level actually names."""
        native = NATIVE_SEGMENTS.get(self.level.currentText())
        if native is not None and self.nseg.value() != native:
            self.nseg.setValue(native)
        self.refresh()

    def _mark_stale(self, *_):
        if self.session.has(Stage.STRAIN):
            self.banner.setText("A setting changed. Compute strain again.")

    def _set_pick(self, which, on):
        self._pick = which if on else None
        other = self.pick_ins if which == "origin" else self.pick_origin
        if on and other.isChecked():
            other.setChecked(False)
        if on:
            self.banner.setText(f"Click the strain map to set the {which}.")

    def _clicked(self, event):
        """Take a point from the strain map, in library pixel coordinates."""
        if self._pick is None or event.inaxes is None or event.xdata is None:
            return
        r0, c0 = self._crop
        # the library grids are 1-based, and the panel shows a crop
        x = c0 + float(event.xdata) + 1.0
        y = r0 + float(event.ydata) + 1.0
        if self._pick == "origin":
            self.origin_x.setValue(x); self.origin_y.setValue(y)
            self.pick_origin.setChecked(False)
        else:
            self.ins_x.setValue(x); self.ins_y.setValue(y)
            self.pick_ins.setChecked(False)
        self.banner.setText(f"{self._pick} set to ({x:.1f}, {y:.1f}). "
                            "Compute strain again.")
        self._pick = None

    def _sync_regions(self):
        """Offer the structures that were actually outlined."""
        keys = sorted(self.session.structure_masks)
        names = ["all outlined tissue"] + [_REGION_NAMES.get(k, k) for k in keys]
        want = self.region.currentText()
        blocked = self.region.blockSignals(True)
        self.region.clear()
        self.region.addItems(names)
        if want in names:
            self.region.setCurrentText(want)
        self.region.blockSignals(blocked)
        self.region.setEnabled(len(names) > 1)

    def _region_mask(self):
        """The mask to analyse, which may be a subset of the segmented tissue."""
        text = self.region.currentText()
        for key, name in _REGION_NAMES.items():
            if name == text and key in self.session.structure_masks:
                return self.session.structure_masks[key]
        if text in self.session.structure_masks:
            return self.session.structure_masks[text]
        return self.session.mask

    def _sync_components(self):
        """Offer only the components the result actually holds."""
        res = self.session.result(Stage.STRAIN)
        st = res["strain"] if res else None
        names = [c for c in self.COMPONENTS
                 if st is None or getattr(st, c, None) is not None]
        want = self.component.currentText()
        blocked = self.component.blockSignals(True)
        self.component.clear()
        self.component.addItems([self.LABELS.get(c, c) for c in names])
        if want in [self.LABELS.get(c, c) for c in names]:
            self.component.setCurrentText(want)
        self.component.blockSignals(blocked)

    # ------------------------------------------------------------ computing

    def _times(self, nfr):
        """Frame times from the acquisition, set on the Data tab."""
        ds = self.session.dataset
        if ds is None:
            return None
        t = ds.times
        return t if t is not None and t.size == nfr else None

    def _origin(self):
        if self.origin_auto.isChecked():
            return None
        return (self.origin_x.value(), self.origin_y.value())

    def run(self):
        motion = self.session.result(Stage.MOTION)
        self._sync_regions()
        mask = self._region_mask()
        if motion is None or mask is None:
            self.banner.setText("Estimate motion first.")
            return
        if not mask[:, :, 0].any():
            self.banner.setText("That region is empty on the first frame.")
            return

        nfr = motion.best.shape[-1]
        st = pixel_strain_2d(motion.best, mask[:, :, 0], origin=self._origin(),
                             times=self._times(nfr))
        kw = {}
        if not self.ins_auto.isChecked():
            kw["position_a"] = self.ins_x.value()
            kw["position_b"] = self.ins_y.value()
        seg = strain_by_segments(
            st, n_segments=self.nseg.value(),
            clockwise=self.clockwise.isChecked(), origin=self._origin(),
            outlier_threshold=self.outlier.value(), **kw,
        )
        self.session.set_result(Stage.STRAIN, {"strain": st, "segments": seg})
        self._sync_components()
        ox, oy = st.origin if st.origin is not None else (float("nan"),) * 2
        self.banner.setText(
            f"{int(st.mask.sum())} pixels of {self.region.currentText()} over "
            f"{seg.n_segments} segments, centre ({ox:.1f}, {oy:.1f})."
        )
        self.refresh()

    # -------------------------------------------------------------- display

    def _level(self):
        text = self.level.currentText()
        return None if text == "none" else text

    def _labels(self, n):
        return aha_labels(self._level(), n)

    def _range(self, field):
        if self.auto_range.isChecked():
            lo, hi = nan_range(field, self.pct_lo.value(), self.pct_hi.value())
            if lo == hi:
                lo, hi = lo - 0.1, hi + 0.1
            return lo, hi
        return self.vmin.value(), self.vmax.value()

    def _fill_table(self, seg):
        try:
            header, rows = segment_table(seg, times=self._times(
                (seg.CC if seg.CC is not None else seg.RR).shape[-1]),
                level=self._level())
        except ValueError:
            self.table.setRowCount(0); self.table.setColumnCount(0)
            return
        self.table.setColumnCount(len(header))
        self.table.setHorizontalHeaderLabels(header)
        self.table.setRowCount(len(rows))
        for r, row in enumerate(rows):
            for c, v in enumerate(row):
                text = v if isinstance(v, str) else (
                    f"{v:.4f}" if isinstance(v, float) else str(v))
                item = QTableWidgetItem(text)
                item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                self.table.setItem(r, c, item)
        self.table.resizeColumnsToContents()

    def refresh(self):
        res = self.session.result(Stage.STRAIN)
        self.plot.clear()
        if res is None:
            self.table.setRowCount(0)
            self.plot.draw()
            return
        st, seg = res["strain"], res["segments"]
        comp = self.component.currentText().split()[0] or "CC"
        field = getattr(st, comp, None)
        if field is None:
            self.plot.draw()
            return

        nfr = field.shape[-1]
        fr = min(self.session.current_frame, nfr - 1)
        r0, r1, c0, c1 = mask_bbox(st.mask)
        self._crop = (r0, c0)
        times = self._times(nfr)
        xs = times if times is not None else np.arange(nfr)
        xlabel = "time (ms)" if times is not None else "frame"
        diverging = comp not in {"p1or"}

        ax = self.plot.axes(2, 3, 1)
        lo, hi = self._range(field)
        im = ax.imshow(field[r0:r1, c0:c1, fr],
                       cmap="RdBu_r" if diverging else "twilight",
                       vmin=lo, vmax=hi)
        ax.set_title(f"{comp} frame {fr + 1}"); ax.axis("off")
        self.plot.fig.colorbar(im, ax=ax, shrink=0.8)
        if st.origin is not None and np.all(np.isfinite(st.origin)):
            ax.plot(st.origin[0] - 1 - c0, st.origin[1] - 1 - r0,
                    "k+", ms=9, mew=1.5)
        if not self.ins_auto.isChecked():
            ax.plot(self.ins_x.value() - 1 - c0, self.ins_y.value() - 1 - r0,
                    "kx", ms=8, mew=1.5)

        ax2 = self.plot.axes(2, 3, 2)
        for name, style in (("CC", "-o"), ("RR", "-s")):
            if getattr(st, name, None) is not None:
                ax2.plot(xs, st.global_curve(name), style, ms=3,
                         label=f"G{name[0]}S")
        ax2.axhline(0, color="#888", lw=0.5)
        ax2.set_xlabel(xlabel); ax2.set_ylabel("strain")
        ax2.legend(fontsize=8); ax2.grid(alpha=0.3)
        ax2.set_title("global curves")

        ax3 = self.plot.axes(2, 3, 3)
        ax3.imshow(seg.segments[r0:r1, c0:c1], cmap="tab20")
        ax3.set_title(f"{seg.n_segments} segments"); ax3.axis("off")

        curves = seg.CC if seg.CC is not None else seg.RR
        which = "CC" if seg.CC is not None else "RR"
        labels = self._labels(seg.n_segments)

        ax4 = self.plot.axes(2, 3, 4)
        if curves is not None:
            for s in range(curves.shape[0]):
                ax4.plot(xs, curves[s], lw=1, label=labels[s])
            ax4.legend(fontsize=5, ncol=2)
        ax4.axhline(0, color="#888", lw=0.5)
        ax4.set_xlabel(xlabel); ax4.grid(alpha=0.3)
        ax4.set_title(f"{which} per segment")

        if curves is not None:
            self._bullseye(seg, curves, which, labels)
            self._peak_bars(curves, which, labels)
        self._fill_table(seg)
        self.plot.draw()

    def _bullseye(self, seg, curves, which, labels):
        """Peak strain per segment on a polar ring."""
        m = segment_metrics(curves, self._times(curves.shape[-1]), which)
        ax = self.plot.fig.add_subplot(2, 3, 5, projection="polar")
        rings = bullseye_rings(m["peak"])
        vals = m["peak"]
        good = vals[np.isfinite(vals)]
        lo, hi = (float(good.min()), float(good.max())) if good.size else (-1.0, 0.0)
        if lo == hi:
            lo, hi = lo - 0.01, hi + 0.01
        sm = ScalarMappable(norm=Normalize(vmin=lo, vmax=hi), cmap="RdBu_r")
        for ring in rings:
            for theta, width, value in ring["wedges"]:
                ax.bar(theta + width / 2, ring["height"], width=width,
                       bottom=ring["radius"],
                       color=sm.to_rgba(value) if np.isfinite(value) else "#dddddd",
                       edgecolor="white", linewidth=1.0)
        step = 2 * np.pi / seg.n_segments
        ax.set_xticks([i * step + step / 2 for i in range(seg.n_segments)])
        ax.set_xticklabels([str(i + 1) for i in range(seg.n_segments)], fontsize=6)
        ax.set_yticks([]); ax.set_ylim(0, 1)
        ax.set_title(f"peak {which}", fontsize=9)
        self.plot.fig.colorbar(sm, ax=ax, shrink=0.7)

    def _peak_bars(self, curves, which, labels):
        m = segment_metrics(curves, self._times(curves.shape[-1]), which)
        ax = self.plot.axes(2, 3, 6)
        pos = np.arange(len(labels))
        ax.barh(pos, m["peak"], color="#4878a8")
        ax.set_yticks(pos)
        ax.set_yticklabels(labels, fontsize=5)
        ax.invert_yaxis()
        ax.axvline(0, color="#888", lw=0.5)
        ax.set_xlabel(f"peak {which}")
        ax.grid(alpha=0.3, axis="x")
        ax.set_title("peak per segment", fontsize=9)

    # ------------------------------------------------------------ snapshots

    def settings(self):
        """Everything that shaped the current result, for the settings diff."""
        motion = self.session.result(Stage.MOTION)
        meta = dict(getattr(motion, "meta", {}) or {}) if motion else {}
        out = {
            "region": self.region.currentText(),
            "segments": self.nseg.value(),
            "level": self.level.currentText(),
            "clockwise": self.clockwise.isChecked(),
            "outlier above": self.outlier.value(),
            "centre": ("mask centroid" if self.origin_auto.isChecked()
                       else (self.origin_x.value(), self.origin_y.value())),
            "insertion": ("default angle" if self.ins_auto.isChecked()
                          else (self.ins_x.value(), self.ins_y.value())),
        }
        for key in ("kernel", "a", "eta", "downsample", "reference", "njumps",
                    "phase_jump_rule", "connectivity", "seed", "temporal_fit",
                    "temporal_order"):
            if key in meta:
                v = meta[key]
                out[key] = tuple(np.ravel(v)) if isinstance(v, np.ndarray) else v
        return out

    def take_snapshot(self):
        res = self.session.result(Stage.STRAIN)
        if res is None:
            self.banner.setText("Compute strain before taking a snapshot.")
            return
        n = len(self.session.snapshots) + 1
        snap = Snapshot(
            label=f"run {n}: {self.region.currentText()}, {self.nseg.value()} seg",
            settings=self.settings(),
            strain=res["strain"], segments=res["segments"],
        )
        self.session.snapshots.add(snap)
        self.session.snapshotsChanged.emit()
        self.banner.setText(f"Kept as \"{snap.label}\". "
                            f"{len(self.session.snapshots)} snapshots held.")

    # --------------------------------------------------------------- export

    def export(self):
        res = self.session.result(Stage.STRAIN)
        if res is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export", "harpi_result.mat",
            "MATLAB (*.mat);;Segment table (*.csv)")
        if not path:
            return
        st, seg = res["strain"], res["segments"]
        if path.endswith(".csv"):
            import csv

            nfr = (seg.CC if seg.CC is not None else seg.RR).shape[-1]
            header, rows = segment_table(seg, times=self._times(nfr),
                                         level=self._level())
            with open(path, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(header)
                w.writerows(rows)
        else:
            from ..io import save_mat

            save_mat(path, motion=self.session.result(Stage.MOTION), strain=st)
        self.banner.setText(f"Wrote {path}")
