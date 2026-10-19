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
from ..analysis import (NATIVE_SEGMENTS, expected_septal_arc,
                        short_labels)
from .session import Stage
from .snapshots import Snapshot
from .widgets import (BackgroundControls, MplCanvas, StageBanner,
                      background_image, mask_bbox, nan_range)
from .worker import run_in_background

__all__ = ["MotionTab", "StrainTab"]


def draw_background(ax, controls, session, frame, box):
    """Put the tagged image under a map, and say whether one was drawn.

    The data layer keeps its NaN outside the mask, which matplotlib renders
    transparent, so the anatomy shows through the holes with no masking work.
    """
    if not controls.wanted():
        return False
    bg = background_image(session.dataset, session.result(Stage.FILTER),
                          controls.choice(), frame, session.current_slice)
    if bg is None:
        return False
    r0, r1, c0, c1 = box
    ax.imshow(bg.array[r0:r1, c0:c1], cmap=bg.cmap or "gray",
              vmin=bg.vmin, vmax=bg.vmax, zorder=0)
    return True

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
        self.jump_rule.addItem("nearest coordinate", "nearest")
        self.jump_rule.addItem("published", "published")
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
        # Only the published rule has a search budget. Set the state as well as
        # connecting the signal: a connection alone leaves the initial state
        # whatever the widget defaults to, which was right only while the
        # published rule happened to be the first item in the box.
        self.jump_rule.currentIndexChanged.connect(lambda _: self._sync_jump_budget())
        self._sync_jump_budget()

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
        self.bg = BackgroundControls()
        self.bg.changed.connect(self.refresh)

        # Arrows are drawn at their true length by default, so what is on
        # screen is the displacement in pixels and can be read against the
        # image under it.
        self.arrow_scale = QDoubleSpinBox()
        self.arrow_scale.setRange(0.1, 20.0)
        self.arrow_scale.setDecimals(1)
        self.arrow_scale.setSingleStep(0.5)
        self.arrow_scale.setValue(1.0)
        self.arrow_scale.setPrefix("x ")
        self.arrow_scale.setToolTip(
            "How much longer than the displacement itself the arrows are "
            "drawn. One is true length.")
        self.arrow_scale.valueChanged.connect(lambda _: self.refresh())
        arrow_form = QFormLayout()
        arrow_form.addRow("arrow length", self.arrow_scale)
        self.arrow_box = QGroupBox("Arrows")
        self.arrow_box.setLayout(arrow_form)

        side.addWidget(self.run_btn)
        side.addWidget(self.cancel_btn)
        side.addWidget(self.progress)
        side.addWidget(self.flip_box)
        side.addWidget(self.arrow_box)
        side.addWidget(self.bg)
        side.addWidget(QLabel("phase jumps (x2π), per frame:"))
        side.addWidget(self.jump_lbl)
        side.addWidget(self.banner)
        side.addStretch(1)

        self.plot = MplCanvas(7, 5)
        panel = QWidget(); panel.setLayout(side); panel.setMaximumWidth(320)
        # the controls are taller than the window, as on the Strain tab
        scroll = QScrollArea()
        scroll.setWidget(panel); scroll.setWidgetResizable(True)
        scroll.setMaximumWidth(340)
        lay = QHBoxLayout(self)
        lay.addWidget(scroll)
        lay.addWidget(self.plot, 1)

    def _mark_edited(self, name):
        """Remember that a value was set by hand, so defaults leave it alone."""
        if not self._syncing:
            self._edited.add(name)

    def _sync_jump_budget(self):
        self.njumps.setEnabled(self.jump_rule.currentData() == "published")

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
        draw_background(ax, self.bg, self.session, fr, (r0, r1, c0, c1))
        im = ax.imshow(crop, cmap="magma", alpha=self.bg.opacity(), zorder=1)
        ax.set_title(f"|u| frame {fr + 1}")
        ax.axis("off")
        # from a mappable rather than the faded image, so the bar stays solid
        bar = ScalarMappable(norm=Normalize(*im.get_clim()), cmap="magma")
        self.plot.fig.colorbar(bar, ax=ax, shrink=0.8, label="px")

        ax2 = self.plot.axes(1, 2, 2)
        ux = u[r0:r1, c0:c1, 0, fr]
        uy = u[r0:r1, c0:c1, 1, fr]
        step = max(crop.shape[0] // 30, 1)
        ys, xs = np.mgrid[0:crop.shape[0]:step, 0:crop.shape[1]:step]
        # u[..., 1] is the displacement along rows, and invert_yaxis already
        # makes the drawn y axis increase downwards to match them. Negating it
        # as well flipped the arrows a second time, so a contraction fanned
        # outwards top and bottom while converging left and right.
        draw_background(ax2, self.bg, self.session, fr, (r0, r1, c0, c1))
        # the arrows stay solid: they are line art covering almost no area,
        # and the image is already visible between them
        # scale is data units per unit of arrow length, so the reciprocal of
        # the exaggeration asked for: 1 draws the displacement itself
        ax2.quiver(xs, ys, ux[::step, ::step], uy[::step, ::step],
                   color="#3a7fe8", scale_units="xy", angles="xy",
                   scale=1.0 / max(self.arrow_scale.value(), 1e-6),
                   zorder=2)
        ax2.set_aspect("equal")
        # Set the limits rather than inverting them. imshow already puts row 0
        # at the top, so a background image plus invert_yaxis would flip twice,
        # and this holds whether or not one was drawn.
        ax2.set_xlim(-0.5, crop.shape[1] - 0.5)
        ax2.set_ylim(crop.shape[0] - 0.5, -0.5)
        ax2.axis("off")
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
        self._pick = None          # which point a click would set, or None
        self._arming = False       # see _set_pick
        self._crop = (0, 0)        # row, col offset of the displayed crop
        self._map_ax = None        # the panel a click is taken from

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

        # The anterior insertion separates the septum from the anterior wall,
        # so it is the boundary between segment 1 and segment 2. The inferior
        # one closes the septum, and the two together say which way round the
        # numbering runs.
        self.ins_ant_x = QDoubleSpinBox(); self.ins_ant_x.setRange(0, 4096)
        self.ins_ant_y = QDoubleSpinBox(); self.ins_ant_y.setRange(0, 4096)
        self.ins_inf_x = QDoubleSpinBox(); self.ins_inf_x.setRange(0, 4096)
        self.ins_inf_y = QDoubleSpinBox(); self.ins_inf_y.setRange(0, 4096)
        for w in (self.ins_ant_x, self.ins_ant_y,
                  self.ins_inf_x, self.ins_inf_y):
            w.setDecimals(1); w.setSingleStep(0.5)
        self.ins_auto = QCheckBox("default reference angle")
        self.ins_auto.setChecked(True)
        self.ins_auto.toggled.connect(self._ins_mode)
        self.pick_ant = QPushButton("pick anterior on the map")
        self.pick_ant.setCheckable(True)
        self.pick_ant.toggled.connect(lambda on: self._set_pick("anterior", on))
        self.pick_inf = QPushButton("pick inferior on the map")
        self.pick_inf.setCheckable(True)
        self.pick_inf.toggled.connect(lambda on: self._set_pick("inferior", on))
        self.septum_lbl = QLabel("-")
        self.septum_lbl.setWordWrap(True)
        self.septum_lbl.setToolTip(
            "What the two insertion points span, against what the model "
            "expects. The wedges stay equal either way.")

        # one entry per point that can be picked, so arming is not a set of
        # hard-coded cases
        self._pickers = {"origin": self.pick_origin,
                         "anterior": self.pick_ant,
                         "inferior": self.pick_inf}

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
                  self.ins_ant_x, self.ins_ant_y,
                  self.ins_inf_x, self.ins_inf_y):
            w.valueChanged.connect(self._mark_stale)
        for w in (self.ins_ant_x, self.ins_ant_y,
                  self.ins_inf_x, self.ins_inf_y):
            w.valueChanged.connect(self._sync_clockwise)
        for w in (self.clockwise, self.origin_auto, self.ins_auto):
            w.toggled.connect(self._mark_stale)

        geo = QFormLayout()
        geo.addRow("", self.origin_auto)
        geo.addRow("centre x", self.origin_x)
        geo.addRow("centre y", self.origin_y)
        geo.addRow("", self.pick_origin)
        geo.addRow("", self.ins_auto)
        geo.addRow("anterior x", self.ins_ant_x)
        geo.addRow("anterior y", self.ins_ant_y)
        geo.addRow("", self.pick_ant)
        geo.addRow("inferior x", self.ins_inf_x)
        geo.addRow("inferior y", self.ins_inf_y)
        geo.addRow("", self.pick_inf)
        geo.addRow("septum", self.septum_lbl)
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

        self.bg = BackgroundControls()
        self.bg.changed.connect(self.refresh)

        side = QVBoxLayout()
        for w in (geo_box, seg_box, disp_box, self.bg, self.run_btn,
                  self.snap_btn, self.export_btn, self.banner):
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
        for w in (self.ins_ant_x, self.ins_ant_y, self.pick_ant,
                  self.ins_inf_x, self.ins_inf_y, self.pick_inf):
            w.setEnabled(not auto)
        self._sync_clockwise()

    def _sync_clockwise(self, *_):
        """The two insertions decide the direction, so the box is then inert."""
        both = self._insertions()[1] is not None
        self.clockwise.setEnabled(not both)
        self.clockwise.setToolTip(
            "Which way the segment numbers run. The two insertion points "
            "decide this when both are set, so it only applies to one point "
            "or none.")

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
        """Arm or disarm a picker, leaving exactly one armed.

        Un-checking the other button emits toggled and lands back here. It
        used to clear the mode that had just been set, so the new button
        stayed lit while the banner invited a click that went nowhere. The
        fix is the order: disarm the others first, then record the mode. The
        guard is belt and braces, so the re-entrant call does nothing at all,
        and a genuine un-toggle still clears.
        """
        if self._arming:
            return
        if on:
            self._arming = True
            try:
                for name, button in self._pickers.items():
                    if name != which and button.isChecked():
                        button.setChecked(False)
            finally:
                self._arming = False
            self._pick = which
            self.banner.setText(f"Click the strain map to set the {which}.")
        elif self._pick == which:
            self._pick = None

    def _clicked(self, event):
        """Take a point from the strain map, in library pixel coordinates."""
        # held, because disarming the button below clears self._pick
        which = self._pick
        if which is None or event.xdata is None:
            return
        # only the map carries pixel coordinates; the curves and the bullseye
        # would hand over a time or an angle
        if self._map_ax is not None and event.inaxes is not self._map_ax:
            return
        if event.inaxes is None:
            return
        r0, c0 = self._crop
        # the library grids are 1-based, and the panel shows a crop
        x = c0 + float(event.xdata) + 1.0
        y = r0 + float(event.ydata) + 1.0
        boxes = {"origin": (self.origin_x, self.origin_y),
                 "anterior": (self.ins_ant_x, self.ins_ant_y),
                 "inferior": (self.ins_inf_x, self.ins_inf_y)}[which]
        boxes[0].setValue(x); boxes[1].setValue(y)
        self._pickers[which].setChecked(False)
        self.banner.setText(f"{which} set to ({x:.1f}, {y:.1f}). "
                            "Compute strain again.")

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

    def _insertions(self):
        """The anterior and inferior points, either of which may be absent.

        A point left at the origin of the spin boxes has not been placed, so
        it is read as absent and the inferior one stays optional.
        """
        if self.ins_auto.isChecked():
            return None, None
        anterior = (self.ins_ant_x.value(), self.ins_ant_y.value())
        inferior = (self.ins_inf_x.value(), self.ins_inf_y.value())
        return (anterior if any(anterior) else None,
                inferior if any(inferior) else None)

    def _describe_septum(self, geometry):
        """What the two points span, against what the model expects."""
        if not geometry:
            return "-"
        arc = geometry.get("septal_arc")
        if arc is None:
            if geometry.get("anterior") is None:
                return "no insertion points, segment 1 starts at the right"
            return "one point, direction from the clockwise box"
        measured = np.degrees(arc)
        ideal = expected_septal_arc(self._level(), self.nseg.value())
        if ideal is None:
            return f"{measured:.0f} deg between the insertions"
        off = abs(measured - np.degrees(ideal))
        text = (f"{measured:.0f} deg between the insertions, "
                f"{off:.0f} from the model's {np.degrees(ideal):.0f}")
        if off > 45:
            text += ". Check which point is which."
        return text

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
        anterior, inferior = self._insertions()
        seg = strain_by_segments(
            st, n_segments=self.nseg.value(),
            clockwise=self.clockwise.isChecked(), origin=self._origin(),
            outlier_threshold=self.outlier.value(),
            anterior=anterior, inferior=inferior,
        )
        septum = self._describe_septum(seg.geometry)
        self.septum_lbl.setText(septum)
        self.session.set_result(Stage.STRAIN, {"strain": st, "segments": seg})
        self._sync_components()
        ox, oy = st.origin if st.origin is not None else (float("nan"),) * 2
        self.banner.setText(
            f"{int(st.mask.sum())} pixels of {self.region.currentText()} over "
            f"{seg.n_segments} segments, centre ({ox:.1f}, {oy:.1f}). "
            f"Septum: {septum}"
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
        # "".split() is empty, so the old "or CC" could never stand in for
        # an empty box; it raised instead
        label = self.component.currentText().split()
        comp = label[0] if label else "CC"
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
        self._map_ax = ax
        lo, hi = self._range(field)
        cmap = "RdBu_r" if diverging else "twilight"
        draw_background(ax, self.bg, self.session, fr, (r0, r1, c0, c1))
        ax.imshow(field[r0:r1, c0:c1, fr], cmap=cmap, vmin=lo, vmax=hi,
                  alpha=self.bg.opacity(), zorder=1)
        ax.set_title(f"{comp} frame {fr + 1}"); ax.axis("off")
        bar = ScalarMappable(norm=Normalize(vmin=lo, vmax=hi), cmap=cmap)
        self.plot.fig.colorbar(bar, ax=ax, shrink=0.8)
        centre = None
        if st.origin is not None and np.all(np.isfinite(st.origin)):
            centre = (st.origin[0] - 1 - c0, st.origin[1] - 1 - r0)
            ax.plot(*centre, "k+", ms=9, mew=1.5)
        anterior, inferior = self._insertions()
        for point, marker, style in ((anterior, "kx", "--"),
                                     (inferior, "k*", ":")):
            if point is None:
                continue
            px, py = point[0] - 1 - c0, point[1] - 1 - r0
            ax.plot(px, py, marker, ms=9, mew=1.5)
            if centre is not None:
                # the ray, so the numbering can be checked by eye
                ax.plot([centre[0], px], [centre[1], py], style,
                        color="k", lw=1.0)

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
        # seeing the wedges against the wall is the point of this panel
        draw_background(ax3, self.bg, self.session, fr, (r0, r1, c0, c1))
        ax3.imshow(seg.segments[r0:r1, c0:c1], cmap="tab20",
                   alpha=self.bg.opacity(), zorder=1)
        arc = seg.geometry.get("septal_arc")
        ax3.set_title(f"{seg.n_segments} segments", fontsize=9)
        ax3.set_xlabel("" if arc is None else f"septum {np.degrees(arc):.0f} deg",
                       fontsize=7)
        ax3.set_xticks([]); ax3.set_yticks([])
        for side in ax3.spines.values():
            side.set_visible(False)

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
                # centred on theta, not theta + width/2, so segment 1 sits
                # at twelve o'clock rather than starting there
                ax.bar(theta, ring["height"], width=width,
                       bottom=ring["radius"],
                       color=sm.to_rgba(value) if np.isfinite(value) else "#dddddd",
                       edgecolor="white", linewidth=1.0)
        step = 2 * np.pi / seg.n_segments
        # segment 1 is the anterior wall, which belongs at the top, and the
        # wedges run the way they do on the image rather than always
        # anticlockwise
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-seg.geometry.get("direction", -1))
        ax.set_xticks([i * step for i in range(seg.n_segments)])
        ax.set_xticklabels(short_labels(self._level(), seg.n_segments),
                           fontsize=6)
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

    # ----------------------------------------------------------- projects

    def geometry_config(self):
        """The segment geometry, for saving with a project."""
        anterior, inferior = self._insertions()
        return {
            "origin": None if self.origin_auto.isChecked()
            else [self.origin_x.value(), self.origin_y.value()],
            "anterior": None if anterior is None else list(anterior),
            "inferior": None if inferior is None else list(inferior),
            "clockwise": self.clockwise.isChecked(),
            "n_segments": self.nseg.value(),
            "level": self.level.currentText(),
            "outlier_threshold": self.outlier.value(),
        }

    def apply_geometry_config(self, config):
        """Put a saved geometry back on the controls.

        An absent key means the project predates this being saved, so the
        control keeps its default rather than being reset to something the
        user never chose.
        """
        if not config:
            return
        widgets = (self.origin_x, self.origin_y, self.ins_ant_x,
                   self.ins_ant_y, self.ins_inf_x, self.ins_inf_y,
                   self.origin_auto, self.ins_auto, self.clockwise,
                   self.nseg, self.level, self.outlier)
        blocked = [w.blockSignals(True) for w in widgets]
        try:
            origin = config.get("origin")
            self.origin_auto.setChecked(origin is None)
            if origin is not None:
                self.origin_x.setValue(float(origin[0]))
                self.origin_y.setValue(float(origin[1]))
            anterior = config.get("anterior")
            inferior = config.get("inferior")
            self.ins_auto.setChecked(anterior is None)
            for point, boxes in ((anterior, (self.ins_ant_x, self.ins_ant_y)),
                                 (inferior, (self.ins_inf_x, self.ins_inf_y))):
                if point is not None:
                    boxes[0].setValue(float(point[0]))
                    boxes[1].setValue(float(point[1]))
            if "clockwise" in config:
                self.clockwise.setChecked(bool(config["clockwise"]))
            if "n_segments" in config:
                self.nseg.setValue(int(config["n_segments"]))
            if config.get("level"):
                self.level.setCurrentText(str(config["level"]))
            if "outlier_threshold" in config:
                self.outlier.setValue(float(config["outlier_threshold"]))
        finally:
            for w, was in zip(widgets, blocked):
                w.blockSignals(was)
        self._origin_mode(self.origin_auto.isChecked())
        self._ins_mode(self.ins_auto.isChecked())

    # ------------------------------------------------------------ snapshots

    def settings(self):
        """Everything that shaped the current result, for the settings diff."""
        motion = self.session.result(Stage.MOTION)
        meta = dict(getattr(motion, "meta", {}) or {}) if motion else {}
        anterior, inferior = self._insertions()
        res = self.session.result(Stage.STRAIN) or {}
        geometry = getattr(res.get("segments"), "geometry", {})
        out = {
            "region": self.region.currentText(),
            "segments": self.nseg.value(),
            "level": self.level.currentText(),
            "clockwise": self.clockwise.isChecked(),
            "outlier above": self.outlier.value(),
            "centre": ("mask centroid" if self.origin_auto.isChecked()
                       else (self.origin_x.value(), self.origin_y.value())),
            "anterior insertion": ("default angle" if anterior is None
                                   else anterior),
            "inferior insertion": ("not set" if inferior is None else inferior),
            "septal arc": self._describe_septum(geometry),
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
            # every other guard here says why nothing happened
            self.banner.setText("Compute strain before exporting.")
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
