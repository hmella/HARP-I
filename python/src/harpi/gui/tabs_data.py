"""Data and Filtering tabs."""

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
    QPushButton,
    QSlider,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from .. import (
    HarpFilter2D,
    centre_to_frequency,
    frequency_to_centre,
    load_mat_2d,
    load_mat_3d,
)
from ..harmonic import harmonic_center
from ..fft import itok
from .canvas import ImageCanvas
from .session import Stage
from .widgets import StageBanner, nan_range

__all__ = ["DataTab", "FilterTab"]


class DataTab(QWidget):
    """Load a dataset and check its acquisition parameters."""

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session

        self.open_mat = QPushButton("Open .mat…")
        self.open_dicom = QPushButton("Open DICOM folder…")
        self.open_dicom_file = QPushButton("Open DICOM file…")
        self.open_dicom_file.setToolTip(
            "A single file, which may hold every frame of a cine.")
        self.open_dicom_file.clicked.connect(self._open_dicom_file)
        self.open_mat.clicked.connect(self._open_mat)
        self.open_dicom.clicked.connect(self._open_dicom)

        self.shape_lbl = QLabel("-")
        self.frames_lbl = QLabel("-")

        # Editable, since the DICOM tag spacing is often wrong and the shape
        # parameter depends on it.
        # One box per axis. Non-square pixels change ke and so the suggested
        # shape parameter, and collapsing them to one number is wrong.
        self.pixel_size = QDoubleSpinBox()       # rows, also the linked value
        self.pixel_size_col = QDoubleSpinBox()
        self.pixel_size_slice = QDoubleSpinBox()
        for w in (self.pixel_size, self.pixel_size_col, self.pixel_size_slice):
            w.setDecimals(6); w.setRange(1e-6, 1e3); w.setSingleStep(0.1)
        self.link_size = QCheckBox("square pixels")
        self.link_size.setChecked(True)
        self.link_size.toggled.connect(self._link_changed)

        self.tag_spacing = QDoubleSpinBox()
        self.tag_spacing.setDecimals(6); self.tag_spacing.setRange(1e-6, 1e3)
        self.tag_spacing.setSingleStep(0.1)

        # Frame times come from the DICOM TriggerTime when the series was read
        # from DICOM. A .mat file carries none, so they can be set by hand as a
        # constant interval; zero leaves every time axis in frames.
        self.frame_interval = QDoubleSpinBox()
        self.frame_interval.setDecimals(2); self.frame_interval.setRange(0.0, 10000.0)
        self.frame_interval.setSingleStep(5.0); self.frame_interval.setSuffix(" ms")
        self.frame_interval.valueChanged.connect(self._interval_edited)
        self.times_lbl = QLabel("-")

        self.pixel_size.valueChanged.connect(self._row_size_changed)
        for w in (self.pixel_size_col, self.pixel_size_slice, self.tag_spacing):
            w.valueChanged.connect(self._params_edited)

        self.ppw_lbl = QLabel("-")
        self.a_lbl = QLabel("-")
        self.override_lbl = QLabel("")
        self.override_lbl.setStyleSheet("color:#d08b26;")

        form = QFormLayout()
        form.addRow("image", self.shape_lbl)
        form.addRow("frames", self.frames_lbl)
        form.addRow("", self.link_size)
        form.addRow("pixel size, rows", self.pixel_size)
        form.addRow("pixel size, columns", self.pixel_size_col)
        # The interface takes 2-D studies only, which _adopt enforces before
        # session.dataset is ever set, so this can never be enabled. It is
        # kept, hidden, because the library and _params_edited still handle a
        # volume and the row would be wanted again the day the tab does.
        self.slice_row = QLabel("slice thickness")
        form.addRow(self.slice_row, self.pixel_size_slice)
        self.slice_row.setVisible(False)
        self.pixel_size_slice.setVisible(False)
        form.addRow("tag spacing", self.tag_spacing)
        form.addRow("frame interval", self.frame_interval)
        form.addRow("frame times", self.times_lbl)
        form.addRow("pixels / wavelength", self.ppw_lbl)
        form.addRow("suggested a", self.a_lbl)
        form.addRow("", self.override_lbl)

        box = QGroupBox("Acquisition")
        box.setLayout(form)

        self.canvas = ImageCanvas(session, "magnitude")
        self.banner = StageBanner("Open a dataset to begin.")

        buttons = QHBoxLayout()
        buttons.addWidget(self.open_mat)
        buttons.addWidget(self.open_dicom)
        buttons.addWidget(self.open_dicom_file)
        buttons.addStretch(1)

        side = QVBoxLayout()
        side.addLayout(buttons)
        side.addWidget(box)
        side.addWidget(self.banner)
        side.addStretch(1)

        lay = QHBoxLayout(self)
        panel = QWidget(); panel.setLayout(side); panel.setMaximumWidth(340)
        lay.addWidget(panel)
        lay.addWidget(self.canvas, 1)

        self._loading = False

    # ------------------------------------------------------------- loading
    def _open_mat(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open dataset", "", "MATLAB (*.mat)")
        if not path:
            return
        try:
            ds = load_mat_2d(path)
        except Exception:
            try:
                ds = load_mat_3d(path)
            except Exception as exc:
                self.banner.setText(f"Could not read {path}:\n{exc}")
                return
        self._adopt(ds, path)

    def _open_dicom(self):
        path = QFileDialog.getExistingDirectory(self, "Open DICOM folder")
        if path:
            self._load_dicom(path)

    def _open_dicom_file(self):
        """A single file, which may hold a whole cine as multiple frames."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Open DICOM file", "", "DICOM (*.dcm *.ima *.IMA);;Any file (*)")
        if path:
            self._load_dicom(path)

    def _load_dicom(self, path):
        from ..io.dicomio import read_dicom_series, scan_dicom

        from .dicom_picker import choose_series

        try:
            catalogue = scan_dicom(path)
        except Exception as exc:
            self.banner.setText(f"Could not scan {path}:\n{exc}")
            return
        chosen = choose_series(catalogue, self)
        if chosen is None:
            self.banner.setText("No series chosen.")
            return
        try:
            ds = read_dicom_series(path, series_uid=chosen.uid, catalogue=catalogue)
        except Exception as exc:
            self.banner.setText(f"Could not read {chosen.description or path}:\n{exc}")
            return
        if self._adopt(ds, path):
            where = ds.meta.get("tag_spacing_source", "unknown")
            name = chosen.description or chosen.protocol or chosen.uid[-12:]
            self.banner.setText(
                f"Loaded {name}: {ds.n_frames} frames, tag spacing "
                f"{ds.tag_spacing:.2f} ({where})."
            )

    def _adopt(self, ds, path):
        # One place to refuse 3-D. Everything downstream assumes a 2-D series:
        # the canvases, the contour tools, the segment wedges and the strain
        # panels. Half-accepting a volume used to fail later in several
        # different ways, each of them confusing. The library handles 3-D, so
        # the message points at it.
        if getattr(ds, "image", None) is not None and ds.image.ndim != 4:
            self.banner.setText(
                f"This looks like a {ds.image.ndim - 2}-D volume "
                f"{ds.image.shape}. The interface works on a 2-D series.\n\n"
                "Use harpi.estimate_motion_3d from a script for volumes, or "
                "load a single slice here."
            )
            return False

        self._loading = True
        self.session.dataset = ds
        self.session.set_result(Stage.DATA, ds)
        self.session.mask = None
        self.session.contours.clear()
        self.session.labels.clear()

        self.shape_lbl.setText(str(ds.image.shape))
        self.frames_lbl.setText(str(ds.n_frames))
        size = ds.pixel_size if hasattr(ds, "pixel_size") else ds.voxel_size
        self.pixel_size.setValue(float(size[0]))
        self.pixel_size_col.setValue(float(size[1] if size.size > 1 else size[0]))
        self.pixel_size_slice.setValue(float(size[2] if size.size > 2 else size[0]))
        self.link_size.setChecked(bool(ds.is_isotropic))
        self.tag_spacing.setValue(float(ds.tag_spacing))
        times = ds.times
        # a constant interval read back from the file, so the box agrees with it
        step = 0.0
        if times is not None and times.size > 1:
            diffs = np.diff(times)
            if np.allclose(diffs, diffs[0]):
                step = float(diffs[0])
        self.frame_interval.setValue(step)
        self.override_lbl.setText("")
        self.banner.setText(f"Loaded {path}")
        self._loading = False
        self._refresh_derived()
        self.refresh()
        return True

    def _params_edited(self):
        if self._loading or self.session.dataset is None:
            return
        ds = self.session.dataset
        row = self.pixel_size.value()
        col = row if self.link_size.isChecked() else self.pixel_size_col.value()
        if hasattr(ds, "pixel_size"):
            ds.pixel_size = np.array([row, col], dtype=float)
        else:
            thick = (row if self.link_size.isChecked()
                     else self.pixel_size_slice.value())
            ds.voxel_size = np.array([row, col, thick], dtype=float)
        ds.tag_spacing = self.tag_spacing.value()
        self.override_lbl.setText("parameters overridden by hand")
        self.session.invalidate_from(Stage.FILTER)
        self._refresh_derived()

    def _interval_edited(self):
        """A non-zero interval replaces the frame times with an even spacing."""
        if self._loading or self.session.dataset is None:
            return
        ds = self.session.dataset
        step = self.frame_interval.value()
        if step > 0:
            ds.trigger_times = np.arange(ds.n_frames, dtype=float) * step
        elif ds.trigger_times is not None and np.allclose(
                np.diff(ds.trigger_times), np.diff(ds.trigger_times)[:1]):
            # an evenly spaced set came from the box, so clearing it removes them
            ds.trigger_times = None
        self.session.invalidate_from(Stage.STRAIN)
        self._refresh_times()

    def _refresh_times(self):
        ds = self.session.dataset
        if ds is None:
            self.times_lbl.setText("-")
            return
        times = ds.times
        if times is None:
            self.times_lbl.setText("none, axes will show frame numbers")
        else:
            self.times_lbl.setText(
                f"{times[0]:.0f} to {times[-1]:.0f} ms over {times.size} frames")

    def _update_size_enablement(self):
        """Enablement only. Kept apart from propagation to avoid a cycle."""
        linked = self.link_size.isChecked()
        self.pixel_size_col.setEnabled(not linked)
        is_3d = self.session.dataset is not None and not hasattr(
            self.session.dataset, "pixel_size")
        self.pixel_size_slice.setEnabled(is_3d and not linked)
        for w in (self.slice_row, self.pixel_size_slice):
            w.setVisible(is_3d)

    def _link_changed(self, linked):
        """Follow the row spacing when the pixels are declared square."""
        self._update_size_enablement()
        if linked:
            self._row_size_changed(self.pixel_size.value())

    def _row_size_changed(self, value):
        if self.link_size.isChecked():
            for w in (self.pixel_size_col, self.pixel_size_slice):
                blocked = w.blockSignals(True)
                w.setValue(value)
                w.blockSignals(blocked)
        self._params_edited()

    def _refresh_derived(self):
        ds = self.session.dataset
        if ds is None:
            return
        per_axis = ds.pixels_per_wavelength_per_axis
        if ds.is_isotropic:
            self.ppw_lbl.setText(f"{ds.pixels_per_wavelength:.3f}")
        else:
            self.ppw_lbl.setText(
                f"{ds.pixels_per_wavelength:.3f}   "
                f"(rows {per_axis[0]:.3f}, columns {per_axis[1]:.3f})")
        self.a_lbl.setText(f"{ds.suggested_a:.3f}")
        self._update_size_enablement()
        self._refresh_times()

    # ------------------------------------------------------------- display
    def refresh(self):
        ds = self.session.dataset
        if ds is None:
            return
        fr = min(self.session.current_frame, ds.n_frames - 1)
        img = ds.image
        frame = img[:, :, 0, fr] if img.ndim == 4 else \
            img[:, :, min(self.session.current_slice, img.shape[2] - 1), 0, fr]
        lo, hi = nan_range(frame)
        self.canvas.set_image(np.real(frame), lo, hi)
        self.canvas.fit(force=False)


def _response_edge(response, level):
    """Where a filter response crosses ``level``, as scene polylines.

    contourpy comes with matplotlib, which is already required, and returns
    ``(x, y)`` points, which is the canvas convention. A raster edge would
    stair-step once the panel is zoomed.
    """
    peak = float(np.nanmax(response))
    if not np.isfinite(peak) or peak <= 0:
        return []
    try:
        from contourpy import contour_generator

        return contour_generator(z=np.asarray(response, dtype=float)).lines(
            level * peak)
    except Exception:                                 # noqa: BLE001
        return []


class FilterTab(QWidget):
    """Choose the harmonic bandpass and show its effect in k-space."""

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session

        self.kind = QComboBox(); self.kind.addItems(["transmission", "butterworth"])
        self.cutoff = QDoubleSpinBox(); self.cutoff.setRange(1, 200); self.cutoff.setValue(17)
        self.order = QSpinBox(); self.order.setRange(1, 40); self.order.setValue(10)
        self.kind.currentTextChanged.connect(self._kind_changed)
        for w in (self.cutoff, self.order):
            w.valueChanged.connect(self._controls_changed)
        self.kind_note = QLabel()
        self.kind_note.setWordWrap(True)

        # How strongly the response covers the k-space it is drawn over.
        self.shading = QSlider(Qt.Horizontal)
        self.shading.setRange(0, 100)
        self.shading.setValue(55)
        self.shading_lbl = QLabel("55%")
        self.shading.valueChanged.connect(self._shading_changed)

        # Where the filter effectively stops. The shading says how much is
        # passed; this says where that has fallen to almost nothing.
        self.edge = QDoubleSpinBox()
        self.edge.setRange(0.1, 50.0)
        self.edge.setDecimals(1)
        self.edge.setSingleStep(1.0)
        self.edge.setValue(5.0)
        self.edge.setSuffix(" %")
        self.edge.setToolTip(
            "A line is drawn where the response falls to this fraction of "
            "its peak, as an amplitude rather than a power.")
        self.edge.valueChanged.connect(lambda _: self._draw_preview())

        # Harmonic centre per encoding, in k-space pixels. The wave vector is
        # derived from it, so moving the centre moves the passband for either
        # filter kind.
        self.centres = np.zeros((2, 2))
        self._syncing = False
        self.encoding = QComboBox()
        self.encoding.addItems(["1", "2"])
        self.encoding.currentIndexChanged.connect(self._encoding_changed)

        self.c_row = QDoubleSpinBox(); self.c_row.setDecimals(2); self.c_row.setRange(-1e4, 1e4)
        self.c_col = QDoubleSpinBox(); self.c_col.setDecimals(2); self.c_col.setRange(-1e4, 1e4)
        self.c_row.setSingleStep(0.5); self.c_col.setSingleStep(0.5)
        for w in (self.c_row, self.c_col):
            w.valueChanged.connect(self._centre_edited)

        self.freq_lbl = QLabel("-")
        self.freq_lbl.setWordWrap(True)
        self.source_lbl = QLabel("-")
        self.source_lbl.setWordWrap(True)

        # the k-space of the raw frame depends on no control, so it is kept.
        # The dataset is held by identity rather than by id(): a freed one can
        # have its address reused, and the panel would then show the k-space
        # of the series before last.
        self._k_cache = None       # (encoding, frame, dataset, array)
        self._mag_range = None
        self._dragging = False

        self.auto_btn = QPushButton("Detect peaks again")
        self.auto_btn.clicked.connect(self.detect_centres)
        self.nominal_btn = QPushButton("Reset to nominal")
        self.nominal_btn.clicked.connect(self.reset_centres)

        self.apply = QPushButton("Apply filter")
        self.apply.clicked.connect(self.run)

        form = QFormLayout()
        form.addRow("kind", self.kind)
        form.addRow("cutoff", self.cutoff)
        form.addRow("order", self.order)
        # spanning the row, so the wrapped text is given its full height
        form.addRow(self.kind_note)
        shade = QHBoxLayout()
        shade.addWidget(self.shading)
        shade.addWidget(self.shading_lbl)
        holder = QWidget(); holder.setLayout(shade)
        form.addRow("shading", holder)
        form.addRow("edge at", self.edge)

        cform = QFormLayout()
        cform.addRow("encoding", self.encoding)
        cform.addRow("centre row", self.c_row)
        cform.addRow("centre col", self.c_col)
        cform.addRow("", self.freq_lbl)
        cform.addRow("", self.auto_btn)
        cform.addRow("", self.nominal_btn)
        cbox = QGroupBox("Harmonic centre"); cbox.setLayout(cform)

        box = QGroupBox("Bandpass"); box.setLayout(form)
        self.banner = StageBanner("Load a dataset first.")

        side = QVBoxLayout()
        side.addWidget(box)
        side.addWidget(cbox)
        side.addWidget(self.apply)
        side.addWidget(self.source_lbl)
        side.addWidget(self.banner)
        side.addStretch(1)

        self.k_before = ImageCanvas(session, "k-space")
        self.k_before.markers_draggable = True
        self.k_before.markerMoved.connect(self._centre_dragged)
        self.k_before.markerDragged.connect(self._centre_dragging)
        self.k_after = ImageCanvas(session, "filtered")
        self.mag = ImageCanvas(session, "harmonic magnitude")
        grid = QHBoxLayout()
        for c, lbl in ((self.k_before, "k-space"), (self.k_after, "after filter"),
                       (self.mag, "harmonic magnitude")):
            col = QVBoxLayout()
            t = QLabel(lbl); t.setAlignment(Qt.AlignCenter)
            col.addWidget(t); col.addWidget(c)
            w = QWidget(); w.setLayout(col)
            grid.addWidget(w, 1)

        panel = QWidget(); panel.setLayout(side); panel.setMaximumWidth(300)
        lay = QHBoxLayout(self)
        lay.addWidget(panel)
        lay.addLayout(grid, 1)

        # set, not merely connected: the constructed default is transmission,
        # whose two bandpass boxes do not apply
        self._update_kind_enablement()
        self.k_before.set_overlay_opacity(self.shading.value() / 100.0)

    # ---------------------------------------------------------- the preview

    def _update_kind_enablement(self):
        """Which bandpass controls apply. Enablement only, no side effects.

        Kept apart from the redraw: a method that both sets enablement and
        redraws, reached from a redraw, loops.
        """
        butterworth = self.kind.currentText() == "butterworth"
        self.cutoff.setEnabled(butterworth)
        self.order.setEnabled(butterworth)
        self.kind_note.setText(
            "" if butterworth else
            "The transmission passband is set by the wave vector, so cutoff "
            "and order do not apply to it.")

    def _kind_changed(self, *_):
        self._update_kind_enablement()
        self._controls_changed()

    def _shading_changed(self, value):
        """Strength only, so this neither recomputes nor touches the stage."""
        self.shading_lbl.setText(f"{value}%")
        self.k_before.set_overlay_opacity(value / 100.0)

    def _controls_changed(self, *_):
        """A bandpass control moved: redraw, and say the result is behind."""
        if self._syncing:
            return
        self._draw_preview()
        self.session.note_preview(Stage.FILTER)

    def _settings(self):
        """What the controls say, as the stored result records it.

        Cutoff and order are left out for the transmission window, which
        takes neither, so moving them cannot mark a transmission result as
        out of date.
        """
        kind = self.kind.currentText()
        butterworth = kind == "butterworth"
        return {
            "kind": kind,
            "cutoff": float(self.cutoff.value()) if butterworth else None,
            "order": int(self.order.value()) if butterworth else None,
            "centres": [[round(float(v), 6) for v in c] for c in self.centres],
        }

    def settings_match_result(self):
        """Whether the applied filter is the one the controls describe."""
        res = self.session.result(Stage.FILTER)
        return res is not None and res.get("settings") == self._settings()

    def _build_filter(self):
        """The filter the controls on screen describe.

        Both the preview and Apply go through this, so what is drawn is what
        would be stored.
        """
        ds = self.session.dataset
        if ds is None:
            raise ValueError("Load a dataset first.")
        if ds.image.ndim != 4:
            raise ValueError("3-D filtering is not wired into the GUI yet.")
        if not self.centres.any():
            self.detect_centres()
        if not self.centres.any():
            # detect_centres draws, which comes back here, so this has to
            # stop rather than ask again
            raise ValueError("No harmonic peak was found. Set a centre by hand.")

        shape = ds.image.shape[:2]
        # Derive the frequency and the direction from the chosen centres, so
        # the transmission filter (which uses the wave vectors) and the
        # Butterworth filter (which uses the centre) both follow it.
        omegas = [centre_to_frequency(self.centres[i], shape) for i in range(2)]
        mags = [float(np.hypot(*w)) for w in omegas]
        if min(mags) < 1e-6:
            raise ValueError("A centre sits on DC. Move it onto the harmonic peak.")
        direction = tuple(float(np.arctan2(w[1], w[0])) for w in omegas)
        return HarpFilter2D(ds.image, mags, direction,
                            kind=self.kind.currentText(),
                            cutoff=self.cutoff.value(), order=self.order.value(),
                            centre=np.round(self.centres).astype(int))

    def preview(self):
        """What the controls would produce, from one frame of one encoding.

        Bitwise what Apply would store for that frame, since filter_plane is
        the same arithmetic on the same slice, and 0.6 ms rather than the 52
        ms of the whole series.
        """
        ds = self.session.dataset
        enc = self._current_encoding()
        fr = min(self.session.current_frame, ds.n_frames - 1)
        flt = self._build_filter()

        cached = self._k_cache
        if cached is not None and cached[:2] == (enc, fr) and cached[2] is ds:
            k = cached[3]
        else:
            k = np.log1p(np.abs(itok(ds.image[:, :, enc, fr])))
            self._k_cache = (enc, fr, ds, k)

        Ih = flt.filter_plane(ds.image[:, :, enc, fr], enc)
        return {"filter": flt, "response": flt.response(enc), "k_before": k,
                "k_after": np.log1p(np.abs(itok(Ih))), "magnitude": np.abs(Ih)}

    def _draw_preview(self, markers=True):
        """Paint the three panels from the controls on screen.

        Never from the stored result: the point is to see what the settings
        do before committing them. Guarded, because a handler reached on
        every mouse move must not raise into the event loop.
        """
        if self.session.dataset is None:
            self.banner.setText("Load a dataset first.")
            return
        try:
            pv = self.preview()
        except Exception as exc:                      # noqa: BLE001
            # leave the panels as they were rather than blanking them
            self.banner.setText(str(exc))
            self.source_lbl.setText("-")
            return

        lo, hi = nan_range(pv["k_before"])
        self.k_before.set_image(pv["k_before"], lo, hi)
        self.k_before.set_mask_overlay(pv["response"])
        self.k_before.set_outlines(_response_edge(pv["response"],
                                                  self.edge.value() / 100.0))
        # the same range as the panel beside it, so the suppression is
        # visible instead of being normalised away, and so the greys do not
        # pump while a marker is dragged
        self.k_after.set_image(pv["k_after"], lo, hi)
        if self._mag_range is None or not self._dragging:
            self._mag_range = nan_range(pv["magnitude"])
        self.mag.set_image(pv["magnitude"], *self._mag_range)

        if markers:
            self._draw_markers()
        for c in (self.k_before, self.k_after, self.mag):
            c.fit(force=False)
        self._say_what_is_shown()

    def _say_what_is_shown(self):
        """Which filter is on screen, and whether it is the applied one."""
        res = self.session.result(Stage.FILTER)
        if res is None:
            self.source_lbl.setText("live preview; nothing applied yet")
            self.banner.setText(
                "The panels show the filter the controls describe. Press "
                "Apply filter to use it for the rest of the pipeline.")
            return
        if self.settings_match_result():
            self.source_lbl.setText("showing the filter that was applied")
            return
        stored = res.get("settings") or {}
        centres = stored.get("centres")
        where = "" if centres is None else f" with centres {centres} (row, col)"
        self.source_lbl.setText(
            "preview only; the applied filter is the earlier one")
        self.banner.setText(
            f"The panels show the controls on screen. The applied filter, and "
            f"everything solved from it, is still "
            f"{stored.get('kind', 'the earlier one')}{where}. Press Apply "
            f"filter to replace it.")

    # ------------------------------------------------------- harmonic centre
    def _current_encoding(self):
        return max(self.encoding.currentIndex(), 0)

    def _set_centre_boxes(self, centre):
        self._syncing = True
        try:
            self.c_row.setValue(float(centre[0]))
            self.c_col.setValue(float(centre[1]))
        finally:
            self._syncing = False
        self._update_freq_label()

    def _update_freq_label(self):
        ds = self.session.dataset
        if ds is None:
            return
        shape = ds.image.shape[:2]
        w = centre_to_frequency(self.centres[self._current_encoding()], shape)
        mag = float(np.hypot(*w))
        ang = float(np.degrees(np.arctan2(w[1], w[0])))
        wl = (2 * np.pi / mag) if mag > 1e-9 else np.inf
        self.freq_lbl.setText(
            f"|k| = {mag:.4f} rad/px at {ang:.1f} deg\n"
            f"wavelength {wl:.2f} px  (nominal {ds.pixels_per_wavelength:.2f})"
        )

    def _encoding_changed(self, *_):
        """Show the centre of the encoding now selected."""
        if self.session.dataset is None:
            return
        self._set_centre_boxes(self.centres[self._current_encoding()])
        self.refresh()

    def _centre_edited(self, *_):
        if self._syncing or self.session.dataset is None:
            return
        i = self._current_encoding()
        self.centres[i] = (self.c_row.value(), self.c_col.value())
        self._update_freq_label()
        self._draw_preview()
        self.session.note_preview(Stage.FILTER)

    def _centre_dragging(self, index, x, y):
        """A marker is moving, with the button still down.

        Deliberately narrow. It must not call _draw_markers, because
        set_markers destroys every marker item while the canvas holds an
        index into the list being dragged, and it must not touch the session,
        because that reaches _update_headers and, through any refresh, the
        same rebuild. The release handler does both.
        """
        if self.session.dataset is None:
            return
        self._dragging = True
        self.centres[index] = (y, x)
        if index == self._current_encoding():
            self._set_centre_boxes(self.centres[index])
        self._draw_preview(markers=False)

    def _centre_dragged(self, index, x, y):
        """A marker was dropped, so x is the column and y the row."""
        if self.session.dataset is None:
            return
        self.centres[index] = (y, x)
        if index == self._current_encoding():
            self._set_centre_boxes(self.centres[index])
        self._dragging = False
        self._draw_preview()
        self.session.note_preview(Stage.FILTER)

    def detect_centres(self):
        """Locate the harmonic peaks from the data."""
        ds = self.session.dataset
        if ds is None or ds.image.ndim != 4:
            return
        shape = ds.image.shape[:2]
        measured = self._measured_centres()
        angles = ((ds.meta or {}).get("tag_angles")
                  if measured is not None else None)
        for i in range(2):
            ang = (float(np.asarray(angles)[i]) if angles is not None
                   else (0.0 if i == 0 else np.pi / 2))
            omega = ds.ke[min(i, len(ds.ke) - 1)] * np.array(
                [np.cos(ang), np.sin(ang)])
            # harmonic_center takes (col, row); the measurement is (row, col)
            centre = (None if measured is None
                      else np.round(measured[i][::-1]).astype(int))
            self.centres[i] = harmonic_center(itok(ds.image[:, :, i, 0]), omega,
                                              (2, 2), shape, centre)
        self._set_centre_boxes(self.centres[self._current_encoding()])
        self._draw_preview()
        self.session.note_preview(Stage.FILTER)
        self.banner.setText(f"Peaks at {self.centres.round(1).tolist()} (row, col).")

    def _measured_centres(self):
        """Harmonic peaks measured when the series was read, if any.

        A DICOM export rarely records the tag spacing and the tags are not
        always along the image axes, so the loader measures both. Seeding from
        that beats the nominal 0 and 90 degrees, which is wrong for a tag grid
        prescribed at 45.
        """
        ds = self.session.dataset
        if ds is None:
            return None
        c = (ds.meta or {}).get("tag_centres")
        if c is None:
            return None
        c = np.asarray(c, dtype=float)
        return c if c.shape == (2, 2) else None

    def reset_centres(self):
        """Put the centres back where the data say they are."""
        ds = self.session.dataset
        if ds is None:
            return
        shape = ds.image.shape[:2]
        measured = self._measured_centres()
        if measured is not None:
            self.centres = measured.copy()
            note = "Centres reset to the peaks measured from the images."
        else:
            for i, ang in enumerate((0.0, np.pi / 2)):
                omega = ds.ke[min(i, len(ds.ke) - 1)] * np.array(
                    [np.cos(ang), np.sin(ang)])
                self.centres[i] = frequency_to_centre(omega, shape)
            note = "Centres reset to the nominal tag frequency."
        self._set_centre_boxes(self.centres[self._current_encoding()])
        self._draw_preview()
        self.session.note_preview(Stage.FILTER)
        self.banner.setText(note)

    def _draw_markers(self):
        """Show both centres, the selected one highlighted."""
        cur = self._current_encoding()
        colours = ["#e8503a", "#3a7fe8"]
        self.k_before.set_markers([
            (self.centres[i][1], self.centres[i][0],
             colours[i] if i == cur else "#777777",
             f"encoding {i + 1}")
            for i in range(2)
        ])

    def run(self):
        """Apply the filter the controls describe to the whole series."""
        try:
            flt = self._build_filter()
            Ih = flt(self.session.dataset.image)
        except ValueError as exc:
            self.banner.setText(str(exc))
            return
        except Exception as exc:                      # noqa: BLE001
            self.banner.setText(f"Filtering failed: {exc}")
            return
        self.session.set_result(
            Stage.FILTER,
            {"filter": flt, "image": Ih, "settings": self._settings()})
        self.refresh()
        self.banner.setText(
            "Filtered with centres at "
            f"{np.round(self.centres, 1).tolist()} (row, col)."
        )

    def refresh(self):
        """Draw the controls on screen, whether or not any filter is applied.

        There is one drawing path, and it is the preview: the panels always
        show what the settings would produce, which is the only way to judge
        them before committing.
        """
        if self.session.dataset is None:
            return
        if not self.centres.any():
            self.detect_centres()
            return
        self._draw_preview()
