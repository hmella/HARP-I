# harpi

A Python implementation of **HARP-I**, which recovers cardiac motion and
strain from tagged MR images by interpolating unwrapped harmonic phase with
radial basis functions.

> H. Mella, J. Mura, H. Sanchez, S. Uribe, C. Bertoglio, J. Sotelo,
> "HARP-I: A Harmonic Phase Interpolation Method for the Estimation of Motion
> From Tagged MR Images", *IEEE TMI* 40(4):1240-1252, 2021.
> [doi:10.1109/TMI.2021.3051092](https://doi.org/10.1109/TMI.2021.3051092)

It accompanies the MATLAB toolbox in the parent directory and needs neither
MATLAB nor DENSEanalysis: the quality-guided phase unwrapper and the contour
handling are implemented here.

## Install

```sh
pip install -e python                # library and CLI
pip install -e 'python[gui]'         # adds the PySide6 interface
pip install -e 'python[dev]'         # adds pytest, ruff, mypy
```

## Use

```python
import harpi

ds  = harpi.load_mat_2d("data/data_2D.mat")
flt = harpi.HarpFilter2D(ds.image, ds.ke, kind="transmission")

res = harpi.estimate_motion_2d(
    flt(ds.image), ds.mask,
    kernel="multiquadric3",
    a=ds.suggested_a,          # eq. 17: 150 / pixels per wavelength
    eta=15.0,                  # Tikhonov weight of eq. 8
    smooth_reference=True, temporal_fit=True,
)

st  = harpi.pixel_strain_2d(res.best, ds.mask[..., 0])
gcs = st.global_curve("CC")    # global circumferential strain
seg = harpi.strain_by_segments(st, n_segments=6)
```

From the shell:

```sh
harpi info    data/data_2D.mat
harpi process data/data_2D.mat --smooth-reference --temporal-fit --out result.mat
harpi gui
```

`--out` writes the field names used by the MATLAB toolbox, so
`examples/strain_curves2D.m` reads the output unchanged.

## The interface

`harpi gui` opens one window with six tabs, in the order the data flows:

**Data, Filtering, Segmentation, Unwrapping, Motion, Strain**

Segmentation comes after filtering because the borders are drawn on the
harmonic phase, which only exists once the filter has run.

Each tab header shows the state of that stage, and changing a setting clears
every stage after it.

On the **Filtering** tab the passband is shaded on the k-space panel and its
centre is a draggable marker, one per encoding direction. Drag it onto the
harmonic peak, or type the row and column, and the readout gives the
resulting frequency, angle and tag wavelength next to the nominal one.
*Detect peaks again* re-runs the search and *Reset to nominal* puts the
centre back on the tag frequency. The wave vector is derived from the
centre, so the transmission and the Butterworth filter both follow it.

The segmentation tab supports undo and redo (Ctrl+Z, Ctrl+Shift+Z), deleting
points with a right click, scroll wheel zoom and drag to pan, keyboard frame
and slice navigation, a live mask overlay, and propagation over a chosen
frame range as a single undoable action.

Each contour is labelled with the wall it outlines, chosen from the **wall**
box: LV or RV, epicardium or endocardium. Contours are grouped per structure
and each contributes the region between its outer and inner boundary, so an
LV and an RV outline give two separate annuli rather than one. Each wall is
drawn in its own colour. Contours left as *Other* keep the simpler
convention, where the first is the outer boundary and the rest are holes.

The **draw on** box selects what the editable panel shows: the raw
magnitude, the harmonic magnitude, or either harmonic phase. Drawing on the
harmonic magnitude is usually easier, since the tag pattern is suppressed
and the wall stands out.

If a mask already exists there is no need to draw at all. **Use mask from
dataset** takes the one stored in the file, and **Load mask from file** reads
one from any `.mat`, `.npy` or `.npz`, picking the first 2-D or 3-D binary
array it finds and preferring anything named `mask`. A 2-D mask is applied to
every frame. The box below the buttons says where the current mask came from;
drawing a contour replaces it with one built from the contours.

Two things the scripts left invisible are shown here: whether unwrapping
completed on each frame, and the phase jumps picked by the consistency step.

The phase jump of a frame is the whole number of periods, one per tag
direction, that unwrapping leaves undetermined. The **jump rule** box on the
Motion tab chooses how it is found. The *published* rule, the default and the
one in the MATLAB toolbox, compares the largest distance between the virtual
coordinates of the frame and those of the reference, searching the number of
steps set in *jump search*. The *nearest coordinate* rule starts from the
difference of the mean coordinates and keeps the offset that brings the
coordinates of the frame closest to reference ones. It is steadier when the
mask of a frame does not hold exactly the tissue of the reference mask. It is
available for 2-D data in the Python port only, from the library as

```python
res = harpi.estimate_motion_2d(flt(ds.image), ds.mask, a=ds.suggested_a,
                               phase_jump_rule="nearest")
```

and from the shell as `harpi process --phase-jump-rule nearest`.

The Motion tab has a **Displacement components** box with *swap x and y*,
*flip x* and *flip y*. The sign and the order of the components follow the
tag directions given to the filter and the row and column convention of the
images; when those do not match the acquisition, the recovered motion comes
out mirrored or with the axes exchanged, and the strain collapses. The boxes
correct it without solving again, and the same transform is available from
the library:

```python
res = harpi.flip_components(res, swap_xy=True, flip_y=True)
```

After an estimate the view moves to the frame with the most motion, since
the reference frame has none by definition.

Long solves run on a worker thread with a progress bar. Cancel takes effect
at the next frame boundary.

### Shortcuts

| key | action |
|---|---|
| `[` `]` | previous / next frame |
| `,` `.` | previous / next slice |
| Space | play / pause |
| `M` | toggle the mask overlay |
| `F` | fit the view |
| Ctrl+Z, Ctrl+Shift+Z | undo, redo |
| Ctrl+S | save the project |
| Enter, Esc | close or cancel a contour |

## Conventions

* Arrays keep the MATLAB layout, `(row, col[, slice], encoding, frame)`.
* Pixel coordinates are 0-based and converted at the I/O boundary. Contour
  rasterisation takes an explicit `origin`, where 1 gives the MATLAB grid.
* `NaN` marks pixels outside the mask or not reached by unwrapping.

## Tests

```sh
cd python && PYTHONPATH=src python3 -m pytest tests -q
```

Set `QT_QPA_PLATFORM=offscreen` to run the GUI tests without a display, and
`HARPI_NO_NUMBA=1` to exercise the plain NumPy unwrapper.

Tests that compare against MATLAB skip unless a reference file is present.
See `tests/matlab/README.md` for how to produce one.

## Licence

MPL-2.0, as for the MATLAB toolbox. `harpi/unwrap/` and
`harpi/segmentation/spline.py` derive from DENSEanalysis, which is also
MPL-2.0. Their module docstrings carry the attribution.
