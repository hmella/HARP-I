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

The segmentation tab supports undo and redo (Ctrl+Z, Ctrl+Shift+Z), deleting
points with a right click, scroll wheel zoom and drag to pan, keyboard frame
and slice navigation, any number of contours with the first as the outer
boundary, a live mask overlay, and propagation over a chosen frame range as
a single undoable action.

Two things the scripts left invisible are shown here: whether unwrapping
completed on each frame, and the phase jumps picked by the consistency step.

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
