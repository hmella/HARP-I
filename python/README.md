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

*Track contours through frames* saves drawing the same borders on every
frame. Propagation copies a contour unchanged, but the tissue moves: on
`data_2D.mat` a node travels up to 4.9 px against a wall about 8 px thick, so
by mid systole a copied contour is half a wall out. Tracking instead follows
each node through the harmonic phase, which labels tissue rather than space.
Draw on whichever frame is clearest and press the button: it runs outwards in
both directions and fills the whole cycle, as one undoable action. A window
reports which contour and frame it is on and cancels the run, which then leaves
every frame as it was rather than writing a half-tracked sequence. Nodes land
within about 0.2 px of the motion solve on this dataset, and no node past half
a tag period.

It needs the filter to have run, and nothing else. In particular it needs no
mask, which is what lets it help with producing one, and no motion solve,
which would be circular. A few limits are worth knowing. The tags fade through
the cycle, so late diastole is where the phase is weakest and a node can lock
onto the neighbouring tag line; a guard compares each node with the rest of
the contour and pulls back the ones that slipped, but it is a heuristic and
the result is a starting point to look over, not a final answer. The
endocardium is the harder border, since it meets the blood pool where the tags
wash out. It is 2-D only, it depends on the harmonic centres on the Filtering
tab being right, and it follows the nodes rather than the curve between them,
so a contour drawn with very few nodes describes the tissue only as well as
those nodes do.

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
Motion tab chooses how it is found.

The *nearest coordinate* rule is the default. It starts from the difference of
the mean virtual coordinates and keeps the offset that brings the coordinates of
the frame closest to reference ones, using all of them.

The *published* rule is the one in the MATLAB toolbox. It reduces each frame to
the largest distance between its virtual coordinates and the reference ones, and
searches the number of steps set in *jump search*. That distance changes by only
a fraction of a period between neighbouring candidates when the mask of a frame
does not hold exactly the tissue of the reference mask, so a whole frame can
land one period out. On the bundled `data_2D.mat` it does: six of eighteen
frames are a period off in x, and the median x displacement jumps by 3.96 px
against a 4 px tag period, where the true excursion is about 1.3 px. The nearest
rule keeps the largest step to 0.69 px.

It is kept because it reproduces the published results and because the
zero-tolerance MATLAB regression is written against it, not because it is
recommended.

Worth knowing before reading too much into the change: strain is a spatial
derivative, so a whole-period offset of one frame largely cancels in it. On
`data_2D.mat` the two rules give global circumferential strain within about six
per cent of peak of each other. The difference that matters is in the
displacement field and in anything built directly on it, such as trajectories
or a through-plane measure.

Either rule is available in 2-D and 3-D, from the library as

```python
res = harpi.estimate_motion_2d(flt(ds.image), ds.mask, a=ds.suggested_a,
                               phase_jump_rule="published")
```

and from the shell as `harpi process --phase-jump-rule published`. The rule used
is recorded in `res.meta["phase_jump_rule"]` and printed by `harpi process`.

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

The Motion and Strain maps are drawn over the tagged image, so the numbers
sit on the anatomy they describe rather than on a blank background. The
**Background** box turns it off, chooses which image goes underneath, and
sets how strongly the data covers it. The arrows stay solid whatever the
opacity, being line art that hides almost nothing.

On the Strain tab the AHA segments are oriented by the two
right-ventricular insertion points, anterior and inferior, each set
numerically or by clicking the map. The anterior one is the boundary
between segment 1, anterior, and segment 2, anteroseptal, which is what the
AHA statement means by using the insertion to separate the septum from the
anterior wall. The pair also fixes which way the numbering runs, the septum
being the shorter arc between them, so the *clockwise* box applies only
when a single point is given.

The wedges stay equiangular, 60 degrees for six and 90 for four, as the
model specifies. The points orient it, they do not stretch it, so what the
two actually span is reported beside them against the ideal. A departure
beyond about 45 degrees usually means the two points are the wrong way
round, and it says so.

Before this the segments were numbered from half a turn away from the
insertion. The wedges and their boundaries were right, but segment 1 sat
nearly opposite the ray, so every name landed on the wrong wall and "mid
anterior" reported the inferolateral segment. Results read off a bullseye
produced by an earlier version need that rotation applied before they are
compared with these.

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
