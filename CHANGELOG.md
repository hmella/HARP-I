# Changelog

All notable changes to HARP-I are recorded here.

## [Unreleased]

### Fixed — blocking

- **Search-path collision with DENSEanalysis.** DENSEanalysis ships its own,
  incompatible `pixelstrain` (and a copy of `parseinputs`). Whichever directory came
  first on the path won, so following the README's install order made `demo_2D.m` fail
  with the misleading `Invalid spline data`. `HARPIsetup` now enforces
  DENSEanalysis-first / HARP-I-last via the new `tools/harpi_paths.m`, which also
  reports any HARP-I function that is still shadowed.
- **The segmentation GUI could not start.** `src/common/getSegmentation.fig` was deleted
  in commit `5621244`; since `getSegmentation.m` is a GUIDE app with
  `'gui_LayoutFcn', []`, `gui_mainfcn` could not build the window at all. The layout has
  been restored from history (byte-identical to the last committed version).
- **The GUI always produced an empty mask.** `getSegmentation` passed
  `'Contours'`/`'Resolution'` while `getMask` declared `'Contour'`/`'ContourRes'`.
  `parseinputs` silently discards unknown names, so `api.Contour` stayed `[]` and
  `getMask` threw on every call. Both call sites now use the declared names.
- **`getMask` discarded the user's contour flags.** `IsClosed`, `IsCurved` and
  `IsCorner` were hardcoded to `true`, all-`true` and `false`, so corner points placed
  in the GUI had no effect on the mask. They are now read from the contour object.
- **Loading saved contours crashed the GUI.** `handles.contours = {varargin{1}.Contours}`
  wrapped the whole `1xNfr` cell in another cell; the subsequent `.Position` access threw,
  and the half-initialised state then took down the opening function.
- **Masks were shifted one pixel.** The GUI drew on a 0-based grid
  (`meshgrid(0:cols-1, 0:rows-1)`) while `getMask` rasterises on a 1-based grid. The GUI
  now uses a 1-based grid. *Contour `.mat` files saved by earlier versions are offset by
  one pixel relative to masks produced now.*
- **3-D strain was wrong on every slice but the first.** `pixelstrain3D` indexed the 3-D
  angle maps with two subscripts (`ct(i,j)`), collapsing dimensions 2 and 3. Now
  `ct(i,j,k)`.
- **Isolated voxels reported a spurious −0.5 strain in 3-D.** The 2-D code guards the
  whole block with `Nneighbor > 1`; the 3-D code did not, so a voxel with no qualifying
  neighbour produced `Fave = 0` and hence `E = -0.5*I`. The guard has been added.
- **`HARPI_3D` silently reused another subject's phases.** Unwrapped phases, phase jumps
  and smoothed phases were cached in `[pwd '/tmp_*.mat']` and reloaded inside a bare
  `try`/`catch`, keyed only on the file name. All three caches have been removed.
- **`HARPI_3D` gave the reference frame a non-zero displacement.** Its interpolation loop
  started at frame 1 and overwrote the identity initialisation; the 2-D path starts at
  frame 2. Now consistent.
- **`getStrainBySegments` was unusable through its documented interface.** `X`/`Y` were
  only built when the caller *omitted* them, so supplying them left both undefined;
  `errid` was never defined, so both error paths themselves errored; `meshgrid`'s
  arguments were transposed. All fixed, and `'Frames'` (documented in
  `examples/README.md`) is now honoured rather than silently ignored.
- **`TemporalPhaseConsistency` ignored its `njumps` argument.** The default guard tested
  `nargin < 5` although `njumps` is the 4th positional argument, so `HARPI`'s request for
  8 candidate jumps was silently clamped to 1. On the bundled dataset the wider search
  selects the same jumps, so this had been costing robustness margin rather than accuracy.
- **`WindowFilter` could never run** (it referenced a non-existent `image_size` property)
  and nothing used it. Removed.
- **`HARPIsetup` destroyed the user's `startup.m`.** It opened the file with `'wt'` and
  wrote a single line. It now appends, and backs up any existing file to `startup.m.bak`.

### Fixed — correctness

- `HARPI` used a leaked loop index when smoothing the reference phase
  (`smooth_f = a(i)/5`, where `i` was left over from the unwrapping loop and equalled
  `Nfr`); it now uses `a(1)`.
- The `RBFFacDist` schedule is extracted into `src/common/rbfShapeSchedule.m`,
  matched case-insensitively, and validated. Previously the documented default
  `'quadratic'` matched none of the case-sensitive branches, leaving `a` at its input
  length so that `a(i)` threw for `i >= 3`. The `'Linear'` schedule ramped
  `max -> 2*max-min` and never visited `min`; it now ramps `min -> max`.
- Default `Mask` values had the wrong rank in `HARPI`, `HARPI_3D` and `TemporalFitting`
  (`true(size(Ih))` is 4-D/5-D), so `mask(:,:,fr)` indexed the encoding channel instead
  of the frames.
- `maskSA` hardcoded exactly two contours while `getMask` looped over all of them.
  It now handles any number, and a contour drawn outside the outer boundary can no
  longer leak pixels into the mask.
- `getStrainBySegments` now intersects each angular wedge with the mask (segment means
  previously averaged over a full-image wedge), wraps `theta` into the rotated frame (a
  wedge of width `|angle0|` used to fall outside every bin), defaults `PositionB` to
  `origin(2)` rather than `origin(1)`, and exposes the previously hardcoded 0.5 outlier
  cut as `'OutlierThreshold'`.
- `HARPI_3D` no longer reports every temporal-fitting failure as "not yet implemented".
  `TemporalFitting` returns only `[dxr, dyr]`, so the 3-D path referenced an undefined
  `dz`; it now raises `HARPI_3D:temporalFittingUnsupported`.

### Changed

- **Unknown parameters are no longer silent.** `src/common/assertNoUnknownArgs.m` is
  called at all nine `parseinputs` sites and warns about unrecognised names; set the
  global `HARPI_STRICT_ARGS` to `true` to make it an error. `parseinputs` itself is
  deliberately unchanged, because it shadows an identical copy in DENSEanalysis whose
  callers pass extra arguments on purpose.
- **Dead options removed**: `'FOV'` and `'PixelSize'` from `HARPI`/`HARPI_3D` (never
  read; `[1 1]` is hardcoded into the `unwrap2` calls), `'ROI'` from `HARPI` (reachable
  only from an unreachable debug block), `'Connectivity'` from `HARPI_3D` (`unwrap3` has
  no such parameter), `'Butterworth_cuttoff'` and `'Gabor_lambda'` from `HARPFilter3D`
  (stored but never used — `Butterworth23D` computes its own `0.8/dGrid` cutoff). The
  demos were updated in the same change.
- **`demo_2D.m` no longer mislabels the strain it displays.** It assigned the principal
  strains `st.p1`/`st.p2` to variables *named* `RR`/`CC`, while `strain_curves2D.m` reads
  the genuine `st.RR`/`st.CC` from the same file — so the animation and the curves showed
  different quantities under the same labels. Both are now computed under their own names
  (`RR`, `CC`, `P1`, `P2`) and a `STRAIN_FOR_DISPLAY` switch selects what is animated.
  It defaults to `'p2'` so the committed GIFs stay reproducible.

  The two conventions are **not** interchangeable. On the bundled 2-D dataset:

  | quantity | peak (global mean) | max difference | relative |
  |---|---|---|---|
  | `E_cc` vs 2nd principal | 0.169 | 0.063 | 37% |
  | `E_rr` vs 1st principal | 0.080 | 0.063 | 78% |

  Pixelwise the difference reaches 0.335. `examples/compare_principal_vs_rotated.m`
  reproduces this. The paper reports the rotated components, so `'CC'` is the setting
  that matches the publication.

### Performance

Measured on the bundled datasets, single-threaded (`maxNumCompThreads(1)`).

| | before | after | |
|---|---|---|---|
| `HARPI` (2-D demo, 160x160x18) | 2.29 s | **1.02 s** | 2.2x |
| — `TemporalFitting` | 1.17 s | 0.18 s | 6.5x |
| — `TemporalPhaseConsistency` (njumps=8) | 0.35 s | 0.018 s | 19x |
| `HARPI_3D` (112³x3, per frame) | 20.0 s | **15.9 s** | 1.3x |
| `pixelstrain3D` (3 frames) | — | 1.9 s | |

Every change below is bitwise identical on the regression reference except the
temporal fit, which is a documented correction.

- **`TemporalPhaseConsistency` reduced to convex hulls.** The jump search built a dense
  `M x N` distance matrix for each of `(2*njumps+2)^2` candidates per frame and reduced
  each to a single `max()`. The maximum distance between two point sets is attained at
  convex-hull vertices, and a candidate jump only translates one set, so the hulls are
  built once per frame. `test_maxPairDistance` checks the reduction against brute force
  over 360 randomised 2-D and 3-D point sets, including collinear and duplicate cases.
- **`TemporalFitting` solved in batch.** Six `lsqlin` calls per masked pixel became one
  null-space least-squares solve for all pixels at once, with `Z`, `pinv(Aeq)` and the QR
  factorisation of `A*Z` computed once. The null-space route avoids KKT normal equations,
  which would square an already poor condition number (`cond(A) = 6.4e7` but
  `cond(A'A) = 4.0e15`). **The Optimization Toolbox is no longer required.**

  *This changes results, and the change is a correction.* Over 200 sampled pixels the new
  solution has a strictly lower residual than lsqlin's on every one (mean 4.065 vs 4.265);
  lsqlin was stopping short of the minimum, and tightening its tolerances does not close
  the gap. Fitted displacements move by up to **1.12 px (rms 0.23 px)** and the global
  curves by up to 0.010 (GCS) and 0.015 (GRS).
- **`itok`/`ktoi` shift only the transformed dimension.** Called without a `dim`
  argument they rotated every axis on every pass; the pairs cancelled, so results were
  correct, but `Dims=[1 2 3]` cost six full-array circshifts instead of two.
- **The k-space filter is broadcast, not replicated.** `repmat` tiled the filter across
  frames purely to line up an elementwise multiply — about 0.7 GB of avoidable
  allocation for the 3-D dataset. This also removes a latent bug: `HARPFilter` sized the
  copy from the *constructor's* frame count, so reusing a filter on a stack of different
  length mis-sized it.
- **Distance matrices no longer build unused components.** `distanceMatrix2dNew`/`3dNew`
  always materialised `rx`, `ry` (and `rz`) although every caller binds only `r`.
- **`RefPhaseSmoothing` builds one matrix instead of four.** Centres and evaluation
  points coincide, so `RBFInterp2D/3D` was constructing the distance matrix and the
  kernel twice each; written directly it is `B*((B + eta*I)\f)`.
- **Strain loops iterate over qualifying pixels.** The frame-independent mask/neighbour
  test was re-evaluated for every pixel of every frame: 460800 visits for ~500 pixels in
  2-D, 1.4 million voxels per frame for ~100 thousand in 3-D.
- `getStrainBySegments` sliced each frame once per segment; `getHarmonicCenter` used a
  hand-written double argmax that left its indices undefined for an all-zero window.

### Added

- `tools/harpi_paths.m` — path setup with collision reporting.
- `tools/regression/capture_reference.m` / `compare_reference.m` — stage-by-stage
  snapshots of the 2-D pipeline and a field-by-field comparator (including an
  integer wrap-count check on unwrapped phase). Single-threaded runs are bitwise
  reproducible, so the tolerance for "should not change" fixes is exactly zero.
- `tests/` — 26 MATLAB unit tests (`runtests('tests')`).
- `src/common/extremePoints.m`, `src/common/maxPairDistance.m` — hull-reduced maximum
  distance between point sets.
- `src/common/rbfShapeSchedule.m`, `src/common/assertNoUnknownArgs.m`.
- `examples/compare_principal_vs_rotated.m`.

### Known issues

- Fixing the `TemporalPhaseConsistency` `njumps` off-by-one multiplies its cost by ~324×
  (the search is `(2*njumps+2)^2` per frame), so the hull-based optimisation lands first.
- Temporal fitting is not implemented for 3-D data.
- `HARPI`'s `UnwrappedPhases` output is at the downsampled resolution while every other
  output is full resolution.
