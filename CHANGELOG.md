# Changelog

## [Unreleased]

### Fixed

- **Search path collision with DENSEanalysis.** DENSEanalysis ships an
  incompatible `pixelstrain` and its own `parseinputs`. Whichever directory
  comes first on the path wins, so the install order in the README could make
  `demo_2D.m` fail with `Invalid spline data`. `HARPIsetup` now adds
  DENSEanalysis first and HARP-I last through the new `tools/harpi_paths.m`,
  which also reports any function that is still shadowed.
- **The segmentation GUI could not start.** `getSegmentation.m` is a GUIDE
  app declared with `'gui_LayoutFcn', []`, so `gui_mainfcn` needs
  `getSegmentation.fig`. That file had been deleted; it is restored.
- **The GUI produced an empty mask.** `getSegmentation` passed
  `'Contours'`/`'Resolution'` while `getMask` declares
  `'Contour'`/`'ContourRes'`. Unknown names are dropped by `parseinputs`, so
  `api.Contour` stayed empty and `getMask` threw on every call.
- **`getMask` ignored the contour flags.** `IsClosed`, `IsCurved` and
  `IsCorner` were hardcoded, so corner points placed in the GUI had no
  effect on the mask. They are now read from the contour object.
- **Loading saved contours crashed the GUI.** The stored cell array was
  wrapped in another cell, and the following `.Position` access threw during
  figure opening.
- **Masks were shifted by one pixel.** The GUI drew on a 0-based grid while
  `getMask` rasterises on a 1-based one. Contour `.mat` files saved by
  earlier versions are offset by one pixel relative to masks produced now.
- **3-D strain was wrong away from the first slice.** `pixelstrain3D`
  indexed the 3-D angle maps with two subscripts, collapsing dimensions 2
  and 3.
- **Isolated voxels reported a strain of -0.5 in 3-D.** The 2-D code guards
  the block with `Nneighbor > 1`; the 3-D code did not, so a voxel with no
  qualifying neighbour gave `Fave = 0` and `E = -0.5*I`. The guard is added.
- **`HARPI_3D` reused phases across datasets.** Unwrapped phases, phase
  jumps and smoothed phases were cached in `[pwd '/tmp_*.mat']` behind bare
  try/catch blocks, keyed only on the file name. The caches are removed.
- **`HARPI_3D` gave the reference frame a non-zero displacement.** Its
  interpolation loop started at frame 1 and overwrote the identity
  initialisation. It now starts at frame 2, as the 2-D path does.
- **`getStrainBySegments` did not work through its documented interface.**
  `X` and `Y` were only built when the caller omitted them, `errid` was
  never defined so both error paths themselves errored, and the `meshgrid`
  arguments were transposed. Segment means now intersect the wedge with the
  mask, `theta` is wrapped into the rotated frame so no wedge falls outside
  every bin, `PositionB` defaults to `origin(2)`, the outlier cut is exposed
  as `'OutlierThreshold'`, and `'Frames'` is honoured.
- **`TemporalPhaseConsistency` ignored `njumps`.** The default guard tested
  `nargin < 5` although `njumps` is the fourth positional argument, so a
  call that omitted the optional `ROI` had its jump budget reset to 1.
- **`HARPI` used a leaked loop index** when smoothing the reference phase.
- **Default `Mask` values had the wrong rank** in `HARPI`, `HARPI_3D` and
  `TemporalFitting`, so `mask(:,:,fr)` indexed the encoding channel instead
  of the frames.
- **`maskSA` handled only two contours** while `getMask` looped over all of
  them. It now takes any number, and a contour outside the outer boundary
  can no longer add pixels to the mask.
- **`HARPIsetup` overwrote `startup.m`.** It now appends, and backs up any
  existing file to `startup.m.bak`.
- **`WindowFilter` could not run** and nothing used it. Removed.
- `HARPI_3D` reported every temporal fitting failure as "not yet
  implemented", because the block referenced an undefined `dz` inside a
  try/catch. It now raises `HARPI_3D:temporalFittingUnsupported`.

### Changed

- The RBF shape schedule moved to `src/common/rbfShapeSchedule.m`, is
  matched case insensitively, and raises an error on an unknown name. The
  documented default `'quadratic'` previously matched none of the case
  sensitive branches. The `'Linear'` schedule now ramps from the lower to
  the upper bound.
- Unrecognised parameters are reported by `src/common/assertNoUnknownArgs.m`
  at all nine `parseinputs` call sites. Set the global `HARPI_STRICT_ARGS`
  to turn the warning into an error. `parseinputs` itself is unchanged,
  since it shadows an identical copy in DENSEanalysis whose callers pass
  extra arguments on purpose.
- Removed options that had no effect: `'FOV'` and `'PixelSize'` from
  `HARPI` and `HARPI_3D`, `'ROI'` from `HARPI`, `'Connectivity'` from
  `HARPI_3D`, and `'Butterworth_cuttoff'` and `'Gabor_lambda'` from
  `HARPFilter3D`. The demos were updated in the same change.
- `demo_2D.m` assigned the principal strains `st.p1` and `st.p2` to
  variables named `RR` and `CC`, while `strain_curves2D.m` reads the rotated
  components from the same file. All four are now computed under their own
  names and a `STRAIN_FOR_DISPLAY` switch selects what the animation shows.
  It defaults to `'p2'`, which is what the committed GIFs show; `'CC'` is
  the quantity reported in the paper. The two differ by a substantial
  margin, which `examples/compare_principal_vs_rotated.m` quantifies.

### Performance

The changes below leave the regression reference unchanged, except for the
temporal fit.

- `TemporalPhaseConsistency` reduces both point sets to their convex hulls.
  The criterion needs only the largest distance between two sets, which is
  reached at a hull vertex, and a candidate jump translates one set, so the
  hulls are built once per frame rather than once per candidate.
- `TemporalFitting` solves all pixels in one pass. Six `lsqlin` calls per
  masked pixel became one null space least squares solve, with `Z`,
  `pinv(Aeq)` and the QR factorisation of `A*Z` computed once. The null
  space route avoids the normal equations, which would square an already
  poor condition number. The Optimization Toolbox is no longer required.

  This changes the result. The new solution has a lower residual than the
  one `lsqlin` returned, which was stopping short of the minimum on this
  problem. Fitted displacements and the global strain curves move slightly.
- `itok` and `ktoi` shift only the dimension being transformed. Without a
  `dim` argument they rotated every axis on every pass.
- The k-space filter is broadcast over frames rather than replicated with
  `repmat`. This also removes a latent bug: `HARPFilter` sized the copy from
  the frame count of the image passed to the constructor.
- `distanceMatrix2dNew` and `distanceMatrix3dNew` build the component
  differences only when they are requested.
- `RefPhaseSmoothing` builds one kernel matrix instead of four, since its
  centres and evaluation points coincide.
- The strain loops evaluate the mask and neighbour test once instead of once
  per frame, and iterate only over the pixels that pass.
- `getStrainBySegments` slices each frame once rather than once per segment.
  `getHarmonicCenter` uses `max` and `ind2sub` instead of a double loop that
  left its indices undefined for an all-zero window.

### Added

- `tools/harpi_paths.m`, for path setup with collision reporting.
- `tools/regression/capture_reference.m` and `compare_reference.m`, which
  snapshot each stage of the 2-D pipeline and compare two snapshots field by
  field, including an integer wrap count check on the unwrapped phase.
- `tests/`, with MATLAB unit tests. Run `runtests('tests')`.
- `src/common/extremePoints.m` and `src/common/maxPairDistance.m`.
- `src/common/rbfShapeSchedule.m` and `src/common/assertNoUnknownArgs.m`.
- `examples/compare_principal_vs_rotated.m`.

### Known issues

- Temporal fitting is not implemented for 3-D data.
- `HARPI` returns `UnwrappedPhases` at the downsampled resolution while the
  other outputs are at full resolution.
