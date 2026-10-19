"""HARP-I motion estimation.

The method (Mella et al., IEEE TMI 40(4):1240-1252, 2021):

1. take the wrapped harmonic phase of each tag direction (eq. 3);
2. unwrap each frame under the mask, giving the *virtual coordinates*
   ``X_V = Omega*X + noise`` (eq. 6);
3. resolve the per-frame 2*pi offsets left by independent unwrapping
   (eq. 10-12);
4. per frame, fit an RBF whose centres are the deformed-frame virtual
   coordinates and whose values are the pixel positions, regularised as
   ``(Psi + eta*I) lam = G`` (eq. 8);
5. evaluate it at the *reference* frame's virtual coordinates to get the
   deformed position of each reference pixel (eq. 7, Appendix B);
6. subtract the reference grid for a Lagrangian displacement (eq. 9);
7. optionally smooth each trajectory with a 10th-order polynomial.
"""

from __future__ import annotations

import numpy as np

from .consistency import nearest_coordinate_consistency, temporal_phase_consistency
from .interp import distance_matrix, rbf_interpolate, solve_tikhonov
from .rbf import get_kernel
from .smoothing import smooth_reference_phase
from .temporal import fit_temporal
from .types import MotionResult
from .unwrap import phase_quality_2d, unwrap2, unwrap3

# Rules for the per-frame 2*pi offsets; the first is the default.
PHASE_JUMP_RULES = ("nearest", "published")

__all__ = ["estimate_motion_2d", "estimate_motion_3d", "shape_schedule",
           "flip_components"]


def flip_components(result, swap_xy=False, flip_x=False, flip_y=False, flip_z=False):
    """Rearrange the displacement components of a result.

    The sign and the order of the components follow the tag directions given
    to the filter and the row and column convention of the images. When those
    do not match the acquisition, the recovered motion comes out mirrored or
    with x and y exchanged. This applies the correction after the fact, so
    there is no need to solve again.

    Parameters
    ----------
    result
        A :class:`~harpi.types.MotionResult`.
    swap_xy
        Exchange the first two components.
    flip_x, flip_y, flip_z
        Negate a component. ``flip_z`` is ignored for 2-D results.

    Returns
    -------
    A new :class:`~harpi.types.MotionResult`. Trajectories are rebuilt from
    the reference grid so they stay consistent with the displacement.
    """
    from copy import copy

    if not (swap_xy or flip_x or flip_y or flip_z):
        return result

    ndim = result.u.shape[-2]
    order = list(range(ndim))
    if swap_xy and ndim >= 2:
        order[0], order[1] = order[1], order[0]
    signs = np.ones(ndim)
    for i, f in enumerate((flip_x, flip_y, flip_z)[:ndim]):
        if f:
            signs[i] = -1.0

    def apply(u):
        if u is None:
            return None
        out = np.take(u, order, axis=-2)
        shape = [1] * u.ndim
        shape[-2] = ndim
        return out * signs.reshape(shape)

    out = copy(result)
    out.u = apply(result.u)
    out.u_fitted = apply(result.u_fitted)

    if result.trajectories is not None:
        isz = result.u.shape[:-2]
        if ndim == 2:
            Xg, Yg = np.meshgrid(np.arange(1, isz[1] + 1), np.arange(1, isz[0] + 1))
            grid = np.stack([Xg, Yg], axis=-1)
        else:
            Xg, Yg, Zg = np.meshgrid(np.arange(1, isz[1] + 1),
                                     np.arange(1, isz[0] + 1),
                                     np.arange(1, isz[2] + 1))
            grid = np.stack([Xg, Yg, Zg], axis=-1)
        out.trajectories = grid[..., np.newaxis] + out.u

    out.meta = dict(result.meta)
    out.meta["flip"] = {"swap_xy": bool(swap_xy), "flip_x": bool(flip_x),
                        "flip_y": bool(flip_y), "flip_z": bool(flip_z)}
    return out


def shape_schedule(a, dist: str, nfr: int) -> np.ndarray:
    """Per-frame RBF shape parameter.

    Eq. 16 and 17 use a single value for every frame. The other schedules
    are an extension.
    """
    a = np.atleast_1d(np.asarray(a, dtype=float))
    if a.size == 0 or not np.all(np.isfinite(a)):
        raise ValueError("the shape parameter 'a' must be finite and non-empty")
    if a.size == 1:
        return np.full(nfr, float(a[0]))

    lo, hi = float(a.min()), float(a.max())
    key = str(dist).lower().replace("_", "").replace("-", "")
    t = np.linspace(0.0, 1.0, nfr)
    if key == "constant":
        return np.full(nfr, lo)
    if key == "decreasinglinear":
        return hi - (hi - lo) * t
    if key == "linear":
        return lo + (hi - lo) * t
    if key == "quadratic":
        return lo + np.abs((hi - lo) ** 0.4 * np.linspace(-1, 1, nfr)) ** 2.5
    raise ValueError(
        f"unknown a_dist {dist!r}; use constant, linear, decreasinglinear or quadratic"
    )


def _downsample(arr, factor, ndim):
    if factor <= 1:
        return arr
    sl = tuple(slice(None, None, factor) for _ in range(ndim))
    return arr[sl + (Ellipsis,)]


def estimate_motion_2d(
    harmonic_image,
    mask,
    kernel="multiquadric3",
    a=None,
    a_dist="constant",
    eta=15.0,
    connectivity=4,
    seed="auto",
    smooth_reference=False,
    temporal_fit=False,
    temporal_order=10,
    downsample=1,
    reference=0,
    njumps=8,
    progress=None,
    phase_jump_rule="nearest",
):
    """Estimate 2-D Lagrangian displacement from a filtered harmonic image.

    Parameters
    ----------
    harmonic_image
        Complex ``(rows, cols, 2, frames)`` from :class:`~harpi.filters.HarpFilter2D`.
    mask
        ``(rows, cols, frames)`` or ``(rows, cols)``.
    kernel, a, eta
        RBF kernel, shape parameter (eq. 17: ``150 / pixels per wavelength``)
        and Tikhonov weight (eq. 8).
    a_dist
        Per-frame schedule for ``a`` when it is given as a range.
    downsample
        Stride for the interpolation centres, to bound the dense solve.
    progress
        Optional ``callable(frame_index, n_frames, stage)``.
    phase_jump_rule
        ``"nearest"``, the default, chooses the 2*pi offset of each frame with
        :func:`~harpi.consistency.nearest_coordinate_consistency` and ignores
        ``njumps``. ``"published"`` reproduces the MATLAB toolbox with
        :func:`~harpi.consistency.temporal_phase_consistency`, searching
        ``njumps`` steps around its estimate; on ``data_2D.mat`` it puts six of
        eighteen frames one tag period out, so it is kept for reproducing the
        paper rather than recommended.

    Returns
    -------
    :class:`~harpi.types.MotionResult`
    """
    if phase_jump_rule not in PHASE_JUMP_RULES:
        raise ValueError(f"phase_jump_rule must be one of {PHASE_JUMP_RULES}, "
                         f"not {phase_jump_rule!r}")
    Ih = np.asarray(harmonic_image)
    if Ih.ndim != 4 or Ih.shape[2] != 2:
        raise ValueError("harmonic_image must be (rows, cols, 2, frames)")
    isz = Ih.shape[:2]
    nfr = Ih.shape[3]

    mask = np.asarray(mask, dtype=bool)
    if mask.ndim == 2:
        mask = np.repeat(mask[:, :, None], nfr, axis=2)
    if mask.shape[:2] != isz:
        raise ValueError("mask and image disagree on spatial shape")

    if a is None:
        raise ValueError(
            "the shape parameter 'a' is required: eq. 17 uses "
            "150 / (pixels per tag wavelength)"
        )
    a_sched = shape_schedule(a, a_dist, nfr)
    kern = get_kernel(kernel)

    xpha = np.angle(Ih[:, :, 0, :])
    ypha = np.angle(Ih[:, :, 1, :])

    dfac = int(downsample)
    dxpha = _downsample(xpha, dfac, 2).copy()
    dypha = _downsample(ypha, dfac, 2).copy()
    dmask = _downsample(mask, dfac, 2)

    quality = np.full(dxpha.shape, np.nan)
    complete = np.zeros(nfr, dtype=bool)
    xref = yref = None

    for i in range(nfr):
        if progress:
            progress(i, nfr, "unwrap")
        qx = phase_quality_2d(dxpha[:, :, i], dmask[:, :, i], (1.0, 1.0), connectivity)
        quality[:, :, i] = qx
        dxpha[:, :, i], cx = unwrap2(dxpha[:, :, i], dmask[:, :, i], (1.0, 1.0),
                                     seed, connectivity, quality=qx, warn=False)
        dypha[:, :, i], cy = unwrap2(dypha[:, :, i], dmask[:, :, i], (1.0, 1.0),
                                     seed, connectivity, warn=False)
        complete[i] = cx and cy

        if i == reference:
            if dfac != 1:
                # Seed the full resolution reference from the downsampled one.
                uvx = np.full(isz, np.nan)
                uvx[::dfac, ::dfac] = dxpha[:, :, i]
                xref, _ = unwrap2(xpha[:, :, i], mask[:, :, i], (1.0, 1.0), seed,
                                  connectivity, unwrapped_image=uvx, warn=False)
                uvy = np.full(isz, np.nan)
                uvy[::dfac, ::dfac] = dypha[:, :, i]
                yref, _ = unwrap2(ypha[:, :, i], mask[:, :, i], (1.0, 1.0), seed,
                                  connectivity, unwrapped_image=uvy, warn=False)
            else:
                xref = dxpha[:, :, i].copy()
                yref = dypha[:, :, i].copy()

    if progress:
        progress(0, nfr, "consistency")
    if phase_jump_rule == "nearest":
        jumps, jump_info = nearest_coordinate_consistency(
            [dxpha, dypha], reference=reference, downsample=2, return_info=True,
        )
    else:
        jumps, jump_info = temporal_phase_consistency(
            [dxpha, dypha], reference=reference, njumps=njumps, downsample=2,
            return_info=True,
        )
    dxpha = dxpha + jumps[0][None, None, :]
    dypha = dypha + jumps[1][None, None, :]

    if smooth_reference:
        if progress:
            progress(0, nfr, "smooth")
        xref, yref = smooth_reference_phase([xref, yref], kern, a_sched[0] / 5, eta)

    # 1-based grids, as in the MATLAB code.
    dXg, dYg = np.meshgrid(np.arange(1, isz[1] + 1, dfac),
                           np.arange(1, isz[0] + 1, dfac))
    Xg, Yg = np.meshgrid(np.arange(1, isz[1] + 1), np.arange(1, isz[0] + 1))

    tf_ref = ~np.isnan(xref) & ~np.isnan(yref)
    pha_ref = np.stack([xref[tf_ref], yref[tf_ref]], axis=0)
    # C order throughout, so that the flattened grid lines up with the
    # boolean indexing below and centres stay paired with their values.
    G = np.stack([dXg.ravel(), dYg.ravel()], axis=0).astype(float)

    traj = np.full(isz + (2, nfr), np.nan)
    traj[:, :, 0, reference] = np.where(tf_ref, Xg, np.nan)
    traj[:, :, 1, reference] = np.where(tf_ref, Yg, np.nan)

    for i in range(nfr):
        if i == reference:
            continue
        if progress:
            progress(i, nfr, "interpolate")
        tf = ~np.isnan(dxpha[:, :, i]) & ~np.isnan(dypha[:, :, i])
        centres = np.stack([dxpha[:, :, i][tf], dypha[:, :, i][tf]], axis=0)
        values = G[:, tf.ravel()].T

        xy = rbf_interpolate(centres, values, pha_ref, kern, a_sched[i], eta)

        xt = np.full(isz, np.nan); xt[tf_ref] = xy[:, 0]
        yt = np.full(isz, np.nan); yt[tf_ref] = xy[:, 1]
        traj[:, :, 0, i] = xt
        traj[:, :, 1, i] = yt

    u = traj - np.stack([Xg, Yg], axis=-1)[:, :, :, None]

    u_fitted = None
    if temporal_fit:
        if progress:
            progress(0, nfr, "temporal fit")
        u_fitted = fit_temporal(u, mask[:, :, reference], order=temporal_order)

    return MotionResult(
        u=u,
        u_fitted=u_fitted,
        trajectories=traj,
        unwrapped_phase=np.stack([dxpha, dypha], axis=2),
        phase_jumps=jumps,
        quality=quality,
        unwrap_complete=complete,
        meta={"kernel": kern.name, "a": a_sched, "eta": eta, "downsample": dfac,
              "njumps": njumps, "jump_search": jump_info,
              "phase_jump_rule": phase_jump_rule,
              "temporal_order": temporal_order if temporal_fit else None},
    )


def estimate_motion_3d(
    harmonic_image,
    mask,
    kernel="multiquadric3",
    a=None,
    a_dist="constant",
    eta=15.0,
    connectivity=6,
    seed="auto",
    smooth_reference=False,
    downsample=2,
    reference=0,
    njumps=2,
    chunks=4,
    progress=None,
    phase_jump_rule="nearest",
):
    """Estimate 3-D Lagrangian displacement. See :func:`estimate_motion_2d`.

    ``chunks`` splits the evaluation of each frame into interleaved slabs, so
    the dense evaluation matrix is never held whole. The factorisation is
    reused across slabs.

    Both phase jump rules are available, as in 2-D. The nearest-coordinate
    rule is dimension generic, so nothing here is specific to three axes.
    """
    if phase_jump_rule not in PHASE_JUMP_RULES:
        raise ValueError(f"phase_jump_rule must be one of {PHASE_JUMP_RULES}, "
                         f"not {phase_jump_rule!r}")

    Ih = np.asarray(harmonic_image)
    if Ih.ndim != 5 or Ih.shape[3] != 3:
        raise ValueError("harmonic_image must be (rows, cols, slices, 3, frames)")
    isz = Ih.shape[:3]
    nfr = Ih.shape[4]

    mask = np.asarray(mask, dtype=bool)
    if mask.ndim == 3:
        mask = np.repeat(mask[..., None], nfr, axis=3)
    if mask.shape[:3] != isz:
        raise ValueError("mask and image disagree on spatial shape")

    if a is None:
        raise ValueError("the shape parameter 'a' is required")
    a_sched = shape_schedule(a, a_dist, nfr)
    kern = get_kernel(kernel)

    phases = [np.angle(Ih[:, :, :, d, :]) for d in range(3)]
    dfac = int(downsample)
    dpha = [_downsample(p, dfac, 3).copy() for p in phases]
    dmask = _downsample(mask, dfac, 3)

    complete = np.zeros(nfr, dtype=bool)
    ref_pha = [None] * 3

    for i in range(nfr):
        if progress:
            progress(i, nfr, "unwrap")
        ok = True
        for d in range(3):
            dpha[d][..., i], c = unwrap3(dpha[d][..., i], dmask[..., i], (1.0, 1.0, 1.0),
                                         seed, connectivity, warn=False)
            ok &= c
        complete[i] = ok
        if i == reference:
            for d in range(3):
                if dfac != 1:
                    uv = np.full(isz, np.nan)
                    uv[::dfac, ::dfac, ::dfac] = dpha[d][..., i]
                    ref_pha[d], _ = unwrap3(phases[d][..., i], mask[..., i],
                                            (1.0, 1.0, 1.0), seed, connectivity,
                                            unwrapped_volume=uv, warn=False)
                else:
                    ref_pha[d] = dpha[d][..., i].copy()

    if progress:
        progress(0, nfr, "consistency")
    if phase_jump_rule == "nearest":
        jumps, jump_info = nearest_coordinate_consistency(
            dpha, reference=reference, downsample=2, return_info=True,
        )
    else:
        jumps, jump_info = temporal_phase_consistency(
            dpha, reference=reference, njumps=njumps, downsample=2,
            return_info=True,
        )
    for d in range(3):
        dpha[d] = dpha[d] + jumps[d][None, None, None, :]

    if smooth_reference:
        # HARPI_3D.m:204-208 smooths the downsampled phases of every frame, the
        # ones that become the interpolation centres, and leaves the
        # full-resolution reference alone. Smoothing the reference instead, as
        # this did, is both the wrong array and a dense system the size of the
        # volume.
        for i in range(nfr):
            if progress:
                progress(i, nfr, "smooth")
            smoothed = smooth_reference_phase(
                [dpha[d][..., i] for d in range(3)], kern, a_sched[0] / 5, eta)
            for d in range(3):
                dpha[d][..., i] = smoothed[d]

    grids = np.meshgrid(np.arange(1, isz[1] + 1), np.arange(1, isz[0] + 1),
                        np.arange(1, isz[2] + 1))
    Xg, Yg, Zg = grids[0], grids[1], grids[2]
    dgrids = [g[::dfac, ::dfac, ::dfac] for g in (Xg, Yg, Zg)]

    tf_ref = np.ones(isz, dtype=bool)
    for d in range(3):
        tf_ref &= ~np.isnan(ref_pha[d])

    G = np.stack([g.ravel() for g in dgrids], axis=0).astype(float)

    traj = np.full(isz + (3, nfr), np.nan)
    for d, g in enumerate((Xg, Yg, Zg)):
        traj[..., d, reference] = np.where(tf_ref, g, np.nan)

    ref_idx = np.argwhere(tf_ref)
    for i in range(nfr):
        if i == reference:
            continue
        if progress:
            progress(i, nfr, "interpolate")
        tf = np.ones(dpha[0].shape[:3], dtype=bool)
        for d in range(3):
            tf &= ~np.isnan(dpha[d][..., i])
        centres = np.stack([dpha[d][..., i][tf] for d in range(3)], axis=0)
        values = G[:, tf.ravel()].T

        # One factorisation, evaluated over interleaved row slabs.
        psi = kern(distance_matrix(centres), a_sched[i])
        w = solve_tikhonov(psi, values, eta)
        del psi

        for j in range(chunks):
            rows = np.arange(j, isz[0], chunks)
            sel = np.isin(ref_idx[:, 0], rows)
            if not sel.any():
                continue
            pts = ref_idx[sel]
            query = np.stack([ref_pha[d][tuple(pts.T)] for d in range(3)], axis=0)
            vals = kern(distance_matrix(query, centres), a_sched[i]) @ w
            for d in range(3):
                traj[pts[:, 0], pts[:, 1], pts[:, 2], d, i] = vals[:, d]

    u = traj - np.stack([Xg, Yg, Zg], axis=-1)[..., None]

    return MotionResult(
        u=u,
        trajectories=traj,
        phase_jumps=jumps,
        unwrap_complete=complete,
        meta={"kernel": kern.name, "a": a_sched, "eta": eta, "downsample": dfac,
              "njumps": njumps, "jump_search": jump_info,
              "phase_jump_rule": phase_jump_rule},
    )
