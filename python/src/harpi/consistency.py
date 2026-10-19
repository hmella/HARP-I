# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

"""Temporal phase consistency.

Each frame is unwrapped on its own, so its phase carries an unknown additive
multiple of 2*pi (the constant ``c_k`` of eq. 10 to 12). These are recovered
by matching the spread of the virtual coordinates against the reference
frame.

The search runs on convex hulls. The criterion only needs the largest
distance between two point sets, which is always reached at hull vertices,
and a candidate jump just translates one of the sets.

:func:`nearest_coordinate_consistency` is an alternative rule that compares
every virtual coordinate of a frame with the reference ones. It is available
in the Python port only, and :func:`temporal_phase_consistency` stays the
default everywhere.
"""

from __future__ import annotations

from itertools import product

import numpy as np
from scipy.spatial import ConvexHull, cKDTree

__all__ = ["temporal_phase_consistency", "nearest_coordinate_consistency",
           "max_pair_distance", "extreme_points"]


def extreme_points(points: np.ndarray) -> np.ndarray:
    """Convex hull vertices of a ``(D, N)`` point set.

    Returns all points when no hull exists, for example when they are
    collinear or duplicated.
    """
    p = np.asarray(points, dtype=float)
    d, n = p.shape
    if n <= d + 1:
        return p
    try:
        hull = ConvexHull(p.T)
        return p[:, np.unique(hull.vertices)]
    except Exception:
        return p


def max_pair_distance(a: np.ndarray, b: np.ndarray) -> float:
    """Largest Euclidean distance between two ``(D, K)`` point sets."""
    diff = a[:, :, None] - b[:, None, :]
    return float(np.sqrt((diff**2).sum(axis=0)).max())


def _downsample(arr: np.ndarray, factor: int, ndim: int) -> np.ndarray:
    sl = tuple(slice(None, None, factor) for _ in range(ndim)) + (slice(None),)
    return arr[sl]


def temporal_phase_consistency(phases, reference=0, njumps=8, downsample=2,
                               return_info=False):
    """Resolve the per-frame 2*pi offsets between unwrapped frames.

    Parameters
    ----------
    phases
        Sequence of ``ndim`` arrays, each ``(..., frames)`` of unwrapped phase
        with NaN outside the mask.
    reference
        Index of the reference frame.
    njumps
        Half width of the search around the coarse estimate.
    downsample
        Stride applied before searching.

    return_info
        Also return a dict describing the search, with ``delta`` (the coarse
        estimate per frame) and ``offset`` (how many steps the chosen jump
        sits away from it). An offset that never reaches ``njumps`` means
        the budget is not binding and raising it will change nothing.

    Returns
    -------
    ``(ndim, frames)`` array of additive phase offsets, each a multiple of
    2*pi, and optionally the search info.
    """
    phases = [np.asarray(p, dtype=float) for p in phases]
    ndim = len(phases)
    spatial = phases[0].ndim - 1
    nfr = phases[0].shape[-1]

    if downsample and downsample > 1:
        phases = [_downsample(p, downsample, spatial) for p in phases]

    valid = np.ones(phases[0].shape, dtype=bool)
    for p in phases:
        valid &= ~np.isnan(p)

    def points(fr):
        sel = valid[..., fr]
        return np.stack([p[..., fr][sel] for p in phases], axis=0)

    ref_pts = points(reference)
    if ref_pts.shape[1] == 0:
        raise ValueError("the reference frame has no valid pixels")
    ref_hull = extreme_points(ref_pts)
    ref_max = max_pair_distance(ref_hull, ref_hull)
    ref_peak = ref_pts.max(axis=1)

    jumps = np.zeros((ndim, nfr))
    coarse = np.zeros((ndim, nfr))
    offset = np.zeros((ndim, nfr))
    for fr in range(nfr):
        if fr == reference:
            continue
        pts = points(fr)
        if pts.shape[1] == 0:
            continue
        hull = extreme_points(pts)

        # Coarse estimate, then a small search around it.
        delta = np.floor((ref_peak - pts.max(axis=1)) / (2 * np.pi))
        candidates = [
            np.concatenate(([0.0], 2 * np.pi * np.arange(d - njumps, d + njumps + 1)))
            for d in delta
        ]

        best = (np.inf, np.zeros(ndim))
        for combo in np.ndindex(*(len(c) for c in candidates)):
            shift = np.array([candidates[d][combo[d]] for d in range(ndim)])
            r = max_pair_distance(hull + shift[:, None], ref_hull)
            diff = abs(r - ref_max)
            if diff < best[0]:
                best = (diff, shift)

        jumps[:, fr] = best[1]
        coarse[:, fr] = delta
        offset[:, fr] = np.round(best[1] / (2 * np.pi)) - delta

    if return_info:
        return jumps, {"delta": coarse, "offset": offset, "njumps": int(njumps)}
    return jumps


def nearest_coordinate_consistency(phases, reference=0, radius=2, downsample=2,
                                   return_info=False):
    """Resolve the per-frame 2*pi offsets with the nearest-coordinate rule.

    :func:`temporal_phase_consistency` reduces each set of virtual coordinates
    to its largest pairwise distance. When the mask of a frame does not hold
    exactly the tissue of the reference mask, that distance changes by only a
    fraction of a period between neighbouring candidates, and a small change
    of the images or of the mask can move a whole frame by one period. This
    rule uses every coordinate instead. It starts from the rounded difference
    between the mean coordinates of the reference and of the frame, tries the
    integer offsets within ``radius`` of it in every direction, and keeps the
    one that minimises the mean distance from each coordinate of the frame to
    the nearest reference coordinate. The distance runs from the frame to the
    reference only, since the mask of a frame may hold less tissue than the
    reference mask.

    Parameters
    ----------
    phases
        Sequence of ``ndim`` arrays, each ``(..., frames)`` of unwrapped phase
        with NaN outside the mask.
    reference
        Index of the reference frame.
    radius
        Half width of the search around the starting estimate, in periods.
    downsample
        Stride applied before searching, as in :func:`temporal_phase_consistency`.
    return_info
        Also return a dict with ``delta`` (the starting estimate per frame),
        ``offset`` (how many steps the chosen jump sits away from it) and
        ``njumps`` (the radius), the keys :func:`temporal_phase_consistency`
        returns.

    Returns
    -------
    ``(ndim, frames)`` array of additive phase offsets, each a multiple of
    2*pi, and optionally the search info.
    """
    phases = [np.asarray(p, dtype=float) for p in phases]
    ndim = len(phases)
    spatial = phases[0].ndim - 1
    nfr = phases[0].shape[-1]

    if downsample and downsample > 1:
        phases = [_downsample(p, downsample, spatial) for p in phases]

    valid = np.ones(phases[0].shape, dtype=bool)
    for p in phases:
        valid &= ~np.isnan(p)

    def points(fr):
        """``(N, ndim)`` coordinates of one frame."""
        sel = valid[..., fr]
        return np.stack([p[..., fr][sel] for p in phases], axis=1)

    ref_pts = points(reference)
    if len(ref_pts) == 0:
        raise ValueError("the reference frame has no valid pixels")
    tree = cKDTree(ref_pts)
    ref_mean = ref_pts.mean(axis=0)
    steps = [np.array(s) for s in product(range(-radius, radius + 1), repeat=ndim)]

    jumps = np.zeros((ndim, nfr))
    start = np.zeros((ndim, nfr))
    offset = np.zeros((ndim, nfr))
    for fr in range(nfr):
        if fr == reference:
            continue
        pts = points(fr)
        if len(pts) == 0:
            continue
        guess = np.round((ref_mean - pts.mean(axis=0)) / (2 * np.pi))
        best, best_c = np.inf, guess
        for step in steps:
            c = guess + step
            d = tree.query(pts + 2 * np.pi * c)[0].mean()
            if d < best:
                best, best_c = d, c
        jumps[:, fr] = 2 * np.pi * best_c
        start[:, fr] = guess
        offset[:, fr] = best_c - guess

    if return_info:
        return jumps, {"delta": start, "offset": offset, "njumps": int(radius)}
    return jumps
