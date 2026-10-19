"""Temporal phase consistency.

Each frame is unwrapped on its own, so its phase carries an unknown additive
multiple of 2*pi (the constant ``c_k`` of eq. 10 to 12). These are recovered
by matching the spread of the virtual coordinates against the reference
frame.

The search runs on convex hulls. The criterion only needs the largest
distance between two point sets, which is always reached at hull vertices,
and a candidate jump just translates one of the sets.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import ConvexHull

__all__ = ["temporal_phase_consistency", "max_pair_distance", "extreme_points"]


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


def temporal_phase_consistency(phases, reference=0, njumps=8, downsample=2):
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

    Returns
    -------
    ``(ndim, frames)`` array of additive phase offsets, each a multiple of 2*pi.
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

    return jumps
