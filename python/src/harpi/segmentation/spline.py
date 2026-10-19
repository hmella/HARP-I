"""Contour resampling.

Derived from clinesegments.m in DENSEanalysis (Copyright (c) 2016
DENSEanalysis Contributors, MPL-2.0).

A contour is a list of control points, each flagged curved or corner. Curved
runs are interpolated with a cubic spline parameterised by cumulative chord
length. Corner points break the spline, and straight segments are sampled
linearly.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

__all__ = ["resample_contour"]


def _chord_length(pts: np.ndarray) -> np.ndarray:
    d = np.sqrt((np.diff(pts, axis=0) ** 2).sum(axis=1))
    return np.concatenate([[0.0], np.cumsum(d)])


def _dedupe(pts: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    """Drop consecutive duplicate points.

    Chord length parameterisation needs strictly increasing knots, and
    duplicates turn up easily, for instance from a double click.
    """
    if len(pts) < 2:
        return pts
    keep = np.ones(len(pts), dtype=bool)
    keep[1:] = (np.abs(np.diff(pts, axis=0)) > tol).any(axis=1)
    return pts[keep]


def resample_contour(points, closed=True, curved=None, corner=None, resolution=0.5):
    """Densely sample a contour.

    Parameters
    ----------
    points
        ``(N, 2)`` control points as ``(x, y)``.
    closed
        Whether the contour closes back on the first point.
    curved
        Per-point flag; ``True`` (default) interpolates, ``False`` draws straight.
    corner
        Per-point flag; a corner breaks the spline there.
    resolution
        Target spacing between output samples, in pixels.

    Returns
    -------
    ``(M, 2)`` sampled points.
    """
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError("points must be (N, 2)")
    if closed and len(pts) > 1 and np.allclose(pts[0], pts[-1]):
        pts = pts[:-1]          # the closing point is implicit
    n = pts.shape[0]
    if n == 0:
        return pts.copy()
    if n == 1:
        return pts.copy()

    curved = np.ones(n, dtype=bool) if curved is None else np.asarray(curved, dtype=bool).ravel()
    corner = np.zeros(n, dtype=bool) if corner is None else np.asarray(corner, dtype=bool).ravel()
    if curved.size != n:
        curved = np.resize(curved, n)
    if corner.size != n:
        corner = np.resize(corner, n)

    if n == 2 or not curved.any():
        return _sample_polyline(pts, closed, resolution)

    # Corner points split the contour into independently splined runs.
    breaks = np.flatnonzero(corner | ~curved)
    if breaks.size == 0:
        return _spline_run(pts, closed, resolution)

    if not closed:
        segments = np.split(np.arange(n), breaks[breaks > 0])
    else:
        roll = breaks[0]
        order = np.roll(np.arange(n), -roll)
        cuts = np.flatnonzero((corner | ~curved)[order])
        segments = [order[a:b + 1] for a, b in zip(cuts, list(cuts[1:]) + [n - 1])]
        if segments and segments[-1][-1] != order[-1]:
            segments.append(np.concatenate([order[cuts[-1]:], order[:1]]))

    out = []
    for seg in segments:
        if len(seg) < 2:
            continue
        sub = pts[seg]
        piece = _spline_run(sub, False, resolution) if len(seg) > 2 else \
            _sample_polyline(sub, False, resolution)
        out.append(piece if not out else piece[1:])
    if not out:
        return _sample_polyline(pts, closed, resolution)
    result = np.vstack(out)
    if closed:
        result = np.vstack([result, result[:1]])
    return result


def _sample_polyline(pts, closed, resolution):
    loop = np.vstack([pts, pts[:1]]) if closed else pts
    out = []
    for k in range(len(loop) - 1):
        a, b = loop[k], loop[k + 1]
        d = np.hypot(*(b - a))
        m = max(int(np.ceil(d / resolution)), 1)
        t = np.linspace(0, 1, m + 1)[:-1]
        out.append(a + t[:, None] * (b - a))
    out.append(loop[-1:])
    return np.vstack(out)


def _spline_run(pts, closed, resolution):
    pts = _dedupe(np.asarray(pts, dtype=float))
    if len(pts) < 3:
        return _sample_polyline(pts, closed and len(pts) > 2, resolution)

    if closed:
        loop = np.vstack([pts, pts[:1]])
        bc = "periodic"
    else:
        loop = pts
        bc = "natural"

    s = _chord_length(loop)
    if s[-1] <= 0 or not np.all(np.diff(s) > 0):
        return _sample_polyline(pts, closed, resolution)

    try:
        cs = CubicSpline(s, loop, bc_type=bc, axis=0)
    except ValueError:
        return _sample_polyline(pts, closed, resolution)

    m = max(int(np.ceil(s[-1] / resolution)), 1)
    return cs(np.linspace(0, s[-1], m + 1))
