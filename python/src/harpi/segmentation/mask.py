"""Turning contours into binary masks.

The mask is everything inside the outermost contour and outside every inner
one. Points lying exactly on a contour count as wall, but only where they
are inside the outer boundary.
"""

from __future__ import annotations

import numpy as np
from matplotlib.path import Path

from .spline import resample_contour

__all__ = ["contours_to_mask", "Contour"]


class Contour:
    """A single closed contour with per-point curve/corner flags."""

    __slots__ = ("points", "closed", "curved", "corner", "label")

    def __init__(self, points, closed=True, curved=None, corner=None, label=""):
        self.points = np.asarray(points, dtype=float)
        self.closed = bool(closed)
        n = len(self.points)
        self.curved = np.ones(n, bool) if curved is None else np.asarray(curved, bool)
        self.corner = np.zeros(n, bool) if corner is None else np.asarray(corner, bool)
        self.label = label

    def sampled(self, resolution=0.5):
        return resample_contour(self.points, self.closed, self.curved,
                                self.corner, resolution)

    def __len__(self):
        return len(self.points)


def contours_to_mask(contours, shape, resolution=0.5, origin=1):
    """Rasterise contours onto a pixel grid.

    Parameters
    ----------
    contours
        Sequence of :class:`Contour` or ``(N, 2)`` arrays. The first is the
        outer boundary; the rest are holes.
    shape
        ``(rows, cols)``.
    resolution
        Contour sampling step.
    origin
        Coordinate of the first pixel centre. Use ``1`` for the 1-based grid
        of the MATLAB toolbox, ``0`` for 0-based contours.

    Returns
    -------
    Boolean ``(rows, cols)``.
    """
    if len(contours) == 0:
        raise ValueError("at least one contour is required")

    rows, cols = int(shape[0]), int(shape[1])
    X, Y = np.meshgrid(np.arange(cols) + origin, np.arange(rows) + origin)
    pts = np.column_stack([X.ravel(), Y.ravel()])

    polys = []
    for c in contours:
        sampled = c.sampled(resolution) if isinstance(c, Contour) else \
            resample_contour(np.asarray(c, dtype=float), True, None, None, resolution)
        polys.append(sampled)

    # A small positive radius includes points lying on the boundary.
    outer = Path(polys[0])
    inside = outer.contains_points(pts, radius=1e-9)
    on = outer.contains_points(pts, radius=1e-9) & ~outer.contains_points(pts, radius=-1e-9)

    mask = inside.copy()
    for poly in polys[1:]:
        p = Path(poly)
        strictly_in = p.contains_points(pts, radius=-1e-9)
        on |= p.contains_points(pts, radius=1e-9) & ~strictly_in
        mask &= ~p.contains_points(pts, radius=1e-9)

    mask |= on & inside
    return mask.reshape(rows, cols)
