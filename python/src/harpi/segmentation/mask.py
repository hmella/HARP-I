# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

"""Turning contours into binary masks.

A contour may be labelled with the wall it outlines. Contours are then
grouped per structure, and each structure contributes the region between
its outer and inner boundary. The mask is the union of those regions.

Without labels the first contour is taken as the outer boundary and the
rest as holes, which is the short axis case with a single wall.

Points lying exactly on a contour count as wall, but only where they are
inside an outer boundary.
"""

from __future__ import annotations

import numpy as np
from matplotlib.path import Path

from .spline import resample_contour

__all__ = ["contours_to_mask", "structure_masks", "Contour", "WALLS",
           "wall_role", "wall_structure"]


# Walls that can be outlined, as (key, label shown to the user).
WALLS = (
    ("lv_epi", "LV epicardium"),
    ("lv_endo", "LV endocardium"),
    ("rv_epi", "RV epicardium"),
    ("rv_endo", "RV endocardium"),
    ("other", "Other"),
)

_WALL_LABELS = dict(WALLS)


def wall_structure(key: str) -> str:
    """Structure a wall belongs to, for example 'lv' for 'lv_endo'."""
    key = (key or "").lower()
    return key.split("_")[0] if "_" in key else key


def wall_role(key: str) -> str:
    """Whether a wall is the outer or the inner boundary of its structure."""
    key = (key or "").lower()
    if key.endswith("_endo"):
        return "inner"
    if key.endswith("_epi"):
        return "outer"
    return "outer"


def wall_label(key: str) -> str:
    """Name shown to the user for a wall key."""
    return _WALL_LABELS.get((key or "").lower(), key or "contour")


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
        Sequence of :class:`Contour` or ``(N, 2)`` arrays. Labelled contours
        are grouped per structure (see :data:`WALLS`); unlabelled ones take
        the first as the outer boundary and the rest as holes.
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
    labels = []
    for c in contours:
        if isinstance(c, Contour):
            polys.append(c.sampled(resolution))
            labels.append(c.label)
        else:
            polys.append(resample_contour(np.asarray(c, dtype=float), True, None,
                                          None, resolution))
            labels.append("")

    groups = _group_by_structure(polys, labels)

    mask = np.zeros(pts.shape[0], dtype=bool)
    for outer_polys, inner_polys in groups:
        mask |= _region(pts, outer_polys, inner_polys)
    return mask.reshape(rows, cols)


def structure_masks(contours, shape, resolution=0.5, origin=1):
    """One mask per labelled structure, rather than their union.

    ``contours_to_mask`` merges every structure, which is what motion estimation
    wants: it should see all the tissue that was outlined. Regional analysis
    wants them apart, because segment wedges spanning both ventricles are
    meaningless. Keys are ``lv``, ``rv`` or ``other``; unlabelled input yields a
    single ``other`` entry, matching the union.

    Returns
    -------
    dict of structure key to boolean ``(rows, cols)``.
    """
    if len(contours) == 0:
        raise ValueError("at least one contour is required")

    rows, cols = int(shape[0]), int(shape[1])
    X, Y = np.meshgrid(np.arange(cols) + origin, np.arange(rows) + origin)
    pts = np.column_stack([X.ravel(), Y.ravel()])

    polys, labels = [], []
    for c in contours:
        if isinstance(c, Contour):
            polys.append(c.sampled(resolution))
            labels.append(c.label)
        else:
            polys.append(resample_contour(np.asarray(c, dtype=float), True, None,
                                          None, resolution))
            labels.append("")

    keys = _structure_keys(labels)
    out = {}
    for key, (outer_polys, inner_polys) in zip(keys, _group_by_structure(polys, labels)):
        region = _region(pts, outer_polys, inner_polys).reshape(rows, cols)
        out[key] = out[key] | region if key in out else region
    return out


def _structure_keys(labels):
    """Structure key per group, in the order _group_by_structure returns them."""
    known = [lb for lb in labels if wall_structure(lb) in {"lv", "rv"}]
    if not known:
        return ["other"]
    order = []
    for lb in labels:
        key = wall_structure(lb) if wall_structure(lb) in {"lv", "rv"} else "other"
        if key not in order:
            order.append(key)
    return order


def _group_by_structure(polys, labels):
    """Split sampled contours into (outer, inner) pairs per structure.

    Unlabelled input keeps the original convention: the first contour is the
    outer boundary and the rest are holes.
    """
    known = [lb for lb in labels if wall_structure(lb) in {"lv", "rv"}]
    if not known:
        return [(polys[:1], polys[1:])]

    order, grouped = [], {}
    for poly, lb in zip(polys, labels):
        key = wall_structure(lb) if wall_structure(lb) in {"lv", "rv"} else "other"
        if key not in grouped:
            grouped[key] = ([], [])
            order.append(key)
        outer, inner = grouped[key]
        (inner if wall_role(lb) == "inner" else outer).append(poly)

    out = []
    for key in order:
        outer, inner = grouped[key]
        if not outer:
            # Only an inner boundary was drawn, so treat it as a solid region
            # rather than dropping the structure.
            outer, inner = inner, []
        out.append((outer, inner))
    return out


def _region(pts, outer_polys, inner_polys):
    """Points inside any outer boundary and outside every inner one."""
    if not outer_polys:
        return np.zeros(pts.shape[0], dtype=bool)

    # A small positive radius includes points lying on the boundary.
    inside = np.zeros(pts.shape[0], dtype=bool)
    on = np.zeros(pts.shape[0], dtype=bool)
    for poly in outer_polys:
        p = Path(poly)
        loose = p.contains_points(pts, radius=1e-9)
        inside |= loose
        on |= loose & ~p.contains_points(pts, radius=-1e-9)

    region = inside.copy()
    for poly in inner_polys:
        p = Path(poly)
        loose = p.contains_points(pts, radius=1e-9)
        on |= loose & ~p.contains_points(pts, radius=-1e-9)
        region &= ~loose

    return region | (on & inside)
