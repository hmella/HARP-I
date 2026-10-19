# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

"""Pixelwise Lagrangian strain.

For each pixel the deformation gradient is estimated by least squares from
the displacement differences to its neighbours inside the mask, then

    E = 0.5 * (F' F - I)

is rotated into the radial and circumferential frame, plus longitudinal in
3-D. The derivatives are not smoothed.

Pixels are grouped by neighbour pattern, of which there are at most 2^4 in
2-D and 2^6 in 3-D. Each group shares one design matrix and is solved in one
go. Rank deficient patterns, such as two collinear neighbours, fall back to
a least norm solve.
"""

from __future__ import annotations

import numpy as np

from .analysis import segment_geometry, segment_index
from .types import SegmentStrain, StrainResult

__all__ = ["pixel_strain_2d", "pixel_strain_3d", "strain_by_segments"]

# Neighbour order: -i, +i, -j, +j.
_OFF2 = ((-1, 0), (1, 0), (0, -1), (0, 1))
_OFF3 = ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1))


def _shift(a, off, fill=False):
    """Shift an array by ``off``, filling vacated entries."""
    out = np.full_like(a, fill, dtype=a.dtype)
    src, dst = [], []
    for o, n in zip(off, a.shape):
        if o < 0:
            src.append(slice(-o, n)); dst.append(slice(0, n + o))
        elif o > 0:
            src.append(slice(0, n - o)); dst.append(slice(o, n))
        else:
            src.append(slice(0, n)); dst.append(slice(0, n))
    out[tuple(dst)] = a[tuple(src)]
    return out


def _valid_neighbours(mask, offsets):
    """Boolean ``(..., K)``: is neighbour k of this pixel inside the mask?"""
    return np.stack([_shift(mask, o, False) & mask for o in offsets], axis=-1)


def _solve_groups(dxn, dXn, avail, ndim):
    """Least-squares deformation gradient per pixel, grouped by neighbour pattern.

    ``dxn`` is ``(P, ndim, K)`` deformed differences, ``dXn`` the reference
    differences, ``avail`` a ``(P, K)`` mask of usable neighbours.
    """
    npix, _, nnb = dxn.shape
    F = np.full((npix, ndim, ndim), np.nan)

    codes = (avail * (1 << np.arange(nnb))).sum(axis=1)
    for code in np.unique(codes):
        if code == 0:
            continue
        sel = codes == code
        cols = [k for k in range(nnb) if code & (1 << k)]
        if len(cols) < 1:
            continue
        # The grid is regular, so the reference offsets are the same for
        # every pixel in the group.
        A = dXn[sel][0][:, cols]
        B = dxn[sel][:, :, cols]                 # (n_group, ndim, n)
        # F @ A = B, solved as F = B @ pinv(A).
        F[sel] = B @ np.linalg.pinv(A)
    return F


def pixel_strain_2d(u, mask, origin=None, orientation=None, times=None):
    """Lagrangian strain from a 2-D displacement field.

    Parameters
    ----------
    u
        ``(rows, cols, 2, frames)`` displacement on the reference grid.
    mask
        ``(rows, cols)`` reference-frame mask.
    origin
        ``(x, y)`` centre for the radial/circumferential frame; defaults to
        the mask centroid.
    orientation
        Per-pixel angle, overriding ``origin``.

    Returns
    -------
    :class:`~harpi.types.StrainResult`
    """
    u = np.asarray(u, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    isz = u.shape[:2]
    nfr = u.shape[-1]
    if times is None:
        times = np.arange(nfr)

    Xg, Yg = np.meshgrid(np.arange(1, isz[1] + 1), np.arange(1, isz[0] + 1))

    # Drop pixels without an orthogonal neighbour pair.
    quad = [((-1, 0), (0, -1)), ((-1, 0), (0, 1)), ((1, 0), (0, -1)), ((1, 0), (0, 1))]
    keep = np.zeros(isz, dtype=bool)
    for a, b in quad:
        keep |= _shift(mask, a, False) & _shift(mask, b, False)
    mask = keep & mask

    nb = _valid_neighbours(mask, _OFF2)
    nneighbour = nb.sum(axis=-1)
    active = mask & (nneighbour > 1)

    if orientation is None:
        if origin is None:
            origin = (float(Xg[mask].mean()), float(Yg[mask].mean()))
        theta = np.arctan2(Yg - origin[1], Xg - origin[0])
    else:
        theta = np.asarray(orientation, dtype=float)
        origin = origin if origin is not None else (np.nan, np.nan)

    idx = np.argwhere(active)
    npix = idx.shape[0]
    out = {k: np.full(isz + (nfr,), np.nan)
           for k in ("XX", "YY", "XY", "RR", "CC", "RC", "p1", "p2", "p1or")}

    if npix == 0:
        return StrainResult(mask=active, **out, origin=np.asarray(origin))

    avail = nb[active]                                       # (P, 4)
    dXn = np.zeros((npix, 2, 4))
    for k, off in enumerate(_OFF2):
        dXn[:, 0, k] = off[1]        # x offset is the column shift
        dXn[:, 1, k] = off[0]

    ct = np.cos(theta[active])
    st = np.sin(theta[active])

    for fr in range(nfr):
        xt = Xg + u[:, :, 0, fr]
        yt = Yg + u[:, :, 1, fr]
        dxn = np.zeros((npix, 2, 4))
        for k, off in enumerate(_OFF2):
            dxn[:, 0, k] = (_shift(xt, off, np.nan) - xt)[active]
            dxn[:, 1, k] = (_shift(yt, off, np.nan) - yt)[active]
        dxn = np.nan_to_num(dxn)

        F = _solve_groups(dxn, dXn, avail, 2)
        E = 0.5 * (np.einsum("pki,pkj->pij", F, F) - np.eye(2))

        R = np.empty((npix, 2, 2))
        R[:, 0, 0] = ct; R[:, 0, 1] = st
        R[:, 1, 0] = -st; R[:, 1, 1] = ct
        Erot = R @ E @ np.transpose(R, (0, 2, 1))

        # Symmetric 2x2 eigenvalues in closed form, in ascending order
        a11, a12, a22 = E[:, 0, 0], E[:, 0, 1], E[:, 1, 1]
        tr = a11 + a22
        disc = np.sqrt(np.maximum(((a11 - a22) / 2) ** 2 + a12**2, 0.0))
        lo, hi = tr / 2 - disc, tr / 2 + disc
        # eigenvector for the larger eigenvalue
        vx, vy = np.where(np.abs(a12) > 0, hi - a22, 1.0), np.where(np.abs(a12) > 0, a12, 0.0)

        for name, v in (("XX", E[:, 0, 0]), ("XY", E[:, 0, 1]), ("YY", E[:, 1, 1]),
                        ("RR", Erot[:, 0, 0]), ("RC", Erot[:, 0, 1]),
                        ("CC", Erot[:, 1, 1]),
                        ("p1", hi), ("p2", lo),
                        ("p1or", np.arctan2(vy, vx))):
            out[name][active, fr] = v

    return StrainResult(mask=active, **out, origin=np.asarray(origin),
                        meta={"times": np.asarray(times)})


def pixel_strain_3d(u, mask, origin=None, times=None):
    """Lagrangian strain from a 3-D displacement field.

    ``u`` is ``(rows, cols, slices, 3, frames)``; ``mask`` is the reference
    volume. Returns radial, circumferential and longitudinal components.
    """
    u = np.asarray(u, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    isz = u.shape[:3]
    nfr = u.shape[-1]
    if times is None:
        times = np.arange(nfr)

    Xg, Yg, Zg = np.meshgrid(np.arange(1, isz[1] + 1), np.arange(1, isz[0] + 1),
                             np.arange(1, isz[2] + 1))

    nb = _valid_neighbours(mask, _OFF3)
    in_plane = nb[..., :4].sum(axis=-1)
    thru = nb[..., 4:].sum(axis=-1)
    active = mask & ((in_plane > 1) | (thru > 1))

    if origin is None:
        origin = (float(Xg[mask].mean()), float(Yg[mask].mean()))
    theta = np.arctan2(Yg - origin[1], Xg - origin[0])

    idx = np.argwhere(active)
    npix = idx.shape[0]
    names = ("XX", "YY", "ZZ", "RR", "CC", "LL", "RC")
    out = {k: np.full(isz + (nfr,), np.nan) for k in names}
    if npix == 0:
        return StrainResult(mask=active, **out, origin=np.asarray(origin))

    avail = nb[active]
    dXn = np.zeros((npix, 3, 6))
    for k, off in enumerate(_OFF3):
        dXn[:, 0, k] = off[1]
        dXn[:, 1, k] = off[0]
        dXn[:, 2, k] = off[2]

    ct, st = np.cos(theta[active]), np.sin(theta[active])
    R = np.zeros((npix, 3, 3))
    R[:, 0, 0] = ct; R[:, 0, 1] = st
    R[:, 1, 0] = -st; R[:, 1, 1] = ct
    R[:, 2, 2] = 1.0

    for fr in range(nfr):
        trj = [Xg + u[..., 0, fr], Yg + u[..., 1, fr], Zg + u[..., 2, fr]]
        dxn = np.zeros((npix, 3, 6))
        for k, off in enumerate(_OFF3):
            for d in range(3):
                dxn[:, d, k] = (_shift(trj[d], off, np.nan) - trj[d])[active]
        dxn = np.nan_to_num(dxn)

        F = _solve_groups(dxn, dXn, avail, 3)
        E = 0.5 * (np.einsum("pki,pkj->pij", F, F) - np.eye(3))
        Erot = R @ E @ np.transpose(R, (0, 2, 1))

        for name, v in (("XX", E[:, 0, 0]), ("YY", E[:, 1, 1]), ("ZZ", E[:, 2, 2]),
                        ("RR", Erot[:, 0, 0]), ("RC", Erot[:, 0, 1]),
                        ("CC", Erot[:, 1, 1]), ("LL", Erot[:, 2, 2])):
            out[name][active, fr] = v

    return StrainResult(mask=active, **out, origin=np.asarray(origin),
                        meta={"times": np.asarray(times)})


def strain_by_segments(strain, n_segments=6, clockwise=False, origin=None, *,
                       anterior=None, inferior=None, outlier_threshold=0.5,
                       position_a=None, position_b=None):
    """Average strain over equiangular segments.

    The angular wedge is intersected with the mask and ``theta`` is wrapped
    into the rotated frame, so every pixel falls in exactly one segment.

    ``anterior`` and ``inferior`` are the right-ventricular insertion points as
    ``(x, y)``. The anterior one is the boundary between segment 1 and segment
    2 and the two together fix the direction; see
    :func:`~harpi.analysis.segment_geometry`, which holds the layout.

    ``position_a`` and ``position_b`` are the MATLAB spelling of the anterior
    point, its x and its y, kept so existing callers keep working.

    The layout no longer reads 1:1 against ``getStrainBySegments.m``. That
    function numbers from half a turn away, which put every AHA name on the
    wrong wall once this package started attaching them. The MATLAB regression
    stops at ``pixelstrain``, so the numbering is not pinned by it.
    """
    if anterior is not None and (position_a is not None
                                 or position_b is not None):
        raise ValueError("pass anterior=(x, y) or the MATLAB "
                         "position_a/position_b, not both")
    mask = np.asarray(strain.mask, dtype=bool)
    isz = mask.shape
    if mask.ndim != 2:
        # The wedge is built from a 2-D meshgrid over the leading axes, so on a
        # volume the angular assignment is wrong for every slice but one. On a
        # cubic volume that used to pass silently and return plausible numbers.
        raise ValueError(
            f"strain_by_segments works on a 2-D slice, got a {mask.ndim}-D mask "
            f"of shape {isz}. Equiangular segments are a short-axis construct; "
            "pass one slice at a time."
        )
    Xg, Yg = np.meshgrid(np.arange(1, isz[1] + 1), np.arange(1, isz[0] + 1))

    if origin is None:
        origin = strain.origin if strain.origin is not None else None
    if origin is None or not np.all(np.isfinite(origin)):
        origin = (float(Xg[mask].mean()), float(Yg[mask].mean()))

    if anterior is None and (position_a is not None or position_b is not None):
        anterior = (float(Xg.max()) if position_a is None else position_a,
                    float(origin[1]) if position_b is None else position_b)

    geometry = segment_geometry(n_segments=n_segments, origin=origin,
                                anterior=anterior, inferior=inferior,
                                clockwise=clockwise)
    theta = np.arctan2(Yg - origin[1], Xg - origin[0])
    index = segment_index(theta, geometry)

    nfr = None
    for comp in ("CC", "RR"):
        v = getattr(strain, comp)
        if v is not None:
            nfr = v.shape[-1]
            break
    if nfr is None:
        raise ValueError("strain has neither CC nor RR")

    segments = np.full(isz, np.nan)
    curves = {"CC": None, "RR": None}
    for comp in ("CC", "RR"):
        if getattr(strain, comp) is not None:
            curves[comp] = np.full((n_segments, nfr), np.nan)

    for seg in range(n_segments):
        pos = index == seg + 1
        segments[pos] = seg + 1
        inseg = pos & mask
        for comp in ("CC", "RR"):
            v = getattr(strain, comp)
            if v is None:
                continue
            for fr in range(nfr):
                frame = v[..., fr]
                sel = inseg & ~(frame > outlier_threshold)
                vals = frame[sel]
                vals = vals[np.isfinite(vals)]
                if vals.size:
                    curves[comp][seg, fr] = vals.mean()

    segments[~mask] = np.nan
    return SegmentStrain(segments=segments, CC=curves["CC"], RR=curves["RR"],
                         n_segments=n_segments, geometry=geometry)
