"""Phase quality maps for quality-guided unwrapping.

Derived from phasequality2_mask.m and phasequality3.m in DENSEanalysis
(Copyright (c) 2016 DENSEanalysis Contributors, MPL-2.0, written by Drew
Gilliam). Only the masked variant is implemented, which is the one HARP-I
uses.

Quality is ``exp(-sqrt(Gvar))``, where ``Gvar`` is the distance weighted
variance of the locally unwrapped phase gradients to each neighbour
(Ghiglia and Pritt 1998, Spottiswoode et al. 2007).
"""

from __future__ import annotations

import numpy as np

__all__ = ["phase_quality_2d", "phase_quality_3d", "neighbour_offsets"]


def neighbour_offsets(ndim: int, connectivity: int) -> np.ndarray:
    """Neighbour offsets for a connectivity, as ``(K, ndim)``.

    2-D accepts 4 or 8, 3-D accepts 6 or 26.
    """
    if ndim == 2:
        if connectivity == 4:
            return np.array([[-1, 0], [0, -1], [0, 1], [1, 0]], dtype=np.intp)
        if connectivity == 8:
            return np.array(
                [[i, j] for i in (-1, 0, 1) for j in (-1, 0, 1) if (i, j) != (0, 0)],
                dtype=np.intp,
            )
        raise ValueError("2-D connectivity must be 4 or 8")
    if ndim == 3:
        if connectivity == 6:
            return np.array(
                [[-1, 0, 0], [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1]],
                dtype=np.intp,
            )
        if connectivity == 26:
            return np.array(
                [
                    [i, j, k]
                    for i in (-1, 0, 1)
                    for j in (-1, 0, 1)
                    for k in (-1, 0, 1)
                    if (i, j, k) != (0, 0, 0)
                ],
                dtype=np.intp,
            )
        raise ValueError("3-D connectivity must be 6 or 26")
    raise ValueError("only 2-D and 3-D are supported")


def _wrap(x: np.ndarray) -> np.ndarray:
    """Wrap to (-pi, pi]."""
    return np.mod(x + np.pi, 2 * np.pi) - np.pi


def _quality(
    phase: np.ndarray,
    mask: np.ndarray,
    offsets: np.ndarray,
    spacing: np.ndarray,
) -> np.ndarray:
    """Shared implementation for 2-D and 3-D."""
    phase = np.asarray(phase, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    shape = phase.shape

    # Distance weight per neighbour, so closer neighbours count for more.
    # The row offset is scaled by the column spacing, as in the MATLAB code.
    scale = np.asarray(spacing, dtype=float)[::-1]
    dist = np.sqrt(((offsets * scale) ** 2).sum(axis=1))
    weight = 1.0 / dist

    idx = np.argwhere(mask)                       # (Npts, ndim)
    npts = idx.shape[0]
    qual = np.full(shape, np.nan)
    if npts == 0:
        return qual

    here = phase[tuple(idx.T)]
    grads = np.zeros((npts, len(offsets)))
    wts = np.zeros((npts, len(offsets)))

    for n, off in enumerate(offsets):
        nb = idx - off
        inside = np.ones(npts, dtype=bool)
        for d in range(idx.shape[1]):
            inside &= (nb[:, d] >= 0) & (nb[:, d] < shape[d])
        sel = np.zeros(npts, dtype=bool)
        sel[inside] = mask[tuple(nb[inside].T)]
        if not sel.any():
            continue
        grads[sel, n] = phase[tuple(nb[sel].T)] - here[sel]
        wts[sel, n] = weight[n]

    grads = _wrap(grads)

    total = wts.sum(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        wnorm = np.where(total > 0, wts / total, 0.0)

    gmean = (wnorm * grads).sum(axis=1, keepdims=True)
    gvar = (wnorm * (grads - gmean) ** 2).sum(axis=1)

    qual[tuple(idx.T)] = np.exp(-np.sqrt(gvar))
    return qual


def phase_quality_2d(
    phase: np.ndarray,
    mask: np.ndarray | None = None,
    pixel_size: tuple[float, float] = (1.0, 1.0),
    connectivity: int = 4,
) -> np.ndarray:
    """Phase quality in [0, 1]; NaN outside the mask."""
    if phase.ndim != 2:
        raise ValueError("phase_quality_2d expects a 2-D array")
    if mask is None:
        mask = np.ones(phase.shape, dtype=bool)
    return _quality(phase, mask, neighbour_offsets(2, connectivity), np.asarray(pixel_size))


def phase_quality_3d(
    phase: np.ndarray,
    mask: np.ndarray | None = None,
    voxel_size: tuple[float, float, float] = (1.0, 1.0, 1.0),
    connectivity: int = 6,
) -> np.ndarray:
    """Phase quality in [0, 1]; NaN outside the mask."""
    if phase.ndim != 3:
        raise ValueError("phase_quality_3d expects a 3-D array")
    if mask is None:
        mask = np.ones(phase.shape, dtype=bool)
    return _quality(phase, mask, neighbour_offsets(3, connectivity), np.asarray(voxel_size))
