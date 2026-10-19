# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

"""Quality-guided path-following phase unwrapping.

Derived from unwrap2.m and unwrap3.m in DENSEanalysis (Copyright (c) 2016
DENSEanalysis Contributors, MPL-2.0, written by Drew Gilliam), following
Ghiglia and Pritt 1998.

Pixels are unwrapped in descending order of phase quality, each against its
best already unwrapped neighbour. Two details keep the result in step with
the MATLAB version:

* Ties on quality are resolved by ascending column-major linear index. The
  sorted order is precomputed once and used as the heap key.
* Where several unwrapped neighbours tie on ``distance_weight * quality``,
  the first in neighbour order wins. That order comes from
  :func:`harpi.unwrap.quality.neighbour_offsets`.

Set ``HARPI_NO_NUMBA=1`` to use the plain NumPy kernel instead of the
compiled one.
"""

from __future__ import annotations

import heapq
import os
import warnings

import numpy as np

from .quality import neighbour_offsets, phase_quality_2d, phase_quality_3d

__all__ = ["unwrap2", "unwrap3", "HAVE_NUMBA"]

_NO_NUMBA = os.environ.get("HARPI_NO_NUMBA", "") not in ("", "0", "false", "False")

try:  # pragma: no cover - import guard
    if _NO_NUMBA:
        raise ImportError
    from numba import njit

    HAVE_NUMBA = True
except ImportError:  # pragma: no cover
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        def wrap(f):
            return f

        return wrap(args[0]) if args and callable(args[0]) else wrap


def _wrap(x):
    """Wrap to (-pi, pi], keeping the sign of the variation at +/-pi."""
    dps = x - np.floor((x + np.pi) / (2 * np.pi)) * (2 * np.pi)
    if np.isscalar(dps):
        if dps == -np.pi and x > 0:
            return np.pi
        return dps
    dps = np.where((dps == -np.pi) & (x > 0), np.pi, dps)
    return dps


@njit(cache=True, nogil=True)
def _flood_numba(wrapped, unwrapped, isunwrap, mask, qual, rank, order,
                 offs, qualfac, shape):  # pragma: no cover - compiled
    ndim = offs.shape[1]
    n = wrapped.size

    heap = np.empty(n + 1, dtype=np.int64)
    hsize = 0
    inheap = np.zeros(n, dtype=np.bool_)

    # strides for column-major (Fortran) linear indices
    strides = np.ones(ndim, dtype=np.int64)
    for d in range(1, ndim):
        strides[d] = strides[d - 1] * shape[d - 1]

    coord = np.empty(ndim, dtype=np.int64)

    def _push(h, size, val):
        h[size] = val
        i = size
        while i > 0:
            p = (i - 1) // 2
            if h[p] > h[i]:
                tmp = h[p]; h[p] = h[i]; h[i] = tmp
                i = p
            else:
                break
        return size + 1

    # seed the heap with masked, not-yet-unwrapped neighbours of unwrapped pixels
    for idx in range(n):
        if not isunwrap[idx]:
            continue
        rem = idx
        for d in range(ndim - 1, -1, -1):
            coord[d] = rem // strides[d]
            rem = rem % strides[d]
        for k in range(offs.shape[0]):
            ok = True
            nb = 0
            for d in range(ndim):
                c = coord[d] + offs[k, d]
                if c < 0 or c >= shape[d]:
                    ok = False
                    break
                nb += c * strides[d]
            if ok and mask[nb] and not isunwrap[nb] and not inheap[nb]:
                hsize = _push(heap, hsize, rank[nb])
                inheap[nb] = True

    while hsize > 0:
        # pop smallest rank == highest quality
        top = heap[0]
        hsize -= 1
        heap[0] = heap[hsize]
        i = 0
        while True:
            l = 2 * i + 1
            r = l + 1
            s = i
            if l < hsize and heap[l] < heap[s]:
                s = l
            if r < hsize and heap[r] < heap[s]:
                s = r
            if s == i:
                break
            tmp = heap[s]; heap[s] = heap[i]; heap[i] = tmp
            i = s

        idx = order[top]
        if isunwrap[idx] or not mask[idx]:
            continue

        rem = idx
        for d in range(ndim - 1, -1, -1):
            coord[d] = rem // strides[d]
            rem = rem % strides[d]

        best_f = -1.0e308
        best_nb = -1
        for k in range(offs.shape[0]):
            ok = True
            nb = 0
            for d in range(ndim):
                c = coord[d] + offs[k, d]
                if c < 0 or c >= shape[d]:
                    ok = False
                    break
                nb += c * strides[d]
            if not ok or not mask[nb] or not isunwrap[nb]:
                continue
            f = qualfac[k] * qual[nb]
            if f > best_f:
                best_f = f
                best_nb = nb

        if best_nb < 0:
            continue

        ref = unwrapped[best_nb]
        dp = wrapped[idx] - ref
        dps = dp - np.floor((dp + np.pi) / (2 * np.pi)) * (2 * np.pi)
        if dps == -np.pi and dp > 0:
            dps = np.pi
        unwrapped[idx] = ref + dps
        isunwrap[idx] = True

        for k in range(offs.shape[0]):
            ok = True
            nb = 0
            for d in range(ndim):
                c = coord[d] + offs[k, d]
                if c < 0 or c >= shape[d]:
                    ok = False
                    break
                nb += c * strides[d]
            if ok and mask[nb] and not isunwrap[nb] and not inheap[nb]:
                hsize = _push(heap, hsize, rank[nb])
                inheap[nb] = True

    return unwrapped, isunwrap


def _flood_python(wrapped, unwrapped, isunwrap, mask, qual, rank, order,
                  offs, qualfac, shape):
    """Plain NumPy version of the kernel above. Same result, much slower."""
    ndim = offs.shape[1]
    n = wrapped.size
    strides = np.ones(ndim, dtype=np.int64)
    for d in range(1, ndim):
        strides[d] = strides[d - 1] * shape[d - 1]

    def coords(idx):
        out = np.empty(ndim, dtype=np.int64)
        rem = idx
        for d in range(ndim - 1, -1, -1):
            out[d] = rem // strides[d]
            rem = rem % strides[d]
        return out

    def neighbours(idx):
        c = coords(idx)
        for k in range(offs.shape[0]):
            cc = c + offs[k]
            if np.any(cc < 0) or np.any(cc >= shape):
                continue
            yield k, int((cc * strides).sum())

    heap: list[int] = []
    inheap = np.zeros(n, dtype=bool)
    for idx in np.flatnonzero(isunwrap):
        for _, nb in neighbours(int(idx)):
            if mask[nb] and not isunwrap[nb] and not inheap[nb]:
                heapq.heappush(heap, int(rank[nb]))
                inheap[nb] = True

    while heap:
        idx = int(order[heapq.heappop(heap)])
        if isunwrap[idx] or not mask[idx]:
            continue
        best_f, best_nb = -np.inf, -1
        for k, nb in neighbours(idx):
            if mask[nb] and isunwrap[nb]:
                f = qualfac[k] * qual[nb]
                if f > best_f:
                    best_f, best_nb = f, nb
        if best_nb < 0:
            continue
        ref = unwrapped[best_nb]
        unwrapped[idx] = ref + _wrap(wrapped[idx] - ref)
        isunwrap[idx] = True
        for _, nb in neighbours(idx):
            if mask[nb] and not isunwrap[nb] and not inheap[nb]:
                heapq.heappush(heap, int(rank[nb]))
                inheap[nb] = True

    return unwrapped, isunwrap


def _unwrap_nd(phase, mask, quality, seed, unwrapped_in, offsets, qualfac, warn):
    shape = np.asarray(phase.shape, dtype=np.int64)

    # Column-major throughout, so linear indices match MATLAB.
    flat = lambda a: np.asarray(a).ravel(order="F")  # noqa: E731

    wrapped_f = flat(phase).astype(float)
    mask_f = flat(mask).astype(bool)
    qual_f = np.nan_to_num(flat(quality), nan=-np.inf)

    unwrapped_f = np.full(wrapped_f.size, np.nan)
    isunwrap_f = np.zeros(wrapped_f.size, dtype=bool)

    if unwrapped_in is not None:
        seeded = flat(unwrapped_in).astype(float)
        known = np.isfinite(seeded) & mask_f
        unwrapped_f[known] = seeded[known]
        isunwrap_f[known] = True
        if not known.any():
            raise ValueError("UnwrappedImage contained no finite values inside the mask")
    else:
        if isinstance(seed, str):
            if seed != "auto":
                raise ValueError(f"unknown seed mode {seed!r}")
            masked_q = np.where(mask_f, qual_f, -np.inf)
            if not np.isfinite(masked_q).any():
                raise ValueError("no valid pixel to seed from")
            s = int(np.argmax(masked_q))
        else:
            s = int(np.ravel_multi_index(tuple(np.asarray(seed)), phase.shape, order="F"))
            if not mask_f[s]:
                raise ValueError("seed point lies outside the mask")
        unwrapped_f[s] = wrapped_f[s]
        isunwrap_f[s] = True

    # Stable descending sort, so ties fall to the lowest linear index.
    order_idx = np.argsort(-qual_f, kind="stable").astype(np.int64)
    rank = np.empty_like(order_idx)
    rank[order_idx] = np.arange(order_idx.size, dtype=np.int64)

    kernel = _flood_numba if HAVE_NUMBA else _flood_python
    unwrapped_f, isunwrap_f = kernel(
        wrapped_f, unwrapped_f, isunwrap_f, mask_f, qual_f,
        rank, order_idx, offsets.astype(np.int64), qualfac.astype(float), shape,
    )

    complete = bool(np.array_equal(isunwrap_f, mask_f))
    if warn and not complete:
        missing = int(np.sum(mask_f & ~isunwrap_f))
        warnings.warn(
            f"{missing} masked pixel(s) could not be unwrapped; the mask is "
            "probably disconnected.",
            RuntimeWarning,
            stacklevel=3,
        )

    out = np.full(wrapped_f.size, np.nan)
    out[isunwrap_f] = unwrapped_f[isunwrap_f]
    return out.reshape(phase.shape, order="F"), complete


def unwrap2(phase, mask=None, pixel_size=(1.0, 1.0), seed="auto", connectivity=4,
            quality=None, unwrapped_image=None, warn=True):
    """Unwrap a 2-D phase image.

    Parameters
    ----------
    phase
        Wrapped phase, ``(rows, cols)``.
    mask
        Pixels to unwrap; everything else stays NaN.
    pixel_size
        ``(row, col)`` spacing, used to weight neighbours.
    seed
        ``'auto'`` (highest-quality masked pixel) or an ``(row, col)`` index.
    connectivity
        4 or 8.
    quality
        Precomputed quality map; computed if omitted.
    unwrapped_image
        Partially unwrapped image whose finite values seed the flood fill.
        Used to unwrap a full resolution frame from a downsampled solution.
    warn
        Warn when the mask is not fully unwrapped.

    Returns
    -------
    (unwrapped, complete)
        ``unwrapped`` is NaN outside the unwrapped region; ``complete`` says
        whether every masked pixel was reached.
    """
    phase = np.asarray(phase, dtype=float)
    if phase.ndim != 2:
        raise ValueError("unwrap2 expects a 2-D array")
    if mask is None:
        mask = np.ones(phase.shape, dtype=bool)
    mask = np.asarray(mask, dtype=bool)

    if quality is None:
        quality = phase_quality_2d(phase, mask, pixel_size, connectivity)

    offs = neighbour_offsets(2, connectivity)
    # The row offset is scaled by the column spacing, as in the MATLAB code.
    px = np.asarray(pixel_size, dtype=float)
    qualfac = 1.0 / np.sqrt((offs[:, 0] * px[1]) ** 2 + (offs[:, 1] * px[0]) ** 2)

    return _unwrap_nd(phase, mask, quality, seed, unwrapped_image, offs, qualfac, warn)


def unwrap3(phase, mask=None, voxel_size=(1.0, 1.0, 1.0), seed="auto", connectivity=6,
            quality=None, unwrapped_volume=None, warn=True):
    """Unwrap a 3-D phase volume. See :func:`unwrap2`."""
    phase = np.asarray(phase, dtype=float)
    if phase.ndim != 3:
        raise ValueError("unwrap3 expects a 3-D array")
    if mask is None:
        mask = np.ones(phase.shape, dtype=bool)
    mask = np.asarray(mask, dtype=bool)

    if quality is None:
        quality = phase_quality_3d(phase, mask, voxel_size, connectivity)

    offs = neighbour_offsets(3, connectivity)
    vx = np.asarray(voxel_size, dtype=float)
    qualfac = 1.0 / np.sqrt(
        (offs[:, 0] * vx[1]) ** 2 + (offs[:, 1] * vx[0]) ** 2 + (offs[:, 2] * vx[2]) ** 2
    )

    return _unwrap_nd(phase, mask, quality, seed, unwrapped_volume, offs, qualfac, warn)
