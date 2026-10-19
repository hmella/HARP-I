# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

"""Centred FFT helpers.

Equivalent of itok.m and ktoi.m. Both apply a centred (fftshift-ed)
transform along the selected axes only.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np

__all__ = ["itok", "ktoi"]


def _as_axes(axes: int | Iterable[int] | None, ndim: int) -> tuple[int, ...]:
    if axes is None:
        return tuple(range(ndim))
    if isinstance(axes, (int, np.integer)):
        return (int(axes),)
    return tuple(int(a) for a in axes)


def itok(data: np.ndarray, axes: int | Sequence[int] | None = None) -> np.ndarray:
    """Image space -> centred k-space.

    Parameters
    ----------
    data
        Input array.
    axes
        Axes to transform. ``None`` transforms every axis (matching MATLAB's
        ``itok(Data)`` with no ``Dims``).
    """
    ax = _as_axes(axes, data.ndim)
    out = np.asarray(data, dtype=complex)
    for a in ax:
        out = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(out, axes=a), axis=a), axes=a)
    return out


def ktoi(data: np.ndarray, axes: int | Sequence[int] | None = None) -> np.ndarray:
    """Centred k-space -> image space (inverse of :func:`itok`)."""
    ax = _as_axes(axes, data.ndim)
    out = np.asarray(data, dtype=complex)
    for a in ax:
        out = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(out, axes=a), axis=a), axes=a)
    return out
