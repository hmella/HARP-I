"""Locating the first harmonic peak in k-space."""

from __future__ import annotations

import numpy as np

__all__ = ["harmonic_center"]


def harmonic_center(kspace, omega, window=(2, 2), fov=None, center=None):
    """Refine the predicted harmonic peak by a local search.

    Parameters
    ----------
    kspace
        Centred k-space of one encoding, ``(rows, cols)``.
    omega
        Predicted angular frequency ``(wx, wy)`` in rad/pixel.
    window
        Half-size ``(col, row)`` of the search window.
    fov
        ``(rows, cols)`` used to convert omega into a bin offset.
    center
        Skip the prediction and search around this ``(col, row)``.

    Returns
    -------
    ``(row, col)`` of the peak.
    """
    k = np.asarray(kspace)
    if fov is None:
        fov = np.asarray(k.shape, dtype=float)
    fov = np.asarray(fov, dtype=float)
    window = np.asarray(window, dtype=int)

    if center is None:
        s = np.round(0.5 * np.asarray(k.shape, dtype=float))
        c = np.round(s + np.asarray(omega, dtype=float) * fov / (2 * np.pi)).astype(int)
    else:
        c = np.asarray(center, dtype=int)

    # c is (col, row) here.
    row0 = c[1] - window[1]
    row1 = c[1] + window[1]
    col0 = c[0] - window[0]
    col1 = c[0] + window[0]

    row0 = max(row0, 0)
    col0 = max(col0, 0)
    row1 = min(row1, k.shape[0] - 1)
    col1 = min(col1, k.shape[1] - 1)
    if row1 < row0 or col1 < col0:
        raise ValueError("harmonic search window fell outside k-space")

    win = np.abs(k[row0 : row1 + 1, col0 : col1 + 1])
    ii, jj = np.unravel_index(int(np.argmax(win)), win.shape)
    return np.array([row0 + ii, col0 + jj], dtype=int)
