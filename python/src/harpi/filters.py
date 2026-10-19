"""Harmonic bandpass filters.

A Butterworth bandpass centred on the first harmonic, with cutoff and order
chosen to suppress the higher harmonics and DC. The raised cosine
``transmission`` window used by the 2-D demo is also available.
"""

from __future__ import annotations

import numpy as np

from .fft import itok, ktoi
from .harmonic import harmonic_center

__all__ = ["HarpFilter2D", "HarpFilter3D", "angular_frequencies", "butterworth_3d"]


def _freq_range(n: int) -> np.ndarray:
    """Normalised frequency axis.

    Odd lengths divide by ``n-1``, even lengths by ``n``, as in HARPFilter.m.
    """
    if n % 2:
        return np.arange(-(n - 1) / 2, (n - 1) / 2 + 1) / (n - 1)
    return np.arange(-n / 2, n / 2) / n


def angular_frequencies(shape, k1, k2):
    """Frequency component along each wave vector, ``(rows, cols, 2)``."""
    rows, cols = shape[0], shape[1]
    X, Y = np.meshgrid(_freq_range(cols), _freq_range(rows))
    return np.stack([k1[0] * X + k1[1] * Y, k2[0] * X + k2[1] * Y], axis=-1)


def _transmission(shape, wave1, wave2):
    """Raised-cosine bandpass in the log-radial coordinate."""
    W = angular_frequencies(shape, wave1, wave2)
    Wx, Wy = W[..., 0], W[..., 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        # log(0) at the DC bin gives -inf, which falls in the stop band.
        Z = np.abs(np.log(Wx + 1j * Wy))
    H = np.zeros_like(Z)
    passband = np.isfinite(Z) & (Z < 1)
    H[passband] = np.cos(np.pi * Z[passband] / 2) ** 2
    return H


def _butterworth_2d(shape, centre, cutoff, order):
    """Butterworth lowpass centred on the harmonic peak."""
    rows, cols = shape[0], shape[1]
    ox, oy = np.meshgrid(np.arange(1, cols + 1), np.arange(1, rows + 1))
    # centre arrives as (row, col)
    r = np.sqrt((ox - centre[1]) ** 2 + (oy - centre[0]) ** 2)
    H = 1.0 / (1.0 + (r / cutoff) ** (2 * order))
    H[H < 0.005] = 0.0
    return H


class HarpFilter2D:
    """Bandpass filter isolating the first harmonic of each tag direction.

    Parameters
    ----------
    image
        ``(rows, cols, encoding, frames)``, used to locate the harmonic peaks.
    central_freq
        ``(2,)`` encoding frequency per direction, rad/pixel.
    direction
        ``(2,)`` tag angles in radians.
    kind
        ``'transmission'``, ``'butterworth'``.
    cutoff, order
        Butterworth parameters.
    centre
        Optional ``(2, 2)`` of ``(row, col)`` peaks, skipping the search.
    """

    def __init__(self, image, central_freq, direction=(0.0, np.pi / 2),
                 kind="transmission", cutoff=25.0, order=7, search_window=(2, 2),
                 centre=None):
        image = np.asarray(image)
        self.shape = image.shape[:2]
        self.kind = str(kind).lower()
        self.central_freq = np.atleast_1d(np.asarray(central_freq, dtype=float))
        if self.central_freq.size == 1:
            self.central_freq = np.repeat(self.central_freq, 2)
        self.direction = np.asarray(direction, dtype=float)

        cc, ss = np.cos(self.direction), np.sin(self.direction)
        k1 = self.central_freq[0] * np.array([cc[0], ss[0]])
        k2 = self.central_freq[1] * np.array([cc[1], ss[1]])
        self.wave_vectors = (
            2 * np.pi / np.sum(k1**2) * k1,
            2 * np.pi / np.sum(k2**2) * k2,
        )

        if centre is None:
            centre = np.vstack([
                harmonic_center(itok(image[:, :, 0, 0]), k1, search_window, self.shape),
                harmonic_center(itok(image[:, :, 1, 0]), k2, search_window, self.shape),
            ])
        self.centre = np.asarray(centre, dtype=int)

        w1, w2 = self.wave_vectors
        if self.kind == "transmission":
            h1 = _transmission(self.shape, w1, w2)
            h2 = _transmission(self.shape, w2, w1)
        elif self.kind == "butterworth":
            h1 = _butterworth_2d(self.shape, self.centre[0], cutoff, order)
            h2 = _butterworth_2d(self.shape, self.centre[1], cutoff, order)
        else:
            raise ValueError(
                f"unknown filter kind {kind!r}; use 'transmission' or 'butterworth'"
            )

        self.kspace_filter = np.stack([h1, h2], axis=-1)

    def __call__(self, image):
        return self.filter(image)

    def filter(self, image):
        """Apply the filter, returning the complex harmonic image."""
        image = np.asarray(image)
        if image.shape[:2] != self.shape:
            raise ValueError(
                f"image is {image.shape[:2]} but this filter was built for {self.shape}"
            )
        # Broadcast over the frame axis.
        H = self.kspace_filter[..., np.newaxis]
        return ktoi(H * itok(image, axes=(0, 1)), axes=(0, 1))


def butterworth_3d(shape, wave_vec, order, cutoff=None):
    """3-D Butterworth bandpass centred on one wave vector.

    Follows the 3-D branch of Butterworth23D.m.
    """
    shape = tuple(int(s) for s in shape)
    wave_vec = np.asarray(wave_vec, dtype=float)

    dgrid = np.sqrt(np.sum(wave_vec**2))
    shift = wave_vec / dgrid**2
    if cutoff is None:
        cutoff = 0.8 * (1.0 / dgrid)

    axes = [_freq_range(n) for n in shape]
    X, Y, Z = np.meshgrid(axes[1], axes[0], axes[2], indexing="xy")
    X = X + 1e-8
    Y = Y + 1e-8
    Z = Z + 1e-8

    r = np.sqrt((X - shift[1]) ** 2 + (Y - shift[0]) ** 2 + (Z - shift[2]) ** 2)
    f = 1.0 / (1.0 + (r / cutoff) ** (2 * order))
    return f


class HarpFilter3D:
    """3-D harmonic bandpass, one Butterworth per encoding direction."""

    def __init__(self, shape, wave_vectors, order=5, cutoff=None):
        self.shape = tuple(int(s) for s in shape[:3])
        self.wave_vectors = np.asarray(wave_vectors, dtype=float)
        if self.wave_vectors.shape != (3, 3):
            raise ValueError("wave_vectors must be 3x3, one row per encoding direction")
        self.kspace_filter = np.stack(
            [butterworth_3d(self.shape, self.wave_vectors[i], order, cutoff) for i in range(3)],
            axis=-1,
        )

    def __call__(self, volume):
        return self.filter(volume)

    def filter(self, volume):
        volume = np.asarray(volume)
        if volume.shape[:3] != self.shape:
            raise ValueError(
                f"volume is {volume.shape[:3]} but this filter was built for {self.shape}"
            )
        H = self.kspace_filter[..., np.newaxis]
        return ktoi(H * itok(volume, axes=(0, 1, 2)), axes=(0, 1, 2))
