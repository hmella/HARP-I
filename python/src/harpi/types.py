"""Data containers shared across the pipeline.

Array layout follows the MATLAB toolbox: ``(row, col[, slice], encoding,
frame)``. Pixel coordinates are 0-based and converted at the I/O boundary.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["Dataset2D", "Dataset3D", "MotionResult", "StrainResult", "SegmentStrain"]


@dataclass
class Dataset2D:
    """A 2-D tagged acquisition."""

    image: np.ndarray
    """``(rows, cols, encoding, frames)``."""
    mask: np.ndarray
    """``(rows, cols, frames)`` boolean."""
    pixel_size: np.ndarray
    """``(2,)`` in metres."""
    tag_spacing: float
    """Tag period in metres."""
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.image = np.asarray(self.image)
        self.mask = np.asarray(self.mask).astype(bool)
        self.pixel_size = np.atleast_1d(np.asarray(self.pixel_size, dtype=float))

    @property
    def shape(self) -> tuple[int, int]:
        return self.image.shape[0], self.image.shape[1]

    @property
    def n_frames(self) -> int:
        return self.image.shape[-1]

    @property
    def ke(self) -> np.ndarray:
        """Encoding frequency in radians per pixel, per axis."""
        return 2 * np.pi / self.tag_spacing * self.pixel_size

    @property
    def pixels_per_wavelength(self) -> float:
        return float(self.tag_spacing / self.pixel_size[0])

    @property
    def suggested_a(self) -> float:
        """Eq. (17): a = 150 / (pixels per wavelength)."""
        return 150.0 / self.pixels_per_wavelength


@dataclass
class Dataset3D:
    """A 3-D tagged acquisition."""

    image: np.ndarray
    """``(rows, cols, slices, encoding, frames)``."""
    mask: np.ndarray
    """``(rows, cols, slices, frames)`` boolean."""
    voxel_size: np.ndarray
    """``(3,)`` in the same units as ``tag_spacing``."""
    tag_spacing: float
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.image = np.asarray(self.image)
        self.mask = np.asarray(self.mask).astype(bool)
        self.voxel_size = np.atleast_1d(np.asarray(self.voxel_size, dtype=float))

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.image.shape[:3]

    @property
    def n_frames(self) -> int:
        return self.image.shape[-1]

    @property
    def pixels_per_wavelength(self) -> float:
        return float(self.tag_spacing / self.voxel_size[0])

    @property
    def suggested_a(self) -> float:
        return 150.0 / self.pixels_per_wavelength


@dataclass
class MotionResult:
    """Output of motion estimation."""

    u: np.ndarray
    """Raw displacement, ``(rows, cols[, slices], ndim, frames)``."""
    u_fitted: np.ndarray | None = None
    """Temporally fitted displacement, same shape as ``u``."""
    trajectories: np.ndarray | None = None
    unwrapped_phase: np.ndarray | None = None
    phase_jumps: np.ndarray | None = None
    """``(ndim, frames)`` multiples of 2*pi, one per frame."""
    quality: np.ndarray | None = None
    unwrap_complete: np.ndarray | None = None
    """Per frame, whether every masked pixel was unwrapped."""
    meta: dict = field(default_factory=dict)

    @property
    def best(self) -> np.ndarray:
        """Fitted motion when available, otherwise raw."""
        return self.u if self.u_fitted is None else self.u_fitted


@dataclass
class StrainResult:
    """Pixelwise strain tensor components, each ``(..., frames)`` with NaN outside the mask."""

    mask: np.ndarray
    XX: np.ndarray | None = None
    YY: np.ndarray | None = None
    XY: np.ndarray | None = None
    ZZ: np.ndarray | None = None
    RR: np.ndarray | None = None
    CC: np.ndarray | None = None
    LL: np.ndarray | None = None
    RC: np.ndarray | None = None
    p1: np.ndarray | None = None
    p2: np.ndarray | None = None
    p1or: np.ndarray | None = None
    origin: np.ndarray | None = None
    meta: dict = field(default_factory=dict)

    def global_curve(self, component: str) -> np.ndarray:
        """Spatial mean of a component over the mask, per frame."""
        v = getattr(self, component)
        if v is None:
            raise ValueError(f"{component} was not computed")
        axes = tuple(range(v.ndim - 1))
        return np.nanmean(v, axis=axes)


@dataclass
class SegmentStrain:
    """Segment-wise strain curves."""

    segments: np.ndarray
    """Segment index per pixel, NaN outside the mask."""
    CC: np.ndarray | None = None
    """``(n_segments, frames)``."""
    RR: np.ndarray | None = None
    n_segments: int = 6
