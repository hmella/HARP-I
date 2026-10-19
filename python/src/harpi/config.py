"""Configuration objects, serialisable to TOML/JSON."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

__all__ = ["FilterConfig", "MotionConfig", "StrainConfig", "HarpiConfig"]


@dataclass
class FilterConfig:
    kind: str = "transmission"
    """'transmission' or 'butterworth'."""
    cutoff: float = 25.0
    order: int = 7
    direction_deg: tuple[float, float] = (0.0, 90.0)
    search_window: tuple[int, int] = (2, 2)


@dataclass
class MotionConfig:
    kernel: str = "multiquadric3"
    a: float | tuple[float, float] | None = None
    """Shape parameter. Eq. 17 uses 150 / pixels per wavelength."""
    a_dist: str = "constant"
    eta: float = 15.0
    """Tikhonov weight of eq. 8."""
    connectivity: int = 4
    seed: str = "auto"
    smooth_reference: bool = False
    temporal_fit: bool = False
    temporal_order: int = 10
    downsample: int = 1
    reference: int = 0
    njumps: int = 8
    phase_jump_rule: str = "nearest"
    """'published' (as in the MATLAB toolbox) or 'nearest' (2-D only)."""


@dataclass
class StrainConfig:
    n_segments: int = 6
    clockwise: bool = False
    outlier_threshold: float = 0.5
    origin: tuple[float, float] | None = None


@dataclass
class HarpiConfig:
    filter: FilterConfig = field(default_factory=FilterConfig)
    motion: MotionConfig = field(default_factory=MotionConfig)
    strain: StrainConfig = field(default_factory=StrainConfig)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "HarpiConfig":
        return cls(
            filter=FilterConfig(**d.get("filter", {})),
            motion=MotionConfig(**d.get("motion", {})),
            strain=StrainConfig(**d.get("strain", {})),
        )
