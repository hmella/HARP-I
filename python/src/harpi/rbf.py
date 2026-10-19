"""Radial basis function kernels.

Each kernel is a function plus the metadata callers need to use it.

The shape parameter ``a`` means different things across the family: a width
for the Gaussian, an inverse width for the inverse multiquadric, a support
radius for the compactly supported kernels, an additive offset for the
multiquadrics, and nothing at all for the thin plate spline. The suggested
values below come from the table in examples/README.md.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

__all__ = ["RBFKernel", "KERNELS", "get_kernel"]


@dataclass(frozen=True)
class RBFKernel:
    """A radial basis function and the metadata needed to use it."""

    name: str
    phi: Callable[[np.ndarray, float], np.ndarray]
    compact: bool = False
    """True when phi(r) == 0 for r >= a."""
    suggested_a_scale: float | None = None
    """Divide by the pixels per tag wavelength for a starting shape parameter."""
    suggested_eta: float = 1e-8
    """Starting Tikhonov weight for eq. (8)."""
    uses_shape: bool = True
    """False when phi ignores a."""

    def __call__(self, r: np.ndarray, a: float) -> np.ndarray:
        return self.phi(r, a)


def _gaussian(r, a):
    return np.exp(-((r / a) ** 2))


def _imq(r, a):
    return 1.0 / np.sqrt(1.0 + (a * r) ** 2)


def _linear(r, a):
    return r / a


def _multiquadric(r, a):
    return np.sqrt(r**2 + a**2)


def _multiquadric3(r, a):
    return (r**2 + a**2) ** 1.5


def _thin_plate(r, a):
    # r^2 log r, with r == 0 set to 0.
    out = np.zeros_like(r)
    nz = r > 0
    out[nz] = r[nz] ** 2 * np.log(r[nz])
    return out


def _wendland_c2(r, a):
    u = np.minimum(r / a, 1.0)
    return (1 - u) ** 4 * (4 * u + 1)


def _wendland_c4(r, a):
    u = np.minimum(r / a, 1.0)
    return (1 - u) ** 6 * (35 * u**2 + 18 * u + 3)


def _wu_c2(r, a):
    u = np.minimum(r / a, 1.0)
    return (1 - u) ** 5 * (8 + 40 * u + 48 * u**2 + 25 * u**3 + 5 * u**4)


def _wu_c4(r, a):
    u = np.minimum(r / a, 1.0)
    return (1 - u) ** 6 * (6 + 36 * u + 82 * u**2 + 72 * u**3 + 30 * u**4 + 5 * u**5)


KERNELS: dict[str, RBFKernel] = {
    k.name: k
    for k in (
        RBFKernel("gaussian", _gaussian, suggested_a_scale=140.0, suggested_eta=2e-3),
        RBFKernel("imq", _imq, suggested_a_scale=0.2, suggested_eta=1e-4),
        RBFKernel("linear", _linear, suggested_a_scale=100.0, suggested_eta=4.0),
        RBFKernel("multiquadric", _multiquadric, suggested_a_scale=100.0, suggested_eta=4.0),
        # Eq. (16), with a = 150 / (pixels per wavelength).
        RBFKernel("multiquadric3", _multiquadric3, suggested_a_scale=150.0, suggested_eta=15.0),
        RBFKernel("thinplate", _thin_plate, suggested_eta=45.0, uses_shape=False),
        RBFKernel("wendlandc2", _wendland_c2, compact=True, suggested_a_scale=900.0,
                  suggested_eta=5e-3),
        RBFKernel("wendlandc4", _wendland_c4, compact=True, suggested_a_scale=600.0,
                  suggested_eta=1e-2),
        RBFKernel("wuc2", _wu_c2, compact=True, suggested_a_scale=1200.0, suggested_eta=1e-2),
        RBFKernel("wuc4", _wu_c4, compact=True, suggested_a_scale=600.0, suggested_eta=1e-2),
    )
}

# MATLAB class names, accepted so existing scripts keep working.
_ALIASES = {
    "Gaussian": "gaussian",
    "imqx": "imq",
    "Linear": "linear",
    "Multiquadric": "multiquadric",
    "Multiquadric3": "multiquadric3",
    "ThinPlate": "thinplate",
    "WendlandC2": "wendlandc2",
    "WendlandC4": "wendlandc4",
    "WuC2": "wuc2",
    "WuC4": "wuc4",
}


def get_kernel(name: str | RBFKernel) -> RBFKernel:
    """Look up a kernel by name. MATLAB class names also work."""
    if isinstance(name, RBFKernel):
        return name
    key = _ALIASES.get(name, str(name).lower().replace("_", "").replace("-", ""))
    try:
        return KERNELS[key]
    except KeyError:
        raise ValueError(
            f"unknown RBF kernel {name!r}; available: {', '.join(sorted(KERNELS))}"
        ) from None
