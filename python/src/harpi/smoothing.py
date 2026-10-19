"""Reference phase smoothing.

An RBF self interpolation: the centres and the evaluation points are the
same pixels, so the smoothing comes entirely from the Tikhonov term. With
``eta = 0`` it is the identity.
"""

from __future__ import annotations

import numpy as np

from .interp import distance_matrix, solve_tikhonov
from .rbf import get_kernel

__all__ = ["smooth_reference_phase"]


def smooth_reference_phase(phases, kernel, a, eta):
    """Smooth a set of co-registered phase maps.

    Parameters
    ----------
    phases
        Sequence of arrays of equal shape, NaN outside the mask.
    kernel, a, eta
        Kernel name or record, shape parameter, Tikhonov weight.

    Returns
    -------
    List of smoothed arrays, NaN preserved.
    """
    k = get_kernel(kernel)
    phases = [np.asarray(p, dtype=float) for p in phases]
    shape = phases[0].shape

    valid = np.ones(shape, dtype=bool)
    for p in phases:
        valid &= ~np.isnan(p)

    grids = np.meshgrid(*[np.arange(1, n + 1) for n in shape], indexing="ij")
    # X varies along the columns.
    coords = np.stack([g[valid].astype(float) for g in grids[::-1]], axis=0)

    # Centres and query points coincide, so one kernel matrix is enough:
    # smoothed = B @ ((B + eta*I) \ f)
    B = k(distance_matrix(coords), a)
    f = np.stack([p[valid] for p in phases], axis=1)
    w = solve_tikhonov(B, f, eta)
    out_vals = B @ w

    out = []
    for i in range(len(phases)):
        o = np.full(shape, np.nan)
        o[valid] = out_vals[:, i]
        out.append(o)
    return out
