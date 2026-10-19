# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

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


def smooth_reference_phase(phases, kernel, a, eta, max_points=12000):
    """Smooth a set of co-registered phase maps.

    Parameters
    ----------
    phases
        Sequence of arrays of equal shape, NaN outside the mask.
    kernel, a, eta
        Kernel name or record, shape parameter, Tikhonov weight.
    max_points
        Refuse rather than attempt a dense system larger than this. One kernel
        matrix is held for all the valid points at once, so the memory grows as
        the square of their number: 12000 points is about 1.1 GiB, and a
        full-resolution volume reaches tens of gibibytes. Raise it deliberately
        if the machine can take it.

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

    n_valid = int(valid.sum())
    if max_points is not None and n_valid > max_points:
        need = n_valid * n_valid * 8 / 2 ** 30
        raise ValueError(
            f"smoothing {n_valid} points needs a dense {n_valid} by {n_valid} "
            f"matrix, about {need:.1f} GiB. Smooth a downsampled set instead, "
            f"or raise max_points past {n_valid} if there is room."
        )

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
