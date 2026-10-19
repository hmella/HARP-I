"""Temporal polynomial fitting of the estimated motion.

A 10th order polynomial is fitted per component to smooth the tissue
trajectories, with the value at t=0 and the slope at both ends held fixed
(the slopes come from local cubic fits):

    min ||A c - f||^2   subject to   Aeq c = beq

A and Aeq are the same for every pixel. The null space method is used rather
than normal equations, which would square the condition number of an already
poorly conditioned Vandermonde matrix.
"""

from __future__ import annotations

import warnings

import numpy as np

__all__ = ["fit_temporal"]


def _design(time: np.ndarray, order: int) -> np.ndarray:
    powers = np.arange(order, -1, -1)
    return time[:, None] ** powers[None, :]


def fit_temporal(u: np.ndarray, mask: np.ndarray, order: int = 10) -> np.ndarray:
    """Fit a constrained polynomial to each masked pixel's trajectory.

    Parameters
    ----------
    u
        Displacement ``(..., ndim, frames)``.
    mask
        ``(...)`` over the spatial dimensions. A trailing frame axis is
        accepted, in which case the first frame is used.
    order
        Polynomial order.

    Returns
    -------
    Fitted displacement, same shape as ``u``, NaN outside the mask.
    """
    u = np.asarray(u, dtype=float)
    nfr = u.shape[-1]
    ndim = u.shape[-2]
    spatial = u.shape[:-2]

    mask = np.asarray(mask, dtype=bool)
    if mask.ndim == len(spatial) + 1:
        mask = mask[..., 0]
    if mask.shape != spatial:
        raise ValueError(f"mask is {mask.shape} but u has spatial shape {spatial}")

    if nfr < 4:
        raise ValueError(f"temporal fitting needs at least 4 frames, got {nfr}")
    if order + 1 > nfr:
        raise ValueError(f"order {order} needs at least {order + 1} frames, got {nfr}")

    time = np.linspace(0.0, 1.0, nfr)
    powers = np.arange(order, -1, -1)
    M = _design(time, order)
    M4 = _design(time, 3)

    # value at t=0, slope at t=0, slope at t=end
    aeq = np.vstack([
        time[0] ** powers,
        powers * time[0] ** np.maximum(0, powers - 1),
        powers * time[-1] ** np.maximum(0, powers - 1),
    ])

    keep = np.ones(nfr, dtype=bool)
    keep[0] = False                      # t == 0 is fixed by the constraint
    A = M[keep, :]

    # Null space basis and a minimum norm particular solution
    _, sv, vt = np.linalg.svd(aeq)
    rank = int((sv > sv.max() * 1e-12).sum())
    Z = vt[rank:].T
    P = np.linalg.pinv(aeq)
    AZ = A @ Z

    cond = np.linalg.cond(AZ)
    if not np.isfinite(cond) or cond > 1e12:
        warnings.warn(
            f"temporal fit is ill-conditioned (cond = {cond:.2e}); "
            f"consider an order below {order} for {nfr} frames",
            RuntimeWarning,
            stacklevel=2,
        )

    out = np.full_like(u, np.nan)
    idx = np.flatnonzero(mask.ravel(order="F"))
    if idx.size == 0:
        return out

    for c in range(ndim):
        comp = u[..., c, :].reshape(-1, nfr, order="F")[idx, :].T   # (nfr, npix)

        c4a = np.linalg.lstsq(M4[:4, :], comp[:4, :], rcond=None)[0]
        c4b = np.linalg.lstsq(M4[-4:, :], comp[-4:, :], rcond=None)[0]
        t1, t2 = time[0], time[-1]
        d1 = 3 * c4a[0] * t1**2 + 2 * c4a[1] * t1 + c4a[2]
        d2 = 3 * c4b[0] * t2**2 + 2 * c4b[1] * t2 + c4b[2]

        beq = np.vstack([comp[0, :], d1, d2])
        c0 = P @ beq
        y = np.linalg.lstsq(AZ, comp[keep, :] - A @ c0, rcond=None)[0]
        coeff = c0 + Z @ y

        vals = (M @ coeff).T                                        # (npix, nfr)
        flat = out[..., c, :].reshape(-1, nfr, order="F")
        flat[idx, :] = vals
        out[..., c, :] = flat.reshape(spatial + (nfr,), order="F")

    return out
