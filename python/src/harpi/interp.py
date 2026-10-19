"""Tikhonov-regularised RBF interpolation.

Eq. (8) replaces the exact interpolation condition with

    (Psi + eta*I) lam = G

where ``Psi[m, n] = phi(||c_m - c_n||)`` over the centres and ``G`` holds
the values being interpolated. ``eta`` is an absolute shift on the diagonal,
so its useful magnitude depends on the kernel and on the scale of the data.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.spatial.distance import cdist

from .rbf import RBFKernel, get_kernel

__all__ = ["distance_matrix", "solve_tikhonov", "rbf_interpolate"]


def distance_matrix(a: np.ndarray, b: np.ndarray | None = None) -> np.ndarray:
    """Euclidean distances between point sets given as ``(D, N)`` arrays.

    Returns ``(N_a, N_b)``; with ``b`` omitted, the self-distance matrix.
    """
    A = np.ascontiguousarray(np.asarray(a, dtype=float).T)
    B = A if b is None else np.ascontiguousarray(np.asarray(b, dtype=float).T)
    return cdist(A, B)


def solve_tikhonov(psi: np.ndarray, g: np.ndarray, eta: float) -> np.ndarray:
    """Solve ``(psi + eta*I) lam = g``.

    Tries Cholesky first and falls back to a general solve. The cubic
    multiquadric and thin plate kernels are conditionally positive definite,
    so Cholesky can fail on them.
    """
    n = psi.shape[0]
    if psi.shape[1] != n:
        raise ValueError(f"psi must be square, got {psi.shape}")
    if g.shape[0] != n:
        raise ValueError(f"g has {g.shape[0]} rows but psi is {n}x{n}")

    a = psi
    if eta:
        a = psi + eta * np.eye(n, dtype=psi.dtype)

    try:
        return cho_solve(cho_factor(a, lower=True, check_finite=False), g, check_finite=False)
    except np.linalg.LinAlgError:
        return np.linalg.solve(a, g)


def rbf_interpolate(
    centres: np.ndarray,
    values: np.ndarray,
    query: np.ndarray,
    kernel: str | RBFKernel,
    a: float,
    eta: float,
) -> np.ndarray:
    """Fit an RBF through ``values`` at ``centres`` and evaluate at ``query``.

    Parameters
    ----------
    centres
        ``(D, N)`` interpolation centres.
    values
        ``(N, C)`` values to reproduce (C components).
    query
        ``(D, M)`` evaluation points.
    kernel, a, eta
        Kernel name or record, its shape parameter, and the Tikhonov weight.

    Returns
    -------
    ``(M, C)`` interpolated values.
    """
    k = get_kernel(kernel)
    values = np.atleast_2d(np.asarray(values, dtype=float))
    if values.shape[0] != centres.shape[1]:
        values = values.T

    # Release the N x N system before allocating the M x N evaluation matrix.
    psi = k(distance_matrix(centres), a)
    w = solve_tikhonov(psi, values, eta)
    del psi

    return k(distance_matrix(query, centres), a) @ w
