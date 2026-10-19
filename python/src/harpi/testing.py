"""Synthetic data with known ground truth.

A tagged ring undergoing radial contraction and twist, with optional tag
fading and noise, following the in-silico design on p. 1244. The
deformation is prescribed, so the recovered displacement and strain can be
checked against it.
"""

from __future__ import annotations

import numpy as np

__all__ = ["make_phantom", "analytic_displacement"]


def analytic_displacement(x, y, centre, contraction, twist):
    """Displacement of a radial contraction combined with a rigid twist.

    A material point at radius r moves to ``(1 - contraction) * r`` and
    rotates by ``twist * r / r_max``.
    """
    dx = x - centre[0]
    dy = y - centre[1]
    r = np.hypot(dx, dy)
    th = np.arctan2(dy, dx)

    rmax = np.nanmax(r) if np.isfinite(r).any() else 1.0
    scale = 1.0 - contraction
    angle = twist * (r / max(rmax, 1e-12))

    nx = centre[0] + scale * r * np.cos(th + angle)
    ny = centre[1] + scale * r * np.sin(th + angle)
    return nx - x, ny - y


def make_phantom(shape=(128, 128), n_frames=12, tag_spacing=8.0,
                 r_inner=18.0, r_outer=32.0, contraction=0.15, twist=0.25,
                 t1=0.72, noise=0.0, seed=0):
    """Build a synthetic tagged ring with a known deformation.

    Returns
    -------
    dict with
        ``image``   ``(rows, cols, 2, frames)`` real tagged images,
        ``mask``    ``(rows, cols, frames)``,
        ``u_true``  ``(rows, cols, 2, frames)`` Lagrangian displacement,
        ``tag_spacing``, ``pixel_size``, ``centre``, ``params``.

    Notes
    -----
    The tag pattern is sampled at the undeformed position of each pixel.
    """
    rng = np.random.default_rng(seed)
    rows, cols = shape
    centre = ((cols - 1) / 2.0, (rows - 1) / 2.0)

    X, Y = np.meshgrid(np.arange(cols, dtype=float), np.arange(rows, dtype=float))
    omega = 2 * np.pi / tag_spacing

    image = np.zeros((rows, cols, 2, n_frames))
    mask = np.zeros((rows, cols, n_frames), dtype=bool)
    u_true = np.zeros((rows, cols, 2, n_frames))

    for k in range(n_frames):
        # Contract and untwist over the first half, relax over the second.
        phase = np.sin(np.pi * k / max(n_frames - 1, 1))
        c = contraction * phase
        w = twist * phase

        ux, uy = analytic_displacement(X, Y, centre, c, w)
        u_true[:, :, 0, k] = ux
        u_true[:, :, 1, k] = uy

        # Material coordinates of the pixel currently at (X, Y).
        rx, ry = analytic_displacement(X, Y, centre,
                                       -c / max(1 - c, 1e-9), -w)
        Xref = X + rx
        Yref = Y + ry

        fade = 1.0 if t1 is None else np.exp(-(k / max(n_frames - 1, 1)) / t1)
        dc = 0.0 if t1 is None else 1.0 - fade

        image[:, :, 0, k] = dc + fade * np.cos(omega * Xref)
        image[:, :, 1, k] = dc + fade * np.cos(omega * Yref)

        r = np.hypot(X - centre[0], Y - centre[1])
        mask[:, :, k] = (r >= r_inner) & (r <= r_outer)

    if noise:
        image = image + rng.normal(0.0, noise, image.shape)

    return {
        "image": image,
        "mask": mask,
        "u_true": u_true,
        "tag_spacing": tag_spacing,
        "pixel_size": np.array([1.0, 1.0]),
        "centre": np.array(centre),
        "params": {"contraction": contraction, "twist": twist,
                   "r_inner": r_inner, "r_outer": r_outer, "noise": noise},
    }
