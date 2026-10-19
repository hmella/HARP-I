"""Following material points through a tagged sequence.

The unwrapped harmonic phase labels tissue rather than space, which is what
lets HARP-I treat it as a material coordinate. The same fact lets a single
point be followed: the place it moved to in the next frame is the place
carrying the phase it had in this one.

Written against the wrapped phase, so nothing here unwraps and nothing here
needs a mask. That matters for segmentation, where the mask is the thing being
produced: :func:`~harpi.motion.estimate_motion_2d` cannot help until a mask
exists, and this can.
"""

from __future__ import annotations

import numpy as np

__all__ = ["track_points_2d"]


def _sample(field, x, y):
    """Bilinear sample of a complex image at ``(x=column, y=row)``."""
    from scipy.ndimage import map_coordinates

    coords = np.vstack([y, x])
    re = map_coordinates(field.real, coords, order=1, mode="nearest")
    im = map_coordinates(field.imag, coords, order=1, mode="nearest")
    return re + 1j * im


def _wrap(angle):
    """Fold an angle into (-pi, pi]."""
    return np.angle(np.exp(1j * np.asarray(angle)))


def _hop(harmonic, x, y, target, frame, iterations, max_step, gradient_step):
    """Newton from one frame to the next, on the wrapped phase.

    The residual is the wrapped difference between the phase here and the phase
    the point carried, so a whole number of periods never enters and no
    unwrapping is required. The Jacobian is measured from the image rather than
    taken from the nominal wave vectors, so it follows the local deformation.
    """
    x = np.asarray(x, dtype=float).copy()
    y = np.asarray(y, dtype=float).copy()
    ndim = harmonic.shape[2]
    h = float(gradient_step)

    for _ in range(int(iterations)):
        jac = np.zeros((x.size, ndim, 2))
        residual = np.zeros((x.size, ndim))
        for d in range(ndim):
            plane = harmonic[:, :, d, frame]
            residual[:, d] = _wrap(np.angle(_sample(plane, x, y)) - target[d])
            for k, (dx, dy) in enumerate(((h, 0.0), (0.0, h))):
                ahead = np.angle(_sample(plane, x + dx, y + dy))
                behind = np.angle(_sample(plane, x - dx, y - dy))
                jac[:, d, k] = _wrap(ahead - behind) / (2 * h)

        det = jac[:, 0, 0] * jac[:, 1, 1] - jac[:, 0, 1] * jac[:, 1, 0]
        usable = np.abs(det) > 1e-6
        step_x = np.zeros_like(x)
        step_y = np.zeros_like(y)
        step_x[usable] = (jac[usable, 1, 1] * residual[usable, 0]
                          - jac[usable, 0, 1] * residual[usable, 1]) / det[usable]
        step_y[usable] = (-jac[usable, 1, 0] * residual[usable, 0]
                          + jac[usable, 0, 0] * residual[usable, 1]) / det[usable]
        x -= np.clip(step_x, -max_step, max_step)
        y -= np.clip(step_y, -max_step, max_step)
    return x, y


def _consensus(step, modes, iterations=8, cut=2.0):
    """Where the contour as a whole says each node went.

    The nodes lie on a closed curve, so their displacements are periodic around
    it and the first few Fourier modes already describe a translation, a
    dilation and the usual bending. Reweighting the fit by its own residual
    keeps the slipped nodes from dragging it along, which the median of three
    immediate neighbours cannot do once two of the three slip together.

    The weights start from the median step rather than from unity. A run of
    adjacent nodes that slipped by the same period looks like a legitimate bump
    the modes can follow, and an even start lets the fit settle on it and report
    the slip as the answer; a median of every node is not moved by a few of
    them, so the first weights already mark them out.
    """
    n = len(step)
    angle = 2 * np.pi * np.arange(n) / n
    columns = [np.ones(n)]
    for k in range(1, modes + 1):
        columns += [np.cos(k * angle), np.sin(k * angle)]
    basis = np.stack(columns, axis=1)

    fit = np.empty_like(step)
    fit[:] = np.median(step, axis=0)
    for _ in range(iterations):
        residual = np.hypot(*(step - fit).T)
        scale = np.median(residual) + 1e-6
        root = np.sqrt(1.0 / (1.0 + (residual / (cut * scale)) ** 2))
        for c in range(2):
            coefficients, *_ = np.linalg.lstsq(basis * root[:, None],
                                               step[:, c] * root, rcond=None)
            fit[:, c] = basis @ coefficients
    return fit


def _regularise(before, after, tolerance, modes=3):
    """Replace a node that disagrees with the contour about where it went.

    The phase repeats every tag period, so a node can settle on the
    neighbouring tag line and land a whole period away. The rest of the contour
    does not make the same mistake, which is what makes it detectable. Only the
    nodes past the tolerance are replaced, so a node that tracked correctly
    keeps the position measured for it.

    Short contours fall back to the median of a node and its two neighbours,
    there being too few nodes to fit the modes. The nodes form a ring, so that
    neighbourhood wraps.

    Returns the corrected positions and how many nodes were replaced.
    """
    step = after - before
    n = len(step)
    if n < 3:
        return after, 0

    if n >= 2 * modes + 2:
        reference = _consensus(step, modes)
    else:
        index = np.arange(n)
        ring = np.stack([step[(index - 1) % n], step, step[(index + 1) % n]])
        reference = np.median(ring, axis=0)

    fixed = step.copy()
    bad = np.abs(step - reference) > tolerance
    fixed[bad] = reference[bad]
    return before + fixed, int(np.any(bad, axis=1).sum())


def track_points_2d(harmonic_image, points, reference=0, period=None,
                    iterations=12, max_step=2.0, gradient_step=0.5,
                    regularise=True, progress=None, return_info=False):
    """Follow points through every frame of a tagged sequence.

    Parameters
    ----------
    harmonic_image
        ``(rows, cols, 2, frames)`` complex, as
        :class:`~harpi.filters.HarpFilter2D` returns. No mask is needed.
    points
        ``(n, 2)`` of ``(x, y)``, that is ``(column, row)``, in the reference
        frame. Treated as one contour: the regularisation assumes they are
        neighbours along a ring.
    reference
        Frame the points were placed on. Tracking runs outwards in both
        directions from it, one frame at a time, because a single hop is small
        enough to stay inside the half-period the phase can resolve. Going
        straight from the reference to a distant frame is measurably worse.
    period
        Tag period in pixels, used only to size the regularisation tolerance,
        which is half of it. Taken as 2*pi over the median phase gradient when
        omitted.
    iterations, max_step, gradient_step
        Newton iterations per hop, the largest step a single iteration may take
        in pixels, and the half-width of the central difference used for the
        phase gradient.
    regularise
        Replace a node that moves a period away from the rest of the
        contour. Rare, and
        worth having: one uncorrected slip corrupts the starting point of every
        later frame.
    progress
        Called as ``progress(done, total, stage)`` after each hop, counting the
        hops of the whole run. Raising from it stops the run, which is how the
        interface cancels.
    return_info
        Also return a dict with ``corrected``, the number of nodes replaced per
        frame.

    Returns
    -------
    ``(frames, n, 2)`` of tracked positions. The reference frame holds the input
    unchanged.
    """
    harmonic = np.asarray(harmonic_image)
    if harmonic.ndim != 4 or harmonic.shape[2] != 2:
        raise ValueError("harmonic_image must be (rows, cols, 2, frames)")
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError(f"points must be (n, 2) of (x, y), got {pts.shape}")

    nfr = harmonic.shape[3]
    reference = int(reference)
    if not 0 <= reference < nfr:
        raise ValueError(f"reference frame {reference} is outside 0..{nfr - 1}")

    if period is None:
        period = _period_from(harmonic, reference)
    tolerance = 0.5 * float(period)

    out = np.repeat(pts[None, :, :], nfr, axis=0)
    corrected = np.zeros(nfr, dtype=int)
    done, total = 0, max(nfr - 1, 1)

    for direction in (+1, -1):
        x, y = pts[:, 0].copy(), pts[:, 1].copy()
        frame = reference
        while True:
            nxt = frame + direction
            if not 0 <= nxt < nfr:
                break
            target = np.stack([np.angle(_sample(harmonic[:, :, d, frame], x, y))
                               for d in range(2)])
            moved_x, moved_y = _hop(harmonic, x, y, target, nxt,
                                    iterations, max_step, gradient_step)
            before = np.column_stack([x, y])
            after = np.column_stack([moved_x, moved_y])
            if regularise:
                after, n_fixed = _regularise(before, after, tolerance)
                corrected[nxt] = n_fixed
            out[nxt] = after
            x, y = after[:, 0], after[:, 1]
            frame = nxt
            done += 1
            if progress is not None:
                progress(done, total, "track")

    if return_info:
        return out, {"corrected": corrected, "period": float(period),
                     "reference": reference}
    return out


def _period_from(harmonic, frame):
    """Tag period in pixels, from the typical phase gradient of a frame."""
    grads = []
    for d in range(harmonic.shape[2]):
        phase = np.angle(harmonic[:, :, d, frame])
        gy, gx = np.gradient(phase)
        # the wrapped jumps are outliers, so a median is the robust summary
        grads.append(np.median(np.abs(_wrap(gx))))
        grads.append(np.median(np.abs(_wrap(gy))))
    strongest = max(grads)
    if strongest < 1e-6:
        return 8.0
    return float(2 * np.pi / strongest)
