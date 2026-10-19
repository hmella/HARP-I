"""Locating the first harmonic peak in k-space."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["harmonic_center", "TagGeometry", "detect_tag_harmonics"]


def harmonic_center(kspace, omega, window=(2, 2), fov=None, center=None):
    """Refine the predicted harmonic peak by a local search.

    Parameters
    ----------
    kspace
        Centred k-space of one encoding, ``(rows, cols)``.
    omega
        Predicted angular frequency ``(wx, wy)`` in rad/pixel.
    window
        Half-size ``(col, row)`` of the search window.
    fov
        ``(rows, cols)`` used to convert omega into a bin offset.
    center
        Skip the prediction and search around this ``(col, row)``.

    Returns
    -------
    ``(row, col)`` of the peak.
    """
    k = np.asarray(kspace)
    if fov is None:
        fov = np.asarray(k.shape, dtype=float)
    fov = np.asarray(fov, dtype=float)
    window = np.asarray(window, dtype=int)

    if center is None:
        s = np.round(0.5 * np.asarray(k.shape, dtype=float))
        c = np.round(s + np.asarray(omega, dtype=float) * fov / (2 * np.pi)).astype(int)
    else:
        c = np.asarray(center, dtype=int)

    # c is (col, row) here.
    row0 = c[1] - window[1]
    row1 = c[1] + window[1]
    col0 = c[0] - window[0]
    col1 = c[0] + window[0]

    row0 = max(row0, 0)
    col0 = max(col0, 0)
    row1 = min(row1, k.shape[0] - 1)
    col1 = min(col1, k.shape[1] - 1)
    if row1 < row0 or col1 < col0:
        raise ValueError("harmonic search window fell outside k-space")

    win = np.abs(k[row0 : row1 + 1, col0 : col1 + 1])
    ii, jj = np.unravel_index(int(np.argmax(win)), win.shape)
    return np.array([row0 + ii, col0 + jj], dtype=int)


@dataclass
class TagGeometry:
    """The tag pattern read off the harmonic peaks of one frame.

    Clinical exports rarely record the tag spacing, and the angle is not always
    along the image axes, so both are measured from the image instead.
    """

    centres: np.ndarray
    """``(2, 2)`` of ``(row, col)`` peaks, one per direction."""
    angles: np.ndarray
    """``(2,)`` tag angles in radians."""
    spacing_px: np.ndarray
    """``(2,)`` tag period along each direction, in pixels."""
    confidence: float
    """Weakest accepted peak over the background median. Above about 5 is clear."""

    @property
    def pixels_per_wavelength(self) -> float:
        return float(np.mean(self.spacing_px))

    @property
    def suggested_a(self) -> float:
        """Eq. (17), from the measured spacing rather than a nominal one."""
        return 150.0 / self.pixels_per_wavelength

    def spacing_mm(self, pixel_spacing) -> float:
        """Tag period in the units of ``pixel_spacing``, usually millimetres."""
        ps = np.atleast_1d(np.asarray(pixel_spacing, dtype=float))
        return float(np.mean(self.spacing_px) * np.mean(ps[:2]))

    @property
    def plausible(self) -> bool:
        """Whether this looks like a tag pattern rather than anatomy.

        A tag period outside roughly 3 to 15 pixels is not a tagging sequence
        anyone would prescribe, and a weak second peak means only one direction
        was found, which rules out a grid.
        """
        return bool(np.all((self.spacing_px > 3.0) & (self.spacing_px < 15.0))
                    and self.confidence > 4.0)


def detect_tag_harmonics(frame, n_peaks=2, min_separation=6, dc_radius=4):
    """Find the tag harmonics of one frame with no prior.

    :func:`harmonic_center` refines a peak that is already predicted. This
    searches with nothing to go on, which is what a clinical export needs: the
    tag spacing is usually absent from the header and the angle is not
    necessarily along the image axes.

    Parameters
    ----------
    frame
        One image, ``(rows, cols)``.
    n_peaks
        How many directions to find. Two for a grid.
    min_separation
        Peaks closer than this many bins are treated as one.
    dc_radius
        Bins around the centre to ignore, where the low-frequency anatomy sits.

    Returns
    -------
    TagGeometry
    """
    img = np.asarray(frame, dtype=float)
    if img.ndim != 2:
        raise ValueError(f"expected one 2-D frame, got shape {img.shape}")
    rows, cols = img.shape

    mag = np.abs(np.fft.fftshift(np.fft.fft2(img)))
    r0, c0 = rows // 2, cols // 2

    rr, cc = np.ogrid[:rows, :cols]
    dist = np.hypot(rr - r0, cc - c0)
    search = mag.copy()
    search[dist <= dc_radius] = 0.0

    background = float(np.median(search[search > 0])) or 1.0

    # Work in cycles per field of view so a non-square matrix does not bias the
    # radius comparison below.
    order = np.argsort(search.ravel())[::-1][: 400 * max(n_peaks, 1)]
    cand = []
    for idx in order:
        r, c = divmod(int(idx), cols)
        dr, dc = r - r0, c - c0
        cand.append((dr, dc, float(search[r, c]),
                     float(np.hypot(dr / rows, dc / cols))))

    if not cand:
        raise ValueError("k-space held no peak outside the centre")

    accepted = [cand[0]]                       # the strongest is the first tag
    radius0 = accepted[0][3]

    # The remaining directions of a grid share the period of the first, so the
    # search is restricted to a similar radius. Without this the second peak is
    # often a strong low-frequency anatomy component instead of the orthogonal
    # tag: on the amrgatlas grid series that mistook 51 pixels for 5.3.
    for dr, dc, val, radius in cand[1:]:
        if len(accepted) >= n_peaks:
            break
        if not 0.6 * radius0 <= radius <= 1.6 * radius0:
            continue
        if any(np.hypot(dr - ar, dc - ac) < min_separation
               or np.hypot(dr + ar, dc + ac) < min_separation
               for ar, ac, _, _ in accepted):
            continue
        # A direction already covered is not a new one. Compare modulo pi,
        # since a peak and its conjugate describe the same direction.
        ang = np.arctan2(dr, dc)
        if any(abs(((ang - np.arctan2(ar, ac)) + np.pi / 2) % np.pi - np.pi / 2)
               < np.radians(25) for ar, ac, _, _ in accepted):
            continue
        accepted.append((dr, dc, val, radius))

    if len(accepted) < n_peaks:
        raise ValueError(
            f"found {len(accepted)} harmonic peaks, needed {n_peaks}; "
            "the image may not be tagged, or only one direction is tagged"
        )
    # Order by angle so direction 0 is always the shallower one. The peak
    # strengths decide which is found first, and that varies between slices of
    # the same study, which would otherwise swap the displacement components
    # from one slice to the next.
    accepted.sort(key=lambda t: np.arctan2(t[0], t[1]) % np.pi)
    accepted = [(a, b, v) for a, b, v, _ in accepted]

    centres = np.array([[ar + r0, ac + c0] for ar, ac, _ in accepted], dtype=float)
    angles = np.array([np.arctan2(ar, ac) for ar, ac, _ in accepted], dtype=float)
    # cycles per pixel along each axis, so the period is the reciprocal norm
    spacing = np.array([1.0 / np.hypot(ar / rows, ac / cols)
                        for ar, ac, _ in accepted], dtype=float)
    confidence = min(v for _, _, v in accepted) / background

    return TagGeometry(centres=centres, angles=angles, spacing_px=spacing,
                       confidence=float(confidence))
