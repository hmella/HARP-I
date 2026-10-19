"""Regional analysis of strain results.

:func:`~harpi.strain.strain_by_segments` already averages strain over
equiangular wedges. What it does not do is name those wedges, or reduce the
curves to the few numbers a report needs. Both are added here, kept free of Qt
so they can be tested and scripted.

Segment numbering follows ``strain_by_segments``: segment 1 begins at the
reference angle fixed by ``position_a`` and ``position_b``, and the rest follow
in the chosen direction. For the labels below to mean what they say, that
reference angle has to be the anterior right-ventricular insertion point.
"""

from __future__ import annotations

import numpy as np

__all__ = ["AHA_LEVELS", "aha_labels", "segment_metrics", "bullseye_rings",
           "segment_table"]


# The AHA 17-segment model, in the order segments are numbered starting from the
# anterior insertion. Circulation 2002;105:539-542.
AHA_LEVELS: dict[str, tuple[str, ...]] = {
    "basal": ("basal anterior", "basal anteroseptal", "basal inferoseptal",
              "basal inferior", "basal inferolateral", "basal anterolateral"),
    "mid": ("mid anterior", "mid anteroseptal", "mid inferoseptal",
            "mid inferior", "mid inferolateral", "mid anterolateral"),
    "apical": ("apical anterior", "apical septal", "apical inferior",
               "apical lateral"),
}

NATIVE_SEGMENTS: dict[str, int] = {"basal": 6, "mid": 6, "apical": 4}


def aha_labels(level: str | None, n_segments: int) -> list[str]:
    """Names for ``n_segments`` wedges at a short-axis ``level``.

    ``level`` is ``basal``, ``mid`` or ``apical``. When it is None, or when
    ``n_segments`` does not match the model for that level, plain numbers are
    returned instead of names, because a label that does not correspond to the
    wedge it sits on is worse than no label.
    """
    if level is None:
        return [str(i + 1) for i in range(n_segments)]
    key = str(level).lower()
    if key not in AHA_LEVELS:
        raise ValueError(
            f"unknown level {level!r}; use one of {', '.join(AHA_LEVELS)}")
    names = AHA_LEVELS[key]
    if n_segments != len(names):
        return [str(i + 1) for i in range(n_segments)]
    return list(names)


def segment_metrics(curves, times=None, component="CC"):
    """Reduce per-segment curves to peak strain and time to peak.

    Parameters
    ----------
    curves
        ``(n_segments, frames)`` as produced by ``strain_by_segments``.
    times
        Frame times. Frame indices starting at 0 are used when omitted, so the
        time to peak is then a frame count rather than a duration.
    component
        Used only to choose the sign convention. Circumferential and radial
        strain shorten, so their peak is the most negative value; the principal
        strains and the direct components are taken at their largest magnitude.

    Returns
    -------
    dict of arrays, each ``(n_segments,)``: ``peak``, ``time_to_peak``,
    ``peak_frame`` and ``final``.
    """
    curves = np.asarray(curves, dtype=float)
    if curves.ndim != 2:
        raise ValueError(f"expected (n_segments, frames), got {curves.shape}")
    n_seg, nfr = curves.shape
    t = np.arange(nfr, dtype=float) if times is None else np.asarray(times, dtype=float)
    if t.size != nfr:
        raise ValueError(f"times has {t.size} entries for {nfr} frames")

    shortening = str(component).upper() in {"CC", "RR", "LL"}
    peak = np.full(n_seg, np.nan)
    frame = np.full(n_seg, -1, dtype=int)
    for s in range(n_seg):
        row = curves[s]
        good = np.isfinite(row)
        if not good.any():
            continue
        idx = np.flatnonzero(good)
        pick = idx[np.argmin(row[idx])] if shortening else idx[np.argmax(np.abs(row[idx]))]
        peak[s] = row[pick]
        frame[s] = pick

    ttp = np.where(frame >= 0, t[np.clip(frame, 0, nfr - 1)], np.nan)
    final = np.array([
        curves[s][np.flatnonzero(np.isfinite(curves[s]))[-1]]
        if np.isfinite(curves[s]).any() else np.nan
        for s in range(n_seg)
    ])
    return {"peak": peak, "time_to_peak": ttp, "peak_frame": frame, "final": final}


def bullseye_rings(values, n_segments=None):
    """Wedge geometry for a polar bullseye plot.

    ``values`` is either a 1-D array of one number per segment, giving a single
    ring, or a list of such arrays, drawn from the outside in so that basal sits
    outermost as convention expects.

    Returns a list of ``(theta_start, width, value)`` tuples per ring, together
    with the radial extent of each ring, ready to hand to ``matplotlib`` bars on
    a polar axis. The caller keeps control of the colour mapping.
    """
    rings = [np.asarray(values, dtype=float)] if np.ndim(values) == 1 else \
        [np.asarray(v, dtype=float) for v in values]
    if n_segments is not None:
        rings = [r for r in rings if r.size == n_segments]
    if not rings:
        raise ValueError("no ring had any values")

    n_rings = len(rings)
    width = 1.0 / n_rings
    out = []
    for i, ring in enumerate(rings):
        # ring 0 is drawn outermost
        r_inner = 1.0 - (i + 1) * width
        step = 2 * np.pi / ring.size
        wedges = [(k * step, step, float(ring[k])) for k in range(ring.size)]
        out.append({"radius": r_inner, "height": width, "wedges": wedges})
    return out


def segment_table(seg, times=None, level=None, components=("CC", "RR")):
    """A row per segment, for display or export.

    ``seg`` is a ``SegmentStrain``. Returns ``(header, rows)`` with plain Python
    values, so a CSV writer or a table widget can consume it directly.
    """
    labels = aha_labels(level, seg.n_segments)
    present = [c for c in components if getattr(seg, c, None) is not None]
    if not present:
        raise ValueError("the segment result holds none of "
                         f"{', '.join(components)}")

    header = ["segment", "name"]
    stats = {}
    for comp in present:
        stats[comp] = segment_metrics(getattr(seg, comp), times, comp)
        header += [f"peak {comp}", f"time to peak {comp}", f"final {comp}"]

    rows = []
    for s in range(seg.n_segments):
        row = [s + 1, labels[s]]
        for comp in present:
            m = stats[comp]
            row += [float(m["peak"][s]), float(m["time_to_peak"][s]),
                    float(m["final"][s])]
        rows.append(row)
    return header, rows
