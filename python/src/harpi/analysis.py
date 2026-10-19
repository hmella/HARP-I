"""Regional analysis of strain results.

:func:`~harpi.strain.strain_by_segments` already averages strain over
equiangular wedges. What it does not do is name those wedges, or reduce the
curves to the few numbers a report needs. Both are added here, kept free of Qt
so they can be tested and scripted.

The angular layout lives here too, beside the names, because the numbering only
means anything through the names: when the two were in separate files they
disagreed, and every label sat on the wrong wall.

The anterior right-ventricular insertion separates the septum from the anterior
wall, so it is the boundary between segment 1 and segment 2, and segment 1
occupies the step before it. The wedges stay equiangular, 60 degrees for six
and 90 for four, as the AHA statement specifies: the insertions orient the
model, they do not stretch it.
"""

from __future__ import annotations

import numpy as np

__all__ = ["AHA_LEVELS", "aha_labels", "segment_metrics", "bullseye_rings",
           "segment_table", "segment_geometry", "segment_index",
           "short_labels",
           "expected_septal_arc"]


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

# How many of those segments the septum occupies, which is what the two
# insertion points should span if the heart matched the model exactly.
AHA_SEPTAL_SEGMENTS: dict[str, int] = {"basal": 2, "mid": 2, "apical": 1}

# Below this the two insertions are too close, or too nearly opposite, to say
# which way round the septum runs.
_DEGENERATE = np.radians(5.0)


def _wrap(angle):
    """Fold an angle into (-pi, pi]."""
    return np.angle(np.exp(1j * np.asarray(angle, dtype=float)))


def expected_septal_arc(level, n_segments):
    """How wide the septum should be, in radians, or None.

    None when the level is unknown or the segment count does not match the
    model for it, since then there is nothing to compare a measurement with.
    """
    if level is None:
        return None
    key = str(level).lower()
    if NATIVE_SEGMENTS.get(key) != n_segments:
        return None
    return AHA_SEPTAL_SEGMENTS[key] * 2 * np.pi / n_segments


def _angle_to(origin, point):
    """Angle of a point seen from the origin, or None if there is no ray."""
    if point is None:
        return None
    x, y = float(point[0]), float(point[1])
    if not (np.isfinite(x) and np.isfinite(y)):
        return None
    dx, dy = x - float(origin[0]), y - float(origin[1])
    if np.hypot(dx, dy) < 1e-9:
        return None
    return float(np.arctan2(dy, dx))


def segment_geometry(n_segments=6, origin=(0.0, 0.0), anterior=None,
                     inferior=None, clockwise=False):
    """Where the segment boundaries fall and which way the numbering runs.

    ``anterior`` and ``inferior`` are the right-ventricular insertion points as
    ``(x, y)``, in the same grid as ``origin``. The anterior one fixes the
    boundary between segment 1 and segment 2. The two together fix the
    direction: the septum is the shorter arc between them, so going from
    anterior to inferior through it is the direction of increasing segment
    number, and ``clockwise`` is only consulted when that cannot be read.

    Returns a dict. ``direction`` is +1 when the segment number increases with
    the angle and -1 when it decreases, ``septal_arc`` is what the two points
    actually span, and ``direction_from`` says which rule decided.
    """
    step = 2 * np.pi / int(n_segments)
    theta_ant = _angle_to(origin, anterior)
    theta_inf = _angle_to(origin, inferior)

    # with no anterior point the ray runs along +x, which is what the MATLAB
    # defaults of position_a and position_b produce
    anterior_angle = 0.0 if theta_ant is None else theta_ant

    septal_arc = None
    if theta_ant is not None and theta_inf is not None:
        gap = _wrap(theta_inf - theta_ant)
        septal_arc = abs(float(gap))
        if _DEGENERATE < septal_arc < np.pi - _DEGENERATE:
            direction = 1 if gap > 0 else -1
            return {"anterior_angle": anterior_angle, "direction": direction,
                    "step": step, "edges": _edges(anterior_angle, direction,
                                                  step, n_segments),
                    "septal_arc": septal_arc, "direction_from": "insertions",
                    "origin": tuple(origin), "anterior": anterior,
                    "inferior": inferior}

    direction = 1 if clockwise else -1
    return {"anterior_angle": anterior_angle, "direction": direction,
            "step": step,
            "edges": _edges(anterior_angle, direction, step, n_segments),
            "septal_arc": septal_arc, "direction_from": "clockwise flag",
            "origin": tuple(origin), "anterior": anterior,
            "inferior": inferior}


def _edges(anterior_angle, direction, step, n_segments):
    """Boundary angles, starting at the one that opens segment 1."""
    first = anterior_angle - direction * step
    return np.array([_wrap(first + direction * step * k)
                     for k in range(int(n_segments) + 1)])


def segment_index(theta, geometry):
    """1-based wedge index for each angle in ``theta``.

    The ``+ step`` is what puts the anterior ray on the boundary between
    segment 1 and segment 2 rather than at the start of segment 1.
    """
    step = geometry["step"]
    phi = np.mod(geometry["direction"]
                 * (np.asarray(theta, dtype=float) - geometry["anterior_angle"])
                 + step, 2 * np.pi)
    index = np.floor(phi / step).astype(int) + 1
    # an angle landing exactly on 2*pi after the modulo would overflow
    return np.clip(index, 1, int(round(2 * np.pi / step)))


# Short forms for a crowded axis. Truncating the names instead turns both
# anteroseptal and anterolateral into "antero".
WALL_ABBREVIATIONS: dict[str, str] = {
    "anterior": "ant", "anteroseptal": "ant-sept",
    "inferoseptal": "inf-sept", "inferior": "inf",
    "inferolateral": "inf-lat", "anterolateral": "ant-lat",
    "septal": "sept", "lateral": "lat",
}


def short_labels(level: str | None, n_segments: int) -> list[str]:
    """Names short enough for a bullseye axis, numbers when there are none."""
    out = []
    for name in aha_labels(level, n_segments):
        wall = name.split()[-1]
        out.append(WALL_ABBREVIATIONS.get(wall, wall))
    return out


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
