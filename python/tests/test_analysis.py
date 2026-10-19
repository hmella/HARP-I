import numpy as np
import pytest

from harpi import (aha_labels, bullseye_rings, pixel_strain_2d, segment_metrics,
                   segment_table, strain_by_segments)
from harpi.analysis import AHA_LEVELS


def test_aha_labels_match_the_model():
    assert len(aha_labels("basal", 6)) == 6
    assert len(aha_labels("apical", 4)) == 4
    assert aha_labels("mid", 6)[0] == "mid anterior"


def test_labels_fall_back_to_numbers_when_the_count_disagrees():
    # a label that does not match its wedge is worse than no label
    assert aha_labels("apical", 6) == [str(i) for i in range(1, 7)]
    assert aha_labels(None, 3) == ["1", "2", "3"]


def test_unknown_level_is_rejected():
    with pytest.raises(ValueError, match="unknown level"):
        aha_labels("septal", 6)


def test_every_level_is_labelled_consistently():
    for level, names in AHA_LEVELS.items():
        assert aha_labels(level, len(names)) == list(names)


def test_metrics_pick_the_most_negative_for_shortening():
    curves = np.array([[0.0, -0.1, -0.25, -0.2],
                       [0.0, -0.3, -0.1, -0.05]])
    m = segment_metrics(curves, component="CC")
    assert m["peak"] == pytest.approx([-0.25, -0.3])
    assert list(m["peak_frame"]) == [2, 1]
    assert m["final"] == pytest.approx([-0.2, -0.05])


def test_metrics_use_magnitude_for_principal_strain():
    curves = np.array([[0.0, 0.4, -0.1]])
    m = segment_metrics(curves, component="p1")
    assert m["peak"] == pytest.approx([0.4])


def test_metrics_report_time_when_times_are_given():
    curves = np.array([[0.0, -0.2, -0.1]])
    m = segment_metrics(curves, times=[0.0, 40.0, 80.0], component="CC")
    assert m["time_to_peak"] == pytest.approx([40.0])


def test_metrics_tolerate_an_all_nan_segment():
    curves = np.array([[np.nan, np.nan], [0.0, -0.1]])
    m = segment_metrics(curves, component="CC")
    assert np.isnan(m["peak"][0]) and m["peak_frame"][0] == -1
    assert m["peak"][1] == pytest.approx(-0.1)


def test_metrics_reject_a_times_length_mismatch():
    with pytest.raises(ValueError, match="times has"):
        segment_metrics(np.zeros((2, 4)), times=[0.0, 1.0])


def test_bullseye_rings_are_ordered_outside_in():
    rings = bullseye_rings([np.zeros(6), np.zeros(6)])
    assert rings[0]["radius"] > rings[1]["radius"]
    assert sum(r["height"] for r in rings) == pytest.approx(1.0)


def test_bullseye_accepts_a_single_ring():
    rings = bullseye_rings(np.zeros(4))
    assert len(rings) == 1
    assert len(rings[0]["wedges"]) == 4
    # wedges tile the circle exactly
    assert sum(w[1] for w in rings[0]["wedges"]) == pytest.approx(2 * np.pi)


def _ring_strain(nfr=4):
    rows = cols = 40
    yy, xx = np.mgrid[0:rows, 0:cols]
    r = np.hypot(xx - 20.0, yy - 20.0)
    mask = (r > 7) & (r < 15)
    u = np.zeros((rows, cols, 2, nfr))
    for k in range(nfr):
        s = -0.08 * k / max(nfr - 1, 1)
        u[:, :, 0, k] = s * (xx - 20.0)
        u[:, :, 1, k] = s * (yy - 20.0)
    return pixel_strain_2d(u, mask)


def test_segment_table_covers_every_segment():
    st = _ring_strain()
    seg = strain_by_segments(st, n_segments=6)
    header, rows = segment_table(seg, level="mid")
    assert len(rows) == 6
    assert rows[0][1] == "mid anterior"
    assert header[:2] == ["segment", "name"]
    # a peak column for each component present
    assert any(h.startswith("peak CC") for h in header)


def test_segment_table_numbers_when_the_level_is_unknown():
    st = _ring_strain()
    seg = strain_by_segments(st, n_segments=5)
    _, rows = segment_table(seg, level="mid")
    assert [r[1] for r in rows] == ["1", "2", "3", "4", "5"]


def test_segment_table_reports_shortening_on_a_contracting_ring():
    st = _ring_strain()
    seg = strain_by_segments(st, n_segments=6)
    _, rows = segment_table(seg, level="mid")
    peaks = [r[2] for r in rows]
    assert all(p < 0 for p in peaks), peaks


def test_the_outlier_cut_is_one_sided():
    """strain_by_segments drops strain ABOVE the threshold, not large magnitudes.

    Shortening is negative, so a tight threshold leaves a contracting segment
    untouched while removing implausible lengthening. Worth pinning down: a
    symmetric reading of the parameter would be wrong.
    """
    st = _ring_strain()
    loose = strain_by_segments(st, n_segments=6, outlier_threshold=100.0)
    tight = strain_by_segments(st, n_segments=6, outlier_threshold=0.01)
    # every value is negative, so nothing is above either threshold
    assert np.allclose(loose.CC, tight.CC, equal_nan=True)

    # now add lengthening that only the loose threshold keeps
    st.CC = st.CC.copy()
    st.CC[st.mask, -1] = 0.5
    loose = strain_by_segments(st, n_segments=6, outlier_threshold=100.0)
    tight = strain_by_segments(st, n_segments=6, outlier_threshold=0.1)
    assert np.nanmax(loose.CC[:, -1]) > 0.4
    assert np.all(np.isnan(tight.CC[:, -1]))
