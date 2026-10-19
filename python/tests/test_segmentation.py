import numpy as np

from harpi.segmentation import Contour, contours_to_mask, resample_contour


def ring(cx, cy, r, n=12):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([cx + r * np.cos(th), cy + r * np.sin(th)])


def test_resampling_is_dense_and_closed():
    pts = ring(40, 40, 20)
    s = resample_contour(pts, closed=True, resolution=0.5)
    assert len(s) > len(pts) * 10
    assert np.allclose(s[0], s[-1], atol=1e-9)
    step = np.hypot(*np.diff(s, axis=0).T)
    assert step.max() < 1.0


def test_annulus_area_and_topology():
    mask = contours_to_mask([Contour(ring(40, 40, 20)), Contour(ring(40, 40, 10))],
                            (80, 80), 0.5, origin=1)
    ideal = np.pi * (20**2 - 10**2)
    assert abs(mask.sum() / ideal - 1) < 0.05
    assert not mask[39, 39]      # cavity (1-based (40,40))
    assert mask[54, 39]          # mid-wall


def test_grid_origin_convention():
    """origin=1 gives the 1-based grid used by the MATLAB toolbox."""
    c = Contour(ring(3, 3, 2.5, 16))
    one = contours_to_mask([c], (8, 8), 0.25, origin=1)
    zero = contours_to_mask([c], (8, 8), 0.25, origin=0)
    assert one.sum() > 0 and zero.sum() > 0
    assert not np.array_equal(one, zero), "the two conventions differ"


def test_more_than_two_contours():
    mask = contours_to_mask(
        [Contour(ring(40, 40, 20)), Contour(ring(34, 40, 5)), Contour(ring(46, 40, 5))],
        (80, 80), 0.5, origin=1,
    )
    assert not mask[39, 33]
    assert not mask[39, 45]
    assert mask[21, 39]


def test_contour_outside_outer_boundary_does_not_leak():
    inside = contours_to_mask([Contour(ring(40, 40, 20)), Contour(ring(40, 40, 10))],
                              (80, 80), 0.5, origin=1)
    with_stray = contours_to_mask(
        [Contour(ring(40, 40, 20)), Contour(ring(40, 40, 10)), Contour(ring(70, 70, 3))],
        (80, 80), 0.5, origin=1,
    )
    assert with_stray.sum() == inside.sum()
