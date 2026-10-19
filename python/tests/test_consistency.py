"""The rules that choose the 2*pi offset of each frame."""

import numpy as np
import pytest

import harpi

TAG = 8.0                  # tag period, px
KE = 2 * np.pi / TAG


def _ring(shape=(96, 96), centre=48.0, r_inner=16, r_outer=30):
    """Virtual coordinates of a ring, contracting and moving over four frames.

    The coordinates of a pixel are those of the material imaged there, so
    frame k holds the reference coordinates sampled at other points, plus an
    unknown whole number of periods from unwrapping.
    """
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]].astype(float)
    scale = (1.0, 0.95, 0.9, 0.9)
    shift = ((0.0, 0.0), (3.5, -2.0), (13.0, 6.0), (13.0, 6.0))
    phx = np.full(shape + (4,), np.nan)
    phy = np.full(shape + (4,), np.nan)
    for k, (s, (dx, dy)) in enumerate(zip(scale, shift)):
        X = centre + (xx - centre - dx) / s
        Y = centre + (yy - centre - dy) / s
        r = np.hypot(X - centre, Y - centre)
        inside = (r > r_inner) & (r < r_outer)
        if k == 3:
            # 40 % of the ring is missing from the mask of the last frame
            angle = np.mod(np.arctan2(Y - centre, X - centre), 2 * np.pi)
            inside &= angle > 0.4 * 2 * np.pi
        phx[..., k] = np.where(inside, KE * X, np.nan)
        phy[..., k] = np.where(inside, KE * Y, np.nan)
    return phx, phy


# Periods added by unwrapping to frames 1, 2 and 3.
OFFSETS = np.array([[0, 3, -1, 2], [0, -2, 4, -3]])


@pytest.fixture(scope="module")
def ring():
    phx, phy = _ring()
    phx = phx + 2 * np.pi * OFFSETS[0]
    phy = phy + 2 * np.pi * OFFSETS[1]
    return phx, phy


def test_the_nearest_rule_recovers_the_offsets(ring):
    jumps = harpi.nearest_coordinate_consistency(list(ring))
    assert np.array_equal(jumps, -2 * np.pi * OFFSETS)


def test_the_published_rule_agrees_on_the_ring(ring):
    """So the expected offsets do not come from the new rule alone."""
    jumps = harpi.temporal_phase_consistency(list(ring))
    assert np.array_equal(jumps, -2 * np.pi * OFFSETS)


def test_the_search_info_adds_up(ring):
    jumps, info = harpi.nearest_coordinate_consistency(list(ring), return_info=True)
    chosen = np.round(jumps / (2 * np.pi))
    assert np.array_equal(chosen, info["delta"] + info["offset"])
    assert info["njumps"] == 2
    assert np.abs(info["offset"]).max() <= 2


def test_the_reference_frame_can_be_any_frame(ring):
    jumps = harpi.nearest_coordinate_consistency(list(ring), reference=2)
    want = -2 * np.pi * (OFFSETS - OFFSETS[:, 2:3])
    assert np.array_equal(jumps, want)


def test_an_empty_frame_keeps_no_jump(ring):
    phx, phy = (p.copy() for p in ring)
    phx[..., 1] = np.nan
    jumps = harpi.nearest_coordinate_consistency([phx, phy])
    assert np.array_equal(jumps[:, 1], [0.0, 0.0])
    assert np.array_equal(jumps[:, 2:], -2 * np.pi * OFFSETS[:, 2:])


def test_an_empty_reference_is_refused(ring):
    phx, phy = (p.copy() for p in ring)
    phx[..., 0] = np.nan
    with pytest.raises(ValueError, match="reference"):
        harpi.nearest_coordinate_consistency([phx, phy])


def test_the_nearest_rule_on_the_matlab_phases(reference):
    """Given the unwrapped phases of the MATLAB toolbox.

    The two rules agree in y. In x the published rule puts frames 6 to 9, 11
    and 17 one period away from the nearest-coordinate rule.
    """
    phases = [reference["uXpha"], reference["uYpha"]]
    jumps = harpi.nearest_coordinate_consistency(phases, downsample=2)
    # whole periods; 2*pi*k / (2*pi) need not give k exactly, so compare radians
    assert np.array_equal(jumps, 2 * np.pi * np.round(jumps / (2 * np.pi)))
    assert np.array_equal(jumps[:, 0], [0.0, 0.0])

    published = reference["phase_jumps_n8"]
    differ = np.flatnonzero(np.any(jumps != published, axis=0))
    assert differ.tolist() == [5, 6, 7, 8, 10, 16]
    assert np.array_equal(jumps[1], published[1])
    assert np.allclose(np.abs(jumps[0, differ] - published[0, differ]), 2 * np.pi)
