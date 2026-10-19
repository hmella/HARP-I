"""The 3-D path, which had no coverage at all.

The three defects pinned here were found by reading `HARPI_3D.m` against the
port: the mask shape was never validated, the reference-phase smoothing was
applied to the wrong array, and `strain_by_segments` accepted a volume and
returned plausible but wrong numbers.
"""

from __future__ import annotations

import numpy as np
import pytest

from harpi import (Dataset3D, HarpFilter3D, estimate_motion_3d, pixel_strain_3d,
                   smooth_reference_phase, strain_by_segments)


def _volume(rows=24, cols=24, slices=6, nfr=3, spacing=6.0, contraction=0.08):
    """A tagged shell contracting in plane, with an analytic displacement."""
    ke = 2 * np.pi / spacing
    cy, cx = (rows - 1) / 2, (cols - 1) / 2

    Y = (np.arange(rows) - cy)[:, None, None]
    X = (np.arange(cols) - cx)[None, :, None]
    Z = (np.arange(slices) - (slices - 1) / 2)[None, None, :]
    Y, X, Z = np.broadcast_arrays(Y, X, Z)

    r = np.hypot(X, Y)
    mask = (r > 4) & (r < 9)                      # (rows, cols, slices)

    image = np.zeros((rows, cols, slices, 3, nfr), dtype=complex)
    for k in range(nfr):
        s = -contraction * k / max(nfr - 1, 1)
        for d, coord in enumerate((X / (1 + s), Y / (1 + s), Z)):
            image[:, :, :, d, k] = np.exp(1j * ke * coord)
    return image, mask, spacing


# ------------------------------------------------------------ mask validation

def test_a_mask_of_the_wrong_shape_is_refused():
    image, _, spacing = _volume()
    bad = np.ones((image.shape[0] + 1, image.shape[1], image.shape[2]), bool)
    with pytest.raises(ValueError, match="disagree on spatial shape"):
        estimate_motion_3d(image, bad, a=150 / spacing)


def test_a_three_dimensional_mask_is_broadcast_over_frames():
    """One mask for the whole sequence is repeated, as in the 2-D path."""
    image, mask, spacing = _volume()
    assert mask.ndim == 3
    res = estimate_motion_3d(image, mask, a=150 / spacing, downsample=2)
    assert res.u.shape == image.shape[:3] + (3, image.shape[-1])


def test_a_per_frame_mask_is_accepted_too():
    image, mask, spacing = _volume()
    per_frame = np.repeat(mask[..., None], image.shape[-1], axis=3)
    res = estimate_motion_3d(image, per_frame, a=150 / spacing, downsample=2)
    assert res.u.shape == image.shape[:3] + (3, image.shape[-1])


def test_a_bad_image_shape_is_refused():
    with pytest.raises(ValueError, match="rows, cols, slices, 3, frames"):
        estimate_motion_3d(np.zeros((4, 4, 2, 3)), np.ones((4, 4, 2), bool), a=10.0)


# ------------------------------------------------------- the smoothing guard

def test_the_smoother_refuses_a_system_it_cannot_hold():
    """A full-resolution volume would need tens of gibibytes."""
    big = [np.zeros((60, 60, 60)) for _ in range(3)]
    with pytest.raises(ValueError, match="GiB"):
        smooth_reference_phase(big, "multiquadric3", 10.0, 15.0)


def test_the_guard_names_the_point_count():
    big = [np.zeros((40, 40, 40)) for _ in range(3)]
    with pytest.raises(ValueError, match="64000 by 64000"):
        smooth_reference_phase(big, "multiquadric3", 10.0, 15.0)


def test_the_limit_is_honoured_in_both_directions():
    small = [np.zeros((6, 6, 6)) for _ in range(3)]        # 216 points
    with pytest.raises(ValueError, match="216 by 216"):
        smooth_reference_phase(small, "multiquadric3", 10.0, 15.0, max_points=100)
    out = smooth_reference_phase(small, "multiquadric3", 10.0, 15.0,
                                 max_points=1000)
    assert len(out) == 3 and np.isfinite(out[0]).all()


def test_the_limit_can_be_switched_off():
    small = [np.zeros((5, 5, 5)) for _ in range(3)]
    out = smooth_reference_phase(small, "multiquadric3", 10.0, 15.0,
                                 max_points=None)
    assert len(out) == 3


def test_two_dimensional_smoothing_is_unaffected_by_the_guard():
    """The 2-D reference is a few hundred points, far under the limit."""
    phases = [np.full((40, 40), np.nan) for _ in range(2)]
    for p in phases:
        p[10:25, 10:25] = 0.0
    out = smooth_reference_phase(phases, "multiquadric3", 10.0, 15.0)
    assert len(out) == 2
    assert np.isfinite(out[0][10:25, 10:25]).all()


def test_smoothing_runs_on_the_downsampled_volume():
    """With downsample 2 the per-frame systems are small enough to solve."""
    image, mask, spacing = _volume(rows=20, cols=20, slices=4, nfr=2)
    res = estimate_motion_3d(image, mask, a=150 / spacing, downsample=2,
                             smooth_reference=True)
    assert np.isfinite(res.u).any()


def test_smoothing_changes_the_result():
    image, mask, spacing = _volume(rows=20, cols=20, slices=4, nfr=2)
    kw = dict(a=150 / spacing, downsample=2)
    plain = estimate_motion_3d(image, mask, **kw)
    smooth = estimate_motion_3d(image, mask, smooth_reference=True, **kw)
    both = np.isfinite(plain.u) & np.isfinite(smooth.u)
    assert both.any()
    assert not np.allclose(plain.u[both], smooth.u[both])


# ------------------------------------------------------- segments refuse 3-D

def _strain_volume():
    n, nfr = 12, 3
    u = np.zeros((n, n, n, 3, nfr))
    return pixel_strain_3d(u, np.ones((n, n, n), bool))


def test_segments_refuse_a_cubic_volume():
    """This used to pass silently and return plausible numbers."""
    with pytest.raises(ValueError, match="2-D slice"):
        strain_by_segments(_strain_volume())


def test_the_refusal_names_the_shape_and_the_remedy():
    with pytest.raises(ValueError) as e:
        strain_by_segments(_strain_volume())
    msg = str(e.value)
    assert "(12, 12, 12)" in msg
    assert "one slice at a time" in msg


def test_a_slice_of_a_volume_is_accepted():
    """The documented way round: take one slice."""
    st = _strain_volume()
    sl = type(st)(mask=st.mask[:, :, 6],
                  CC=st.CC[:, :, 6, :], RR=st.RR[:, :, 6, :])
    seg = strain_by_segments(sl, n_segments=6)
    assert seg.segments.shape == (12, 12)
    assert seg.CC.shape == (6, 3)


def test_two_dimensional_segments_still_work():
    rows = cols = 32
    yy, xx = np.mgrid[0:rows, 0:cols]
    r = np.hypot(xx - 16.0, yy - 16.0)
    ring = (r > 5) & (r < 12)
    nfr = 4
    u = np.zeros((rows, cols, 2, nfr))
    for k in range(nfr):
        s = -0.1 * k / (nfr - 1)
        u[:, :, 0, k] = s * (xx - 16.0)
        u[:, :, 1, k] = s * (yy - 16.0)
    from harpi import pixel_strain_2d

    seg = strain_by_segments(pixel_strain_2d(u, ring), n_segments=6)
    assert seg.CC.shape == (6, nfr)


# ------------------------------------------------- the interface refuses 3-D

pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication  # noqa: E402

from harpi import Dataset2D  # noqa: E402
from harpi.gui.session import Session, Stage  # noqa: E402
from harpi.gui.tabs_data import DataTab  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def tab(qapp):
    return DataTab(Session())


def _ds3(rows=8, cols=8, slices=4, nfr=3):
    return Dataset3D(image=np.zeros((rows, cols, slices, 3, nfr)),
                     mask=np.ones((rows, cols, slices, nfr)),
                     voxel_size=[2.0, 2.0, 8.0], tag_spacing=8.0)


def _ds2(rows=8, cols=8, nfr=3):
    return Dataset2D(image=np.zeros((rows, cols, 2, nfr)),
                     mask=np.ones((rows, cols, nfr)),
                     pixel_size=[2.0, 2.0], tag_spacing=8.0)


def test_a_volume_is_refused(tab):
    assert tab._adopt(_ds3(), "volume") is False


def test_a_refused_volume_does_not_become_the_session_dataset(tab):
    tab._adopt(_ds3(), "volume")
    assert tab.session.dataset is None
    assert not tab.session.has(Stage.DATA)


def test_the_refusal_says_what_to_do_instead(tab):
    tab._adopt(_ds3(), "volume")
    text = tab.banner.text()
    assert "2-D series" in text
    assert "estimate_motion_3d" in text


def test_a_two_dimensional_series_is_still_accepted(tab):
    assert tab._adopt(_ds2(), "slice") is True
    assert tab.session.dataset is not None
    assert tab.session.has(Stage.DATA)


def test_a_refusal_leaves_an_earlier_dataset_alone(tab):
    tab._adopt(_ds2(), "slice")
    good = tab.session.dataset
    tab._adopt(_ds3(), "volume")
    assert tab.session.dataset is good


# -------------------------------------------------- the jump rule in 3-D

def test_both_jump_rules_are_accepted_in_3d():
    """The nearest rule is dimension generic, so 3-D is not a special case."""
    image, mask, spacing = _volume(rows=20, cols=20, slices=4, nfr=3)
    kw = dict(a=150 / spacing, downsample=2)
    for rule in ("nearest", "published"):
        res = estimate_motion_3d(image, mask, phase_jump_rule=rule, **kw)
        assert res.meta["phase_jump_rule"] == rule


def test_the_nearest_rule_is_the_default_in_3d():
    """2-D and 3-D must not disagree about what the default means."""
    image, mask, spacing = _volume(rows=20, cols=20, slices=4, nfr=3)
    res = estimate_motion_3d(image, mask, a=150 / spacing, downsample=2)
    assert res.meta["phase_jump_rule"] == "nearest"


def test_the_rules_agree_when_no_frame_is_a_period_out():
    """Small motion leaves nothing to disagree about, which is the sanity check
    available without a 3-D reference to compare against."""
    image, mask, spacing = _volume(rows=20, cols=20, slices=4, nfr=3,
                                   contraction=0.02)
    kw = dict(a=150 / spacing, downsample=2)
    a = estimate_motion_3d(image, mask, phase_jump_rule="nearest", **kw)
    b = estimate_motion_3d(image, mask, phase_jump_rule="published", **kw)
    assert np.array_equal(a.phase_jumps, b.phase_jumps)
    both = np.isfinite(a.u) & np.isfinite(b.u)
    assert np.allclose(a.u[both], b.u[both])


def test_an_unknown_rule_is_still_refused_in_3d():
    image, mask, spacing = _volume(rows=12, cols=12, slices=3, nfr=2)
    with pytest.raises(ValueError, match="phase_jump_rule"):
        estimate_motion_3d(image, mask, a=150 / spacing, phase_jump_rule="centroid")


def test_the_search_budget_is_ignored_by_the_nearest_rule_in_3d():
    """njumps is accepted and has no effect, as in 2-D."""
    image, mask, spacing = _volume(rows=20, cols=20, slices=4, nfr=3)
    kw = dict(a=150 / spacing, downsample=2, phase_jump_rule="nearest")
    a = estimate_motion_3d(image, mask, njumps=1, **kw)
    b = estimate_motion_3d(image, mask, njumps=12, **kw)
    assert np.array_equal(a.phase_jumps, b.phase_jumps)
