"""End to end behaviour, against MATLAB and against a synthetic phantom."""

import numpy as np
import pytest

import harpi
from conftest import masked_columns


def _solve(dataset, **kw):
    flt = harpi.HarpFilter2D(dataset.image, dataset.ke, (0.0, np.pi / 2),
                             kind="transmission")
    return harpi.estimate_motion_2d(
        flt(dataset.image), dataset.mask,
        kernel="multiquadric3", a=dataset.suggested_a, eta=15.0,
        connectivity=4, smooth_reference=True, temporal_fit=True,
        temporal_order=10, downsample=1, **kw,
    )


@pytest.fixture(scope="module")
def motion(dataset):
    """The default rule, whichever that is."""
    return _solve(dataset)


@pytest.fixture(scope="module")
def published(dataset):
    """The MATLAB rule, named rather than relied on as the default.

    The equivalence tests take this. Leaving them on the default would move the
    bit-exact guard onto whatever rule happens to be default, where it would
    fail, and loosening a tolerance to make it pass would lose the equivalence
    without anyone noticing.
    """
    return _solve(dataset, phase_jump_rule="published")


@pytest.fixture(scope="module")
def nearest(dataset):
    """The nearest-coordinate rule, named for the tests that contrast them."""
    return _solve(dataset, phase_jump_rule="nearest")


def test_phase_jumps_match_matlab_exactly(reference):
    jumps = harpi.temporal_phase_consistency(
        [reference["uXpha"], reference["uYpha"]], reference=0, njumps=8, downsample=2
    )
    assert np.array_equal(jumps, reference["phase_jumps_n8"])


def test_temporal_fit_matches_matlab(reference, dataset):
    got = harpi.fit_temporal(reference["u_raw"], dataset.mask[:, :, 0], order=10)
    want = reference["u_fit"]
    assert np.array_equal(np.isfinite(got), np.isfinite(want))
    fin = np.isfinite(got)
    assert np.abs(got[fin] - want[fin]).max() < 1e-6


def test_motion_matches_matlab(reference, published):
    for attr, key, tol in (("u", "u_raw", 1e-6), ("u_fitted", "u_fit", 1e-6)):
        got, want = getattr(published, attr), reference[key]
        assert np.array_equal(np.isfinite(got), np.isfinite(want))
        fin = np.isfinite(got)
        assert np.abs(got[fin] - want[fin]).max() < tol, attr


def test_strain_matches_matlab(reference, dataset):
    """Given the MATLAB fitted motion, so only the strain stage is tested."""
    st = harpi.pixel_strain_2d(reference["u_fit"], dataset.mask[:, :, 0])
    ref_mask = reference["st/maskimage"].astype(bool)
    assert np.array_equal(st.mask, ref_mask)
    for comp in ("XX", "YY", "XY", "RR", "CC", "p1", "p2"):
        got = masked_columns(getattr(st, comp), st.mask)
        want = reference[f"st/{comp}"]
        assert np.abs(got - want).max() < 1e-12, comp


def test_strain_end_to_end(reference, dataset, published):
    """Whole pipeline in Python. The tolerance leaves room for the small
    difference in the fitted motion to grow through the derivatives."""
    st = harpi.pixel_strain_2d(published.u_fitted, dataset.mask[:, :, 0])
    for comp in ("RR", "CC"):
        got = masked_columns(getattr(st, comp), st.mask)
        want = reference[f"st/{comp}"]
        assert np.abs(got - want).max() < 1e-7, comp

    gcs = st.global_curve("CC")
    want_gcs = reference["GCS"].ravel()
    assert np.abs(gcs - want_gcs).max() < 1e-8


def test_unwrapping_completes_on_real_data(motion):
    assert motion.unwrap_complete.all()


# ------------------------------------------------------------ phantom tests

def test_phantom_displacement_accuracy(phantom):
    """Against the prescribed deformation."""
    ppw = phantom["tag_spacing"] / phantom["pixel_size"][0]
    ke = 2 * np.pi / phantom["tag_spacing"] * phantom["pixel_size"]
    flt = harpi.HarpFilter2D(phantom["image"], ke, (0.0, np.pi / 2), kind="transmission")
    res = harpi.estimate_motion_2d(flt(phantom["image"]), phantom["mask"],
                                   a=150 / ppw, eta=15.0)

    m = phantom["mask"][:, :, 0]
    errs = []
    for k in range(phantom["image"].shape[-1]):
        sel = m & np.isfinite(res.u[:, :, 0, k])
        e = np.hypot(res.u[:, :, 0, k] - phantom["u_true"][:, :, 0, k],
                     res.u[:, :, 1, k] - phantom["u_true"][:, :, 1, k])[sel]
        errs.append(np.sqrt((e ** 2).mean()))
    assert np.mean(errs) < 0.15, f"RMS displacement error {np.mean(errs):.3f} px"


def test_phantom_reference_frame_is_zero(phantom):
    ppw = phantom["tag_spacing"] / phantom["pixel_size"][0]
    ke = 2 * np.pi / phantom["tag_spacing"] * phantom["pixel_size"]
    flt = harpi.HarpFilter2D(phantom["image"], ke, (0.0, np.pi / 2), kind="transmission")
    res = harpi.estimate_motion_2d(flt(phantom["image"]), phantom["mask"],
                                   a=150 / ppw, eta=15.0)
    u0 = res.u[..., 0]
    assert np.nanmax(np.abs(u0)) == 0.0


def test_phantom_contraction_shows_negative_circumferential_strain(phantom):
    ppw = phantom["tag_spacing"] / phantom["pixel_size"][0]
    ke = 2 * np.pi / phantom["tag_spacing"] * phantom["pixel_size"]
    flt = harpi.HarpFilter2D(phantom["image"], ke, (0.0, np.pi / 2), kind="transmission")
    res = harpi.estimate_motion_2d(flt(phantom["image"]), phantom["mask"],
                                   a=150 / ppw, eta=15.0)
    st = harpi.pixel_strain_2d(res.u, phantom["mask"][:, :, 0],
                               origin=tuple(phantom["centre"] + 1))
    gcs = st.global_curve("CC")
    # A contracting ring shortens circumferentially at peak systole.
    assert gcs[len(gcs) // 2] < -0.05
    assert abs(gcs[0]) < 1e-9


def test_jump_search_records_its_outcome(dataset):
    """The search info lets a caller see whether njumps is binding."""
    flt = harpi.HarpFilter2D(dataset.image, dataset.ke, (0.0, np.pi / 2),
                             kind="transmission")
    Ih = flt(dataset.image)
    res = harpi.estimate_motion_2d(Ih, dataset.mask, a=dataset.suggested_a,
                                   eta=15.0, njumps=8,
                                   phase_jump_rule="published")
    info = res.meta["jump_search"]
    assert info["njumps"] == 8
    assert info["offset"].shape == res.phase_jumps.shape
    # the chosen jump is the coarse estimate plus the offset
    chosen = np.round(res.phase_jumps / (2 * np.pi))
    assert np.allclose(chosen[:, 1:], (info["delta"] + info["offset"])[:, 1:])


def test_widening_the_jump_search_is_safe(dataset):
    """A wider search may change nothing, but must never make things worse."""
    flt = harpi.HarpFilter2D(dataset.image, dataset.ke, (0.0, np.pi / 2),
                             kind="transmission")
    Ih = flt(dataset.image)
    kw = dict(a=dataset.suggested_a, eta=15.0, phase_jump_rule="published")
    narrow = harpi.estimate_motion_2d(Ih, dataset.mask, njumps=1, **kw)
    wide = harpi.estimate_motion_2d(Ih, dataset.mask, njumps=12, **kw)
    used = int(np.abs(narrow.meta["jump_search"]["offset"]).max())
    if used < 1:
        assert np.array_equal(narrow.phase_jumps, wide.phase_jumps)


# ------------------------------------------------------ phase jump rules

def test_the_nearest_rule_is_the_default(motion):
    """It keeps the motion continuous where the published rule does not."""
    assert motion.meta["phase_jump_rule"] == "nearest"


def test_the_published_rule_is_still_reachable(published):
    """Kept so the paper can be reproduced, and so the MATLAB guard has a rule
    to guard."""
    assert published.meta["phase_jump_rule"] == "published"


def test_the_nearest_rule_is_recorded(nearest):
    assert nearest.meta["phase_jump_rule"] == "nearest"


def test_the_rules_differ_only_in_the_jumps(published, nearest):
    """Frames given the same jumps are solved identically, the others differ
    by whole periods."""
    same = np.all(published.phase_jumps == nearest.phase_jumps, axis=0)
    assert same[0]
    assert np.array_equal(published.u[..., same], nearest.u[..., same],
                          equal_nan=True)
    periods = (nearest.phase_jumps - published.phase_jumps) / (2 * np.pi)
    assert np.allclose(periods, np.round(periods), rtol=0, atol=1e-9)


def _largest_step(res):
    centre = np.nanmedian(res.u, axis=(0, 1))           # (2, frames)
    return np.abs(np.diff(centre, axis=1)).max()


def test_the_nearest_rule_keeps_the_motion_continuous(dataset, nearest, published):
    """No frame moves by a tag period relative to the frame before.

    The contrast is the point: on this dataset the published rule puts frames 6
    to 9, 11 and 17 one period off in x, so it fails the same check. That is
    why nearest is the default.
    """
    half = 0.5 * dataset.pixels_per_wavelength
    assert _largest_step(nearest) < half
    assert _largest_step(published) > half


def test_the_default_is_the_continuous_one(dataset, motion):
    """Whatever the default is, it must not leave a period of discontinuity."""
    assert _largest_step(motion) < 0.5 * dataset.pixels_per_wavelength


def test_an_unknown_jump_rule_is_refused():
    with pytest.raises(ValueError, match="phase_jump_rule"):
        harpi.estimate_motion_2d(np.zeros((8, 8, 2, 3), complex), np.ones((8, 8)),
                                 a=1.0, phase_jump_rule="centroid")


def test_an_unknown_jump_rule_is_refused_in_3d():
    with pytest.raises(ValueError, match="phase_jump_rule"):
        harpi.estimate_motion_3d(np.zeros((8, 8, 8, 3, 3), complex),
                                 np.ones((8, 8, 8)), a=1.0,
                                 phase_jump_rule="centroid")
