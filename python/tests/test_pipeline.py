"""End to end behaviour, against MATLAB and against a synthetic phantom."""

import numpy as np
import pytest

import harpi
from conftest import masked_columns


@pytest.fixture(scope="module")
def motion(dataset):
    flt = harpi.HarpFilter2D(dataset.image, dataset.ke, (0.0, np.pi / 2),
                             kind="transmission")
    return harpi.estimate_motion_2d(
        flt(dataset.image), dataset.mask,
        kernel="multiquadric3", a=dataset.suggested_a, eta=15.0,
        connectivity=4, smooth_reference=True, temporal_fit=True,
        temporal_order=10, downsample=1,
    )


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


def test_motion_matches_matlab(reference, motion):
    for attr, key, tol in (("u", "u_raw", 1e-6), ("u_fitted", "u_fit", 1e-6)):
        got, want = getattr(motion, attr), reference[key]
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


def test_strain_end_to_end(reference, dataset, motion):
    """Whole pipeline in Python. The tolerance leaves room for the small
    difference in the fitted motion to grow through the derivatives."""
    st = harpi.pixel_strain_2d(motion.u_fitted, dataset.mask[:, :, 0])
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
