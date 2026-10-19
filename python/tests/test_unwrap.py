import numpy as np
import pytest

from harpi import phase_quality_2d, unwrap2, unwrap3
from harpi.unwrap import pathfollow


def _wrap(x):
    return np.mod(x + np.pi, 2 * np.pi) - np.pi


def test_linear_ramp_is_recovered():
    """A ramp spanning several cycles comes back with the wraps removed."""
    y, x = np.mgrid[0:40, 0:40]
    truth = 0.7 * x + 0.3 * y
    got, complete = unwrap2(_wrap(truth), np.ones((40, 40), bool))
    assert complete
    # Unwrapping recovers the field up to a global 2*pi offset.
    offset = np.round((got - truth).mean() / (2 * np.pi)) * 2 * np.pi
    assert np.abs(got - truth - offset).max() < 1e-9


def test_mask_is_respected():
    y, x = np.mgrid[0:30, 0:30]
    mask = (x - 15) ** 2 + (y - 15) ** 2 < 80
    got, complete = unwrap2(_wrap(0.5 * x), mask)
    assert complete
    assert np.isnan(got[~mask]).all()
    assert np.isfinite(got[mask]).all()


def test_disconnected_mask_reports_incomplete():
    y, x = np.mgrid[0:20, 0:20]
    mask = np.zeros((20, 20), bool)
    mask[2:6, 2:6] = True
    mask[14:18, 14:18] = True      # a second island, not reachable
    with pytest.warns(RuntimeWarning, match="could not be unwrapped"):
        got, complete = unwrap2(_wrap(0.4 * x), mask)
    assert not complete
    assert np.isnan(got).any()


def test_seeded_from_a_partial_solution():
    """A full resolution frame seeded from a downsampled solution."""
    y, x = np.mgrid[0:24, 0:24]
    truth = 0.6 * x - 0.2 * y
    mask = np.ones((24, 24), bool)
    partial = np.full((24, 24), np.nan)
    partial[::4, ::4] = truth[::4, ::4]
    got, complete = unwrap2(_wrap(truth), mask, unwrapped_image=partial)
    assert complete
    assert np.abs(got - truth).max() < 1e-9


def test_3d_ramp():
    z, y, x = np.mgrid[0:12, 0:12, 0:12]
    truth = 0.4 * x + 0.2 * y + 0.1 * z
    got, complete = unwrap3(_wrap(truth), np.ones((12, 12, 12), bool))
    assert complete
    offset = np.round((got - truth).mean() / (2 * np.pi)) * 2 * np.pi
    assert np.abs(got - truth - offset).max() < 1e-9


@pytest.mark.skipif(not pathfollow.HAVE_NUMBA, reason="numba not installed")
def test_numba_and_python_kernels_agree():
    """The compiled and the plain NumPy kernel give the same result."""
    rng = np.random.default_rng(0)
    y, x = np.mgrid[0:28, 0:28]
    phase = _wrap(0.5 * x + 0.3 * y + rng.normal(0, 0.05, (28, 28)))
    mask = (x - 14) ** 2 + (y - 14) ** 2 < 100
    qual = phase_quality_2d(phase, mask)

    fast, _ = unwrap2(phase, mask, quality=qual)

    saved = pathfollow.HAVE_NUMBA
    try:
        pathfollow.HAVE_NUMBA = False
        slow, _ = unwrap2(phase, mask, quality=qual)
    finally:
        pathfollow.HAVE_NUMBA = saved

    assert np.array_equal(np.isnan(fast), np.isnan(slow))
    assert np.array_equal(fast[~np.isnan(fast)], slow[~np.isnan(slow)])


def test_quality_matches_matlab(reference, dataset):
    """Phase quality against the MATLAB reference, frame by frame."""
    Xpha = reference["Xpha"]
    qref = reference["qual_x"]
    for i in range(Xpha.shape[2]):
        q = phase_quality_2d(Xpha[:, :, i], dataset.mask[:, :, i], (1.0, 1.0), 4)
        assert np.array_equal(np.isfinite(q), np.isfinite(qref[:, :, i]))
        fin = np.isfinite(q)
        assert np.abs(q[fin] - qref[:, :, i][fin]).max() < 1e-12


def test_unwrapping_matches_matlab_exactly(reference, dataset):
    """Wrap counts have to agree: a 2*pi shift is a different answer."""
    for key_p, key_u in (("Xpha", "uXpha"), ("Ypha", "uYpha")):
        pha, ref = reference[key_p], reference[key_u]
        for i in range(pha.shape[2]):
            got, _ = unwrap2(pha[:, :, i], dataset.mask[:, :, i], (1.0, 1.0),
                             "auto", 4, warn=False)
            assert np.array_equal(np.isfinite(got), np.isfinite(ref[:, :, i]))
            fin = np.isfinite(ref[:, :, i])
            d = got[fin] - ref[:, :, i][fin]
            assert np.abs(np.round(d / (2 * np.pi))).max() == 0, "wrap count differs"
            assert np.abs(d).max() < 1e-9
