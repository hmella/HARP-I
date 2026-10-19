import numpy as np
import pytest

from harpi import KERNELS, get_kernel
from harpi.interp import distance_matrix, rbf_interpolate, solve_tikhonov


def test_matlab_aliases_resolve():
    assert get_kernel("Multiquadric3").name == "multiquadric3"
    assert get_kernel("WendlandC2").name == "wendlandc2"
    assert get_kernel("imqx").name == "imq"


def test_unknown_kernel_lists_options():
    with pytest.raises(ValueError, match="multiquadric3"):
        get_kernel("banana")


def test_multiquadric3_matches_formula():
    r = np.linspace(0, 5, 11)
    assert np.allclose(get_kernel("multiquadric3")(r, 2.0), (r**2 + 4) ** 1.5)


def test_thin_plate_is_finite_at_zero():
    r = np.array([0.0, 0.5, 1.0, 2.0])
    v = get_kernel("thinplate")(r, np.nan)
    assert np.isfinite(v).all()
    assert v[0] == 0.0


@pytest.mark.parametrize("name", [k for k, v in KERNELS.items() if v.compact])
def test_compact_kernels_vanish_outside_support(name):
    k = get_kernel(name)
    r = np.array([1.5, 2.0, 10.0])
    assert np.allclose(k(r, 1.0), 0.0)


def test_distance_matrix_matches_bruteforce():
    rng = np.random.default_rng(0)
    a, b = rng.normal(size=(2, 7)), rng.normal(size=(2, 5))
    d = distance_matrix(a, b)
    for i in range(7):
        for j in range(5):
            assert abs(d[i, j] - np.hypot(*(a[:, i] - b[:, j]))) < 1e-12


def test_tikhonov_solves_the_regularised_system():
    rng = np.random.default_rng(1)
    n = 20
    pts = rng.normal(size=(2, n))
    psi = get_kernel("multiquadric3")(distance_matrix(pts), 3.0)
    g = rng.normal(size=(n, 2))
    eta = 0.5
    w = solve_tikhonov(psi, g, eta)
    assert np.allclose((psi + eta * np.eye(n)) @ w, g, atol=1e-8)


def test_interpolation_reproduces_a_linear_field():
    """With a small eta the interpolant passes through its data."""
    rng = np.random.default_rng(2)
    c = rng.uniform(-5, 5, size=(2, 40))
    vals = np.stack([2 * c[0] - c[1], c[0] + 3 * c[1]], axis=1)
    got = rbf_interpolate(c, vals, c, "multiquadric3", 5.0, 1e-10)
    assert np.abs(got - vals).max() < 1e-4
