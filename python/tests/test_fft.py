import numpy as np
import pytest

from harpi import itok, ktoi


@pytest.mark.parametrize("shape", [(8, 8), (7, 9), (5, 6, 7), (4, 5, 6, 3)])
@pytest.mark.parametrize("axes", [0, (0, 1), (0, 1, 2), (1, 2)])
def test_round_trip(shape, axes):
    ax = (axes,) if isinstance(axes, int) else axes
    if max(ax) >= len(shape):
        pytest.skip("axis out of range for this shape")
    rng = np.random.default_rng(0)
    x = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    assert np.allclose(ktoi(itok(x, axes), axes), x, atol=1e-10)


def test_matches_explicit_transform():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(6, 7, 3)) + 1j * rng.normal(size=(6, 7, 3))
    want = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(x, axes=1), axis=1), axes=1)
    assert np.allclose(itok(x, 1), want, atol=1e-12)


def test_odd_sizes():
    """fftshift and ifftshift are different permutations at odd lengths."""
    rng = np.random.default_rng(2)
    x = rng.normal(size=(5, 7, 9)) + 1j * rng.normal(size=(5, 7, 9))
    assert np.allclose(ktoi(itok(x, (0, 1, 2)), (0, 1, 2)), x, atol=1e-10)


def test_matches_matlab(reference, dataset):
    from harpi import HarpFilter2D

    flt = HarpFilter2D(dataset.image, dataset.ke, (0.0, np.pi / 2), kind="transmission")
    got = flt(dataset.image)
    want = reference["Ih"]
    assert got.shape == want.shape
    assert np.abs(got - want).max() / np.abs(want).max() < 1e-12
