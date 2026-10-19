"""Moving the harmonic centre of the bandpass."""

import numpy as np
import pytest

from harpi import HarpFilter2D, centre_to_frequency, frequency_to_centre


@pytest.mark.parametrize("centre", [[80, 120], [120, 80], [80, 80], [95, 110.5]])
def test_centre_and_frequency_round_trip(centre):
    shape = (160, 160)
    back = frequency_to_centre(centre_to_frequency(centre, shape), shape)
    assert np.allclose(back, centre)


def test_centre_matches_the_nominal_frequency(dataset):
    """The detected peak of the demo data is the nominal tag frequency."""
    shape = dataset.image.shape[:2]
    omega = centre_to_frequency([80, 120], shape)
    assert np.hypot(*omega) == pytest.approx(dataset.ke[0], rel=1e-9)
    assert np.arctan2(omega[1], omega[0]) == pytest.approx(0.0, abs=1e-9)


def test_non_square_images():
    shape = (64, 128)
    omega = centre_to_frequency([48, 96], shape)
    assert omega[0] == pytest.approx(2 * np.pi * 0.25)
    assert omega[1] == pytest.approx(2 * np.pi * 0.25)
    assert np.allclose(frequency_to_centre(omega, shape), [48, 96])


@pytest.mark.parametrize("kind", ["transmission", "butterworth"])
def test_both_kinds_follow_the_centre(dataset, kind):
    """Moving the centre has to change the filtered image for either kind."""
    shape = dataset.image.shape[:2]

    def build(centres):
        omegas = [centre_to_frequency(c, shape) for c in centres]
        mags = [float(np.hypot(*w)) for w in omegas]
        direction = tuple(float(np.arctan2(w[1], w[0])) for w in omegas)
        return HarpFilter2D(dataset.image, mags, direction, kind=kind,
                            cutoff=17, order=10,
                            centre=np.round(centres).astype(int))

    a = build([[80, 120], [120, 80]])
    b = build([[80, 126], [120, 80]])

    assert not np.allclose(a.kspace_filter[:, :, 0], b.kspace_filter[:, :, 0]), \
        "the filter did not move"
    ia = a(dataset.image)[:, :, 0, 5]
    ib = b(dataset.image)[:, :, 0, 5]
    assert not np.allclose(ia, ib), "the filtered image did not change"


def test_moving_off_the_peak_loses_signal(dataset):
    """A centre away from the harmonic should pass much less energy."""
    shape = dataset.image.shape[:2]

    def energy(centre):
        omega = centre_to_frequency(centre, shape)
        mag = float(np.hypot(*omega))
        ang = float(np.arctan2(omega[1], omega[0]))
        flt = HarpFilter2D(dataset.image, (mag, dataset.ke[1]), (ang, np.pi / 2),
                           kind="butterworth", cutoff=10, order=10,
                           centre=np.array([np.round(centre), [120, 80]], dtype=int))
        return float(np.abs(flt(dataset.image)[:, :, 0, 5]).sum())

    on_peak = energy([80, 120])
    off_peak = energy([60, 100])
    assert off_peak < on_peak
