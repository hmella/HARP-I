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


# ------------------------------------------- where the response actually is

def _peak(H):
    return np.unravel_index(np.argmax(H), H.shape)


@pytest.mark.parametrize("centre", [[80, 120], [120, 80], [40, 40], [5, 150]])
def test_the_butterworth_response_peaks_on_the_centre(centre):
    """Its grid counted from 1 while the centre is a NumPy index, so the
    response sat one pixel up and to the left of the centre placed.

    A sharp filter, because at the interface defaults the top is flat over
    twenty-one pixels and argmax would only report the first of them.
    """
    from harpi.filters import _butterworth_2d

    H = _butterworth_2d((160, 160), centre, 2, 1)
    assert _peak(H) == tuple(centre)


def test_both_kinds_agree_where_the_centre_is():
    """The regression the one-pixel offset was: the two filter kinds put the
    same centre in different places."""
    from harpi.filters import _butterworth_2d, _transmission

    shape = (160, 160)
    centre = [80, 120]
    omega = centre_to_frequency(centre, shape)
    k = np.array(omega)
    w1 = 2 * np.pi / np.sum(k**2) * k
    w2 = 2 * np.pi / np.sum(k**2) * np.array([-k[1], k[0]])

    assert _peak(_transmission(shape, w1, w2)) == tuple(centre)
    assert _peak(_butterworth_2d(shape, centre, 2, 1)) == tuple(centre)


def test_the_butterworth_response_is_symmetric_about_its_centre():
    """Which holds at the interface defaults, where the top is flat."""
    from harpi.filters import _butterworth_2d

    r, c = 80, 120
    H = _butterworth_2d((160, 160), [r, c], 17, 10)
    for d in (1, 5, 13):
        assert H[r, c - d] == pytest.approx(H[r, c + d])
        assert H[r - d, c] == pytest.approx(H[r + d, c])


# ----------------------------------------------------- one frame at a time

def test_filtering_one_plane_matches_the_whole_series(dataset):
    """Bitwise, not approximately. It is what lets the interface preview a
    setting honestly without paying for the full solve."""
    flt = HarpFilter2D(dataset.image, dataset.ke, (0.0, np.pi / 2),
                       kind="transmission")
    whole = flt(dataset.image)
    for enc in (0, 1):
        for fr in (0, 5, dataset.n_frames - 1):
            one = flt.filter_plane(dataset.image[:, :, enc, fr], enc)
            assert np.array_equal(one, whole[:, :, enc, fr]), (enc, fr)


def test_filtering_a_stack_of_planes_keeps_the_frame_axis(dataset):
    flt = HarpFilter2D(dataset.image, dataset.ke, (0.0, np.pi / 2),
                       kind="transmission")
    out = flt.filter_plane(dataset.image[:, :, 0, :], 0)
    assert out.shape == dataset.image.shape[:2] + (dataset.n_frames,)
    assert np.array_equal(out, flt(dataset.image)[:, :, 0, :])


def test_a_plane_of_the_wrong_size_is_refused(dataset):
    flt = HarpFilter2D(dataset.image, dataset.ke, (0.0, np.pi / 2),
                       kind="transmission")
    with pytest.raises(ValueError, match="built for"):
        flt.filter_plane(np.zeros((8, 8)), 0)


def test_the_response_is_one_encoding_of_the_filter(dataset):
    flt = HarpFilter2D(dataset.image, dataset.ke, (0.0, np.pi / 2),
                       kind="transmission")
    for enc in (0, 1):
        assert np.array_equal(flt.response(enc), flt.kspace_filter[:, :, enc])


def test_the_transmission_window_ignores_cutoff_and_order(dataset):
    """Two of the three bandpass controls do nothing on this kind, which is
    why the interface disables them."""
    kw = dict(kind="transmission", centre=[[80, 120], [120, 80]])
    a = HarpFilter2D(dataset.image, dataset.ke, (0.0, np.pi / 2),
                     cutoff=5, order=2, **kw)
    b = HarpFilter2D(dataset.image, dataset.ke, (0.0, np.pi / 2),
                     cutoff=120, order=40, **kw)
    assert np.array_equal(a.kspace_filter, b.kspace_filter)
