"""Following contour nodes through a sequence with the tag pattern."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.ndimage import map_coordinates

import harpi
from harpi.testing import make_phantom


def _at(field, pts):
    """Sample a displacement component at (x, y) points."""
    return map_coordinates(np.nan_to_num(field), np.vstack([pts[:, 1], pts[:, 0]]),
                           order=1)


def _ring(centre, radius, n=24):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([centre[0] + radius * np.cos(t),
                            centre[1] + radius * np.sin(t)])


@pytest.fixture(scope="module")
def phantom_track():
    """A phantom with a known displacement, and nodes tracked through it."""
    p = make_phantom(shape=(128, 128), n_frames=12, tag_spacing=8.0, r_inner=18,
                     r_outer=32, contraction=0.15, twist=0.25, t1=None, noise=0.0)
    ke = 2 * np.pi / p["tag_spacing"] * p["pixel_size"]
    Ih = harpi.HarpFilter2D(p["image"], ke, (0.0, np.pi / 2),
                            kind="transmission")(p["image"])
    pts = _ring(p["centre"], 25.0)
    tracked, info = harpi.track_points_2d(Ih, pts, reference=0, return_info=True)
    return p, Ih, pts, tracked, info


# ----------------------------------------------------- against known motion

def test_tracking_follows_a_known_displacement(phantom_track):
    """The phantom carries u_true, so this is checked against truth."""
    p, _, pts, tracked, _ = phantom_track
    worst = 0.0
    for fr in range(1, p["image"].shape[-1]):
        want = pts + np.column_stack([_at(p["u_true"][:, :, 0, fr], pts),
                                      _at(p["u_true"][:, :, 1, fr], pts)])
        worst = max(worst, np.hypot(*(tracked[fr] - want).T).max())
    assert worst < 0.5, f"worst node {worst:.3f} px"


def test_the_reference_frame_is_returned_unchanged(phantom_track):
    _, _, pts, tracked, _ = phantom_track
    assert np.array_equal(tracked[0], pts)


def test_the_shape_is_one_set_of_points_per_frame(phantom_track):
    p, _, pts, tracked, _ = phantom_track
    assert tracked.shape == (p["image"].shape[-1], len(pts), 2)


def test_the_tag_period_is_measured_from_the_image(phantom_track):
    _, _, _, _, info = phantom_track
    assert info["period"] == pytest.approx(8.0, abs=0.3)


def test_clean_data_needs_no_correction(phantom_track):
    _, _, _, _, info = phantom_track
    assert info["corrected"].sum() == 0


def test_tracking_works_from_a_middle_frame(phantom_track):
    """Both directions, since the user draws on whichever frame is clearest."""
    p, Ih, pts, _, _ = phantom_track
    mid = p["image"].shape[-1] // 2
    here = pts + np.column_stack([_at(p["u_true"][:, :, 0, mid], pts),
                                  _at(p["u_true"][:, :, 1, mid], pts)])
    tracked = harpi.track_points_2d(Ih, here, reference=mid)
    for fr in (0, mid - 2, mid + 2, p["image"].shape[-1] - 1):
        want = pts + np.column_stack([_at(p["u_true"][:, :, 0, fr], pts),
                                      _at(p["u_true"][:, :, 1, fr], pts)])
        assert np.hypot(*(tracked[fr] - want).T).max() < 1.0, f"frame {fr}"


# ------------------------------------------------------------- the guard

@pytest.fixture(scope="module")
def real_track():
    import pathlib

    path = pathlib.Path(__file__).resolve().parents[2] / "data" / "data_2D.mat"
    if not path.exists():
        pytest.skip(f"{path} not present")
    ds = harpi.load_mat_2d(path)
    Ih = harpi.HarpFilter2D(ds.image, ds.ke, (0.0, np.pi / 2),
                            kind="transmission")(ds.image)
    yy, xx = np.nonzero(ds.mask[:, :, 0])
    centre = (xx.mean(), yy.mean())
    radius = float(np.median(np.hypot(yy - centre[1], xx - centre[0])))
    return ds, Ih, _ring(centre, radius)


@pytest.fixture(scope="module")
def reference_field(real_track):
    """The RBF solve, which the tracked nodes are compared against."""
    ds, Ih, _ = real_track
    return harpi.estimate_motion_2d(Ih, ds.mask, a=ds.suggested_a, eta=15.0).u


def _wall_contour(mask, outer, n=20, inset=0.8):
    """A contour on the myocardial boundary, as the user would draw it.

    The wall is not circular: on this dataset the mask reaches from 11 to 20 px
    of the centre. A circle at the outer radius leaves most of its nodes off the
    tissue, where there is no tag signal to track, so the boundary is taken
    sector by sector instead.
    """
    ys, xs = np.nonzero(mask)
    cy, cx = ys.mean(), xs.mean()
    angle = np.arctan2(ys - cy, xs - cx) % (2 * np.pi)
    radius = np.hypot(ys - cy, xs - cx)
    pts = []
    for k in range(n):
        lo, hi = 2 * np.pi * k / n, 2 * np.pi * (k + 1) / n
        here = (angle >= lo) & (angle < hi)
        if not here.any():
            continue
        r = radius[here].max() - inset if outer else radius[here].min() + inset
        mid = 0.5 * (lo + hi)
        pts.append([cx + r * np.cos(mid), cy + r * np.sin(mid)])
    return np.array(pts)


def _worst(field, pts, tracked):
    """Largest node error against a displacement field, over every frame."""
    out = 0.0
    for fr in range(1, field.shape[-1]):
        want = pts + np.column_stack([_at(field[:, :, 0, fr], pts),
                                      _at(field[:, :, 1, fr], pts)])
        out = max(out, np.hypot(*(tracked[fr] - want).T).max())
    return out


def test_the_guard_prevents_a_period_slip(real_track, reference_field):
    """Late in the cycle the tags fade and a node can lock onto the next line.

    One slip is a whole tag period and corrupts the start of every later frame,
    which is why this is worth guarding rather than leaving to the user.
    """
    ds, Ih, pts = real_track
    loose = _worst(reference_field, pts,
                   harpi.track_points_2d(Ih, pts, regularise=False))
    tight = _worst(reference_field, pts,
                   harpi.track_points_2d(Ih, pts, regularise=True))
    period = ds.pixels_per_wavelength
    assert loose > 0.8 * period, "expected a slip of about one period"
    assert tight < 0.5 * period
    assert tight < loose / 2


def test_tracking_beats_copying_the_contour(real_track, reference_field):
    """The behaviour it replaces: the contour stays put while the tissue moves."""
    ds, Ih, pts = real_track
    res = type("R", (), {"u": reference_field})
    tracked = harpi.track_points_2d(Ih, pts)

    copied, followed = [], []
    for fr in range(1, Ih.shape[-1]):
        u = np.column_stack([_at(res.u[:, :, 0, fr], pts),
                             _at(res.u[:, :, 1, fr], pts)])
        copied.append(np.hypot(*u.T).max())           # copying leaves u behind
        followed.append(np.hypot(*(tracked[fr] - (pts + u)).T).max())
    assert max(followed) < max(copied) / 2


def test_the_guard_fires_rarely(real_track):
    """It is a safety net, not a smoother: it must not be reshaping the result."""
    _, Ih, pts = real_track
    _, info = harpi.track_points_2d(Ih, pts, return_info=True)
    node_frames = len(pts) * (Ih.shape[-1] - 1)
    assert info["corrected"].sum() <= 0.02 * node_frames


# ------------------------------------------------------------- the inputs

def test_a_bad_image_shape_is_refused():
    with pytest.raises(ValueError, match="rows, cols, 2, frames"):
        harpi.track_points_2d(np.zeros((8, 8, 3, 4), complex), np.zeros((3, 2)))


def test_bad_points_are_refused():
    with pytest.raises(ValueError, match=r"\(n, 2\)"):
        harpi.track_points_2d(np.zeros((8, 8, 2, 4), complex), np.zeros((3, 3)))


def test_a_reference_outside_the_sequence_is_refused():
    with pytest.raises(ValueError, match="outside"):
        harpi.track_points_2d(np.ones((8, 8, 2, 4), complex), np.zeros((3, 2)),
                              reference=9)


def test_too_few_nodes_to_regularise_is_not_an_error(phantom_track):
    """Two nodes have no neighbourhood to compare against."""
    _, Ih, _, _, _ = phantom_track
    out = harpi.track_points_2d(Ih, np.array([[60.0, 60.0], [64.0, 60.0]]))
    assert out.shape == (Ih.shape[-1], 2, 2)
    assert np.isfinite(out).all()


# ----------------------------------------------- the guard, on its own terms

PERIOD = 4.0


def _uniform(n, step=(1.0, 0.5)):
    """A contour moving as one, which is what a correct hop looks like."""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    before = np.column_stack([60 + 15 * np.cos(t), 60 + 15 * np.sin(t)])
    return before, before + np.asarray(step)


def test_a_slipped_node_is_pulled_back():
    from harpi.tracking import _regularise

    before, after = _uniform(20)
    after[7, 0] += PERIOD
    fixed, count = _regularise(before, after, 0.5 * PERIOD)
    assert count == 1
    assert np.hypot(*(fixed[7] - (before[7] + (1.0, 0.5)))) < 0.1


def test_two_slipped_neighbours_are_both_pulled_back():
    """Where a median of three neighbours fails: two of the three have slipped,
    so the median is the wrong answer and the node is left where it landed."""
    from harpi.tracking import _regularise

    before, after = _uniform(20)
    after[7:9, 0] += PERIOD
    fixed, count = _regularise(before, after, 0.5 * PERIOD)
    assert count == 2
    for node in (7, 8):
        assert np.hypot(*(fixed[node] - (before[node] + (1.0, 0.5)))) < 0.1


def test_a_run_of_slipped_nodes_is_pulled_back():
    """A quarter of the contour slipping together is a bump the modes can
    follow, so the fit has to be kept away from it rather than merely robust."""
    from harpi.tracking import _regularise

    before, after = _uniform(20)
    after[7:12, 0] += PERIOD
    fixed, count = _regularise(before, after, 0.5 * PERIOD)
    assert count == 5
    assert np.abs(fixed[7:12] - (before[7:12] + (1.0, 0.5))).max() < 0.1


def test_a_contour_that_tracked_correctly_is_left_alone():
    """It is a guard, not a smoother."""
    from harpi.tracking import _regularise

    before, after = _uniform(20)
    rng = np.random.default_rng(0)
    after = after + rng.normal(0.0, 0.05, after.shape)
    fixed, count = _regularise(before, after, 0.5 * PERIOD)
    assert count == 0
    assert np.array_equal(fixed, after)


def test_a_real_deformation_is_not_flattened():
    """A dilation with twist is within the modes the fit carries, so none of its
    nodes should look like outliers."""
    from harpi.tracking import _regularise

    t = np.linspace(0, 2 * np.pi, 20, endpoint=False)
    before = np.column_stack([60 + 15 * np.cos(t), 60 + 15 * np.sin(t)])
    radial, twist = 0.12, 0.2
    after = np.column_stack([
        60 + 15 * (1 + radial) * np.cos(t + twist),
        60 + 15 * (1 + radial) * np.sin(t + twist)])
    fixed, count = _regularise(before, after, 0.5 * PERIOD)
    assert count == 0
    assert np.array_equal(fixed, after)


def test_a_short_contour_falls_back_to_its_neighbours():
    """Five nodes cannot carry the modes, and must still be guarded."""
    from harpi.tracking import _regularise

    before, after = _uniform(5)
    after[2, 1] += PERIOD
    fixed, count = _regularise(before, after, 0.5 * PERIOD)
    assert count == 1
    assert np.hypot(*(fixed[2] - (before[2] + (1.0, 0.5)))) < 0.1


# -------------------------------------------- on the contours a user draws

def test_both_wall_boundaries_are_tracked_within_half_a_period(
        real_track, reference_field):
    """The endocardium borders the blood pool, where the tags wash out, and is
    the harder of the two. Half a period is the limit that matters: past it a
    node has locked onto the neighbouring tag line.
    """
    ds, Ih, _ = real_track
    limit = 0.5 * ds.pixels_per_wavelength
    for outer in (True, False):
        pts = _wall_contour(ds.mask[:, :, 0], outer)
        worst = _worst(reference_field, pts, harpi.track_points_2d(Ih, pts))
        assert worst < limit, (
            f"{'epi' if outer else 'endo'}cardial worst {worst:.3f} px")


def test_the_guard_matters_most_on_the_endocardium(real_track, reference_field):
    ds, Ih, _ = real_track
    pts = _wall_contour(ds.mask[:, :, 0], outer=False)
    loose = _worst(reference_field, pts,
                   harpi.track_points_2d(Ih, pts, regularise=False))
    tight = _worst(reference_field, pts, harpi.track_points_2d(Ih, pts))
    assert loose > ds.pixels_per_wavelength, "expected a slip to correct"
    assert tight < loose / 3


# ------------------------------------------------------------- progress

def test_progress_is_reported_once_per_hop(phantom_track):
    _, Ih, pts, _, _ = phantom_track
    seen = []
    harpi.track_points_2d(Ih, pts, reference=3,
                          progress=lambda d, t, s: seen.append((d, t, s)))
    nfr = Ih.shape[-1]
    assert len(seen) == nfr - 1
    assert [d for d, _, _ in seen] == list(range(1, nfr))
    assert {t for _, t, _ in seen} == {nfr - 1}


def test_raising_from_progress_stops_the_run(phantom_track):
    """How the interface cancels."""
    _, Ih, pts, _, _ = phantom_track

    class Stop(Exception):
        pass

    def report(done, total, stage):
        if done == 2:
            raise Stop()

    with pytest.raises(Stop):
        harpi.track_points_2d(Ih, pts, progress=report)


# --------------------------------------- both directions from a middle frame

def test_every_frame_either_side_is_filled(real_track):
    """Drawing on a middle frame must fill the frames before it as well as
    after. An untouched frame comes back holding the input unchanged, which is
    what a one-directional run would leave behind."""
    _, Ih, pts = real_track
    nfr = Ih.shape[-1]
    for ref in (0, 1, nfr // 2, nfr - 2, nfr - 1):
        tracked = harpi.track_points_2d(Ih, pts, reference=ref)
        for fr in range(nfr):
            moved = not np.array_equal(tracked[fr], pts)
            assert moved == (fr != ref), f"reference {ref}, frame {fr}"


def test_the_two_directions_are_a_round_trip(real_track):
    """Tracking out to a frame and back has to return to where it began, or one
    direction is doing something the other is not."""
    _, Ih, pts = real_track
    out = harpi.track_points_2d(Ih, pts, reference=0)
    for k in (5, 9, Ih.shape[-1] - 1):
        back = harpi.track_points_2d(Ih, out[k], reference=k)
        err = np.hypot(*(back[0] - pts).T)
        assert np.median(err) < 0.01, f"0 to {k} and back, median {np.median(err)}"
        assert err.max() < 2.0, f"0 to {k} and back, worst {err.max()}"


def test_a_middle_reference_costs_no_accuracy(real_track, reference_field):
    """Half the frames are then reached backwards, which must be no worse than
    reaching them all forwards."""
    ds, Ih, pts = real_track
    nfr = Ih.shape[-1]
    mid = nfr // 2
    forward = harpi.track_points_2d(Ih, pts, reference=0)
    from_mid = harpi.track_points_2d(Ih, forward[mid], reference=mid)

    worst = 0.0
    for fr in range(nfr):
        want = pts + np.column_stack([_at(reference_field[:, :, 0, fr], pts),
                                      _at(reference_field[:, :, 1, fr], pts)])
        worst = max(worst, np.hypot(*(from_mid[fr] - want).T).max())
    assert worst < ds.pixels_per_wavelength
