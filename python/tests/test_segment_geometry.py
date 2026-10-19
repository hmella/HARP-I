"""Where the AHA segments fall, and which way they are numbered.

The layout used to be numbered from half a turn away from the insertion, so
every name sat on the wrong wall: "mid anterior" landed inferolaterally. The
geometric test below is the one that would have caught it, and it is written
in anatomical terms rather than in angles so it stays readable.

Image angles throughout: y grows downwards, so -90 degrees is straight up and
the septum of a conventionally displayed short axis slice is on the left.
"""

from __future__ import annotations

import numpy as np
import pytest

import harpi
from harpi.analysis import segment_geometry, segment_index
from harpi.types import StrainResult

CENTRE = 41.0
N = 81


def _at(degrees, radius=19.0):
    """A point on a ray from the centre, in the 1-based library grid."""
    a = np.radians(degrees)
    return (CENTRE + radius * np.cos(a), CENTRE + radius * np.sin(a))


@pytest.fixture(scope="module")
def ring():
    """A uniformly contracting annulus, so only the geometry varies."""
    Y, X = np.mgrid[1:N + 1, 1:N + 1]
    r = np.hypot(X - CENTRE, Y - CENTRE)
    mask = (r > 10) & (r < 28)
    cc = np.where(mask[..., None], -0.1, np.nan) * np.ones(2)
    return StrainResult(CC=cc, RR=cc.copy(), mask=mask,
                        origin=(CENTRE, CENTRE))


def _segment_at(seg, degrees):
    x, y = _at(degrees)
    return int(seg.segments[int(round(y)) - 1, int(round(x)) - 1])


def _cut(ring, **kw):
    kw.setdefault("origin", (CENTRE, CENTRE))
    return harpi.strain_by_segments(ring, **kw)


# ------------------------------------------------ the layout, anatomically

#: a septum 120 degrees wide on the left, which is the model exactly
ANTERIOR, INFERIOR = _at(-120), _at(120)

WALLS = [(-90, "mid anterior"), (-150, "mid anteroseptal"),
         (150, "mid inferoseptal"), (90, "mid inferior"),
         (30, "mid inferolateral"), (-30, "mid anterolateral")]


def test_the_aha_names_land_on_the_right_walls(ring):
    """With the insertions placed, every name has to sit where anatomy says.

    Anterior is up, the septum is on the left between the two insertions, and
    the lateral wall is on the right.
    """
    seg = _cut(ring, n_segments=6, anterior=ANTERIOR, inferior=INFERIOR)
    names = harpi.aha_labels("mid", 6)
    for degrees, want in WALLS:
        got = names[_segment_at(seg, degrees) - 1]
        assert got == want, f"at {degrees} deg the wedge is {got}, want {want}"


def test_the_anterior_ray_divides_segment_one_from_segment_two(ring):
    """The insertion separates the septum from the anterior wall, so it is a
    boundary, and segment 1 lies on the side away from the septum."""
    seg = _cut(ring, n_segments=6, anterior=ANTERIOR, inferior=INFERIOR)
    # the septum runs from the anterior ray towards 180, so inside it is the
    # more negative side
    assert _segment_at(seg, -120 - 4) == 2, "just inside the septum"
    assert _segment_at(seg, -120 + 4) == 1, "just outside it"


def test_the_inferior_ray_closes_the_septum(ring):
    seg = _cut(ring, n_segments=6, anterior=ANTERIOR, inferior=INFERIOR)
    assert _segment_at(seg, 120 + 4) == 3, "last septal wedge"
    assert _segment_at(seg, 120 - 4) == 4, "first wedge past it"


def test_apical_puts_the_septum_in_one_wedge(ring):
    """Four segments of 90 degrees, of which the septum is one."""
    ant, inf = _at(-135), _at(135)
    seg = _cut(ring, n_segments=4, anterior=ant, inferior=inf)
    names = harpi.aha_labels("apical", 4)
    assert names[_segment_at(seg, -90) - 1] == "apical anterior"
    assert names[_segment_at(seg, 180) - 1] == "apical septal"
    assert names[_segment_at(seg, 90) - 1] == "apical inferior"
    assert names[_segment_at(seg, 0) - 1] == "apical lateral"


# --------------------------------------------------------- the direction

def test_the_two_points_set_the_direction(ring):
    """Swapping which point is anterior reverses the numbering."""
    forward = _cut(ring, n_segments=6, anterior=ANTERIOR, inferior=INFERIOR)
    back = _cut(ring, n_segments=6, anterior=INFERIOR, inferior=ANTERIOR)
    assert forward.geometry["direction"] == -back.geometry["direction"]


def test_the_clockwise_flag_is_ignored_when_both_points_are_given(ring):
    a = _cut(ring, n_segments=6, anterior=ANTERIOR, inferior=INFERIOR,
             clockwise=False)
    b = _cut(ring, n_segments=6, anterior=ANTERIOR, inferior=INFERIOR,
             clockwise=True)
    assert a.geometry["direction_from"] == "insertions"
    assert np.array_equal(np.nan_to_num(a.segments, nan=-1),
                          np.nan_to_num(b.segments, nan=-1))


def test_one_point_falls_back_to_the_flag(ring):
    """The ray stays the boundary between segments 1 and 2; which side each
    falls on is what the flag decides."""
    seen = []
    for clockwise in (False, True):
        seg = _cut(ring, n_segments=6, anterior=ANTERIOR, clockwise=clockwise)
        assert seg.geometry["direction_from"] == "clockwise flag"
        assert seg.geometry["septal_arc"] is None
        pair = (_segment_at(seg, -120 - 4), _segment_at(seg, -120 + 4))
        assert set(pair) == {1, 2}, f"the ray is not a 1/2 boundary: {pair}"
        seen.append(pair)
    assert seen[0] == seen[1][::-1], "the flag has to reverse them"


def test_degenerate_insertions_fall_back_without_raising(ring):
    """Coincident, antipodal, or sitting on the centre: none of these say
    which way round the septum runs."""
    for anterior, inferior in ((ANTERIOR, ANTERIOR),
                               (_at(-120), _at(60)),
                               (ANTERIOR, (CENTRE, CENTRE))):
        seg = _cut(ring, n_segments=6, anterior=anterior, inferior=inferior)
        assert seg.geometry["direction_from"] == "clockwise flag"


# ------------------------------------------------- equal wedges, measured

def test_the_septum_is_not_stretched_to_fit_the_points(ring):
    """The AHA model is equiangular. Insertions 100 degrees apart still give
    six wedges of 60, which is the whole reason the measured arc is worth
    reporting back."""
    seg = _cut(ring, n_segments=6, anterior=_at(-130), inferior=_at(130))
    counts = np.array([(seg.segments == s + 1).sum() for s in range(6)])
    assert counts.min() > 0
    spread = counts.max() - counts.min()
    assert spread <= 0.1 * counts.mean(), f"wedges differ by {spread} pixels"


def test_the_measured_septal_arc_is_reported(ring):
    seg = _cut(ring, n_segments=6, anterior=_at(-130), inferior=_at(130))
    assert np.degrees(seg.geometry["septal_arc"]) == pytest.approx(100.0, abs=0.5)


def test_the_ideal_septal_arc_is_known_only_for_the_model():
    assert np.degrees(harpi.expected_septal_arc("mid", 6)) == pytest.approx(120)
    assert np.degrees(harpi.expected_septal_arc("basal", 6)) == pytest.approx(120)
    assert np.degrees(harpi.expected_septal_arc("apical", 4)) == pytest.approx(90)
    assert harpi.expected_septal_arc("mid", 5) is None, "not the model"
    assert harpi.expected_septal_arc(None, 6) is None


# ----------------------------------------------------------- the inputs

def test_the_matlab_spelling_still_works(ring):
    """position_a and position_b are the x and y of the anterior point."""
    named = _cut(ring, n_segments=6, anterior=ANTERIOR)
    matlab = _cut(ring, n_segments=6, position_a=ANTERIOR[0],
                  position_b=ANTERIOR[1])
    assert np.array_equal(np.nan_to_num(named.segments, nan=-1),
                          np.nan_to_num(matlab.segments, nan=-1))


def test_giving_both_spellings_is_refused(ring):
    with pytest.raises(ValueError, match="not both"):
        _cut(ring, anterior=ANTERIOR, position_a=3.0)


def test_every_masked_pixel_gets_exactly_one_segment(ring):
    seg = _cut(ring, n_segments=6, anterior=ANTERIOR, inferior=INFERIOR)
    inside = seg.segments[ring.mask]
    assert np.isfinite(inside).all()
    assert set(np.unique(inside)) == {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}
    assert np.isnan(seg.segments[~ring.mask]).all()


def test_the_index_helper_covers_the_whole_circle():
    geom = segment_geometry(n_segments=6, origin=(0.0, 0.0),
                            anterior=(1.0, 0.0))
    theta = np.linspace(-np.pi, np.pi, 721)
    index = segment_index(theta, geom)
    assert index.min() == 1 and index.max() == 6
    assert set(np.unique(index)) == {1, 2, 3, 4, 5, 6}


def test_the_short_labels_stay_distinct():
    """Truncating the names instead turned anteroseptal and anterolateral
    both into "antero"."""
    from harpi.analysis import short_labels

    short = short_labels("mid", 6)
    assert len(set(short)) == 6, short
    assert all(len(s) <= 8 for s in short), short
    assert short[0] == "ant" and short[1] == "ant-sept"
    assert short_labels("apical", 4) == ["ant", "sept", "inf", "lat"]
    assert short_labels(None, 3) == ["1", "2", "3"], "numbers when unnamed"
