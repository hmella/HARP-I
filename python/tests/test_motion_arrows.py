"""The displacement arrows on the Motion tab.

`u[..., 0]` is the displacement along columns and `u[..., 1]` along rows, which
the library pins elsewhere. The drawn axis is inverted so that rows increase
downwards, so the row component goes to `quiver` unchanged. Negating it as well
flipped the arrows a second time: a contraction converged left and right while
fanning outwards top and bottom, which reads as x and y being swapped.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication  # noqa: E402

from harpi import MotionResult  # noqa: E402
from harpi.gui.session import Session, Stage  # noqa: E402
from harpi.gui.tabs_motion import MotionTab  # noqa: E402

ROWS = COLS = 48
CY = CX = (ROWS - 1) / 2.0


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


def _tab_with(u, qapp):
    """A Motion tab holding one displacement field, refreshed."""
    nfr = u.shape[-1]
    sess = Session()
    sess.mask = np.repeat(np.isfinite(u[:, :, 0, -1])[:, :, None], nfr, 2)
    sess.current_frame = nfr - 1
    sess.set_result(Stage.MOTION, MotionResult(
        u=u, phase_jumps=np.zeros((2, nfr)),
        meta={"jump_search": {"njumps": 8,
                              "offset": np.zeros((2, nfr), dtype=int)},
              "phase_jump_rule": "nearest"}))
    tab = MotionTab(sess)
    tab.refresh()
    return tab


def _quiver(tab):
    """The Quiver artist, which is not on the first axes: that is the image,
    and the second is its colour bar."""
    for ax in tab.plot.fig.axes:
        for coll in ax.collections:
            if hasattr(coll, "U") and hasattr(coll, "V"):
                return coll, ax
    raise AssertionError("no quiver was drawn")


def _arrows(tab):
    """U, V and the positions, with the masked entries back as NaN.

    matplotlib masks the entries that were NaN and leaves a fill value in the
    data, so reading .U directly reports arrows where there are none.
    """
    coll, ax = _quiver(tab)
    # set_UVC fills the invalid entries with 1 and keeps the mask in Umask, so
    # reading U and V alone reports unit arrows where there is no data.
    U = np.asarray(coll.U, dtype=float).copy()
    V = np.asarray(coll.V, dtype=float).copy()
    flat = np.asarray(coll.Umask, dtype=bool).ravel()
    U[flat] = np.nan
    V[flat] = np.nan
    off = np.asarray(coll.get_offsets())
    return U, V, off, ax


def _radial(scale=-0.25, nfr=3):
    """A pure contraction towards the centre, with a hole and a rim."""
    Y, X = np.mgrid[0:ROWS, 0:COLS].astype(float)
    r = np.hypot(X - CX, Y - CY)
    mask = (r > 6) & (r < 20)
    u = np.zeros((ROWS, COLS, 2, nfr))
    for k in range(nfr):
        s = scale * k / max(nfr - 1, 1)
        u[:, :, 0, k] = np.where(mask, s * (X - CX), np.nan)
        u[:, :, 1, k] = np.where(mask, s * (Y - CY), np.nan)
    return u


def _translation(dx, dy, nfr=3):
    Y, X = np.mgrid[0:ROWS, 0:COLS].astype(float)
    r = np.hypot(X - CX, Y - CY)
    mask = (r > 6) & (r < 20)
    u = np.zeros((ROWS, COLS, 2, nfr))
    for k in range(nfr):
        f = k / max(nfr - 1, 1)
        u[:, :, 0, k] = np.where(mask, dx * f, np.nan)
        u[:, :, 1, k] = np.where(mask, dy * f, np.nan)
    return u


def test_the_row_component_is_drawn_unchanged(qapp):
    """The axis inversion already orients it; negating as well flips it twice."""
    u = _translation(0.0, 3.0)
    U, V, _, _ = _arrows(_tab_with(u, qapp))
    good = np.isfinite(V) & (np.abs(V) > 1e-9)
    assert good.any()
    assert np.all(V[good] > 0), "a positive row displacement must stay positive"


def test_the_column_component_is_drawn_unchanged(qapp):
    u = _translation(3.0, 0.0)
    U, V, _, _ = _arrows(_tab_with(u, qapp))
    good = np.isfinite(U) & (np.abs(U) > 1e-9)
    assert good.any()
    assert np.all(U[good] > 0)


def test_the_drawn_axis_increases_downwards(qapp):
    """Which is why the row component needs no sign change."""
    _, _, _, ax = _arrows(_tab_with(_radial(), qapp))
    lo, hi = ax.get_ylim()
    assert lo > hi, "the y axis should be inverted to match image rows"


def test_a_contraction_points_inwards_everywhere(qapp):
    """The check that would have caught this: every arrow towards the centre."""
    u, v, off, _ = _arrows(_tab_with(_radial(scale=-0.25), qapp))
    good = np.isfinite(u) & np.isfinite(v) & (np.hypot(u, v) > 1e-9)
    assert good.sum() > 20

    # the crop is centred on the mask, so the mean position is the ring centre
    px, py = off[:, 0], off[:, 1]
    cx, cy = px[good].mean(), py[good].mean()
    inward = -((px - cx) * u + (py - cy) * v)
    assert np.all(inward[good] > 0), "every arrow of a contraction points inwards"


def test_an_expansion_points_outwards_everywhere(qapp):
    u, v, off, _ = _arrows(_tab_with(_radial(scale=+0.25), qapp))
    good = np.isfinite(u) & np.isfinite(v) & (np.hypot(u, v) > 1e-9)
    px, py = off[:, 0], off[:, 1]
    cx, cy = px[good].mean(), py[good].mean()
    outward = (px - cx) * u + (py - cy) * v
    assert np.all(outward[good] > 0)


def test_the_two_components_are_not_interchanged(qapp):
    """A field along columns only must draw nothing along rows."""
    U, V, _, _ = _arrows(_tab_with(_translation(3.0, 0.0), qapp))
    assert np.nanmax(np.abs(V)) < 1e-9
    assert np.nanmax(np.abs(U)) > 1.0


# ------------------------------------------- the tagged image underneath

def _with_dataset(u, qapp, on=True):
    """A Motion tab whose session also holds an image to draw under."""
    from harpi import Dataset2D

    tab = _tab_with(u, qapp)
    nfr = u.shape[-1]
    rng = np.random.default_rng(0)
    tab.session.dataset = Dataset2D(
        image=rng.normal(size=(ROWS, COLS, 2, nfr)),
        mask=np.ones((ROWS, COLS, nfr), bool),
        pixel_size=[2.0, 2.0], tag_spacing=8.0)
    tab.bg.active.setChecked(on)
    tab.refresh()
    return tab


def _images(tab):
    """The background layers only. The data is drawn with imshow too, so the
    two are told apart by the z order the tab gives them."""
    from matplotlib.image import AxesImage

    return [(ax, im) for ax in tab.plot.fig.axes
            for im in ax.images
            if isinstance(im, AxesImage) and im.get_zorder() == 0]


def _panels(tab):
    """The two data panels, leaving out the colour bar."""
    _, arrows = _quiver(tab)
    maps = [ax for ax in tab.plot.fig.axes
            if ax.get_title().startswith("|u|")]
    return maps + [arrows]


def test_the_overlay_is_on_by_default(qapp):
    tab = _tab_with(_radial(), qapp)
    assert tab.bg.active.isChecked()


def test_the_image_is_drawn_under_both_panels(qapp):
    tab = _with_dataset(_radial(), qapp)
    axes = {ax for ax, _ in _images(tab)}
    _, ax_arrows = _quiver(tab)
    assert ax_arrows in axes, "nothing under the arrows"
    assert len(axes) >= 2, "the magnitude panel too"


def test_turning_it_off_leaves_the_arrows_alone(qapp):
    tab = _with_dataset(_radial(), qapp, on=False)
    _, ax_arrows = _quiver(tab)
    assert not [im for ax, im in _images(tab) if ax is ax_arrows]


def test_the_arrows_keep_their_direction_over_an_image(qapp):
    """imshow orients rows downwards by itself, so inverting as well would
    flip them a second time."""
    u = _translation(0.0, 3.0)
    U, V, _, _ = _arrows(_with_dataset(u, qapp))
    good = np.isfinite(V) & (np.abs(V) > 1e-9)
    assert good.any()
    assert np.all(V[good] > 0)


def test_a_contraction_still_points_inwards_over_an_image(qapp):
    u, v, off, _ = _arrows(_with_dataset(_radial(scale=-0.25), qapp))
    good = np.isfinite(u) & np.isfinite(v) & (np.hypot(u, v) > 1e-9)
    px, py = off[:, 0], off[:, 1]
    cx, cy = px[good].mean(), py[good].mean()
    inward = -((px - cx) * u + (py - cy) * v)
    assert np.all(inward[good] > 0)


def test_both_panels_keep_rows_increasing_downwards(qapp):
    for on in (True, False):
        tab = _with_dataset(_radial(), qapp, on=on)
        for ax in _panels(tab):
            lo, hi = ax.get_ylim()
            assert lo > hi, f"overlay={on}: row 0 is not at the top"


def test_the_image_shares_the_crop_with_the_data(qapp):
    """Otherwise the overlay sits off the anatomy it is meant to describe."""
    tab = _with_dataset(_radial(), qapp)
    _, ax_arrows = _quiver(tab)
    under = [im for ax, im in _images(tab) if ax is ax_arrows][0]
    x0, x1, y0, y1 = under.get_extent()
    assert (x1 - x0, y0 - y1) == pytest.approx(
        (ax_arrows.get_xlim()[1] - ax_arrows.get_xlim()[0],
         ax_arrows.get_ylim()[0] - ax_arrows.get_ylim()[1]))


def test_no_dataset_means_no_background(qapp):
    """The tab has to keep working before any image is loaded."""
    tab = _tab_with(_radial(), qapp)
    assert tab.session.dataset is None
    assert not _images(tab)


# ------------------------------------------------------- how long the arrows are

def _lengths(tab):
    u, v, _, _ = _arrows(tab)
    good = np.isfinite(u) & np.isfinite(v)
    return np.hypot(u[good], v[good])


def test_the_arrows_are_drawn_at_their_true_length(qapp):
    """So an arrow on screen is the displacement in pixels, readable against
    the image under it."""
    tab = _tab_with(_translation(3.0, 0.0), qapp)
    assert tab.arrow_scale.value() == pytest.approx(1.0)
    coll, _ = _quiver(tab)
    assert coll.scale == pytest.approx(1.0)
    assert np.nanmax(_lengths(tab)) == pytest.approx(3.0, rel=1e-6)


def test_the_scale_control_lengthens_them(qapp):
    tab = _tab_with(_translation(3.0, 0.0), qapp)
    tab.arrow_scale.setValue(2.5)
    coll, _ = _quiver(tab)
    assert coll.scale == pytest.approx(1 / 2.5)


def test_the_scale_control_does_not_change_the_data(qapp):
    """It is how they are drawn, not what they are."""
    tab = _tab_with(_translation(3.0, 0.0), qapp)
    before = _lengths(tab).copy()
    tab.arrow_scale.setValue(4.0)
    assert np.allclose(before, _lengths(tab))
