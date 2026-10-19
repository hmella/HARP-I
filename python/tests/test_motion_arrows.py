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
