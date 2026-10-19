"""Rearranging the displacement components after a solve."""

import numpy as np
import pytest

from harpi import flip_components
from harpi.types import MotionResult


@pytest.fixture
def result():
    rows, cols, nfr = 4, 5, 3
    u = np.zeros((rows, cols, 2, nfr))
    u[..., 0, :] = 1.0
    u[..., 1, :] = 2.0
    Xg, Yg = np.meshgrid(np.arange(1, cols + 1), np.arange(1, rows + 1))
    traj = np.stack([Xg, Yg], axis=-1)[..., None] + u
    return MotionResult(u=u, u_fitted=u * 10, trajectories=traj)


def test_no_options_returns_the_same_object(result):
    assert flip_components(result) is result


def test_swap_exchanges_the_components(result):
    out = flip_components(result, swap_xy=True)
    assert out.u[0, 0, 0, 0] == 2.0
    assert out.u[0, 0, 1, 0] == 1.0


def test_flip_negates_one_component(result):
    out = flip_components(result, flip_y=True)
    assert out.u[0, 0, 0, 0] == 1.0
    assert out.u[0, 0, 1, 0] == -2.0


def test_flip_applies_after_the_swap(result):
    """'flip x' negates whatever ends up in x."""
    out = flip_components(result, swap_xy=True, flip_x=True)
    assert out.u[0, 0, 0, 0] == -2.0
    assert out.u[0, 0, 1, 0] == 1.0


def test_fitted_motion_is_transformed_too(result):
    out = flip_components(result, swap_xy=True)
    assert out.u_fitted[0, 0, 0, 0] == 20.0
    assert out.u_fitted[0, 0, 1, 0] == 10.0


def test_the_original_is_not_modified(result):
    before = result.u.copy()
    flip_components(result, swap_xy=True, flip_x=True, flip_y=True)
    assert np.array_equal(result.u, before)


def test_trajectories_stay_consistent(result):
    out = flip_components(result, swap_xy=True, flip_y=True)
    rows, cols = result.u.shape[:2]
    Xg, Yg = np.meshgrid(np.arange(1, cols + 1), np.arange(1, rows + 1))
    grid = np.stack([Xg, Yg], axis=-1)
    assert np.allclose(out.trajectories, grid[..., None] + out.u)


def test_applying_twice_restores_the_original(result):
    once = flip_components(result, flip_x=True, flip_y=True)
    twice = flip_components(once, flip_x=True, flip_y=True)
    assert np.allclose(twice.u, result.u)


def test_swapping_twice_restores_the_original(result):
    once = flip_components(result, swap_xy=True)
    twice = flip_components(once, swap_xy=True)
    assert np.allclose(twice.u, result.u)


def test_three_dimensional_results(result):
    u = np.zeros((2, 3, 4, 3, 2))
    for c in range(3):
        u[..., c, :] = c + 1
    out = flip_components(MotionResult(u=u), swap_xy=True, flip_z=True)
    assert list(out.u[0, 0, 0, :, 0]) == [2.0, 1.0, -3.0]


def test_the_options_are_recorded(result):
    out = flip_components(result, flip_x=True)
    assert out.meta["flip"]["flip_x"] is True
    assert out.meta["flip"]["swap_xy"] is False
