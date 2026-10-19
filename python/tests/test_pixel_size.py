"""Anisotropic pixel size, in the data types and in the Data tab."""

from __future__ import annotations

import numpy as np
import pytest

from harpi import Dataset2D, Dataset3D


def _ds2(row=0.002, col=0.002, tag=0.008):
    return Dataset2D(image=np.zeros((8, 8, 2, 3)), mask=np.ones((8, 8, 3)),
                     pixel_size=[row, col], tag_spacing=tag)


def _ds3(row=0.002, col=0.002, thick=0.008, tag=0.008):
    return Dataset3D(image=np.zeros((8, 8, 4, 3, 3)), mask=np.ones((8, 8, 4, 3)),
                     voxel_size=[row, col, thick], tag_spacing=tag)


def test_square_pixels_are_unchanged():
    ds = _ds2()
    assert ds.is_isotropic
    assert ds.pixels_per_wavelength == pytest.approx(4.0)
    assert ds.suggested_a == pytest.approx(37.5)


def test_ke_follows_each_axis():
    ds = _ds2(col=0.004)
    assert ds.ke[1] == pytest.approx(2 * ds.ke[0])


def test_anisotropy_is_detected():
    assert not _ds2(col=0.004).is_isotropic
    assert _ds2().is_isotropic


def test_pixels_per_wavelength_is_reported_per_axis():
    ds = _ds2(col=0.004)
    assert ds.pixels_per_wavelength_per_axis == pytest.approx([4.0, 2.0])


def test_the_single_figure_averages_the_in_plane_axes():
    ds = _ds2(col=0.004)
    assert ds.pixels_per_wavelength == pytest.approx(3.0)
    assert ds.suggested_a == pytest.approx(50.0)


def test_slice_thickness_is_left_out_of_the_average():
    """A thick slice must not drag the suggested shape parameter with it."""
    thin = _ds3(thick=0.002)
    thick = _ds3(thick=0.020)
    assert thin.pixels_per_wavelength == pytest.approx(thick.pixels_per_wavelength)
    assert thin.suggested_a == pytest.approx(thick.suggested_a)


def test_three_dimensional_ke_covers_every_axis():
    ds = _ds3(thick=0.008)
    assert ds.ke.size == 3
    assert ds.ke[2] == pytest.approx(4 * ds.ke[0])


def test_wave_vectors_are_the_diagonal_of_ke():
    ds = _ds3()
    assert ds.wave_vectors.shape == (3, 3)
    assert np.allclose(np.diag(ds.wave_vectors), ds.ke)


# --------------------------------------------------------------------- the tab

pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication  # noqa: E402

from harpi.gui.session import Session, Stage  # noqa: E402
from harpi.gui.tabs_data import DataTab  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def tab(qapp):
    sess = Session()
    sess.dataset = _ds2()
    t = DataTab(sess)
    t.pixel_size.setValue(0.002)
    t.pixel_size_col.setValue(0.002)
    t.tag_spacing.setValue(0.008)
    return t


def test_unlinking_lets_the_axes_differ(tab):
    tab.link_size.setChecked(False)
    tab.pixel_size_col.setValue(0.004)
    assert tab.session.dataset.pixel_size == pytest.approx([0.002, 0.004])
    assert not tab.session.dataset.is_isotropic


def test_linked_axes_follow_the_row_value(tab):
    tab.link_size.setChecked(True)
    tab.pixel_size.setValue(0.003)
    assert tab.session.dataset.pixel_size == pytest.approx([0.003, 0.003])


def test_the_column_box_is_disabled_while_linked(tab):
    tab.link_size.setChecked(True)
    assert not tab.pixel_size_col.isEnabled()
    tab.link_size.setChecked(False)
    assert tab.pixel_size_col.isEnabled()


def test_slice_thickness_stays_disabled_for_two_dimensional_data(tab):
    tab.link_size.setChecked(False)
    assert not tab.pixel_size_slice.isEnabled()


def test_the_suggested_shape_parameter_tracks_anisotropy(tab):
    before = tab.a_lbl.text()
    tab.link_size.setChecked(False)
    tab.pixel_size_col.setValue(0.004)
    assert tab.a_lbl.text() != before
    assert float(tab.a_lbl.text()) == pytest.approx(50.0, rel=1e-3)


def test_the_label_shows_both_axes_when_they_differ(tab):
    tab.link_size.setChecked(False)
    tab.pixel_size_col.setValue(0.004)
    assert "rows" in tab.ppw_lbl.text() and "columns" in tab.ppw_lbl.text()


def test_the_label_stays_short_for_square_pixels(tab):
    tab.link_size.setChecked(True)
    tab.pixel_size.setValue(0.002)
    assert "rows" not in tab.ppw_lbl.text()


def test_editing_the_size_invalidates_from_filtering(tab):
    tab.session.set_result(Stage.FILTER, object())
    tab.link_size.setChecked(False)
    tab.pixel_size_col.setValue(0.005)
    assert not tab.session.has(Stage.FILTER)


# ------------------------------------------------------------- frame times

def test_a_dataset_without_times_reports_none():
    assert _ds2().times is None


def test_times_are_trimmed_to_the_frame_count():
    ds = _ds2()
    ds.trigger_times = np.arange(10) * 40.0
    assert ds.times == pytest.approx([0.0, 40.0, 80.0])


def test_too_few_times_are_refused_rather_than_padded():
    ds = _ds2()
    ds.trigger_times = np.array([0.0, 40.0])
    assert ds.times is None


def test_setting_an_interval_creates_even_frame_times(tab):
    tab.frame_interval.setValue(40.0)
    assert tab.session.dataset.times == pytest.approx([0.0, 40.0, 80.0])
    assert "ms" in tab.times_lbl.text()


def test_clearing_the_interval_removes_even_frame_times(tab):
    tab.frame_interval.setValue(40.0)
    tab.frame_interval.setValue(0.0)
    assert tab.session.dataset.times is None
    assert "frame numbers" in tab.times_lbl.text()


def test_the_interval_invalidates_strain(tab):
    tab.session.set_result(Stage.STRAIN, object())
    tab.frame_interval.setValue(40.0)
    assert not tab.session.has(Stage.STRAIN)


def test_uneven_times_from_a_file_are_kept_and_shown(tab):
    """A real trigger-time list is rarely evenly spaced, and must survive."""
    tab.session.dataset.trigger_times = np.array([0.0, 37.0, 81.0])
    tab._refresh_times()
    assert "0 to 81 ms" in tab.times_lbl.text()
    assert tab.session.dataset.times == pytest.approx([0.0, 37.0, 81.0])
