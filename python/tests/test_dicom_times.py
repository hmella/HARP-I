"""Synthetic DICOM: frame times, multiframe files and line-tagged pairs.

The real examples under ``data/Tagging dataset`` are all grid tagged with one
file per frame, so the cases they do not cover are built here.
"""

from __future__ import annotations

import numpy as np
import pytest

pydicom = pytest.importorskip("pydicom")

from pydicom.dataset import Dataset, FileMetaDataset  # noqa: E402
from pydicom.uid import ExplicitVRLittleEndian, generate_uid  # noqa: E402

from harpi.io.dicomio import read_dicom_series, scan_dicom  # noqa: E402


def _base(rows, cols):
    ds = Dataset()
    ds.file_meta = FileMetaDataset()
    ds.file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    ds.file_meta.MediaStorageSOPClassUID = generate_uid()
    ds.file_meta.MediaStorageSOPInstanceUID = generate_uid()
    ds.SOPInstanceUID = ds.file_meta.MediaStorageSOPInstanceUID
    ds.PixelSpacing = [2.0, 2.0]
    ds.Rows, ds.Columns = rows, cols
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 0
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    return ds


def _grid(rows, cols, period=6.0, angle=0.0, phase=0.0):
    """An image with two orthogonal tag harmonics, so detection has something."""
    yy, xx = np.mgrid[0:rows, 0:cols]
    k = 2 * np.pi / period
    u = xx * np.cos(angle) + yy * np.sin(angle)
    v = -xx * np.sin(angle) + yy * np.cos(angle)
    img = 2000 + 700 * np.cos(k * u + phase) + 700 * np.cos(k * v + phase)
    return img.astype(np.uint16)


def write_series(directory, trigger_times, description="tagged", rows=64, cols=64,
                 angle=0.0, series_uid=None, prefix="s"):
    """One series, one file per frame."""
    series_uid = series_uid or generate_uid()
    for frame, t in enumerate(trigger_times):
        ds = _base(rows, cols)
        ds.SeriesInstanceUID = series_uid
        ds.SeriesDescription = description
        ds.SeriesNumber = 1
        ds.InstanceNumber = frame + 1
        if t is not None:
            ds.TriggerTime = float(t)
        ds.PixelData = _grid(rows, cols, angle=angle,
                             phase=0.3 * frame).tobytes()
        ds.save_as(directory / f"{prefix}_{frame:03d}.dcm", enforce_file_format=True)
    return series_uid


def write_multiframe(directory, n_frames, rows=64, cols=64, trigger=None):
    """One file holding a whole cine."""
    ds = _base(rows, cols)
    ds.SeriesInstanceUID = generate_uid()
    ds.SeriesDescription = "multiframe cine"
    ds.SeriesNumber = 1
    ds.InstanceNumber = 1
    ds.NumberOfFrames = n_frames
    if trigger is not None:
        ds.TriggerTime = float(trigger)
    stack = np.stack([_grid(rows, cols, phase=0.3 * i) for i in range(n_frames)])
    ds.PixelData = stack.tobytes()
    ds.save_as(directory / "cine.dcm", enforce_file_format=True)
    return ds.SeriesInstanceUID


# ------------------------------------------------------------- frame times

def test_trigger_times_are_kept(tmp_path):
    write_series(tmp_path, [0.0, 40.0, 80.0, 120.0])
    ds = read_dicom_series(tmp_path)
    assert ds.trigger_times is not None
    assert ds.times == pytest.approx([0.0, 40.0, 80.0, 120.0])


def test_uneven_trigger_times_are_kept(tmp_path):
    """Real trigger times are rarely evenly spaced."""
    write_series(tmp_path, [0.0, 37.0, 81.0, 118.0])
    ds = read_dicom_series(tmp_path)
    assert ds.times == pytest.approx([0.0, 37.0, 81.0, 118.0])


def test_a_missing_trigger_time_discards_the_whole_set(tmp_path):
    """A partial set would put the wrong numbers on every time axis."""
    write_series(tmp_path, [0.0, 40.0, None, 120.0])
    ds = read_dicom_series(tmp_path)
    assert ds.trigger_times is None


def test_non_increasing_trigger_times_are_discarded(tmp_path):
    write_series(tmp_path, [0.0, 80.0, 80.0, 120.0])
    ds = read_dicom_series(tmp_path)
    assert ds.trigger_times is None


def test_the_frame_count_matches_the_image(tmp_path):
    write_series(tmp_path, [0.0, 40.0, 80.0])
    ds = read_dicom_series(tmp_path)
    assert ds.n_frames == 3
    assert ds.times.size == ds.n_frames


# --------------------------------------------------------------- multiframe

def test_a_multiframe_file_is_one_series_of_many_frames(tmp_path):
    write_multiframe(tmp_path, 12)
    cat = scan_dicom(tmp_path)
    assert len(cat) == 1
    assert cat[0].multiframe
    assert cat[0].n_files == 1
    assert cat[0].n_frames == 12


def test_a_multiframe_file_loads(tmp_path):
    write_multiframe(tmp_path, 10)
    ds = read_dicom_series(tmp_path)
    assert ds.image.shape == (64, 64, 2, 10)
    assert ds.n_frames == 10


def test_a_multiframe_file_can_be_opened_directly(tmp_path):
    """Pointing at the file, not the folder, which the folder dialog cannot do."""
    write_multiframe(tmp_path, 8)
    ds = read_dicom_series(tmp_path / "cine.dcm")
    assert ds.image.shape == (64, 64, 2, 8)


def test_one_trigger_time_for_a_stack_is_not_used_per_frame(tmp_path):
    """A multiframe file carries one time, which cannot time ten frames."""
    write_multiframe(tmp_path, 10, trigger=0.0)
    ds = read_dicom_series(tmp_path)
    assert ds.trigger_times is None


# ------------------------------------------------------------- grid tagging

def test_a_grid_series_fills_both_encoding_channels(tmp_path):
    """One image holds both directions, so both channels carry it."""
    write_series(tmp_path, [0.0, 40.0, 80.0])
    ds = read_dicom_series(tmp_path)
    assert ds.image.shape[2] == 2
    assert np.array_equal(ds.image[:, :, 0, :], ds.image[:, :, 1, :])
    assert ds.meta["tagging"] == "grid"


def test_the_filter_can_be_built_from_a_grid_series(tmp_path):
    """The step that used to raise IndexError on a single encoding channel."""
    from harpi import HarpFilter2D

    write_series(tmp_path, [0.0, 40.0, 80.0])
    ds = read_dicom_series(tmp_path)
    flt = HarpFilter2D(ds.image, ds.ke, tuple(ds.meta["tag_angles"]),
                       centre=np.round(ds.meta["tag_centres"]).astype(int))
    assert flt(ds.image).shape == ds.image.shape


def test_the_tag_spacing_is_measured_when_the_header_lacks_it(tmp_path):
    write_series(tmp_path, [0.0, 40.0], rows=64, cols=64)
    ds = read_dicom_series(tmp_path)
    assert "measured" in ds.meta["tag_spacing_source"]
    # the synthetic period is 6 px at 2 mm per pixel
    assert ds.tag_spacing == pytest.approx(12.0, rel=0.15)


def test_an_explicit_tag_spacing_wins(tmp_path):
    write_series(tmp_path, [0.0, 40.0])
    ds = read_dicom_series(tmp_path, tag_spacing=9.0)
    assert ds.tag_spacing == pytest.approx(9.0)
    assert ds.meta["tag_spacing_source"] == "given"


# ------------------------------------------------------------- line tagging

def test_two_series_can_be_paired_as_two_encodings(tmp_path):
    """Line tagging: one series per direction. No real example to check against."""
    a = write_series(tmp_path, [0.0, 40.0, 80.0], description="tag horizontal",
                     angle=0.0, prefix="h")
    b = write_series(tmp_path, [0.0, 40.0, 80.0], description="tag vertical",
                     angle=np.pi / 2, prefix="v")
    ds = read_dicom_series(tmp_path, series_uid=(a, b))
    assert ds.image.shape == (64, 64, 2, 3)
    assert ds.meta["tagging"] == "line"


def test_pairing_more_than_two_series_is_refused(tmp_path):
    a = write_series(tmp_path, [0.0], prefix="a")
    b = write_series(tmp_path, [0.0], prefix="b")
    c = write_series(tmp_path, [0.0], prefix="c")
    with pytest.raises(ValueError, match="one series for grid"):
        read_dicom_series(tmp_path, series_uid=(a, b, c))


# ----------------------------------------------------------------- choosing

def test_several_series_without_a_choice_is_refused(tmp_path):
    write_series(tmp_path, [0.0, 40.0], description="one", prefix="a")
    write_series(tmp_path, [0.0, 40.0], description="two", prefix="b")
    with pytest.raises(ValueError, match="scan_dicom"):
        read_dicom_series(tmp_path)


def test_a_named_series_is_loaded_from_a_folder_of_several(tmp_path):
    a = write_series(tmp_path, [0.0, 40.0, 80.0], description="wanted", prefix="a")
    write_series(tmp_path, [0.0, 40.0], description="other", prefix="b")
    ds = read_dicom_series(tmp_path, series_uid=a)
    assert ds.n_frames == 3
    assert ds.meta["series_description"] == ["wanted"]


def test_a_non_dicom_file_in_the_folder_is_ignored(tmp_path):
    write_series(tmp_path, [0.0, 40.0])
    (tmp_path / "manifest.cvs").write_text("not dicom at all")
    ds = read_dicom_series(tmp_path)
    assert ds.n_frames == 2


def test_the_ignored_encoding_key_warns_rather_than_misleading(tmp_path):
    write_series(tmp_path, [0.0, 40.0])
    with pytest.warns(RuntimeWarning, match="encoding_key is ignored"):
        read_dicom_series(tmp_path, encoding_key=lambda ds: "x")
