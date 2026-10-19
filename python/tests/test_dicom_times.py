"""Frame times taken from DICOM, which the reader previously discarded."""

from __future__ import annotations

import numpy as np
import pytest

pydicom = pytest.importorskip("pydicom")

from pydicom.dataset import Dataset, FileMetaDataset  # noqa: E402
from pydicom.uid import ExplicitVRLittleEndian, generate_uid  # noqa: E402

from harpi.io.dicomio import read_dicom_series  # noqa: E402


def _write_series(directory, trigger_times, encodings=("tag1", "tag2"),
                  rows=8, cols=8):
    """A minimal two-encoding tagged series, one file per frame per encoding."""
    series_uid = generate_uid()
    for enc_index, name in enumerate(encodings):
        for frame, t in enumerate(trigger_times):
            ds = Dataset()
            ds.file_meta = FileMetaDataset()
            ds.file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
            ds.file_meta.MediaStorageSOPClassUID = generate_uid()
            ds.file_meta.MediaStorageSOPInstanceUID = generate_uid()
            ds.SOPInstanceUID = ds.file_meta.MediaStorageSOPInstanceUID
            ds.SeriesInstanceUID = series_uid
            ds.SequenceName = name
            ds.InstanceNumber = frame + 1
            if t is not None:
                ds.TriggerTime = float(t)
            ds.PixelSpacing = [2.0, 2.0]
            ds.Rows, ds.Columns = rows, cols
            ds.BitsAllocated = 16
            ds.BitsStored = 16
            ds.HighBit = 15
            ds.PixelRepresentation = 0
            ds.SamplesPerPixel = 1
            ds.PhotometricInterpretation = "MONOCHROME2"
            arr = np.full((rows, cols), frame + enc_index, dtype=np.uint16)
            ds.PixelData = arr.tobytes()
            ds.save_as(directory / f"{name}_{frame:03d}.dcm",
                       enforce_file_format=True)
    return series_uid


def test_trigger_times_are_kept(tmp_path):
    _write_series(tmp_path, [0.0, 40.0, 80.0, 120.0])
    ds = read_dicom_series(tmp_path, tag_spacing=14.0)
    assert ds.trigger_times is not None
    assert ds.times == pytest.approx([0.0, 40.0, 80.0, 120.0])


def test_uneven_trigger_times_are_kept(tmp_path):
    """Real trigger times are rarely evenly spaced."""
    _write_series(tmp_path, [0.0, 37.0, 81.0, 118.0])
    ds = read_dicom_series(tmp_path, tag_spacing=14.0)
    assert ds.times == pytest.approx([0.0, 37.0, 81.0, 118.0])


def test_a_missing_trigger_time_discards_the_whole_set(tmp_path):
    """A partial set would put the wrong numbers on every time axis."""
    _write_series(tmp_path, [0.0, 40.0, None, 120.0])
    ds = read_dicom_series(tmp_path, tag_spacing=14.0)
    assert ds.trigger_times is None
    assert ds.times is None


def test_non_increasing_trigger_times_are_discarded(tmp_path):
    _write_series(tmp_path, [0.0, 80.0, 80.0, 120.0])
    ds = read_dicom_series(tmp_path, tag_spacing=14.0)
    assert ds.trigger_times is None


def test_the_frame_count_still_matches_the_image(tmp_path):
    _write_series(tmp_path, [0.0, 40.0, 80.0])
    ds = read_dicom_series(tmp_path, tag_spacing=14.0)
    assert ds.n_frames == 3
    assert ds.times.size == ds.n_frames
