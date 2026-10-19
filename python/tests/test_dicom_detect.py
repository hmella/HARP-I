"""Detection and loading checked against the real clinical exports.

These run against `data/Tagging dataset` and skip when it is absent. The
expected peaks and spacings were measured from the files; the dataset names
carry the nominal tag spacing, which is what makes them worth asserting:
`CSPAMM_GRID_6_0` is a 6 mm grid, `TAG_GRID_7MM_45` a 7 mm grid at 45 degrees.
"""

from __future__ import annotations

import numpy as np
import pytest

pydicom = pytest.importorskip("pydicom")

from pathlib import Path  # noqa: E402

from harpi.harmonic import detect_tag_harmonics  # noqa: E402
from harpi.io.dicomio import read_dicom_series, scan_dicom  # noqa: E402

ROOT = Path(__file__).resolve().parents[2] / "data" / "Tagging dataset"

# name, nominal tag spacing in mm, tag angles in degrees
CASES = [
    ("CSPAMM_GRID_6_0_HCM", 6.0, (0.0, 90.0)),
    ("CSPAMM_GRID_MULTI_Normal", 6.0, (0.0, 90.0)),
    ("TAG_GRID_7MM_45_Normal", 7.0, (45.0, 135.0)),
    ("20211_2_0", 6.0, (47.0, 134.0)),
]


def _folder(name):
    d = ROOT / name
    if not d.is_dir():
        pytest.skip(f"{d} not present")
    return d


def _first_frame(folder):
    p = next(q for q in sorted(folder.rglob("*")) if q.is_file()
             and q.suffix.lower() in (".dcm", ".ima", ""))
    ds = pydicom.dcmread(str(p), force=True)
    return ds.pixel_array.astype(float), np.asarray(ds.PixelSpacing, dtype=float)


@pytest.mark.parametrize("name,nominal,angles", CASES)
def test_the_tag_spacing_is_recovered_from_the_image(name, nominal, angles):
    """No example carries a tag spacing in its header, so it is measured."""
    frame, ps = _first_frame(_folder(name))
    g = detect_tag_harmonics(frame)
    assert g.spacing_mm(ps) == pytest.approx(nominal, abs=0.4)


@pytest.mark.parametrize("name,nominal,angles", CASES)
def test_the_tag_angles_are_recovered(name, nominal, angles):
    """Two of the four are tagged at 45 degrees, not along the image axes."""
    frame, _ = _first_frame(_folder(name))
    g = detect_tag_harmonics(frame)
    got = np.sort(np.degrees(g.angles) % 180)
    assert got == pytest.approx(np.sort(np.asarray(angles)), abs=4.0)


@pytest.mark.parametrize("name,nominal,angles", CASES)
def test_the_pattern_is_called_plausible(name, nominal, angles):
    frame, _ = _first_frame(_folder(name))
    assert detect_tag_harmonics(frame).plausible


@pytest.mark.parametrize("name,nominal,angles", CASES)
def test_the_directions_are_ordered_the_same_for_every_series(name, nominal, angles):
    """Peak strength varies between slices; the component order must not."""
    cat = [s for s in scan_dicom(_folder(name)) if s.looks_tagged]
    assert cat, "no tagged series found"
    first = np.degrees(cat[0].tag.angles) % 180
    for s in cat[1:]:
        assert np.degrees(s.tag.angles) % 180 == pytest.approx(first, abs=4.0)


@pytest.mark.parametrize("name,nominal,angles", CASES)
def test_the_catalogue_describes_every_series(name, nominal, angles):
    cat = scan_dicom(_folder(name))
    assert cat
    for s in cat:
        assert s.n_frames > 1, "a tagged cine has more than one frame"
        assert s.shape[0] > 0 and s.shape[1] > 0
        assert s.label()


@pytest.mark.parametrize("name,nominal,angles", CASES)
def test_the_series_loads_with_two_encoding_channels(name, nominal, angles):
    folder = _folder(name)
    cat = scan_dicom(folder)
    tagged = [s for s in cat if s.looks_tagged]
    ds = read_dicom_series(folder, series_uid=tagged[0].uid, catalogue=cat)
    assert ds.image.ndim == 4
    assert ds.image.shape[2] == 2, "HARP-I needs two encodings"
    assert ds.n_frames == tagged[0].n_frames
    assert ds.tag_spacing == pytest.approx(nominal, abs=0.4)


@pytest.mark.parametrize("name,nominal,angles", CASES)
def test_the_harmonic_filter_can_be_built(name, nominal, angles):
    """This is the step that raised IndexError before: one encoding channel."""
    from harpi import HarpFilter2D

    folder = _folder(name)
    cat = scan_dicom(folder)
    tagged = [s for s in cat if s.looks_tagged]
    ds = read_dicom_series(folder, series_uid=tagged[0].uid, catalogue=cat)
    flt = HarpFilter2D(ds.image, ds.ke, tuple(ds.meta["tag_angles"]),
                       centre=np.round(ds.meta["tag_centres"]).astype(int))
    Ih = flt(ds.image)
    assert Ih.shape == ds.image.shape
    assert np.iscomplexobj(Ih)
    assert np.isfinite(Ih).all()


def test_a_mixed_archive_separates_the_tagged_series():
    """amrgatlas holds 201 series of which 26 are tagged, the rest cine and T1."""
    folder = _folder("amrgatlas")
    cat = scan_dicom(folder)
    assert len(cat) > 100, f"expected a large archive, got {len(cat)}"

    grid = [s for s in cat if "grid" in s.description.lower()]
    other = [s for s in cat if "grid" not in s.description.lower()]
    assert grid, "no grid series found by description"

    missed = [s.description for s in grid if not s.looks_tagged]
    assert not missed, f"tagged series not detected: {missed}"

    false_rate = sum(s.looks_tagged for s in other) / max(len(other), 1)
    assert false_rate < 0.05, f"too many non-tagged series flagged: {false_rate:.2%}"


def test_the_tagged_series_are_offered_first():
    cat = scan_dicom(_folder("amrgatlas"))
    flags = [s.looks_tagged for s in cat]
    # every tagged one comes before every untagged one
    assert flags == sorted(flags, reverse=True)


def test_a_cine_series_is_not_mistaken_for_tagging():
    cat = scan_dicom(_folder("amrgatlas"))
    cines = [s for s in cat if s.description.startswith("tf2d15")]
    assert cines, "no cine series found"
    assert not any(s.looks_tagged for s in cines)


# ------------------------------------------------- through the interface

pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication  # noqa: E402

from harpi.gui.session import Stage  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.mark.parametrize("name,nominal,angles", CASES)
def test_the_data_tab_accepts_a_real_series(qapp, name, nominal, angles):
    from harpi.gui.session import Session
    from harpi.gui.tabs_data import DataTab

    folder = _folder(name)
    cat = scan_dicom(folder)
    tagged = [s for s in cat if s.looks_tagged]
    ds = read_dicom_series(folder, series_uid=tagged[0].uid, catalogue=cat)

    tab = DataTab(Session())
    assert tab._adopt(ds, str(folder)) is True
    assert tab.session.has(Stage.DATA)
    # the tag spacing measured from the images reaches the box
    assert tab.tag_spacing.value() == pytest.approx(ds.tag_spacing, rel=1e-6)


def test_the_filter_tab_seeds_its_centres_from_the_measurement(qapp):
    """A 45 degree grid needs the measured peaks, not the nominal 0 and 90."""
    from harpi.gui.session import Session
    from harpi.gui.tabs_data import DataTab, FilterTab

    folder = _folder("TAG_GRID_7MM_45_Normal")
    cat = scan_dicom(folder)
    tagged = [s for s in cat if s.looks_tagged]
    ds = read_dicom_series(folder, series_uid=tagged[0].uid, catalogue=cat)

    sess = Session()
    DataTab(sess)._adopt(ds, str(folder))
    ft = FilterTab(sess)
    ft.reset_centres()

    measured = np.asarray(ds.meta["tag_centres"], dtype=float)
    assert np.allclose(ft.centres, measured, atol=1.0)
    # the nominal 0 and 90 would put both peaks on the axes through the centre
    rows, cols = ds.shape
    assert not np.isclose(ft.centres[0][0], rows / 2, atol=2.0)


def test_filtering_a_real_series_succeeds(qapp):
    """End to end in the interface, the step that used to raise IndexError."""
    from harpi.gui.session import Session
    from harpi.gui.tabs_data import DataTab, FilterTab

    folder = _folder("CSPAMM_GRID_6_0_HCM")
    cat = scan_dicom(folder)
    tagged = [s for s in cat if s.looks_tagged]
    ds = read_dicom_series(folder, series_uid=tagged[0].uid, catalogue=cat)

    sess = Session()
    DataTab(sess)._adopt(ds, str(folder))
    ft = FilterTab(sess)
    ft.reset_centres()
    ft.run()

    res = sess.result(Stage.FILTER)
    assert res is not None, ft.banner.text()
    assert res["image"].shape == ds.image.shape
    assert np.iscomplexobj(res["image"])
