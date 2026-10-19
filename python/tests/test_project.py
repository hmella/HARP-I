import numpy as np

from harpi.io import HarpiProject


def test_round_trip(tmp_path):
    p = HarpiProject(
        config={"motion": {"a": 37.5, "eta": 15.0}},
        contours={"0,0": [{"points": [[1, 2], [3, 4]], "closed": True}]},
        arrays={"mask": np.array([[True, False], [False, True]])},
        notes="hello",
    )
    f = p.save(tmp_path / "s.harpi")
    q = HarpiProject.load(f)
    assert q.config["motion"]["a"] == 37.5
    assert q.contours["0,0"][0]["closed"] is True
    assert np.array_equal(q.arrays["mask"], p.arrays["mask"])
    assert q.notes == "hello"


def test_source_checksum_detects_change(tmp_path):
    src = tmp_path / "data.bin"
    src.write_bytes(b"original")
    p = HarpiProject()
    p.set_source(src)
    assert p.source_unchanged() is True
    src.write_bytes(b"modified")
    assert p.source_unchanged() is False


def test_missing_source_reports_none(tmp_path):
    src = tmp_path / "gone.bin"
    src.write_bytes(b"x")
    p = HarpiProject()
    p.set_source(src)
    src.unlink()
    assert p.source_unchanged() is None
