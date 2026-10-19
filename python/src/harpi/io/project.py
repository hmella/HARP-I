"""Project files.

A ``.harpi`` file is a zip holding ``meta.json`` and ``arrays/*.npy``. Image
data is referenced by path and checksum rather than embedded, unless
``embed=True``.
"""

from __future__ import annotations

import hashlib
import io
import json
import zipfile
from dataclasses import asdict, dataclass, field, is_dataclass
from pathlib import Path

import numpy as np

__all__ = ["HarpiProject"]

SCHEMA_VERSION = 1


def _sha256(path, chunk=1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


@dataclass
class HarpiProject:
    """Everything needed to reproduce a session."""

    source: str = ""
    source_sha256: str = ""
    config: dict = field(default_factory=dict)
    contours: dict = field(default_factory=dict)
    """``{"slice,frame": [ {points, closed, curved, corner, label}, ... ]}``."""
    arrays: dict = field(default_factory=dict)
    notes: str = ""

    # ---------------------------------------------------------------- save
    def save(self, path, embed=False):
        path = Path(path)
        meta = {
            "schema_version": SCHEMA_VERSION,
            "source": self.source,
            "source_sha256": self.source_sha256,
            "config": _jsonable(self.config),
            "contours": _jsonable(self.contours),
            "notes": self.notes,
            "arrays": sorted(self.arrays),
        }
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("meta.json", json.dumps(meta, indent=2))
            for name, arr in self.arrays.items():
                if not embed and name.startswith("image"):
                    continue
                buf = io.BytesIO()
                np.save(buf, np.asarray(arr), allow_pickle=False)
                z.writestr(f"arrays/{name}.npy", buf.getvalue())
        return path

    # ---------------------------------------------------------------- load
    @classmethod
    def load(cls, path):
        path = Path(path)
        with zipfile.ZipFile(path) as z:
            meta = json.loads(z.read("meta.json"))
            if meta.get("schema_version", 1) > SCHEMA_VERSION:
                raise ValueError(
                    f"{path} was written by a newer version of harpi "
                    f"(schema {meta['schema_version']} > {SCHEMA_VERSION})"
                )
            arrays = {}
            for item in z.namelist():
                if item.startswith("arrays/") and item.endswith(".npy"):
                    arrays[item[len("arrays/"):-4]] = np.load(
                        io.BytesIO(z.read(item)), allow_pickle=False
                    )
        return cls(
            source=meta.get("source", ""),
            source_sha256=meta.get("source_sha256", ""),
            config=meta.get("config", {}),
            contours=meta.get("contours", {}),
            arrays=arrays,
            notes=meta.get("notes", ""),
        )

    # ------------------------------------------------------------- helpers
    def set_source(self, path):
        path = Path(path)
        self.source = str(path)
        self.source_sha256 = _sha256(path)

    def source_unchanged(self) -> bool | None:
        """True or False, or None when the source file is missing."""
        if not self.source or not self.source_sha256:
            return None
        p = Path(self.source)
        if not p.exists():
            return None
        return _sha256(p) == self.source_sha256


def _jsonable(obj):
    if is_dataclass(obj) and not isinstance(obj, type):
        return _jsonable(asdict(obj))
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj
