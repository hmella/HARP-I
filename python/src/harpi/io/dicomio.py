"""Reading tagged MR series from DICOM.

The tag spacing lives in a vendor private field and is often missing, so it
can always be passed in. A warning is raised when it has to be guessed.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np

from ..types import Dataset2D

__all__ = ["read_dicom_series", "list_series"]


def _read_files(directory):
    import pydicom

    paths = sorted(Path(directory).rglob("*"))
    out = []
    for p in paths:
        if not p.is_file():
            continue
        try:
            out.append((p, pydicom.dcmread(str(p), stop_before_pixels=False)))
        except Exception:
            continue
    if not out:
        raise FileNotFoundError(f"no readable DICOM files under {directory}")
    return out


def list_series(directory) -> dict[str, int]:
    """Map SeriesInstanceUID -> number of images, for choosing a series."""
    counts: dict[str, int] = {}
    for _, ds in _read_files(directory):
        uid = getattr(ds, "SeriesInstanceUID", "unknown")
        counts[uid] = counts.get(uid, 0) + 1
    return counts


def _tag_spacing_from(ds) -> float | None:
    for name in ("TagSpacing", "TagThickness"):
        if hasattr(ds, name):
            try:
                return float(getattr(ds, name))
            except (TypeError, ValueError):
                pass
    # Siemens private group
    for elem in ds.iterall() if hasattr(ds, "iterall") else []:
        if elem.name and "tag spacing" in str(elem.name).lower():
            try:
                return float(elem.value)
            except (TypeError, ValueError):
                pass
    return None


def read_dicom_series(directory, series_uid=None, tag_spacing=None,
                      encoding_key=None) -> Dataset2D:
    """Assemble a 2-D tagged series into a :class:`~harpi.types.Dataset2D`.

    Parameters
    ----------
    directory
        Folder to scan recursively.
    series_uid
        Which series to load; required when the folder holds more than one.
    tag_spacing
        Tag period in the same units as PixelSpacing, usually mm. Read from
        the file when omitted, with a warning if it cannot be found.
    encoding_key
        Callable mapping a dataset to its tag-direction index. Defaults to
        grouping by ``SequenceName``.

    Returns
    -------
    Dataset2D with ``image`` of shape ``(rows, cols, encodings, frames)`` and
    an all-true mask (segment it afterwards).
    """
    files = _read_files(directory)

    if series_uid is None:
        uids = {getattr(ds, "SeriesInstanceUID", "unknown") for _, ds in files}
        if len(uids) > 1:
            raise ValueError(
                f"{len(uids)} series found; pass series_uid=. Use list_series() to see them."
            )
        series_uid = uids.pop()
    files = [(p, ds) for p, ds in files
             if getattr(ds, "SeriesInstanceUID", "unknown") == series_uid]
    if not files:
        raise ValueError(f"no images for series {series_uid}")

    if encoding_key is None:
        def encoding_key(ds):
            return getattr(ds, "SequenceName", "") or getattr(ds, "ImageComments", "")

    groups: dict[object, list] = {}
    for p, ds in files:
        groups.setdefault(encoding_key(ds), []).append((p, ds))

    keys = sorted(groups, key=str)
    per_encoding = []
    for k in keys:
        items = groups[k]
        items.sort(key=lambda t: (
            float(getattr(t[1], "TriggerTime", 0) or 0),
            int(getattr(t[1], "InstanceNumber", 0) or 0),
        ))
        frames = []
        for _, ds in items:
            arr = ds.pixel_array.astype(float)
            slope = float(getattr(ds, "RescaleSlope", 1) or 1)
            inter = float(getattr(ds, "RescaleIntercept", 0) or 0)
            frames.append(arr * slope + inter)
        per_encoding.append(np.stack(frames, axis=-1))

    n = min(a.shape[-1] for a in per_encoding)
    image = np.stack([a[..., :n] for a in per_encoding], axis=2)

    ref = files[0][1]
    pixel_spacing = np.asarray(getattr(ref, "PixelSpacing", [1.0, 1.0]), dtype=float)

    if tag_spacing is None:
        tag_spacing = _tag_spacing_from(ref)
        if tag_spacing is None:
            tag_spacing = float(pixel_spacing[0]) * 7.0
            warnings.warn(
                "tag spacing not found in the DICOM header, assuming 7 pixels "
                "per wavelength. Pass tag_spacing= explicitly, since the shape "
                "parameter depends on it.",
                RuntimeWarning,
                stacklevel=2,
            )

    return Dataset2D(
        image=image,
        mask=np.ones(image.shape[:2] + (image.shape[-1],), dtype=bool),
        pixel_size=pixel_spacing,
        tag_spacing=float(tag_spacing),
        meta={"source": str(directory), "series_uid": series_uid,
              "encodings": [str(k) for k in keys]},
    )
