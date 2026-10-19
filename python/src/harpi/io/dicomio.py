# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

"""Reading tagged MR series from DICOM.

Clinical exports arrive in whatever shape the scanner produced. A folder may
hold one series or two hundred, of which only a few are tagged; a cine may be
one file per frame or a single multiframe file; and the tag spacing is usually
absent from the header, with the tag direction not necessarily along the image
axes.

So the reading is in two steps. :func:`scan_dicom` builds a catalogue of what a
path contains, measuring the tag geometry of each series from its harmonic
peaks, and :func:`read_dicom_series` then assembles one chosen series.

Grid tagging, which is what most of these exports use, puts both tag directions
in a single image as two orthogonal harmonic pairs. There is nothing to group by
file in that case: both encodings come from filtering the same image at two
different peaks, so the frames are emitted on both encoding channels and
:class:`~harpi.filters.HarpFilter2D` separates them.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..harmonic import TagGeometry, detect_tag_harmonics
from ..types import Dataset2D

__all__ = ["read_dicom_series", "list_series", "scan_dicom", "SeriesInfo"]


# --------------------------------------------------------------------- reading

def _headers(path):
    """Every readable DICOM header under a path, without the pixel data.

    Reading headers only matters: a folder of a few thousand files would
    otherwise decode every pixel array just to list what is there.
    """
    import pydicom

    path = Path(path)
    paths = [path] if path.is_file() else sorted(p for p in path.rglob("*") if p.is_file())
    out = []
    for p in paths:
        try:
            ds = pydicom.dcmread(str(p), stop_before_pixels=True, force=True)
        except Exception:
            continue
        if "Rows" in ds and "Columns" in ds:
            out.append((p, ds))
    if not out:
        raise FileNotFoundError(f"no readable DICOM images under {path}")
    return out


def _read_files(directory):
    """Headers and pixel data together. Kept for callers that need both."""
    import pydicom

    out = []
    for p, _ in _headers(directory):
        try:
            out.append((p, pydicom.dcmread(str(p), stop_before_pixels=False)))
        except Exception:
            continue
    if not out:
        raise FileNotFoundError(f"no readable DICOM files under {directory}")
    return out


def _frames_of(ds):
    """Pixel frames held by one file, as a list of 2-D arrays."""
    arr = ds.pixel_array.astype(float)
    slope = float(getattr(ds, "RescaleSlope", 1) or 1)
    inter = float(getattr(ds, "RescaleIntercept", 0) or 0)
    arr = arr * slope + inter
    if arr.ndim == 2:
        return [arr]
    # A multiframe file stacks frames on the leading axis.
    return [arr[i] for i in range(arr.shape[0])]


def _sort_key(ds):
    return (float(getattr(ds, "TriggerTime", 0) or 0),
            int(getattr(ds, "InstanceNumber", 0) or 0))


# ------------------------------------------------------------------- catalogue

@dataclass
class SeriesInfo:
    """One series in a DICOM path, and what could be measured about it."""

    uid: str
    description: str = ""
    protocol: str = ""
    series_number: int = 0
    slice_location: float | None = None
    shape: tuple[int, int] = (0, 0)
    pixel_spacing: tuple[float, float] = (1.0, 1.0)
    n_files: int = 0
    n_frames: int = 0
    multiframe: bool = False
    trigger_times: np.ndarray | None = None
    first_path: Path | None = None
    paths: list = field(default_factory=list)
    tag: TagGeometry | None = None
    """Measured tag geometry, or None when detection was skipped or failed."""

    @property
    def looks_tagged(self) -> bool:
        return self.tag is not None and self.tag.plausible

    @property
    def is_cine(self) -> bool:
        return self.n_frames > 1

    def label(self) -> str:
        """A one-line description for a picker or a listing."""
        bits = [self.description or self.protocol or self.uid[-12:]]
        if self.slice_location is not None:
            bits.append(f"slice {self.slice_location:+.1f}")
        bits.append(f"{self.n_frames} frames")
        bits.append(f"{self.shape[0]}x{self.shape[1]}")
        if self.tag is not None and self.tag.plausible:
            ang = np.degrees(self.tag.angles) % 180
            bits.append(f"tags {self.tag.spacing_mm(self.pixel_spacing):.1f} mm "
                        f"at {ang[0]:.0f} and {ang[1]:.0f} deg")
        return ", ".join(bits)


def scan_dicom(path, detect=True) -> list[SeriesInfo]:
    """Catalogue the series under a path, tagged ones first.

    Parameters
    ----------
    path
        A folder, scanned recursively, or a single file.
    detect
        Measure the tag geometry of each series. One frame per series is
        decoded, not one per file, so a folder of thousands stays quick.

    Returns
    -------
    list of :class:`SeriesInfo`, sorted with the likely-tagged series first and
    then by series number, so a caller can offer the useful ones without hiding
    the rest.
    """
    import pydicom

    groups: dict[str, list] = {}
    for p, ds in _headers(path):
        groups.setdefault(str(getattr(ds, "SeriesInstanceUID", "unknown")), []).append((p, ds))

    out = []
    for uid, items in groups.items():
        items.sort(key=lambda t: _sort_key(t[1]))
        first = items[0][1]
        per_file = int(getattr(first, "NumberOfFrames", 1) or 1)
        ps = np.atleast_1d(np.asarray(
            getattr(first, "PixelSpacing", [1.0, 1.0]), dtype=float)).ravel()
        times = np.array([float(getattr(ds, "TriggerTime", np.nan) or np.nan)
                          for _, ds in items], dtype=float)

        info = SeriesInfo(
            uid=uid,
            description=str(getattr(first, "SeriesDescription", "") or ""),
            protocol=str(getattr(first, "ProtocolName", "") or ""),
            series_number=int(getattr(first, "SeriesNumber", 0) or 0),
            slice_location=(float(first.SliceLocation)
                            if getattr(first, "SliceLocation", None) is not None else None),
            shape=(int(first.Rows), int(first.Columns)),
            pixel_spacing=(float(ps[0]), float(ps[1] if ps.size > 1 else ps[0])),
            n_files=len(items),
            n_frames=len(items) * per_file,
            multiframe=per_file > 1,
            trigger_times=times if np.all(np.isfinite(times)) else None,
            first_path=items[0][0],
            paths=[p for p, _ in items],
        )

        if detect:
            try:
                full = pydicom.dcmread(str(info.first_path), stop_before_pixels=False)
                info.tag = detect_tag_harmonics(_frames_of(full)[0])
            except Exception:
                info.tag = None
        out.append(info)

    out.sort(key=lambda s: (not s.looks_tagged, s.series_number, s.description))
    return out


def list_series(directory) -> dict[str, int]:
    """Map SeriesInstanceUID to the number of images.

    Kept for callers that only want the counts. :func:`scan_dicom` says what
    each series actually is.
    """
    return {s.uid: s.n_files for s in scan_dicom(directory, detect=False)}


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


# ---------------------------------------------------------------- assembling

def _stack_series(info, pydicom):
    """Frames of one series in time order, as ``(rows, cols, frames)``."""
    frames = []
    times = []
    for p in info.paths:
        ds = pydicom.dcmread(str(p), stop_before_pixels=False)
        got = _frames_of(ds)
        frames.extend(got)
        t = getattr(ds, "TriggerTime", None)
        base = float(t) if t is not None else np.nan
        # A multiframe file carries one trigger time for the whole stack, so
        # only a single-frame-per-file series can give a time per frame.
        times.extend([base] * len(got))
    return np.stack(frames, axis=-1), np.asarray(times, dtype=float)


def read_dicom_series(directory, series_uid=None, tag_spacing=None,
                      encoding_key=None, catalogue=None) -> Dataset2D:
    """Assemble one tagged series into a :class:`~harpi.types.Dataset2D`.

    Parameters
    ----------
    directory
        Folder to scan recursively, or a single file.
    series_uid
        Which series to load. A single UID for grid tagging, where one image
        holds both directions. A pair of UIDs for line tagging, where each
        series carries one direction. Omitted, the only series is taken, or the
        only tagged one when there are several.
    tag_spacing
        Tag period in the units of PixelSpacing, usually millimetres. Taken
        from the header when present, otherwise measured from the harmonic
        peaks, otherwise guessed with a warning.
    encoding_key
        Accepted and ignored. Grouping files into encodings by header value does
        not work on the exports this reads: every file of a grid series carries
        the same SequenceName, so everything landed in one group and the image
        came out with a single encoding channel.
    catalogue
        A previous :func:`scan_dicom` result, to avoid scanning twice.

    Returns
    -------
    Dataset2D with ``image`` of shape ``(rows, cols, 2, frames)`` and an
    all-true mask. ``meta`` records the tag geometry and where the spacing came
    from.
    """
    import pydicom

    if encoding_key is not None:
        warnings.warn(
            "encoding_key is ignored: grouping files into encodings by header "
            "value does not work on tagged exports, where every file of a grid "
            "series shares one SequenceName.",
            RuntimeWarning, stacklevel=2,
        )

    series = catalogue if catalogue is not None else scan_dicom(directory)
    by_uid = {s.uid: s for s in series}

    if series_uid is None:
        tagged = [s for s in series if s.looks_tagged]
        if len(series) == 1:
            chosen = [series[0]]
        elif len(tagged) == 1:
            chosen = [tagged[0]]
        else:
            raise ValueError(
                f"{len(series)} series found, {len(tagged)} of them tagged. "
                "Pass series_uid=; scan_dicom() lists what each one is."
            )
    elif isinstance(series_uid, (str, bytes)):
        chosen = [by_uid[str(series_uid)]]
    else:
        chosen = [by_uid[str(u)] for u in series_uid]
        if len(chosen) != 2:
            raise ValueError("pass one series for grid tagging, or two for line "
                             f"tagging, not {len(chosen)}")

    stacks, times = [], None
    for info in chosen:
        arr, t = _stack_series(info, pydicom)
        stacks.append(arr)
        if times is None:
            times = t

    if len(stacks) == 1:
        # Grid tagging: both directions live in this one image, and the filter
        # picks them apart at two different harmonic peaks.
        n = stacks[0].shape[-1]
        image = np.stack([stacks[0], stacks[0]], axis=2)
        geometry = chosen[0].tag
    else:
        n = min(a.shape[-1] for a in stacks)
        image = np.stack([a[..., :n] for a in stacks], axis=2)
        geometry = chosen[0].tag
    image = image[:, :, :, :n]

    ref = pydicom.dcmread(str(chosen[0].first_path), stop_before_pixels=True,
                          force=True)
    pixel_spacing = np.asarray(chosen[0].pixel_spacing, dtype=float)

    source = "given"
    if tag_spacing is None:
        tag_spacing = _tag_spacing_from(ref)
        source = "header"
    if tag_spacing is None and geometry is not None and geometry.plausible:
        tag_spacing = geometry.spacing_mm(pixel_spacing)
        source = "measured from the harmonic peaks"
    if tag_spacing is None:
        tag_spacing = float(pixel_spacing[0]) * 7.0
        source = "guessed at 7 pixels per wavelength"
        warnings.warn(
            "tag spacing is not in the header and could not be measured, "
            "assuming 7 pixels per wavelength. Pass tag_spacing= explicitly, "
            "since the shape parameter depends on it.",
            RuntimeWarning, stacklevel=2,
        )

    trigger = (times[:n] if times is not None and times.size >= n
               and np.all(np.isfinite(times[:n])) and np.all(np.diff(times[:n]) > 0)
               else None)

    meta = {
        "source": str(directory),
        "series_uid": [s.uid for s in chosen],
        "series_description": [s.description for s in chosen],
        "tag_spacing_source": source,
        "tagging": "grid" if len(stacks) == 1 else "line",
    }
    if geometry is not None:
        meta["tag_centres"] = geometry.centres
        meta["tag_angles"] = geometry.angles
        meta["tag_confidence"] = geometry.confidence

    return Dataset2D(
        image=image,
        mask=np.ones(image.shape[:2] + (image.shape[-1],), dtype=bool),
        pixel_size=pixel_spacing,
        tag_spacing=float(tag_spacing),
        trigger_times=trigger,
        meta=meta,
    )
