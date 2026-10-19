# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

"""Reading and writing the toolbox's MATLAB files."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scipy.io

from ..types import Dataset2D, Dataset3D

__all__ = ["load_mat_2d", "load_mat_3d", "save_mat", "load_reference"]


def _is_hdf5(path) -> bool:
    with open(path, "rb") as fh:
        return fh.read(8) == b"\x89HDF\r\n\x1a\n"


def _h5(path):
    import h5py

    return h5py.File(path, "r")


def load_mat_2d(path, variable="data") -> Dataset2D:
    """Load a 2-D dataset in the layout of ``data/data_2D.mat``.

    The mask is stored as uint8 and cast to bool here.
    """
    path = Path(path)
    if _is_hdf5(path):
        with _h5(path) as f:
            g = f[variable]
            image = np.array(g["image"]).T
            mask = np.array(g["mask"]).T
            pixel_size = np.array(g["pixelsize"]).ravel()
            tag_spacing = float(np.array(g["tagspacing"]).ravel()[0])
    else:
        m = scipy.io.loadmat(path, struct_as_record=False, squeeze_me=True)
        d = m[variable]
        image = np.asarray(d.image)
        mask = np.asarray(d.mask)
        pixel_size = np.atleast_1d(np.asarray(d.pixelsize, dtype=float)).ravel()
        tag_spacing = float(np.atleast_1d(d.tagspacing).ravel()[0])

    return Dataset2D(
        image=image,
        mask=mask.astype(bool),
        pixel_size=pixel_size,
        tag_spacing=tag_spacing,
        meta={"source": str(path)},
    )


def load_mat_3d(path, reorder=(2, 0, 1)) -> Dataset3D:
    """Load a 3-D dataset in the layout of ``data/data_3D.mat``.

    Parameters
    ----------
    reorder
        Permutation applied to the encoding axis so the tag directions come
        out as (x, y, z). Pass ``None`` to keep the order in the file.

    Notes
    -----
    The z voxel size is taken from ``PixelSpacing[0]``, as in demo_3D.m. The
    ``SliceThickness`` field describes the acquisition slab rather than the
    reconstructed isotropic grid.
    """
    path = Path(path)
    if _is_hdf5(path):
        with _h5(path) as f:
            image = np.array(f["I"]).T
            mask = np.array(f["mask"]).T
            info = {k: np.array(v).ravel() for k, v in f["info"].items()}
            pixel_spacing = np.asarray(info["PixelSpacing"], dtype=float).ravel()
            tag_spacing = float(np.asarray(info["TagSpacing"]).ravel()[0])
    else:
        m = scipy.io.loadmat(path, struct_as_record=False, squeeze_me=True)
        image = np.asarray(m["I"])
        mask = np.asarray(m["mask"])
        info = m["info"]
        pixel_spacing = np.atleast_1d(np.asarray(info.PixelSpacing, dtype=float)).ravel()
        tag_spacing = float(np.atleast_1d(info.TagSpacing).ravel()[0])

    if reorder is not None:
        image = image[:, :, :, list(reorder), :]

    voxel_size = np.array([pixel_spacing[0], pixel_spacing[1], pixel_spacing[0]])
    return Dataset3D(
        image=image,
        mask=mask.astype(bool),
        voxel_size=voxel_size,
        tag_spacing=tag_spacing,
        meta={"source": str(path), "reorder": reorder},
    )


def save_mat(path, motion=None, strain=None, **extra):
    """Write results using the field names of the MATLAB toolbox.

    ``examples/strain_curves2D.m`` reads this output unchanged.
    """
    out = dict(extra)
    if motion is not None:
        out["harpi"] = {
            "RawMotion": motion.u,
            "FittedMotion": np.array([]) if motion.u_fitted is None else motion.u_fitted,
            "RawTrajectories": np.array([]) if motion.trajectories is None
            else motion.trajectories,
            "UnwrappedPhases": np.array([]) if motion.unwrapped_phase is None
            else motion.unwrapped_phase,
        }
    if strain is not None:
        nfr = None
        st = {"maskimage": strain.mask}
        for name in ("XX", "YY", "XY", "RR", "CC", "RC", "LL", "p1", "p2", "p1or"):
            v = getattr(strain, name, None)
            if v is None:
                continue
            nfr = v.shape[-1]
            # (npix, nframes) over maskimage, in column-major order.
            mF = strain.mask.ravel(order="F")
            st[name] = np.stack(
                [v[..., fr].ravel(order="F")[mF] for fr in range(nfr)], axis=1
            )
        out["st"] = st
    scipy.io.savemat(str(path), out)
    return path


def load_reference(path) -> dict:
    """Load a ``capture_reference`` snapshot (v7.3) as plain arrays."""
    import h5py

    out = {}
    with h5py.File(path, "r") as f:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                a = np.array(obj)
                if a.dtype.names and set(a.dtype.names) >= {"real", "imag"}:
                    a = a["real"] + 1j * a["imag"]
                out[name] = a.T
        f.visititems(visit)
    return out
