"""Shared fixtures.

There are two kinds of test here. Reference tests compare against arrays
produced by the MATLAB toolbox with tools/regression/capture_reference.m.
Each stage is given the MATLAB input for that stage, so a difference points
at one step instead of accumulating. The remaining tests use synthetic
phantoms with a known deformation.

The reference file is not committed. Generate it with::

    matlab -batch "maxNumCompThreads(1); \
      addpath(genpath('<repo>/src')); addpath('<repo>/tools'); harpi_paths; \
      addpath('<repo>/tools/regression'); \
      capture_reference('Output','<repo>/python/tests/data/ref_2D.mat')"
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
DATA = Path(__file__).parent / "data"

# Where a reference may live, in order of preference.
_REF_CANDIDATES = [
    DATA / "ref_2D.mat",
    REPO / "tools" / "regression" / "refs" / "py_ref_2D.mat",
]
_DATASET = REPO / "data" / "data_2D.mat"


def _find_reference():
    for p in _REF_CANDIDATES:
        if p.exists():
            return p
    return None


@pytest.fixture(scope="session")
def reference():
    """MATLAB reference arrays, or skip if none is present."""
    path = _find_reference()
    if path is None:
        pytest.skip(
            "no MATLAB reference found; see the note in conftest.py for how to "
            "generate one"
        )
    from harpi.io import load_reference

    return load_reference(path)


@pytest.fixture(scope="session")
def dataset():
    """The 2-D dataset, or skip if it is not downloaded."""
    if not _DATASET.exists():
        pytest.skip(f"{_DATASET} not present (see data/README.md)")
    from harpi import load_mat_2d

    return load_mat_2d(_DATASET)


@pytest.fixture(scope="session")
def phantom():
    from harpi.testing import make_phantom

    return make_phantom(shape=(96, 96), n_frames=10, tag_spacing=8.0,
                        r_inner=16, r_outer=30, contraction=0.15, twist=0.25,
                        t1=None, noise=0.0)


def masked_columns(arr: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Extract ``(npix, nframes)`` in column-major mask order."""
    nfr = arr.shape[-1]
    flat = mask.ravel(order="F")
    return np.stack([arr[..., f].ravel(order="F")[flat] for f in range(nfr)], axis=1)
