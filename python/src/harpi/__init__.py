"""HARP-I: motion and strain from tagged MR images.

A Python implementation of the harmonic phase interpolation method of

    H. Mella, J. Mura, H. Sanchez, S. Uribe, C. Bertoglio, J. Sotelo,
    "HARP-I: A Harmonic Phase Interpolation Method for the Estimation of
    Motion From Tagged MR Images", IEEE Transactions on Medical Imaging
    40(4):1240-1252, 2021. doi:10.1109/TMI.2021.3051092

Typical use::

    import harpi

    ds  = harpi.load_mat_2d("data/data_2D.mat")
    flt = harpi.HarpFilter2D(ds.image, ds.ke, kind="transmission")
    res = harpi.estimate_motion_2d(flt(ds.image), ds.mask,
                                   a=ds.suggested_a, eta=15.0,
                                   smooth_reference=True, temporal_fit=True)
    st  = harpi.pixel_strain_2d(res.best, ds.mask[..., 0])
    gcs = st.global_curve("CC")

This package accompanies the MATLAB toolbox at
https://github.com/hmella/HARP-I and is released under the same licence
(MPL-2.0). The phase unwrapping and contour modules also derive from
DENSEanalysis, as noted in their module docstrings.
"""

from .config import FilterConfig, HarpiConfig, MotionConfig, StrainConfig
from .consistency import temporal_phase_consistency
from .fft import itok, ktoi
from .filters import (
    HarpFilter2D,
    HarpFilter3D,
    centre_to_frequency,
    frequency_to_centre,
)
from .harmonic import harmonic_center
from .interp import rbf_interpolate, solve_tikhonov
from .io import HarpiProject, load_mat_2d, load_mat_3d, save_mat
from .motion import estimate_motion_2d, estimate_motion_3d, flip_components
from .rbf import KERNELS, RBFKernel, get_kernel
from .segmentation import Contour, contours_to_mask
from .smoothing import smooth_reference_phase
from .analysis import (aha_labels, bullseye_rings, segment_metrics,
                       segment_table)
from .strain import pixel_strain_2d, pixel_strain_3d, strain_by_segments
from .temporal import fit_temporal
from .types import Dataset2D, Dataset3D, MotionResult, SegmentStrain, StrainResult
from .unwrap import phase_quality_2d, phase_quality_3d, unwrap2, unwrap3

__version__ = "0.1.0"

__all__ = [
    "__version__",
    # data
    "Dataset2D", "Dataset3D", "MotionResult", "StrainResult", "SegmentStrain",
    "load_mat_2d", "load_mat_3d", "save_mat", "HarpiProject",
    # config
    "HarpiConfig", "FilterConfig", "MotionConfig", "StrainConfig",
    # pipeline
    "itok", "ktoi",
    "HarpFilter2D", "HarpFilter3D", "harmonic_center",
    "centre_to_frequency", "frequency_to_centre",
    "phase_quality_2d", "phase_quality_3d", "unwrap2", "unwrap3",
    "temporal_phase_consistency", "smooth_reference_phase",
    "rbf_interpolate", "solve_tikhonov", "get_kernel", "KERNELS", "RBFKernel",
    "estimate_motion_2d", "estimate_motion_3d", "fit_temporal",
    "flip_components",
    "pixel_strain_2d", "pixel_strain_3d", "strain_by_segments",
    # regional analysis
    "aha_labels", "segment_metrics", "bullseye_rings", "segment_table",
    "Contour", "contours_to_mask",
]
