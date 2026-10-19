# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy was not distributed with this file, You can
# obtain one at https://mozilla.org/MPL/2.0/.

"""Contour handling and mask generation."""

from .mask import (WALLS, Contour, contours_to_mask, structure_masks,
                   wall_label, wall_role, wall_structure)
from .spline import resample_contour

__all__ = [
    "Contour",
    "contours_to_mask", "structure_masks",
    "resample_contour",
    "WALLS",
    "wall_label",
    "wall_role",
    "wall_structure",
]
