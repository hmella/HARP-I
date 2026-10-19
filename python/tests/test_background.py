"""The tagged image that the tabs draw their own layers over.

The choice used to be written out inside the Segmentation tab. Motion and
Strain needed the same thing, and a third copy is how the three drift apart,
so it lives in one place now.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("PySide6")

from harpi import Dataset2D  # noqa: E402
from harpi.gui.widgets import (  # noqa: E402
    BACKGROUND_CHOICES,
    background_image,
)

ROWS = COLS = 24
FRAMES = 4


@pytest.fixture
def dataset():
    rng = np.random.default_rng(0)
    image = rng.normal(size=(ROWS, COLS, 2, FRAMES))
    return Dataset2D(image=image, mask=np.ones((ROWS, COLS, FRAMES), bool),
                     pixel_size=[2.0, 2.0], tag_spacing=8.0)


@pytest.fixture
def filtered(dataset):
    rng = np.random.default_rng(1)
    shape = (ROWS, COLS, 2, FRAMES)
    return {"image": rng.normal(size=shape) + 1j * rng.normal(size=shape)}


def test_every_choice_has_its_own_title(dataset, filtered):
    titles = {background_image(dataset, filtered, key, 0).title
              for _, key in BACKGROUND_CHOICES}
    assert len(titles) == len(BACKGROUND_CHOICES), titles


def test_each_choice_is_reported_back(dataset, filtered):
    for _, key in BACKGROUND_CHOICES:
        assert background_image(dataset, filtered, key, 0).choice == key


def test_a_harmonic_choice_falls_back_before_the_filter_has_run(dataset):
    """And says that it did, so the tab can explain itself."""
    bg = background_image(dataset, None, "harmonic", 0)
    assert bg.choice == "magnitude"
    assert bg.title == "magnitude"
    assert np.array_equal(bg.array, np.real(dataset.image[:, :, 0, 0]))


def test_the_raw_magnitude_needs_no_filter(dataset):
    for key in ("magnitude", "magnitude2"):
        assert background_image(dataset, None, key, 0).choice == key


def test_phase_is_drawn_on_a_cyclic_map_over_a_full_turn(dataset, filtered):
    bg = background_image(dataset, filtered, "phase2", 0)
    assert bg.cmap == "twilight"
    assert (bg.vmin, bg.vmax) == (-np.pi, np.pi)
    assert np.allclose(bg.array, np.angle(filtered["image"][:, :, 1, 0]))


def test_the_grey_images_carry_no_colour_map(dataset, filtered):
    """ImageCanvas wants None there; a matplotlib caller substitutes grey."""
    for key in ("magnitude", "magnitude2", "harmonic"):
        assert background_image(dataset, filtered, key, 0).cmap is None


def test_the_second_encoding_clamps_on_single_encoding_data(filtered):
    """Asking for encoding 2 of a one-encoding dataset used to index past the
    end."""
    one = Dataset2D(image=np.zeros((ROWS, COLS, 1, FRAMES)),
                    mask=np.ones((ROWS, COLS, FRAMES), bool),
                    pixel_size=[2.0, 2.0], tag_spacing=8.0)
    bg = background_image(one, None, "magnitude2", 0)
    assert bg.array.shape == (ROWS, COLS)


def test_the_frame_clamps(dataset):
    assert background_image(dataset, None, "magnitude", 99).array.shape \
        == (ROWS, COLS)


def test_a_volume_takes_the_slice_it_is_given():
    """A 5-D image has a slice axis, which the inlined version ignored."""
    image = np.zeros((ROWS, COLS, 3, 2, FRAMES))
    image[:, :, 2, 0, 0] = 7.0
    ds = Dataset2D(image=image, mask=np.ones((ROWS, COLS, FRAMES), bool),
                   pixel_size=[2.0, 2.0], tag_spacing=8.0)
    bg = background_image(ds, None, "magnitude", 0, slice_index=2)
    assert np.allclose(bg.array, 7.0)
    assert np.allclose(background_image(ds, None, "magnitude", 0, 0).array, 0.0)


def test_no_dataset_gives_nothing(filtered):
    """Which is how the Strain tab stays usable with only a motion result."""
    assert background_image(None, filtered, "magnitude", 0) is None
