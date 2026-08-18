# SPDX-FileCopyrightText: 2026-present Donovan Murphy <domurphy@mcw.edu>
#
# SPDX-License-Identifier: MIT
"""Shared pytest fixtures."""

import numpy as np
import pytest
from PIL import Image


@pytest.fixture
def rng():
    """Seeded generator so the synthetic slide is identical every run."""
    return np.random.default_rng(20260101)


@pytest.fixture
def synthetic_slide(rng):
    """A small fake H&E image: pink background, purple nuclei, red vessels.

    Includes a bright band and a black border strip so the ECF threshold and
    the black-edge fill both have something to act on. The top-left pixel is
    deliberately NOT black, because the fill value is read from it.
    """
    height = width = 600
    image = np.zeros((height, width, 3), np.uint8)
    image[:] = (222, 160, 200)

    for _ in range(400):
        row = int(rng.integers(20, height - 20))
        col = int(rng.integers(20, width - 20))
        image[row - 3 : row + 4, col - 3 : col + 4] = (95, 60, 140)

    for _ in range(10):
        row = int(rng.integers(60, height - 60))
        col = int(rng.integers(60, width - 60))
        image[row - 8 : row + 9, col - 8 : col + 9] = (205, 30, 30)

    image[:80, :] = 250
    image[:, -40:] = 0
    return image


@pytest.fixture
def slide_file(tmp_path, synthetic_slide):
    """The synthetic slide written to disk as a PNG."""
    path = tmp_path / "slide.png"
    Image.fromarray(synthetic_slide).save(path)
    return path
