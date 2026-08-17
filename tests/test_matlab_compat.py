# SPDX-FileCopyrightText: 2026-present Donovan Murphy <domurphy@mcw.edu>
#
# SPDX-License-Identifier: MIT
"""Tests for the helpers that stand in for MATLAB built-ins.

These are the highest value tests here. Every one of these functions caused a
real mismatch against MATLAB at some point during the port, so they are the
places a regression would actually hurt.
"""

import numpy as np
import pytest

from hist2mri.gimmeh2m import (
    _cell_count,
    _graythresh,
    _imadjust,
    _imbinarize,
    _imresize,
    _rgb_hue_sat,
    _stretchlim,
)


class TestGraythresh:
    def test_bimodal_split(self):
        image = np.concatenate(
            [np.full(5000, 40, np.uint8), np.full(5000, 200, np.uint8)]
        ).reshape(100, 100)
        assert 40 / 255 < _graythresh(image) < 200 / 255

    def test_constant_image(self):
        """No valid split exists; MATLAB returns 0."""
        assert _graythresh(np.full((10, 10), 128, np.uint8)) == 0.0

    def test_level_is_a_multiple_of_one_over_255(self):
        """MATLAB returns (idx-1)/255, so graythresh IS on 256 bins."""
        rng = np.random.default_rng(0)
        image = rng.integers(0, 256, (200, 200), dtype=np.uint8)
        level = _graythresh(image)
        assert abs(level * 255 - round(level * 255)) < 1e-9


class TestImbinarize:
    def test_threshold_is_strictly_greater(self):
        image = np.concatenate(
            [np.full(500, 10, np.uint8), np.full(500, 250, np.uint8)]
        ).reshape(20, 50)
        assert _imbinarize(image).sum() == 500


class TestStretchlim:
    def test_uses_65536_bins_not_256(self):
        """The bug that took longest to find.

        MATLAB uses 256 bins only for uint8 input; a double image gets 65536.
        So the limits must come back as multiples of 1/65535. Using 256 here
        shifts the limits by about 0.0008, which mis-thresholds roughly 0.02%
        of nuclei pixels on a real slide.
        """
        rng = np.random.default_rng(1)
        channel = np.clip(rng.normal(0.5, 0.15, (500, 500)), 0, 1)
        for value in _stretchlim(channel):
            assert abs(value * 65535 - round(value * 65535)) < 1e-6

    def test_clips_the_requested_fraction(self):
        channel = np.linspace(0, 1, 10000).reshape(100, 100)
        low, high = _stretchlim(channel)
        assert 0.005 < low < 0.02
        assert 0.98 < high < 0.995

    def test_degenerate_channel_returns_full_range(self):
        assert _stretchlim(np.full((50, 50), 0.5)) == (0.0, 1.0)


class TestImadjust:
    def test_rescales_and_clips(self):
        out = _imadjust(np.array([[0.0, 0.25, 0.5, 0.75, 1.0]]), 0.25, 0.75)
        assert out[0, 0] == 0.0
        assert out[0, 2] == pytest.approx(0.5)
        assert out[0, 4] == 1.0

    def test_degenerate_window_is_passthrough(self):
        channel = np.array([[0.3, 0.7]])
        assert np.array_equal(_imadjust(channel, 0.5, 0.5), channel)


class TestRgbHueSat:
    def test_pure_red(self):
        hue, sat = _rgb_hue_sat(np.array([[[255, 0, 0]]], np.uint8))
        assert hue[0, 0] == pytest.approx(0.0)
        assert sat[0, 0] == pytest.approx(1.0)

    def test_grey_has_no_saturation(self):
        _, sat = _rgb_hue_sat(np.array([[[128, 128, 128]]], np.uint8))
        assert sat[0, 0] == 0.0

    def test_black_does_not_divide_by_zero(self):
        _, sat = _rgb_hue_sat(np.zeros((2, 2, 3), np.uint8))
        assert np.all(sat == 0.0)

    def test_outputs_stay_in_unit_range(self):
        rng = np.random.default_rng(3)
        image = rng.integers(0, 256, (50, 50, 3), dtype=np.uint8)
        hue, sat = _rgb_hue_sat(image)
        assert hue.min() >= 0.0 and hue.max() <= 1.0
        assert sat.min() >= 0.0 and sat.max() <= 1.0


class TestCellCount:
    def test_empty_mask(self):
        assert _cell_count(np.zeros((10, 10), np.uint8)) == 0

    def test_separate_blobs(self):
        mask = np.zeros((10, 10), np.uint8)
        mask[1, 1] = 1
        mask[8, 8] = 1
        assert _cell_count(mask) == 2

    def test_diagonal_touch_counts_as_one_blob(self):
        """MATLAB's bwconncomp defaults to 8-connectivity.

        scipy defaults to 4-connectivity, which would report 2 here and
        silently inflate every cell count in the pipeline.
        """
        mask = np.zeros((5, 5), np.uint8)
        mask[1, 1] = 1
        mask[2, 2] = 1
        assert _cell_count(mask) == 1


class TestImresize:
    def test_target_shape(self):
        assert _imresize(np.zeros((100, 200, 3), np.uint8), (10, 20)).shape == (
            10,
            20,
            3,
        )
