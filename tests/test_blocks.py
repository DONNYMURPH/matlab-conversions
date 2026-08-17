# SPDX-FileCopyrightText: 2026-present Donovan Murphy <domurphy@mcw.edu>
#
# SPDX-License-Identifier: MIT
"""Tests for the blockproc replacements."""

import numpy as np

from hist2mri.gimmeh2m import _block_count, _block_sum, _grid_shape


class TestGridShape:
    def test_exact_multiple(self):
        assert _grid_shape(500, 250) == (10, 5)

    def test_rounds_up_for_partial_tiles(self):
        """blockproc keeps a partial edge tile rather than dropping pixels."""
        assert _grid_shape(501, 251) == (11, 6)

    def test_matches_a_real_slide(self):
        """Dimensions of a production slide, confirmed against MATLAB."""
        assert _grid_shape(28054, 21184) == (562, 424)


class TestBlockSum:
    def test_counts_per_tile(self):
        mask = np.zeros((100, 100), np.uint8)
        mask[:50, :50] = 1
        result = _block_sum(mask)
        assert result.shape == (2, 2)
        assert result[0, 0] == 2500
        assert result[0, 1] == 0

    def test_partial_edge_tiles(self):
        """Edge tiles are smaller; zero padding gives the same sum."""
        result = _block_sum(np.ones((60, 60), np.uint8))
        assert result[0, 0] == 2500
        assert result[1, 1] == 100

    def test_total_is_preserved(self):
        rng = np.random.default_rng(4)
        mask = (rng.random((237, 189)) > 0.5).astype(np.uint8)
        assert _block_sum(mask).sum() == mask.sum()


class TestBlockCount:
    def test_counts_blobs_not_pixels(self):
        mask = np.zeros((50, 50), np.uint8)
        mask[10:20, 10:20] = 1
        assert _block_count(mask)[0, 0] == 1

    def test_blob_on_a_boundary_is_counted_twice(self):
        """Documented original behaviour, kept on purpose.

        Each tile is labelled independently, so a nucleus straddling a tile
        edge is counted once in each tile.
        """
        mask = np.zeros((100, 50), np.uint8)
        mask[45:55, 20:30] = 1
        assert _block_count(mask).sum() == 2
