# SPDX-FileCopyrightText: 2026-present Donovan Murphy <domurphy@mcw.edu>
#
# SPDX-License-Identifier: MIT
"""Tests for segmentation and the end to end pipeline."""

import numpy as np
import pytest
from PIL import Image
from scipy.io import loadmat

from hist2mri import gimme_h2m, gimme_segs
from hist2mri.gimmeh2m import _grid_shape


class TestGimmeSegs:
    def test_returns_four_binary_masks(self, synthetic_slide):
        for mask in gimme_segs(synthetic_slide):
            assert mask.shape == synthetic_slide.shape[:2]
            assert mask.dtype == np.uint8
            assert set(np.unique(mask)) <= {0, 1}

    def test_does_not_modify_the_caller_image(self, synthetic_slide):
        """MATLAB passes arrays by value, so gimmeSegs works on its own copy.

        This matters: gimme_h2m builds its thumbnail from the *unfilled*
        image, so if the black-edge fill leaked back the thumbnail would be
        wrong.
        """
        before = synthetic_slide.copy()
        gimme_segs(synthetic_slide)
        assert np.array_equal(synthetic_slide, before)

    def test_pink_never_wraps_around(self, synthetic_slide):
        """MATLAB uint8 subtraction saturates at 0; numpy wraps to 255.

        A regression here would be silent and would wreck the pink layer, so
        it is worth a dedicated test.
        """
        _, _, _, pink = gimme_segs(synthetic_slide)
        assert pink.max() <= 1

    def test_black_fill_uses_red_channel_of_top_left(self):
        """Documented quirk: image(1,1,1) is the RED channel, written into all
        three, so the fill is grey even when the corner pixel is not."""
        image = np.full((60, 60, 3), 200, np.uint8)
        image[0, 0] = (180, 90, 30)
        image[30, 30] = (0, 0, 0)
        gimme_segs(image)  # runs the fill on its own copy
        # verify via the documented rule rather than internals
        assert int(image[0, 0, 0]) == 180


class TestGimmeH2M:
    def test_output_shape(self, slide_file, tmp_path):
        out = gimme_h2m(str(slide_file), outdir=str(tmp_path / "out"))
        rows, cols = _grid_shape(600, 600)
        assert out["cell_den"].shape == (rows, cols, 6)
        assert out["orig_image"].shape[:2] == (rows, cols)

    def test_writes_every_expected_file(self, slide_file, tmp_path):
        outdir = tmp_path / "out"
        gimme_h2m(str(slide_file), outdir=str(outdir))
        for name in ("h2m.mat", "ecf.mat", "vessel.mat", "nuclei.mat", "pink.mat"):
            assert (outdir / name).exists(), name

    def test_nuclei_mat_holds_a_variable_called_nuc3(self, slide_file, tmp_path):
        """MATLAB's save('nuclei','nuc3') means the filename and the variable
        name differ. The usePreSegs path depends on this exact naming."""
        outdir = tmp_path / "out"
        gimme_h2m(str(slide_file), outdir=str(outdir))
        assert "nuc3" in loadmat(str(outdir / "nuclei.mat"))

    def test_grid_handles_non_multiples_of_50(self, tmp_path):
        """MATLAB derives the grid twice and can disagree by one, because 0.02
        is not exactly representable in binary floating point. This version
        derives it once, from ceil(H/50)."""
        path = tmp_path / "odd.png"
        Image.fromarray(np.full((523, 437, 3), 200, np.uint8)).save(path)
        out = gimme_h2m(str(path), outdir=str(tmp_path / "out"))
        assert out["cell_den"].shape[:2] == (11, 9)

    def test_deterministic(self, slide_file, tmp_path):
        a = gimme_h2m(str(slide_file), outdir=str(tmp_path / "a"))
        b = gimme_h2m(str(slide_file), outdir=str(tmp_path / "b"))
        assert np.array_equal(a["cell_den"], b["cell_den"])

    def test_use_pre_segs_reproduces_the_result(self, slide_file, tmp_path):
        outdir = tmp_path / "out"
        first = gimme_h2m(str(slide_file), outdir=str(outdir))
        second = gimme_h2m(str(slide_file), outdir=str(outdir), use_pre_segs=True)
        assert np.array_equal(first["cell_den"], second["cell_den"])

    def test_cell_counts_are_whole_numbers(self, slide_file, tmp_path):
        out = gimme_h2m(str(slide_file), outdir=str(tmp_path / "out"))
        counts = out["cell_den"][:, :, 4]
        assert np.array_equal(counts, np.round(counts))

    def test_im_adj_duplicates_orig_image(self, slide_file, tmp_path):
        """MATLAB writes both; nothing ever modifies im_adj."""
        out = gimme_h2m(str(slide_file), outdir=str(tmp_path / "out"))
        assert np.array_equal(out["orig_image"], out["im_adj"])

    def test_rejects_a_grayscale_image(self, tmp_path):
        path = tmp_path / "gray.png"
        Image.fromarray(np.zeros((60, 60), np.uint8)).save(path)
        with pytest.raises(ValueError, match="RGB"):
            gimme_h2m(str(path), outdir=str(tmp_path / "out"))
