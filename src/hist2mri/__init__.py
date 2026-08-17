# SPDX-FileCopyrightText: 2026-present Donovan Murphy <domurphy@mcw.edu>
#
# SPDX-License-Identifier: MIT
"""
hist2mri -- histology to MRI tissue density mapping.

A python port of the gimmeh2m.m file. A pathologist scans a
stained slice of tissue and you get a gigapixel image where you can see
individual cells. The same patient's MRI is very low resolution, where one
voxel covers thousands of cells. This pipeline bridges the two: it boils the
huge microscope image down into six coarse density maps at 1/50 scale, so you
can ask what the tissue actually looks like inside one MRI voxel.

Each value in the output summarises one 50x50 pixel tile of the slide:

===== ==================== =========================================
Layer Name                 What it counts
===== ==================== =========================================
1     ECF                  empty space / extracellular fluid pixels
2     Vessel               blood vessel pixels
3     Nuclei               nuclei pixels (area)
4     Pink                 cytoplasm / connective tissue pixels
5     Cell Count           separate nuclei (count, not area)
6     Smoothed Cell Count  layer 5, gaussian smoothed
===== ==================== =========================================

Quick start::

    from hist2mri import gimme_h2m

    out = gimme_h2m("slide.tif", outdir="results/")
    nuclei_density = out["cell_den"][:, :, 2]

Or from the shell::

    hist2mri run slide.tif --outdir results/
"""

from hist2mri.__about__ import __version__
from hist2mri.gimmeh2m import BLOCK, gimme_h2m, gimme_segs

__all__ = ["BLOCK", "__version__", "gimme_h2m", "gimme_segs"]
