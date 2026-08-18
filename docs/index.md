# hist2mri

Histology to MRI tissue density mapping. A python port of the lab's MATLAB
gimmeh2m file.

## Overview

A pathologist scans a stained slice of tissue and you get a gigapixel image
where individual cells are visible. The same patient's MRI is very low
resolution -- one voxel covers thousands of cells. This pipeline bridges the
two: it reduces the huge microscope image into six coarse density maps at 1/50
scale, so you can ask what the tissue actually looks like inside one MRI voxel.

**Input:** one slide image.

**Output:** `h2m.mat` holding a thumbnail plus a `rows x cols x 6` stack, where
every value summarises one 50x50 pixel tile of the slide.

| Layer | Name | What it counts |
|-------|------|----------------|
| 1 | ECF | empty space / extracellular fluid pixels |
| 2 | Vessel | blood vessel pixels |
| 3 | Nuclei | nuclei pixels (area) |
| 4 | Pink | cytoplasm / connective tissue pixels |
| 5 | Cell Count | separate nuclei (count, not area) |
| 6 | Smoothed Cell Count | layer 5, gaussian smoothed |

Layer 5 is the only one counting *objects* rather than pixels. A tile with a
few large nuclei and a tile with many small ones can have identical pixel area
but very different cell counts.

## Quick Start

### Installation

```console
git clone https://github.com/DONNYMURPH/matlab-conversions.git
cd matlab-conversions
pip install -e .
```

Add the plotting extra if you want `--show`:

```console
pip install -e ".[viz]"
```

### From the shell

```console
hist2mri run slide.tif --outdir results/
hist2mri run slide.tif --outdir results/ -vv     # intermediate statistics
hist2mri run slide.tif --outdir results/ --show  # needs [viz]
```

### From python

```python
from hist2mri import gimme_h2m

out = gimme_h2m("slide.tif", outdir="results/")
nuclei_density = out["cell_den"][:, :, 2]
cell_counts = out["cell_den"][:, :, 4]
```

Progress goes through `logging`, so whoever calls the library decides what gets
shown:

```python
import logging

logging.basicConfig(level=logging.INFO)  # stage messages
logging.basicConfig(level=logging.DEBUG)  # plus thresholds and coverages
```

## Which MATLAB file is where

Five MATLAB files were merged into `gimmeh2m.py`, each as its own function.

| MATLAB | Python |
|--------|--------|
| `gimmeH2M.m` | `gimme_h2m()` |
| `gimmeSegs.m` | `gimme_segs()` |
| `cellcount.m` | `_cell_count()` |
| `SeparateStains.m` | `_separate_stains_h()` |
| `normalizeImage.m` | `_normalize_image()` |

## Fidelity

Validated against MATLAB on 4000x4000 crops from five slides:

| Layer | Result |
|-------|--------|
| ECF | 0 of 16,000,000 pixels differ |
| Nuclei | 0 of 16,000,000 pixels differ |
| Cell count | identical, max block difference 0 |
| Smoothed count | agrees to ~1e-14 (float64 rounding) |
| Vessel | 3 pixels on one slide |

Those 3 pixels all have saturation of exactly 0.6, where the vessel rule tests
`> 0.6` and the last bit of a float64 division rounds differently in MATLAB
than in numpy. The two disagree in *both* directions, so neither is more
correct, and all three are near-black background rather than real vessels.

## Memory

The deconvolution is chunked and only the hematoxylin channel is computed, so a
594-megapixel slide runs in roughly 16 GB. The MATLAB original builds the whole
float64 array at once and needs about 40 GB, which is why it runs out of memory
on full-resolution slides where this does not.

## Development

Clone the repository and install [Hatch](https://hatch.pypa.io/latest/install/):

```console
git clone https://github.com/DONNYMURPH/matlab-conversions.git
cd matlab-conversions
pip install hatch
```

### Common Commands

| Task | Command |
|------|---------|
| Run tests | `hatch run test:test` |
| Run tests with coverage | `hatch run test:cov` |
| Lint code | `hatch run lint:check` |
| Format code | `hatch run lint:format` |
| Auto-fix lint issues | `hatch run lint:fix` |
| Type check | `hatch run types:check` |
| Build docs | `hatch run docs:build-docs` |
| Serve docs locally | `hatch run docs:serve-docs` |
| Build wheel | `hatch build` |

### Project Structure

```text
hist2mri/
├── src/
│   └── hist2mri/
│       ├── __init__.py        # public API
│       ├── __about__.py       # version
│       ├── cli.py             # command line wrapper
│       ├── gimmeh2m.py        # the pipeline (all five MATLAB files)
│       └── py.typed           # PEP 561 type marker
├── tests/
│   ├── conftest.py            # synthetic slide fixtures
│   ├── test_matlab_compat.py  # the MATLAB built-in reimplementations
│   ├── test_blocks.py         # blockproc replacements
│   ├── test_gimmeh2m.py       # segmentation + end to end
│   └── test_cli.py
├── docs/
├── pyproject.toml
├── Dockerfile
├── Makefile
└── mkdocs.yml
```

## API Reference

::: hist2mri