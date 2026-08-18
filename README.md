# hist2mri

Histology to MRI tissue density mapping. A python port of the lab's MATLAB
gimmeh2m file.

## What it does

A pathologist scans a stained slice of tissue and you get a gigapixel image
where individual cells are visible. The same patient's MRI is very low
resolution -- one voxel covers thousands of cells. This pipeline bridges the
two: it reduces the huge microscope image to six coarse density maps at 1/50
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

Layers 1-4 sum to 2500 per tile (50 x 50), so they read as proportions. Layer 5
is the only one counting objects rather than pixels -- two tiles can have
identical nuclei area but very different cell counts.

The output format is unchanged from the MATLAB version, so anything downstream
that already reads `h2m.mat` keeps working.

## Install

```console
git clone https://github.com/DONNYMURPH/matlab-conversions.git
cd matlab-conversions
pip install -e .
```

Add the plotting extra if you want `--show`:

```console
pip install -e ".[viz]"
```

```console
python -m venv ~/h2m-venv
source ~/h2m-venv/bin/activate
pip install -e .
```

## Use

From the shell:

```console
hist2mri run slide.tif --outdir results/
hist2mri run slide.tif --outdir results/ -v      # stage messages
hist2mri run slide.tif --outdir results/ -vv     # plus thresholds, coverages, peak memory
hist2mri run slide.tif --outdir results/ --show  # needs [viz]
```

**Never point `--outdir` at a folder containing slide data.** The pipeline
writes `h2m.mat`, `ecf.mat`, `vessel.mat`, `nuclei.mat` and `pink.mat` -- the
same names the MATLAB version uses -- so it would overwrite whatever is
already there.

From python:

```python
from hist2mri import gimme_h2m

out = gimme_h2m("slide.tif", outdir="results/")
nuclei_density = out["cell_den"][:, :, 2]
cell_counts = out["cell_den"][:, :, 4]
```

Progress goes through `logging`, so a caller decides what gets shown:

```python
import logging

logging.basicConfig(level=logging.INFO)   # stage messages
logging.basicConfig(level=logging.DEBUG)  # plus thresholds and coverages
```

Reading the result back in MATLAB works exactly as before:

```matlab
load('h2m.mat')
imagesc(out.cell_den(:,:,3))
```

See [`docs/usage.md`](docs/usage.md) for a fuller guide aimed at people who
do not use Python.

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

Validated against fresh MATLAB runs on 4000x4000 crops from five slides,
compared at full pixel resolution:

| Layer | Result |
|-------|--------|
| ECF | 0 of 16,000,000 pixels differ, on all five slides |
| Nuclei | 0 of 16,000,000 pixels differ, on all five slides |
| Cell count | identical, max block difference 0 |
| Smoothed count | agrees to ~1e-14 (float64 rounding) |
| Pink | identical on four of five slides |
| Vessel | identical on four of five; 3 pixels differ on the fifth |

Those 3 pixels all have saturation of exactly 0.6, where the vessel rule tests
`> 0.6` and the last bit of a float64 division rounds differently in MATLAB
than in numpy. The two disagree in *both* directions, so neither is more
correct, and all three are near-black background rather than real vessels.

## Memory

Measured on a 594-megapixel slide (21184 x 28054): **11.5 GB peak**, with the
whole deconvolution accounting for only about 2 GB of that. Loading the image
dominates.

Rough rule: allow about **20 bytes of RAM per pixel** of slide.

The hue and saturation computation, the Otsu binarization and the deconvolution
are all chunked, and the normalize, complement and adjust steps work in place.
A direct translation of the MATLAB needs roughly 40 GB for the same slide,
which is why the original runs out of memory on full-resolution slides where
this does not.

## Development

With [hatch](https://hatch.pypa.io/):

```console
pip install hatch
hatch run test:test      # tests
hatch run test:cov       # tests with coverage
hatch run lint:check     # ruff
hatch run lint:all       # format + autofix + check
hatch run types:check    # mypy
hatch run docs:serve-docs
hatch build
```

If hatch cannot create its environments -- which happens on shared machines
where the system Python is read-only -- use a plain virtualenv and call the
tools directly:

```console
python -m venv ~/h2m-venv
source ~/h2m-venv/bin/activate
pip install -e . ruff pytest

ruff check src tests
ruff format src tests --check
python -m pytest tests -q
```

Those are the same checks CI runs.

### Testing notes

The suite covers the MATLAB built-in reimplementations most heavily, since
that is where every real discrepancy during the port turned out to be.