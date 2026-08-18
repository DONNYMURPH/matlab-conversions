# Using hist2mri

A guide for anyone in the lab who wants to run the pipeline. No Python
experience assumed.

---

## What it does

You give it one histology slide image. It gives you back six density maps,
each 1/50th the width and height of the slide, plus a small thumbnail.

Every value in those maps summarises one 50×50-pixel tile of the original
slide, which is roughly MRI voxel scale. That's the whole point: it turns a
gigapixel image where you can see individual cells into something coarse
enough to line up against MRI.

| Layer | Name | What the number means |
|-------|------|----------------------|
| 1 | ECF | how many empty-space pixels are in this tile |
| 2 | Vessel | how many blood vessel pixels |
| 3 | Nuclei | how many nuclei pixels (area) |
| 4 | Pink | how many cytoplasm / connective tissue pixels |
| 5 | Cell Count | how many **separate nuclei** (count, not area) |
| 6 | Smoothed Cell Count | layer 5 with a light Gaussian blur |

Layers 1–4 add up to 2500 per tile (50 × 50), so you can read them as
percentages. Layer 5 is the only one counting objects rather than pixels — two
tiles can have identical nuclei *area* but very different cell *counts*.

This is a port of the lab's MATLAB pipeline (`gimmeH2M.m` and friends). It
produces the same numbers; it just runs in Python and needs far less memory.

---

## One-time setup

You need Python 3.9 or newer. Check with `python --version`.

### Option A — you just want to run it

```bash
pip install hist2mri
```

### Option B — from the repository

```bash
git clone https://github.com/LavLabInfrastructure/hist2mri.git
cd hist2mri
pip install -e .
```

### If pip complains about permissions

You're probably on a shared machine where you can't write to the system
Python. Make your own environment:

```bash
python -m venv ~/h2m-venv
source ~/h2m-venv/bin/activate
pip install hist2mri
```

You'll need to run that `source` line once in every new terminal. When it's
active your prompt shows `(h2m-venv)`.

### Check it worked

```bash
hist2mri --version
```

If you get "command not found" but the install succeeded, the venv isn't
active — run the `source` line again.

---

## Running it

### The basic command

```bash
hist2mri run /path/to/slide.tif --outdir /path/to/where/output/goes
```

That's it. Two things: the slide, and where to put the results.

**Tip for long paths:** in most terminals you can drag a file from the file
browser into the terminal window and it pastes the full path for you.

### A worked example

```bash
mkdir -p ~/h2m_results
hist2mri run /Volumes/Siren/Brain_data/1.PatientDirectory/113/Histology/Processed/S14_Large/HE/large_recon_10_nowhite_HE.tiff --outdir ~/h2m_results -v
```

You'll see it work through the stages:

```
Loading image
Segmenting image:
Chopping off black edges
1. Segmenting ECF
2. Segmenting Vessels
3. Segmenting Nuclei
4. Segmenting Pink
segmentation complete
Writing segmentations
Writing to h2m
hist2mri 3.0 complete
wrote /home/you/h2m_results/h2m.mat  (grid 562 x 424)
```

> **⚠ Never point `--outdir` at a folder containing slide data.**
> The pipeline writes files called `h2m.mat`, `ecf.mat`, `vessel.mat`,
> `nuclei.mat` and `pink.mat`. Those are the same names the MATLAB version
> uses, so pointing it at a slide folder will overwrite whatever is already
> there. Always write to your own directory.

---

## Options

| Option | What it does |
|--------|--------------|
| `--outdir DIR` | Where results go. Defaults to the current folder — set it explicitly. |
| `-v` | Show stage messages as it runs. Recommended. |
| `-vv` | Also show thresholds, coverage percentages, and peak memory. Useful when something looks wrong. |
| `-q` | Silent except the final line. |
| `--small` | Halve the image before processing. **Changes the meaning of the output** — see below. |
| `--use-pre-segs` | Skip segmentation and reuse mask files already in `--outdir`. |
| `--show` | Pop up a window with the six maps. Needs `pip install "hist2mri[viz]"` and a display. |
| `--version` | Print the version. |
| `--help` | List all options. |

### About `--small`

This halves the image first, so it runs about four times faster. But the tile
size stays at 50 pixels, which means each output value now covers a 100×100
region of the *original* slide instead of 50×50.

So `--small` output is not comparable to normal output. Use it for a quick
look, not for anything you'll analyse.

### About `--use-pre-segs`

Segmentation is the slow part. If you've already run a slide and only want to
redo the tile counting, this reuses the saved masks. Only useful if you're
changing something downstream.

---

## What you get

In your output folder:

| File | Contents |
|------|----------|
| `h2m.mat` | **The main output.** Thumbnail plus the six-layer stack. |
| `ecf.mat` | Full-resolution empty-space mask |
| `vessel.mat` | Full-resolution vessel mask |
| `nuclei.mat` | Full-resolution nuclei mask (variable inside is called `nuc3`) |
| `pink.mat` | Full-resolution pink mask |

The four mask files are large — tens of megabytes each — because they're at
full slide resolution. Delete them if you only need `h2m.mat`.

### Reading the results in MATLAB

Exactly as before. Nothing changed about the format.

```matlab
load('h2m.mat')
imagesc(out.cell_den(:,:,3))   % nuclei
title('Nuclei density')

nuclei_fraction = out.cell_den(:,:,3) / 2500;   % as a proportion
cell_counts = out.cell_den(:,:,5);
```

### Reading the results in Python

```python
from scipy.io import loadmat

m = loadmat("h2m.mat", squeeze_me=True, struct_as_record=False)["out"]
cell_den = m.cell_den

nuclei = cell_den[:, :, 2]        # note: Python counts from 0
cell_counts = cell_den[:, :, 4]
```

Careful with the index shift — MATLAB layer 3 is Python index 2.

### Using it directly in Python

If you're writing a script or working in a notebook, skip the command line:

```python
from hist2mri import gimme_h2m

out = gimme_h2m("slide.tif", outdir="results/")
nuclei = out["cell_den"][:, :, 2]
```

To see progress in a notebook:

```python
import logging
logging.basicConfig(level=logging.INFO)
```

---

## How long, and how much memory

Depends mostly on slide size. For reference, a 594-megapixel slide
(21184 × 28054):

- **Loading:** a few minutes if the slide is on a network drive. This is
  usually the slowest part, and nothing prints while it happens.
- **Processing:** a few minutes.
- **Peak memory:** about 11.5 GB.

Rough rule: **allow about 20 bytes of RAM per pixel of slide**. A 600 MP slide
wants ~12 GB, a 200 MP slide wants ~4 GB.

If your machine has less than that, either process a crop or use `--small`.

There are two long silences where nothing prints — after `Loading image`, and
during `3. Segmenting Nuclei`. Both are normal. Use `-vv` if you want to watch
memory as it goes.

---

## Troubleshooting

**"command not found: hist2mri"**
The virtual environment isn't active. Run `source ~/h2m-venv/bin/activate`.

**"No such file or directory"**
Check the slide path. A common mistake is putting `~/` in front of an absolute
path — `~/Volumes/...` means `/home/you/Volumes/...`, which isn't the same as
`/Volumes/...`.

**Killed, or the terminal returns with no message**
Out of memory. Check available RAM with `free -g` (Linux). Process a crop or
use `--small`.

**"Expected an RGB slide, got a single-channel image"**
The file is greyscale. This pipeline needs colour — it works by separating the
two H&E stains.

**A wall of TIFF tag output**
You're on an old version. Update; `-vv` used to switch on debug logging for
every library in the process.

**Results look wrong**
Run with `-vv` and check the reported numbers. `ecf coverage` far from 30–50%,
or `nuclei coverage` far from 1–3%, suggests something unusual about the slide.

---

## Making a crop

If a slide is too big for your machine, or you just want to test quickly:

```python
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

im = Image.open("/path/to/slide.tif")
w, h = im.size
size = 4000
x, y = (w - size) // 2, (h - size) // 2
im.crop((x, y, x + size, y + size)).convert("RGB").save("crop.png")
```

Then run the pipeline on `crop.png` as normal.

Note that a 4000×4000 crop is about 2% of a large slide, so it may miss
features that are sparse or clustered — vessels in particular.

---

## Things worth knowing about the results

These are properties of the original MATLAB pipeline, faithfully reproduced.
Worth understanding before drawing conclusions.

**The nuclei threshold is relative to each slide.** Nuclei are found by
measuring hematoxylin stain, then rescaling to that slide's own range and
cutting at 0.7. So the same real stain intensity can land on either side of
the cut on two different slides, depending on how each one was stained. Nuclei
counts are more comparable *within* a slide than *between* slides.

**Nuclei on tile boundaries are counted twice.** Layer 5 counts blobs within
each tile independently, so a nucleus straddling an edge is counted in both
tiles. Total counts are therefore slightly inflated.

**The black-edge removal often does nothing.** It only catches pixels that are
*exactly* black. On one production slide it fixed 828k pixels and missed 1.44M
near-black ones.

**The vessel layer detects very little.** Across five slides tested, vessel
coverage ran between 0.0001% and 0.001%. The rule looks for strongly saturated
red, i.e. actual blood — vessel lumens without blood in them read as empty
space instead. Treat near-zero vessel numbers as expected rather than as an
error, and check with Pete or Sam before using that layer for anything.

---

## Getting help

- `hist2mri --help` lists every option
- Issues: https://github.com/LavLabInfrastructure/hist2mri/issues
- The repository README covers how the port maps onto the original MATLAB
  files, and how it was validated