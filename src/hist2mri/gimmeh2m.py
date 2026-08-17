#!/usr/bin/env python3
"""
Welcome to the new gimmeh2m python file.
This should be as close to a one to one replication of the matlab equivalent.
gimmeh2m uses the following files.
    -gimmeh2m.m
    -cellcounts.m
    -gimmeSegs.m
    -SeparateStains.m
    -normalizeImage.m

I have put all these files and merged them all into this single python file as their own methods.
I am going to try and comment a lot of what this is doing
and compare to the original to help folks unfamiliar with python.

Usage:
    python hist2mri.py slide.tif
    python hist2mri.py slide.tif --show --small
    python hist2mri.py slide.tif --use-pre-segs --outdir results/
    python hist2mri.py slide.tif --outdir results/ --debug     (prints stats)

Requires: numpy, scipy, pillow  (matplotlib only for --show)
"""

import argparse
import os

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.io import loadmat, savemat

# Whole-slides are massive in size and Pillow has a limit that is supposed to stop you
# from using such big files so there isn't a decompression bomb. But that wont happen with us.
# so we are turning that off.
Image.MAX_IMAGE_PIXELS = None


# 50 is the tile sized used and we need 0.02 of it.
BLOCK = 50            # blockproc block size
SCALE = 1.0 / BLOCK   # the 0.02 in the MATLAB source

DEBUG = False         # set by --debug; prints intermediate statistics


def _dbg(msg):
    if DEBUG:
        print(f"    [debug] {msg}")


# the following are some matlab methods that do not have an exact 1-to-1 that we need to use.
# python has a mehtod that could do this but does not have the same tie breaking as matlab and bins.
def _graythresh(ch):
    """MATLAB graythresh: Otsu's method on a 256-bin histogram.
    

    Returns a level in [0, 1].
    """
    counts = np.bincount(ch.ravel(), minlength=256).astype(np.float64)
    p = counts / counts.sum()

    
    omega = np.cumsum(p)
    mu = np.cumsum(p * np.arange(1, 257))
    mu_t = mu[-1]

    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_b_sq = (mu_t * omega - mu) ** 2 / (omega * (1.0 - omega))

    finite = np.isfinite(sigma_b_sq)
    if not finite.any():
        return 0.0

    maxval = sigma_b_sq[finite].max()
    # MATLAB: idx = mean(find(sigma_b_squared == maxval))  -- 1-based
    idx = np.mean(np.flatnonzero(sigma_b_sq == maxval) + 1)
    return (idx - 1.0) / 255.0


def _imbinarize(ch):
    """MATLAB imbinarize(I) for a uint8 channel: global Otsu, strictly greater.
    
    every pixel greater than threshold is True, below is False."""
    level = _graythresh(ch)
    return ch.astype(np.float64) / 255.0 > level


def _stretchlim(ch, tol_low=0.01, tol_high=0.99, nbins=65536):
    """MATLAB stretchlim for a double image already in [0, 1].

    MATLAB uses 256 bins ONLY when the input is uint8; every other class,
    including double, gets 65536. normalizeImage passes a double, so 65536
    is the correct count here. (graythresh above is a different function and
    genuinely does use 256 on the uint8 green channel.)

    We need to find the top 1% brightest pixels and the top 1% 
    darkest pixels so it does not throw off the image.
    """
    counts = np.zeros(nbins, dtype=np.int64)
    flat = ch.reshape(-1)
    step = 1 << 24
    for i in range(0, flat.size, step):
        idx = np.rint(flat[i:i + step] * (nbins - 1)).astype(np.int32)
        np.clip(idx, 0, nbins - 1, out=idx)
        counts += np.bincount(idx, minlength=nbins)

    cdf = np.cumsum(counts) / counts.sum()

    lo = np.flatnonzero(cdf > tol_low)
    hi = np.flatnonzero(cdf >= tol_high)
    ilow = lo[0] if lo.size else 0
    ihigh = hi[0] if hi.size else nbins - 1
    if ilow == ihigh:
        return 0.0, 1.0
    return ilow / (nbins - 1.0), ihigh / (nbins - 1.0)


def _imadjust(ch, low, high):
    """MATLAB imadjust with default gamma=1: linear rescale then clip.
    This will clip everything back betwen 0 - 1 that was clipped from Stretchlim
    """
    if high <= low:
        return ch.copy()
    return np.clip((ch - low) / (high - low), 0.0, 1.0)


def _rgb_hue_sat(im):
    """Hue and saturation channels of MATLAB rgb2hsv, both in [0, 1].
    Takes every pixel and converts it to a Hue, saturation, brightness metric.
    This makes colors more easily chosen if its a vessel or other parts.

    Value is never used downstream, so it isn't computed or stored.
    """
    a = im.astype(np.float64) / 255.0
    r, g, b = a[..., 0], a[..., 1], a[..., 2]

    v = np.maximum(np.maximum(r, g), b)
    mn = np.minimum(np.minimum(r, g), b)
    d = v - mn

    s = np.zeros_like(v)
    nz = v > 0
    s[nz] = d[nz] / v[nz]

    h = np.zeros_like(v)
    m = d > 0
    is_r = m & (v == r)
    is_g = m & (v == g) & ~is_r
    is_b = m & ~is_r & ~is_g

    h[is_r] = ((g[is_r] - b[is_r]) / d[is_r]) % 6.0
    h[is_g] = (b[is_g] - r[is_g]) / d[is_g] + 2.0
    h[is_b] = (r[is_b] - g[is_b]) / d[is_b] + 4.0

    return h / 6.0, s


def _imresize(im, size_hw):
    """MATLAB imresize (bicubic, antialiased) via Pillow.

    Pillow's BICUBIC downscale scales the kernel support by the reduction
    factor, which is what MATLAB's antialiasing does. Both use a = -0.5.
    Results are close but not bit-identical.
    """
    h, w = size_hw
    return np.asarray(Image.fromarray(im).resize((w, h), Image.BICUBIC))

# This is the cellcont.m file as a method
def _cell_count(mask):
    """cellcount.m: bwconncomp with default 8-connectivity, return the count."""
    if not mask.any():
        return 0
    _, n = ndimage.label(mask, structure=np.ones((3, 3), dtype=bool))
    return n


# scikit-image's H&E&DAB reference vectors, as hardcoded in gimmeSegs.m
_HE = np.array([0.6443186, 0.7166757, 0.26688856])
_EO = np.array([0.09283128, 0.9545457, 0.28324])
_RES = np.array([0.63595444, 0.001, 0.7717266])


def _deconv_vector():
    hdab_to_rgb = np.stack(
        [_HE / np.linalg.norm(_HE),
         _EO / np.linalg.norm(_EO),
         _RES / np.linalg.norm(_RES)],
        axis=0,
    )
    return np.linalg.inv(hdab_to_rgb)[:, 0]

# Method version fo the NormalizeImage.m file.
def _normalize_image(ch):
    """normalizeImage(x, 'stretch') applied to a single channel.

    Three steps, in order: min-max to [0,1], invert, then a percentile
    stretch.
    """
    lo, hi = ch.min(), ch.max()
    _dbg(f"raw hematoxylin: min={lo:.6f} max={hi:.6f} range={hi-lo:.6f}")
    _dbg(f"raw percentiles: 1%={np.percentile(ch, 1):.6f} "
         f"50%={np.percentile(ch, 50):.6f} 99%={np.percentile(ch, 99):.6f}")
    if hi <= lo:
        # MATLAB would produce NaN here (0/0). Degenerate input.
        return np.zeros_like(ch)

    ch = (ch - lo) / (hi - lo)
    ch = 1.0 - ch
    slo, shi = _stretchlim(ch)
    _dbg(f"stretchlim: low={slo:.6f} high={shi:.6f}")
    out = _imadjust(ch, slo, shi)
    _dbg(f"post-stretch: mean={out.mean():.6f} "
         f"frac<=0.3 (becomes nuclei)={float((out <= 0.3).mean()):.4%}")
    return out

# method version of SeperateStains.m file
def _separate_stains_h(im, chunk_rows=2048):
    """SeparateStains, returning only the normalized hematoxylin channel.

    Channels 2 and 3 are computed and discarded by gimmeSegs.m.
    """
    w = _deconv_vector()
    h, wid = im.shape[:2]
    od = np.empty((h, wid), dtype=np.float64)

    for i in range(0, h, chunk_rows):
        j = min(i + chunk_rows, h)
        blk = im[i:j].astype(np.float64)
        blk += 2.0                    # MATLAB: +2 to avoid log(0)
        np.log(blk, out=blk)
        np.negative(blk, out=blk)
        od[i:j] = blk @ w

    return _normalize_image(od)



# The method version of gimmeSegs.m file.

def gimme_segs(image, show_me=False, save_segs=False, outdir="."):
    """Segment an RGB slide into ecf / vessel / nuclei / pink uint8 masks.
    """
    image = np.array(image, dtype=np.uint8, copy=True)
    sx, sy = image.shape[:2]

    print("Chopping off black edges")
    black = (image[..., 0] == 0) & (image[..., 1] == 0) & (image[..., 2] == 0)
    fill = int(image[0, 0, 0])
    _dbg(f"top-left pixel = {tuple(int(v) for v in image[0, 0])}, fill value = {fill}")
    _dbg(f"pure-black pixels: {int(black.sum())} ({float(black.mean()):.4%})")
    if DEBUG:
        dark = (image.max(axis=2) <= 10) & ~black
        _dbg(f"near-black but NOT pure black (max channel <=10): "
             f"{int(dark.sum())} ({float(dark.mean()):.4%}) -- these are NOT filled")
    image[black] = fill

    hue, sat = _rgb_hue_sat(image)

    print("1. Segmenting ECF")
    # Otsu on the green channel only.
    if DEBUG:
        _dbg(f"otsu level (green) = {_graythresh(image[..., 1]):.8f}")
    ecf = _imbinarize(image[..., 1]).astype(np.uint8)
    _dbg(f"ecf coverage = {float(ecf.mean()):.4%}")

    print("2. Segmenting Vessels")
    vessel = ((hue < 0.1) & (sat > 0.6)).astype(np.uint8)
    _dbg(f"vessel coverage = {float(vessel.mean()):.6%}")
    del hue, sat

    print("3. Segmenting Nuclei")
    stain_h = _separate_stains_h(image)
    nuc3 = (1.0 - stain_h)
    nuc3 = (nuc3 >= 0.7).astype(np.uint8)
    _dbg(f"nuclei coverage = {float(nuc3.mean()):.4%}")
    if DEBUG and black.any():
        _dbg(f"of the originally-black pixels: "
             f"{float(nuc3[black].mean()):.2%} classified as nuclei, "
             f"{float(ecf[black].mean()):.2%} as ecf")
    del stain_h

    print("4. Segmenting Pink")
    pink = np.clip(
        1 - ecf.astype(np.int16) - vessel.astype(np.int16) - nuc3.astype(np.int16),
        0, 255,
    ).astype(np.uint8)

    print("segmentation complete")

    if show_me:
        _show_segs(image, ecf, vessel, nuc3, pink)

    if save_segs:
        print("Writing segmentations")
        # Variable names match the MATLAB save() calls exactly, including
        # nuclei.mat holding a variable called nuc3.
        savemat(os.path.join(outdir, "ecf.mat"), {"ecf": ecf})
        savemat(os.path.join(outdir, "vessel.mat"), {"vessel": vessel})
        savemat(os.path.join(outdir, "nuclei.mat"), {"nuc3": nuc3})
        savemat(os.path.join(outdir, "pink.mat"), {"pink": pink})

    return ecf, vessel, nuc3, pink


def _grid_shape(h, w, bs=BLOCK):
    return -(-h // bs), -(-w // bs)


def _block_sum(mask, bs=BLOCK):
    """blockproc(mask, [50 50], summify).

    blockproc does not pad partial blocks by default, so edge
    blocks are smaller.
    """
    h, w = mask.shape
    gh, gw = _grid_shape(h, w, bs)
    padded = np.zeros((gh * bs, gw * bs), dtype=np.int64)
    padded[:h, :w] = mask
    return padded.reshape(gh, bs, gw, bs).sum(axis=(1, 3)).astype(np.float64)


def _block_count(mask, bs=BLOCK):
    """blockproc(mask, [50 50], countify).

    MATLAB QUIRK: cellcount runs per block, so a nucleus straddling a block
    boundary is counted once in each block. Zero-padding creates no new
    components, so partial edge blocks behave identically.
    """
    h, w = mask.shape
    gh, gw = _grid_shape(h, w, bs)
    out = np.zeros((gh, gw), dtype=np.float64)

    for i in range(gh):
        r0, r1 = i * bs, min((i + 1) * bs, h)
        for j in range(gw):
            c0, c1 = j * bs, min((j + 1) * bs, w)
            out[i, j] = _cell_count(mask[r0:r1, c0:c1])
    return out



# The main man gimmeh2m.m
def gimme_h2m(slide, show_yn=False, small_yn=False, use_pre_segs=False,
              outdir="."):
    """Port of gimmeH2M. Returns the dict written to h2m.mat."""
    os.makedirs(outdir, exist_ok=True)

    print("Loading image")
    im = np.asarray(Image.open(slide))

    if im.ndim == 2:
        raise ValueError("Expected an RGB slide, got a single-channel image.")
    if im.shape[2] > 3:
        im = im[:, :, :3]
    im = np.ascontiguousarray(im, dtype=np.uint8)

    if small_yn:
        h, w = im.shape[:2]
        im = _imresize(im, (int(round(h * 0.5)), int(round(w * 0.5))))

    if use_pre_segs:
        print("Loading Segmentations")
        ecf = loadmat(os.path.join(outdir, "ecf.mat"))["ecf"]
        vessel = loadmat(os.path.join(outdir, "vessel.mat"))["vessel"]
        nuc = loadmat(os.path.join(outdir, "nuclei.mat"))["nuc3"]
        pink = loadmat(os.path.join(outdir, "pink.mat"))["pink"]
    else:
        print("Segmenting image:")
        ecf, vessel, nuc, pink = gimme_segs(im, show_me=False, save_segs=True,
                                            outdir=outdir)

    print("Writing to h2m")
    gh, gw = _grid_shape(*ecf.shape)
    orig_image = _imresize(im, (gh, gw))

    ecf_ds = _block_sum(ecf)
    vessel_ds = _block_sum(vessel)
    nuc_ds = _block_sum(nuc)
    pink_ds = _block_sum(pink)
    count_nuc = _block_count(nuc)

    # imgaussfilt(x, 2): 9x9 kernel (2*ceil(2*sigma)+1), replicate padding.
    smoothnuc = ndimage.gaussian_filter(count_nuc, sigma=2, mode="nearest",
                                        truncate=2.0)

    cell_den = np.zeros((gh, gw, 6), dtype=np.float64)
    cell_den[:, :, 0] = ecf_ds
    cell_den[:, :, 1] = vessel_ds
    cell_den[:, :, 2] = nuc_ds
    cell_den[:, :, 3] = pink_ds
    cell_den[:, :, 4] = count_nuc
    cell_den[:, :, 5] = smoothnuc

    out = {
        "orig_image": orig_image,
        "im_adj": orig_image,
        "cell_den": cell_den,
    }

    savemat(os.path.join(outdir, "h2m.mat"), {"out": out})
    print("hist2mri 3.0 complete")

    if show_yn:
        _show_h2m(cell_den)

    return out


# matplotlib needed to show.
def _show_segs(image, ecf, vessel, nuc3, pink):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 5, figsize=(16, 6))
    for ax, data, title in zip(
        axes,
        [image, ecf, vessel, nuc3, pink],
        ["orig", "ecf", "vessel", "nuclei", "pink"],
    ):
        ax.imshow(data)
        ax.set_title(title)
        ax.set_axis_off()
    plt.tight_layout()
    plt.show()


def _show_h2m(cell_den):
    import matplotlib.pyplot as plt

    titles = ["ECF", "Vessel", "Nuclei", "Pink", "Cell Count",
              "Smoothed Cell Count"]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    for k, (ax, title) in enumerate(zip(axes.ravel(), titles)):
        ax.imshow(cell_den[:, :, k])
        ax.set_title(title)
        ax.set_axis_off()
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="hist2mri 3.0")
    ap.add_argument("slide", help="path to the slide image")
    ap.add_argument("--show", action="store_true",
                    help="display the six density maps (showYn)")
    ap.add_argument("--small", action="store_true",
                    help="halve the image before processing (smallYn)")
    ap.add_argument("--use-pre-segs", action="store_true",
                    help="reuse saved segmentation .mat files (usePreSegs)")
    ap.add_argument("--outdir", default=".",
                    help="where the .mat files are written (default: cwd)")
    ap.add_argument("--debug", action="store_true",
                    help="print intermediate statistics")
    args = ap.parse_args()

    global DEBUG
    DEBUG = args.debug

    gimme_h2m(args.slide, show_yn=args.show, small_yn=args.small,
              use_pre_segs=args.use_pre_segs, outdir=args.outdir)


if __name__ == "__main__":
    main()