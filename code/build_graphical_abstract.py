"""The graphical abstract."""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_figures as bf

OUT = bf.FIGS
INK, GREY = bf.INK, bf.GREY
CONVERTED, CLASS = bf.CONVERTED, bf.CLASS
BG = "#F8F9FA"
EDGE = "#DEE2E6"
LIPID = "#E7E2D4"
LIPID_EDGE = "#B2AA96"

W, H = 10.6, 4.25
PANELS = (("ga_native.png", 0.022, 0.414), ("ga_converted.png", 0.586, 0.978))
TOP, BOT = 0.865, 0.065
IMG_BOTTOM = 0.105
BAND_TOP, BAND_HALF = 0.260, 0.105
BAND_X0, BAND_X1 = 0.575, 0.990
TITLE = ("Reverse-QTY conversion reaches the obligate 14-3-3 dimer interface "
         "and leaves it intact")


def fatal(msg):
    raise SystemExit("STOP: " + msg)


def gravy(x):
    return "%+.2f" % x if x >= 0 else "%.2f" % x


def load(name):
    p = os.path.join(OUT, name)
    if not os.path.isfile(p):
        fatal("no figures/%s. Run code/render_ga_structures.py with the "
              "open-source PyMOL interpreter first" % name)
    im = plt.imread(p)
    if im.ndim != 3 or im.shape[2] != 4:
        fatal("figures/%s carries no alpha, so it would paint a box over the "
              "panel" % name)
    return im


def crop(images):
    """One box around the content of both, so the two keep one camera."""
    lo = np.array([1 << 30, 1 << 30])
    hi = np.array([0, 0])
    for im in images:
        ys, xs = np.where(im[..., 3] > 0.03)
        if not len(ys):
            fatal("one render carries no opaque pixel")
        lo = np.minimum(lo, [ys.min(), xs.min()])
        hi = np.maximum(hi, [ys.max(), xs.max()])
    pad = 8
    lo = np.maximum(lo - pad, 0)
    hi = np.minimum(hi + pad + 1, [images[0].shape[0], images[0].shape[1]])
    return [im[lo[0]:hi[0], lo[1]:hi[1]] for im in images]


def place(x0, x1, aspect):
    """The largest rect of that aspect that fits the panel, sitting low."""
    wmax, hmax = (x1 - x0) * 0.95, (TOP - BOT) * 0.70
    w = min(wmax, hmax * aspect * (H / W))
    h = w * (W / H) / aspect
    return [x0 + (x1 - x0 - w) / 2.0, IMG_BOTTOM, w, h]


def text(fig, x, y, s, size=8.4, color=INK, weight="normal", ha="left", z=8):
    fig.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha,
             va="center", zorder=z)


def bilayer(fig, x0, x1, top, half, n=86, z=-6):
    ycen = top - half
    fig.patches.append(FancyBboxPatch(
        (x0, ycen - half), x1 - x0, 2 * half,
        boxstyle="round,pad=0,rounding_size=0.003", transform=fig.transFigure,
        facecolor=LIPID, edgecolor="none", zorder=z))
    for i in range(n):
        x = x0 + (x1 - x0) * (i + 0.5) / n
        fig.lines.append(plt.Line2D(
            [x, x], [ycen - half + 0.005, ycen + half - 0.005],
            transform=fig.transFigure, color=LIPID_EDGE, linewidth=0.45,
            zorder=z + 1))
    for yy in (ycen - half, ycen + half):
        fig.lines.append(plt.Line2D(
            [x0, x1], [yy, yy], transform=fig.transFigure, color=LIPID_EDGE,
            linewidth=1.2, zorder=z + 2))


def main():
    V = bf.load("derived_values.json")
    D, I, M = V["design"], V["interface"], V["md"]
    iface = sorted(I["lead_iface"]["conv_in_iface"])
    if len(iface) != len(D["lead_iface_names"]):
        fatal("the derived file names %d interface conversions and lists %d"
              % (len(D["lead_iface_names"]), len(iface)))
    if D["lead_gravy_native"] >= 0 or D["lead_gravy_conv"] <= 0:
        fatal("the lead hydropathicity does not cross zero, so the membrane "
              "panel must go")

    raw = [load(n) for n, _, _ in PANELS]
    if len({im.shape for im in raw}) != 1:
        fatal("the two renders differ in size, so they were not drawn from "
              "one camera")
    ims = crop(raw)
    aspect = ims[0].shape[1] / float(ims[0].shape[0])

    fig = plt.figure(figsize=(W, H), facecolor="white")
    text(fig, 0.022, 0.940, TITLE, size=13.2, weight="bold")

    rects = [place(x0, x1, aspect) for _, x0, x1 in PANELS]

    if rects[1][1] <= BAND_TOP - 2 * BAND_HALF:
        fatal("the fold reaches below the bilayer, which would draw it as "
              "spanning the membrane instead of resting on it")
    bilayer(fig, BAND_X0, BAND_X1, BAND_TOP, BAND_HALF)

    for (name, x0, x1), im, r in zip(PANELS, ims, rects):
        ax = fig.add_axes(r, zorder=3)
        ax.imshow(im, interpolation="antialiased")
        ax.set_axis_off()
        ax.set_facecolor("none")

    text(fig, 0.040, 0.812, "SOLUBLE OBLIGATE DIMER", size=8.8, weight="bold",
         color=GREY)
    text(fig, 0.040, 0.765, "%s, GRAVY %s"
         % (D["lead"], gravy(D["lead_gravy_native"])), size=8.6)

    text(fig, 0.500, 0.765, "reverse-QTY", size=9.8, weight="bold", ha="center")
    fig.patches.append(FancyArrowPatch(
        (0.432, 0.700), (0.568, 0.700), transform=fig.transFigure,
        arrowstyle="-|>", mutation_scale=17, linewidth=1.8, color=INK,
        zorder=8))
    for i, (a, b, k) in enumerate([("Gln", "Leu", "Q"), ("Thr", "Val", "T"),
                                   ("Tyr", "Phe", "Y")]):
        text(fig, 0.500, 0.612 - 0.062 * i, "%s to %s" % (a, b), size=9.6,
             weight="bold", color=CLASS[k], ha="center")
    text(fig, 0.500, 0.376, "%d positions" % D["lead_conv_n"], size=8.4,
         color=GREY, ha="center")
    text(fig, 0.500, 0.334, "%d in the interface" % len(iface), size=8.4,
         color=GREY, ha="center")

    text(fig, 0.604, 0.812, "MEMBRANE-COMPATIBLE TARGET", size=8.8,
         weight="bold", color=CONVERTED)
    text(fig, 0.604, 0.765, "GRAVY %s, interface holds over %g ns"
         % (gravy(D["lead_gravy_conv"]), M["ns"]["native"]), size=8.6)
    text(fig, 0.604, 0.722, "membrane behaviour predicted, the runs are in "
         "water", size=7.4, color=GREY)

    for ext in ("png", "svg"):
        p = os.path.join(OUT, "graphical_abstract." + ext)
        fig.savefig(p, dpi=bf.DPI if ext == "png" else None, facecolor="white")
        print("wrote figures/graphical_abstract." + ext)
    plt.close(fig)

    toc(ims, aspect, D, I)
    return 0


TOC_W, TOC_H = 3.25, 1.75


def toc(ims, aspect, D, I):
    fig = plt.figure(figsize=(TOC_W, TOC_H), facecolor="white")
    gap, top, bottom = 0.015, 0.80, 0.175
    w = (0.5 - gap) * 0.86
    h = w * (TOC_W / TOC_H) / aspect
    if h > top - bottom:
        h = top - bottom
        w = h * aspect * (TOC_H / TOC_W)
    band_top = bottom + h * 0.30

    fig.patches.append(FancyBboxPatch(
        (0.515, band_top - h * 0.30), 0.47, h * 0.30,
        boxstyle="round,pad=0,rounding_size=0.004", transform=fig.transFigure,
        facecolor=LIPID, edgecolor=LIPID_EDGE, linewidth=0.7, zorder=-6))
    for i in range(44):
        x = 0.515 + 0.47 * (i + 0.5) / 44
        fig.lines.append(plt.Line2D(
            [x, x], [band_top - h * 0.30 + 0.012, band_top - 0.012],
            transform=fig.transFigure, color=LIPID_EDGE, linewidth=0.3,
            zorder=-5))

    for im, x0 in zip(ims, (0.5 - gap - w, 0.5 + gap)):
        ax = fig.add_axes([x0, bottom, w, h], zorder=3)
        ax.imshow(im, interpolation="antialiased")
        ax.set_axis_off()

    fig.text(0.5, 0.945, "reverse-QTY", ha="center", va="center", fontsize=7.4,
             fontweight="bold", color=INK)
    fig.text(0.5, 0.862, "Gln to Leu, Thr to Val, Tyr to Phe", ha="center",
             va="center", fontsize=6.0, color=GREY)
    fig.text(0.5, 0.072, "%d of %d converted positions lie in the dimer "
             "interface" % (len(I["lead_iface"]["conv_in_iface"]),
                            D["lead_conv_n"]),
             ha="center", va="center", fontsize=6.0, color=INK)
    fig.patches.append(FancyArrowPatch(
        (0.455, 0.52), (0.545, 0.52), transform=fig.transFigure,
        arrowstyle="-|>", mutation_scale=7, linewidth=1.0, color=INK,
        zorder=9))

    p = os.path.join(OUT, "graphical_abstract_toc.png")
    fig.savefig(p, dpi=600, facecolor="white")
    plt.close(fig)
    print("wrote figures/graphical_abstract_toc.png  %.2f x %.2f in"
          % (TOC_W, TOC_H))


if __name__ == "__main__":
    sys.exit(main())
