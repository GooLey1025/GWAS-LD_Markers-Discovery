#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.path import Path
from matplotlib.colors import PowerNorm
from collections import defaultdict

# =========================
# Config
# =========================
BIN_SIZE = 400_000

FIG_SIZE = (15, 8)
N_COLS = 12
COL_GAP = 2.0
ROW_GAP = 2.4

CHR_W = 0.24
CEN_W = 0.10
BIN_H = 0.08
HOURGLASS_POWER = 1.6

MAX_PER_ROW = 10
DOT_SIZE = 48
DOT_EDGE_LW = 0.4
DOT_X_STEP = 0.10
DOT_Y_STEP = 0.10
DOT_X_OFFSET = 0.60

LINK_LW = 0.30
LINK_COLOR = (0, 0, 0, 0.5)

CHR_LABEL_PAD = 0.55
BOTTOM_SAFE = 0.55

DENSITY_CMAP = plt.cm.Greys
DENSITY_GAMMA = 0.45   # 越小，中低值越容易变黑；可试 0.35–0.6

OUT_PNG = "karyotype.png"
OUT_PDF = "karyotype.pdf"
DPI = 300
TOP_BAND = 0.90

# 左侧 Mb 标尺
SCALE_MB_MAX = 50
SCALE_MB_STEP = 10
SCALE_X = -0.75


# =========================
# Helpers
# =========================
def parse_marker(m):
    s = str(m)
    p = s.split("-")
    if len(p) < 3:
        return None, None, None
    try:
        return p[0], int(p[1]), int(p[2])
    except Exception:
        return None, None, None


def build_colors(cats):
    base = ['#7CB87C', '#6B7CB3', '#E6B84D', '#C75B5B', '#9B7CB6',
            '#8B7CB6', '#D4E4C4', '#4A7CB8', '#7CB8B8']
    cmap = plt.get_cmap("tab20")
    d = {}
    for i, c in enumerate(cats):
        d[c] = base[i] if i < len(base) else cmap(i % 20)
    return d


def chrom_path_hourglass(xc, y0, y1, cy0, cy1, wf, wc,
                         power=1.6, n_body=220, n_arc=60):
    wf = float(wf)
    wc = float(wc)
    r = wf / 2.0

    if (y1 - y0) < 2 * r:
        y1 = y0 + 2 * r + 1e-6

    y_bot = y0
    y_top = y1 - r

    def width_at(y):
        if cy1 <= cy0 or (y < cy0) or (y > cy1):
            return wf
        t = (y - cy0) / (cy1 - cy0)
        s = abs(2 * t - 1)
        return wc + (wf - wc) * (s ** power)

    ys = np.linspace(y_bot, y_top, n_body)

    right = [(xc + width_at(y) / 2.0, y) for y in ys]
    left = [(xc - width_at(y) / 2.0, y) for y in ys[::-1]]

    verts = []
    codes = []

    theta = np.linspace(np.pi, 0, n_arc)
    for i, t in enumerate(theta):
        x = xc + r * np.cos(t)
        y = y0 - r * np.sin(t)
        verts.append((x, y))
        codes.append(Path.MOVETO if i == 0 else Path.LINETO)

    for x, y in right:
        verts.append((x, y))
        codes.append(Path.LINETO)

    cap_cy = y1 - r
    theta = np.linspace(0, np.pi, n_arc)
    for t in theta:
        x = xc + r * np.cos(t)
        y = cap_cy + r * np.sin(t)
        verts.append((x, y))
        codes.append(Path.LINETO)

    for x, y in left:
        verts.append((x, y))
        codes.append(Path.LINETO)

    verts.append(verts[0])
    codes.append(Path.CLOSEPOLY)

    return Path(verts, codes)


def panel_xy(i, panel_h):
    r = i // N_COLS
    c = i % N_COLS
    return c * COL_GAP, -r * (panel_h + ROW_GAP)


def avoid_overlap(ranges, y, half, step):
    for k in range(40):
        shift = ((k + 1) // 2) * step
        if k == 0:
            yy = y
        else:
            yy = y + shift if (k % 2 == 1) else y - shift

        ok = True
        for a, b in ranges:
            if not (yy + half < a or yy - half > b):
                ok = False
                break
        if ok:
            ranges.append((yy - half, yy + half))
            return yy

    ranges.append((y - half, y + half))
    return y


def draw_mb_scale(ax, x, y0, mb_max=50, mb_step=10):
    tick_len = 0.12

    y_start = y0
    y_end = y0 + (mb_max * 1_000_000 / BIN_SIZE) * BIN_H

    ax.plot([x, x], [y_start, y_end], color="black", lw=0.8, zorder=10)

    for mb in range(0, mb_max + 1, mb_step):
        y = y0 + (mb * 1_000_000 / BIN_SIZE) * BIN_H
        ax.plot([x - tick_len / 2, x + tick_len / 2],
                [y, y], color="black", lw=0.8, zorder=10)

        ax.text(x - 0.12, y, f"{mb}",
                ha="right", va="center", fontsize=8)

    ax.text(x, y_end + 0.35, "Mb",
            ha="center", va="bottom", fontsize=9)


# =========================
# Load data
# =========================
chr_df = pd.read_csv("chr_length_cen_annotate.tsv", sep=r"\s+|\t+", engine="python")
chr_df["chr"] = chr_df["chr"].astype(int)
chr_df = chr_df.sort_values("chr")

df = pd.read_csv("all_lead_markers_final.annotated.noredu.tsv", sep="\t")
df["vtype"], df["chr"], df["pos"] = zip(*df["Marker"].map(parse_marker))
df = df.dropna(subset=["chr", "pos"])

df["chr"] = df["chr"].astype(int)
df["pos"] = df["pos"].astype(int)
df["bin"] = df["pos"] // BIN_SIZE

bin_counts = df.groupby(["chr", "bin"]).size().reset_index(name="count")

vmax = int(bin_counts["count"].max()) if not bin_counts.empty else 1
norm = PowerNorm(gamma=DENSITY_GAMMA, vmin=0, vmax=vmax)

cats = sorted(df["Category"].unique())
cat_color = build_colors(cats)

chr_bins = {r.chr: math.ceil(r.length / BIN_SIZE) for _, r in chr_df.iterrows()}
max_bins = max(chr_bins.values())
panel_h = max_bins * BIN_H + 1.2


# =========================
# Figure
# =========================
fig, ax = plt.subplots(figsize=FIG_SIZE, dpi=DPI)
ax.axis("off")

occupied = defaultdict(list)


# =========================
# Draw Mb scale
# =========================
draw_mb_scale(
    ax,
    x=SCALE_X,
    y0=0,
    mb_max=SCALE_MB_MAX,
    mb_step=SCALE_MB_STEP
)


# =========================
# Draw chromosomes
# =========================
for i, row in chr_df.iterrows():
    chrom = int(row.chr)
    L = int(row.length)
    nb = math.ceil(L / BIN_SIZE)

    x0, y0 = panel_xy(i, panel_h)
    y1 = y0 + nb * BIN_H
    xc = x0 + CHR_W / 2

    cy0 = y0 + (min(row.centromere_start, row.centromere_end) / BIN_SIZE) * BIN_H
    cy1 = y0 + (max(row.centromere_start, row.centromere_end) / BIN_SIZE) * BIN_H

    p = chrom_path_hourglass(
        xc, y0, y1, cy0, cy1,
        wf=CHR_W,
        wc=CEN_W,
        power=HOURGLASS_POWER
    )

    patch = PathPatch(p, fc="white", ec="black", lw=0.9, zorder=1)
    ax.add_patch(patch)

    cen_band = Rectangle(
        (x0 - 2.0, cy0),
        6.0,
        max(0.02, cy1 - cy0),
        fc="#E6E6E6",
        ec="none",
        zorder=2
    )
    cen_band.set_clip_path(patch)
    ax.add_patch(cen_band)

    ax.text(
        xc,
        y0 - CHR_LABEL_PAD,
        f"Chr{chrom}",
        ha="center",
        va="top",
        fontsize=9
    )

    sub = bin_counts[bin_counts.chr == chrom]

    for _, r in sub.iterrows():
        b = int(r.bin)
        by = y0 + b * BIN_H
        yc = by + BIN_H / 2
        col = DENSITY_CMAP(norm(r["count"]))

        rect = Rectangle(
            (x0, by),
            CHR_W,
            BIN_H,
            fc=col,
            ec="none",
            zorder=3
        )
        rect.set_clip_path(patch)
        ax.add_patch(rect)

        dots = df[(df.chr == chrom) & (df.bin == b)]
        n = len(dots)

        if n == 0:
            continue

        rows = math.ceil(n / MAX_PER_ROW)
        half = max(0.06, (rows - 1) / 2 * DOT_Y_STEP) + 0.05

        y_center = yc
        if (y_center - y0) < BOTTOM_SAFE:
            y_center += BOTTOM_SAFE - (y_center - y0) + 0.12

        y_center = avoid_overlap(
            occupied[chrom],
            y_center,
            half,
            step=DOT_Y_STEP * 1.3
        )

        xs = x0 + CHR_W
        xb = xs + DOT_X_OFFSET

        ax.plot(
            [xs, xb - 0.06],
            [yc, y_center],
            lw=LINK_LW,
            color=LINK_COLOR,
            zorder=4
        )

        for j, cat in enumerate(dots.Category.tolist()):
            rr = j // MAX_PER_ROW
            cc = j % MAX_PER_ROW

            dx = xb + cc * DOT_X_STEP
            dy = y_center + (rr - (rows - 1) / 2) * DOT_Y_STEP

            ax.scatter(
                dx,
                dy,
                s=DOT_SIZE,
                c=cat_color[cat],
                edgecolors="black",
                linewidths=DOT_EDGE_LW,
                zorder=5
            )


# =========================
# Legends
# =========================
cat_counts_map = df["Category"].value_counts().to_dict()
labels = [f"{c} ({cat_counts_map.get(c, 0)})" for c in cats]

handles = [
    plt.Line2D(
        [0],
        [0],
        marker="o",
        linestyle="",
        markerfacecolor=cat_color[c],
        markeredgecolor="black",
        markersize=4.8
    )
    for c in cats
]

fig.legend(
    handles,
    labels,
    loc="upper left",
    bbox_to_anchor=(0.01, 1.02),
    ncol=math.ceil(len(cats) / 2),
    frameon=False,
    fontsize=9,
    handletextpad=0.4,
    columnspacing=0.8
)

sm = plt.cm.ScalarMappable(cmap=DENSITY_CMAP, norm=norm)
sm.set_array([])

cax = fig.add_axes([0.82, 0.935, 0.16, 0.018])
cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
cbar.set_label("Lead variants per 0.4Mb bin", fontsize=8, labelpad=2)
cbar.ax.tick_params(labelsize=7, length=2)


# =========================
# Layout and save
# =========================
ax.set_xlim(SCALE_X - 0.45, (N_COLS - 1) * COL_GAP + 2.0)

plt.tight_layout(rect=[0, 0, 1, TOP_BAND])

fig.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight")
fig.savefig(OUT_PDF, bbox_inches="tight")

print(f"Done. Saved: {OUT_PNG}, {OUT_PDF}")