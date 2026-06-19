#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# =========================
# 参数
# =========================
BASE_DIR = "pca_results"

COMBINED_GROUPS = [
    "404rice",
    "1439rice",
    "176rice",
    "378rice",
    "3023rice",
    "532rice"
]

SINGLE_GROUPS = [
    "1171rice",
    "705rice",
    "6048rice"
]

COMBINED_FIGSIZE = (8, 12)   # 3x2
SINGLE_FIGSIZE = (4,3)

POINT_SIZE = 18
POINT_ALPHA = 0.8
DPI = 400

# =========================
# 字体
# =========================
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["axes.linewidth"] = 0.8

# =========================
# 读取 PCA
# =========================

def read_pca(group):

    file = f"{BASE_DIR}/{group}/{group}.pca.eigenvec"

    df = pd.read_csv(file, sep=r"\s+", header=None)

    n_pc = df.shape[1] - 2
    df.columns = ["FID","IID"] + [f"PC{i}" for i in range(1,n_pc+1)]

    return df


def read_eigenval(group):

    file = f"{BASE_DIR}/{group}/{group}.pca.eigenval"

    if not os.path.exists(file):
        return None

    eig = pd.read_csv(file, header=None)[0].values

    var = eig / eig.sum() * 100

    return var


# =========================
# PCA绘图函数
# =========================
def plot_pca(ax, df, group):

    ax.scatter(
        df["PC1"],
        df["PC2"],
        s=18,
        facecolor="black",     # 点填充
        edgecolor="black",     # 黑色描边
        linewidth=0.25,
        alpha=0.35             # 透明度形成密度
    )

    var = read_eigenval(group)

    if var is not None:
        ax.set_xlabel(f"PC1 ({var[0]:.2f}%)")
        ax.set_ylabel(f"PC2 ({var[1]:.2f}%)")
    else:
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")

    ax.set_title(group, fontsize=11)

    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)

    # 四面框
    for spine in ax.spines.values():
        spine.set_visible(True)

# =========================
# 6个群体合并图
# =========================

fig, axes = plt.subplots(3,2, figsize=COMBINED_FIGSIZE)

axes = axes.flatten()

for i, group in enumerate(COMBINED_GROUPS):

    df = read_pca(group)

    plot_pca(axes[i], df, group)

for j in range(len(COMBINED_GROUPS), len(axes)):
    axes[j].axis("off")

plt.tight_layout()

plt.savefig("PCA_6populations.pdf", bbox_inches="tight")
plt.savefig("PCA_6populations.png", dpi=DPI, bbox_inches="tight")

plt.close()


# =========================
# 1171rice 和 705rice 单独图
# =========================

for group in SINGLE_GROUPS:

    df = read_pca(group)

    fig, ax = plt.subplots(figsize=SINGLE_FIGSIZE)

    plot_pca(ax, df, group)

    plt.tight_layout()

    plt.savefig(f"PCA_{group}.pdf", bbox_inches="tight")
    plt.savefig(f"PCA_{group}.png", dpi=DPI, bbox_inches="tight")

    plt.close()


print("All PCA plots finished.")