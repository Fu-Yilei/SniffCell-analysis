#!/usr/bin/env python3
"""Render Figure 2e-g from Supplementary Table 5 aggregate tables."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SAMPLES = ["pbmc1", "pbmc2", "bcontrol1", "bcontrol2", "bmsa", "bparkinsons"]
DEPTHS = ["30x", "60x", "full"]


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save(fig, outdir, name):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(outdir / f"{name}.{ext}", dpi=300, bbox_inches="tight", pad_inches=.1)
    plt.close(fig)


def render(data, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 18,
        "font.weight": "normal", "axes.labelweight": "normal", "axes.labelsize": 22,
        "xtick.labelsize": 18, "ytick.labelsize": 21,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    rows = table(data / "figure_2e_counts.tsv")
    lookup = {(r["sample"], r["coverage"]): r for r in rows}
    fig, axes = plt.subplots(1, 6, figsize=(11.26, 5.76), sharey=True)
    for ax, sample in zip(axes, SAMPLES):
        bottom = [0, 0, 0]
        for kind, color in (("SV", "#ff7f0e"), ("TR", "#4c78a8")):
            values = [int(lookup[sample, cov][kind]) for cov in DEPTHS]
            ax.bar(range(3), values, bottom=bottom, width=.72, color=color, label=kind + "s")
            bottom = [a + b for a, b in zip(bottom, values)]
        for i, total in enumerate(bottom):
            if total:
                ax.text(i, total + .45, str(total), ha="center", va="bottom", fontsize=14)
        # Sample facet labels identify data, not a global figure title.
        ax.set_title(sample, fontsize=17, pad=5, fontweight="normal")
        ax.set_xticks(range(3), DEPTHS, rotation=45, ha="right", fontsize=13)
        ax.tick_params(axis="y", labelsize=13)
        ax.set_ylim(0, 36)
        ax.set_yticks(range(0, 36, 5))
        ax.grid(axis="y", color="#e6e6e6")
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Variant count", fontsize=16)
    axes[0].legend(title="Variant class", fontsize=12, title_fontsize=13, loc="upper left")
    fig.subplots_adjust(left=.068, right=.985, bottom=.105, top=.94, wspace=.16)
    save(fig, outdir, "figure_2e_variant_counts")

    for panel in ("f", "g"):
        is_f = panel == "f"
        rows = table(data / ("figure_2f_validation.tsv" if is_f else "figure_2g_composition.tsv"))
        fig, ax = plt.subplots(figsize=(12, 4.8) if is_f else (12, 4.5))
        groups = ["TR", "SV"] if is_f else ["pbmc", "brain"]
        group_key = "variant_class" if is_f else "tissue"
        ticks = []
        for y, group in zip((1, 0), groups):
            parts = [r for r in rows if r[group_key] == group]
            left = 0
            for r in parts:
                percent, count = float(r["percent"]), int(r["count"])
                if is_f:
                    observed = r["status"] == "Yes"
                    label = "Observed" if observed else "N.O." if percent < 20 else "Not observed"
                    color, ink = ("#1b9e77", "white") if observed else ("#bdbdbd", "black")
                else:
                    label = {"T-cell+NK-cell+B-cell": "Lymphocyte group", "Monocyte": "Mono"}.get(r["cell_type"], r["cell_type"])
                    color = {"Lymphocyte group": "#d95f02", "Mono": "#7b3294",
                             "Neuron": "#1b9e77", "Oligodendrocyte": "#0072b2"}[label]
                    ink = "white"
                ax.barh(y, percent, left=left, height=.62, color=color, edgecolor="white", linewidth=1.5)
                if count:
                    ax.text(left + percent / 2, y, f"{label}\n{count} ({percent:.1f}%)",
                            ha="center", va="center", color=ink, fontweight="semibold",
                            fontsize=(16 if percent < 20 else 19) if is_f else
                                     (14 if percent < 10 else 19))
                left += percent
            if abs(left - 100) > 1e-7:
                raise ValueError("Composition does not sum to 100%")
            display = group if is_f else "PBMC" if group == "pbmc" else "Brain"
            ticks.append(f"{display}\nn={parts[0]['total']}")
        ax.set_yticks((1, 0), ticks)
        ax.set_xlim(0, 100)
        ax.set_xticks(range(0, 101, 20))
        ax.set_xlabel("Variant calls (%)")
        ax.grid(axis="x", color="#e6e6e6")
        ax.set_axisbelow(True)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0, pad=14, labelsize=22 if is_f else 21)
        fig.subplots_adjust(left=.16 if is_f else .15, right=.985, bottom=.24, top=.96)
        save(fig, outdir, "figure_2" + panel + ("_fans_validation" if is_f else "_celltype_composition"))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    a = p.parse_args()
    render(a.data_dir, a.outdir)
