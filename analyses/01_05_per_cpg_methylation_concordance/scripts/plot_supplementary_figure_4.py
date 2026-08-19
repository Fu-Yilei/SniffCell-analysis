#!/usr/bin/env python3
"""Recompute or render pooled per-CpG SniffCell-versus-FANS concordance."""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path
from typing import Iterator

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


NBINS = 120
COVERAGE_THRESHOLD = 5
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
CHROM_ORDER = {chrom: index for index, chrom in enumerate(CHROMS)}


class Pooled:
    def __init__(self) -> None:
        self.hist = np.zeros((NBINS, NBINS), dtype=np.int64)
        self.n = 0
        self.sx = self.sy = self.sxy = self.sxx = self.syy = self.sabs = 0.0

    def add(self, x: float, y: float) -> None:
        self.hist[min(int(x * NBINS), NBINS - 1), min(int(y * NBINS), NBINS - 1)] += 1
        self.n += 1
        self.sx += x
        self.sy += y
        self.sxy += x * y
        self.sxx += x * x
        self.syy += y * y
        self.sabs += abs(x - y)

    @property
    def pearson(self) -> float:
        denominator = math.sqrt(
            (self.n * self.sxx - self.sx**2) * (self.n * self.syy - self.sy**2)
        )
        return (self.n * self.sxy - self.sx * self.sy) / denominator if denominator else math.nan

    @property
    def mae(self) -> float:
        return self.sabs / self.n if self.n else math.nan

    def save(self, path: Path) -> None:
        np.savez(
            path,
            hist=self.hist,
            n=self.n,
            sx=self.sx,
            sy=self.sy,
            sxy=self.sxy,
            sxx=self.sxx,
            syy=self.syy,
            sabs=self.sabs,
        )

    @classmethod
    def load(cls, path: Path) -> "Pooled":
        result = cls()
        data = np.load(path)
        result.hist = data["hist"]
        result.n = int(data["n"])
        result.sx, result.sy, result.sxy = float(data["sx"]), float(data["sy"]), float(data["sxy"])
        result.sxx, result.syy, result.sabs = float(data["sxx"]), float(data["syy"]), float(data["sabs"])
        return result


def records(path: Path) -> Iterator[tuple[tuple[int, int, str], int, int]]:
    with gzip.open(path, "rt") as handle:
        for line in handle:
            chrom, position, coverage, modified = line.rstrip("\n").split("\t")
            key = (CHROM_ORDER.get(chrom, len(CHROM_ORDER)), int(position), chrom)
            yield key, int(coverage), int(modified)


def combined_records(paths: list[Path]) -> Iterator[tuple[tuple[int, int, str], int, int]]:
    iterators = [iter(records(path)) for path in paths]
    current = [next(iterator, None) for iterator in iterators]
    while any(record is not None for record in current):
        key = min(record[0] for record in current if record is not None)
        coverage = modified = 0
        for index, record in enumerate(current):
            if record is not None and record[0] == key:
                coverage += record[1]
                modified += record[2]
                current[index] = next(iterators[index], None)
        yield key, coverage, modified


def add_comparison(pool: Pooled, bulk: Path, fans: list[Path]) -> None:
    bulk_iter = iter(records(bulk))
    fans_iter = iter(combined_records(fans))
    left, right = next(bulk_iter, None), next(fans_iter, None)
    while left is not None and right is not None:
        if left[0] == right[0]:
            if left[1] > COVERAGE_THRESHOLD and right[1] > COVERAGE_THRESHOLD:
                pool.add(left[2] / left[1], right[2] / right[1])
            left, right = next(bulk_iter, None), next(fans_iter, None)
        elif left[0] < right[0]:
            left = next(bulk_iter, None)
        else:
            right = next(fans_iter, None)


def comparison_paths(root: Path, dataset: str) -> list[tuple[Path, list[Path]]]:
    comparisons: list[tuple[Path, list[Path]]] = []
    if dataset == "brain":
        for sample in ("b1", "b2", "b3", "b4"):
            for group in ("Neuron", "Oligodendrocyte"):
                comparisons.append(
                    (root / f"brain/{sample}/bulk/{group}.cpg_simple.tsv.gz", [root / f"brain/{sample}/fans/{group}.cpg_simple.tsv.gz"])
                )
    else:
        for sample in ("pbmc1", "pbmc2"):
            comparisons.append(
                (root / f"pbmc/{sample}/bulk/Monocyte.cpg_simple.tsv.gz", [root / f"pbmc/{sample}/fans/Monocyte.cpg_simple.tsv.gz"])
            )
            comparisons.append(
                (
                    root / f"pbmc/{sample}/bulk/Lymphocyte.cpg_simple.tsv.gz",
                    [root / f"pbmc/{sample}/fans/{group}.cpg_simple.tsv.gz" for group in ("T-cell", "NK-cell", "B-cell")],
                )
            )
    return comparisons


def render(pool: Pooled, output: Path, title: str) -> None:
    centers = (np.arange(NBINS) + 0.5) / NBINS
    grid_x, grid_y = np.meshgrid(centers, centers, indexing="ij")
    pct_gt04 = 100.0 * pool.hist[np.abs(grid_x - grid_y) > 0.4].sum() / pool.hist.sum()
    fig, axis = plt.subplots(figsize=(6.2, 5.8))
    image = axis.imshow(
        pool.hist.T,
        origin="lower",
        extent=[0, 1, 0, 1],
        aspect="equal",
        cmap="viridis",
        norm=LogNorm(vmin=1, vmax=int(pool.hist.max())),
        interpolation="nearest",
    )
    axis.plot([0, 1], [0, 1], "--", color="white", lw=1.0, alpha=0.7)
    axis.plot([0.0, 0.6], [0.4, 1.0], "--", color="#d62728", lw=1.6, label="ctDMR threshold (|Δβ| = 0.4)")
    axis.plot([0.4, 1.0], [0.0, 0.6], "--", color="#d62728", lw=1.6)
    axis.legend(loc="lower right", fontsize=8, framealpha=0.85)
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_xlabel("Bulk SniffCell-deconvolved methylation beta")
    axis.set_ylabel("Cell-sorted reference methylation beta")
    text = (
        f"shared CpGs (>5x): {pool.n:,}\n"
        f"Pearson r: {pool.pearson:.3f}\n"
        f"MAE: {pool.mae:.3f}\n"
        f"|Δβ| > 0.4: {pct_gt04:.2f}%"
    )
    axis.text(0.04, 0.96, text, transform=axis.transAxes, va="top", ha="left", fontsize=10, bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})
    colorbar = fig.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
    colorbar.set_label("Shared CpG density (log scale)")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf"):
        fig.savefig(output.with_suffix(f".{extension}"), dpi=300)
    plt.close(fig)


def main() -> int:
    analysis_dir = Path(__file__).resolve().parents[1]
    default_shared = analysis_dir.parent / "shared_data/modkit_cpg_pileups_01_04_to_01_06/data/sniffcell_modkit_cpg_pileups"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pileup-root", type=Path, default=default_shared)
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()

    titles = {
        "brain": "Brain (b1–b4 pooled): deconvolved vs cell-sorted CpG methylation",
        "pbmc": "PBMC (pbmc1–pbmc2 pooled): deconvolved vs cell-sorted CpG methylation",
    }
    for dataset in ("brain", "pbmc"):
        cache = analysis_dir / f"data/{dataset}_concordance_pooled.npz"
        if args.rebuild:
            pool = Pooled()
            comparisons = comparison_paths(args.pileup_root, dataset)
            for index, (bulk, fans) in enumerate(comparisons, start=1):
                print(f"{dataset}: comparison {index}/{len(comparisons)} ({bulk.parent.parent.name}/{bulk.stem})", flush=True)
                add_comparison(pool, bulk, fans)
            pool.save(cache)
        else:
            pool = Pooled.load(cache)
        print(f"{dataset}: n={pool.n:,}, r={pool.pearson:.4f}, MAE={pool.mae:.4f}")
        render(pool, analysis_dir / f"results/supplementary_figure_4_{dataset}", titles[dataset])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
