#!/usr/bin/env python3
"""Plot bulk-free V(D)J deletion evidence from the public Zenodo tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


HERE = Path(__file__).resolve().parents[1]
SITE_ORDER = ["TRB", "TRA_TRD", "TRG", "IGH", "IGK", "IGL"]
COLORS = {"T cell": "#D55E00", "B cell": "#0072B2", "Other": "#9CA3AF"}
MARKERS = {"SniffCell": "o", "FANS": "^"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def jitter(*parts: object) -> float:
    digest = hashlib.sha256("|".join(map(str, parts)).encode()).digest()
    return (int.from_bytes(digest[:4], "big") / 2**32 - 0.5) * 0.14


def load_events(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(path)
    required = {
        "site", "source_group", "expected_celltype", "expected_status",
        "deletion_len_bp", "read_name", "read_length", "mapq", "evidence",
        "split_is_strict",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Missing required columns in {path}")
    if "bam" in rows[0]:
        raise ValueError(f"Public input unexpectedly contains a BAM-path column: {path}")
    kept = []
    for row in rows:
        if row["source_group"] not in MARKERS:
            continue
        if int(row["deletion_len_bp"]) < 1000 or int(row["mapq"]) < 20:
            continue
        if int(row["read_length"]) < 1000:
            continue
        if row["evidence"].lower().startswith("split") and row["split_is_strict"].lower() != "true":
            continue
        kept.append(row)
    return kept


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--events-dir", type=Path, required=True)
    parser.add_argument("--out-prefix", type=Path, default=HERE / "results" / "figure_3ab_vdj_rearrangements")
    args = parser.parse_args()

    ranges = {row["site"]: row for row in read_tsv(HERE / "data" / "expected_deletion_ranges.tsv")}
    samples = []
    for sample in ("pbmc1", "pbmc2"):
        path = args.events_dir / f"{sample}_vdj_events.tsv"
        samples.append((sample.upper(), load_events(path)))
    if not all(rows for _, rows in samples):
        raise ValueError("Both PBMC event tables must contain plotted events")

    fig, axes = plt.subplots(2, 1, figsize=(12.4, 7.2), sharex=True, sharey=True)
    offsets = {"SniffCell": -0.18, "FANS": 0.18}
    for ax, (sample, rows) in zip(axes, samples):
        ax.axvspan(-0.5, 2.5, color="#E0F2FE", alpha=0.18)
        ax.axvspan(2.5, 5.5, color="#FCE7F3", alpha=0.14)
        ax.axvline(2.5, color="#94A3B8", ls=":", lw=0.9)
        for x, site in enumerate(SITE_ORDER):
            low, high = int(ranges[site]["minimum_bp"]), int(ranges[site]["maximum_bp"])
            ax.fill_between([x - 0.42, x + 0.42], low, high, color="#111827", alpha=0.075)
            ax.hlines([low, high], x - 0.42, x + 0.42, color="#111827", alpha=0.34, lw=0.7)
        for row in rows:
            source = row["source_group"]
            expected = row["expected_status"] == "Expected cell type"
            color = COLORS.get(row["expected_celltype"], COLORS["Other"]) if expected else COLORS["Other"]
            baseline = row["evidence"] == "SniffCell_other_local_read"
            x = SITE_ORDER.index(row["site"]) + offsets[source] + jitter(
                sample, source, row["site"], row["read_name"], row["deletion_len_bp"]
            )
            ax.scatter(x, int(row["deletion_len_bp"]), s=8 if baseline else 28,
                       marker=MARKERS[source], color=color, alpha=0.24 if baseline else 0.82,
                       edgecolor="black" if expected and not baseline else "none", linewidth=0.25)
        ax.set_yscale("log")
        ax.set_ylim(800, 1_400_000)
        ax.set_xlim(-0.5, 5.5)
        ax.grid(axis="y", which="major", alpha=0.22)
        ax.set_title(sample, fontweight="bold")
    axes[-1].set_xticks(range(6), SITE_ORDER)
    axes[-1].set_xlabel("Recombination site")
    fig.supylabel("Deletion size (bp, log scale)")
    handles = [
        Line2D([0], [0], marker="o", linestyle="", color=COLORS["T cell"], label="T cell"),
        Line2D([0], [0], marker="o", linestyle="", color=COLORS["B cell"], label="B cell"),
        Line2D([0], [0], marker="o", linestyle="", color=COLORS["Other"], label="Other cell types"),
        Line2D([0], [0], marker="o", linestyle="", color="black", label="SniffCell"),
        Line2D([0], [0], marker="^", linestyle="", color="black", label="FACS"),
        Patch(facecolor="#111827", alpha=0.16, label="Expected range"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.94), ncol=6)
    fig.subplots_adjust(top=0.82, bottom=0.10, left=0.09, right=0.98, hspace=0.34)
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(args.out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    counts = Counter((sample, row["site"], row["source_group"], row["expected_status"]) for sample, rows in samples for row in rows)
    with args.out_prefix.with_name(args.out_prefix.name + "_counts.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "site", "source", "expected_status", "event_count"])
        for key, count in sorted(counts.items()):
            writer.writerow([*key, count])


if __name__ == "__main__":
    main()
