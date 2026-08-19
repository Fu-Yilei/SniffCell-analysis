#!/usr/bin/env python3
"""Plot Monocyte-like ctDMR signal among conflicting brain reads."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
INPUT = ANALYSIS_DIR / "data" / "brain_conflicting_reads_monocyte_like_summary.tsv"
OUTBASE = ANALYSIS_DIR / "results" / "supplementary_figure_6"

SAMPLE_LABELS = {
    "b1": "bcontrol1",
    "b4": "bcontrol2",
    "b2": "bmsa",
    "b3": "bparkinsons",
}
SAMPLE_ORDER = ["b1", "b4", "b2", "b3"]


def load_rows() -> list[dict[str, str]]:
    with INPUT.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_sample = {row["sample"]: row for row in rows}
    missing = set(SAMPLE_ORDER) - set(by_sample)
    if missing:
        raise ValueError(f"Missing expected samples: {sorted(missing)}")
    return [by_sample[sample] for sample in SAMPLE_ORDER]


def main() -> None:
    rows = load_rows()
    labels = [SAMPLE_LABELS[row["sample"]] for row in rows]
    values = [float(row["pct_conflicting_any_monocyte_like"]) for row in rows]
    numerators = [int(row["conflicting_reads_any_monocyte_like"]) for row in rows]
    denominators = [int(row["conflicting_reads"]) for row in rows]

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    bars = ax.bar(labels, values, color="#6F8F72", edgecolor="#263326", linewidth=0.8)
    ax.set_ylabel("Conflicting reads with any\nMonocyte-like ctDMR row (%)")
    ax.set_xlabel("Brain sample")
    ax.set_ylim(0, 72)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.8)
    ax.grid(axis="x", visible=False)
    ax.spines[["top", "right"]].set_visible(False)
    for bar, value, numerator, denominator in zip(bars, values, numerators, denominators):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 1.2,
            f"{value:.1f}%\n{numerator:,}/{denominator:,}",
            ha="center",
            va="bottom",
            fontsize=8.2,
        )

    OUTBASE.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUTBASE.with_suffix(".png"), dpi=300)
    fig.savefig(OUTBASE.with_suffix(".pdf"))
    plt.close(fig)
    print(f"Wrote {OUTBASE.with_suffix('.png')}")
    print(f"Wrote {OUTBASE.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
