#!/usr/bin/env python3
"""Build public PBMC aggregate variant-count tables and Figure 2f component."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = ANALYSIS_DIR / "data"
RESULTS_DIR = ANALYSIS_DIR / "results"

GERMLINE_INPUT = DATA_DIR / "pbmc_bulk_germline_counts.tsv"
CELLTYPE_INPUT = DATA_DIR / "pbmc_celltype_associated_variant_counts.tsv"
MANUSCRIPT_INPUT = DATA_DIR / "manuscript_reported_counts.tsv"

EXPECTED_SAMPLES = {"pbmc1", "pbmc2"}
EXPECTED_CLASSES = {"TR", "SV"}
EXPECTED_ASSIGNMENTS = {"Lymphocyte group", "Monocyte"}
COLORS = {"Lymphocyte group": "#D55E00", "Monocyte": "#0072B2"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_nonnegative_int(value: str, label: str) -> int:
    number = int(value)
    if number < 0:
        raise ValueError(f"{label} must be non-negative: {value}")
    return number


def load_germline() -> list[dict[str, object]]:
    rows = read_tsv(GERMLINE_INPUT)
    if {row["sample_id"] for row in rows} != EXPECTED_SAMPLES:
        raise ValueError("Germline input must contain exactly pbmc1 and pbmc2")
    if len(rows) != len(EXPECTED_SAMPLES):
        raise ValueError("Germline input contains duplicate sample rows")

    output: list[dict[str, object]] = []
    for row in sorted(rows, key=lambda item: item["sample_id"]):
        output.append(
            {
                "sample_id": row["sample_id"],
                "tissue": row["tissue"],
                "small_variant_count": parse_nonnegative_int(
                    row["small_variant_count"], f"{row['sample_id']} small_variant_count"
                ),
                "sv_count": parse_nonnegative_int(row["sv_count"], f"{row['sample_id']} sv_count"),
                "tr_loci_different_from_reference": parse_nonnegative_int(
                    row["tr_loci_different_from_reference"],
                    f"{row['sample_id']} tr_loci_different_from_reference",
                ),
            }
        )
    return output


def load_celltype() -> list[dict[str, object]]:
    rows = read_tsv(CELLTYPE_INPUT)
    observed_keys: set[tuple[str, str, str]] = set()
    output: list[dict[str, object]] = []
    for row in rows:
        key = (row["sample_id"], row["variant_class"], row["cell_type_assignment"])
        if key in observed_keys:
            raise ValueError(f"Duplicate cell-type count row: {key}")
        observed_keys.add(key)
        output.append(
            {
                "sample_id": row["sample_id"],
                "variant_class": row["variant_class"],
                "cell_type_assignment": row["cell_type_assignment"],
                "count": parse_nonnegative_int(row["count"], "cell-type count"),
            }
        )

    expected_keys = {
        (sample, variant_class, assignment)
        for sample in EXPECTED_SAMPLES
        for variant_class in EXPECTED_CLASSES
        for assignment in EXPECTED_ASSIGNMENTS
    }
    if observed_keys != expected_keys:
        missing = sorted(expected_keys - observed_keys)
        extra = sorted(observed_keys - expected_keys)
        raise ValueError(f"Unexpected cell-type table keys; missing={missing}, extra={extra}")
    return output


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_celltype_table(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output = sorted(
        rows,
        key=lambda row: (
            str(row["sample_id"]),
            str(row["variant_class"]),
            str(row["cell_type_assignment"]),
        ),
    )
    combined: defaultdict[tuple[str, str], int] = defaultdict(int)
    for row in rows:
        combined[(str(row["variant_class"]), str(row["cell_type_assignment"]))] += int(row["count"])
    for variant_class in ("TR", "SV"):
        for assignment in ("Lymphocyte group", "Monocyte"):
            output.append(
                {
                    "sample_id": "combined",
                    "variant_class": variant_class,
                    "cell_type_assignment": assignment,
                    "count": combined[(variant_class, assignment)],
                }
            )
    return output


def verified_metric_counts(
    germline: list[dict[str, object]], celltype: list[dict[str, object]]
) -> dict[tuple[str, str, str], int]:
    verified: dict[tuple[str, str, str], int] = {}
    for row in germline:
        sample = str(row["sample_id"])
        verified[("bulk_germline", sample, "small_variant")] = int(row["small_variant_count"])
        verified[("bulk_germline", sample, "sv")] = int(row["sv_count"])
        verified[("bulk_germline", sample, "tr")] = int(row["tr_loci_different_from_reference"])

    combined: defaultdict[tuple[str, str], int] = defaultdict(int)
    for row in celltype:
        combined[(str(row["variant_class"]), str(row["cell_type_assignment"]))] += int(row["count"])
    verified[("celltype_associated", "combined", "tr")] = sum(
        combined[("TR", assignment)] for assignment in EXPECTED_ASSIGNMENTS
    )
    verified[("celltype_associated", "combined", "lymphocyte_tr")] = combined[
        ("TR", "Lymphocyte group")
    ]
    verified[("celltype_associated", "combined", "monocyte_tr")] = combined[("TR", "Monocyte")]
    return verified


def build_reconciliation(
    germline: list[dict[str, object]], celltype: list[dict[str, object]]
) -> list[dict[str, object]]:
    verified = verified_metric_counts(germline, celltype)
    output: list[dict[str, object]] = []
    for row in read_tsv(MANUSCRIPT_INPUT):
        key = (row["scope"], row["sample_id"], row["metric"])
        if key not in verified:
            raise ValueError(f"No verified count for manuscript metric {key}")
        reported = parse_nonnegative_int(row["reported_count"], f"reported count {key}")
        current = verified[key]
        output.append(
            {
                "scope": row["scope"],
                "sample_id": row["sample_id"],
                "metric": row["metric"],
                "manuscript_reported_count": reported,
                "validated_source_count": current,
                "difference": current - reported,
                "status": "match" if current == reported else "discrepancy",
            }
        )
    return output


def plot_pbmc_component(celltype: list[dict[str, object]]) -> None:
    counts = {assignment: 0 for assignment in EXPECTED_ASSIGNMENTS}
    for row in celltype:
        counts[str(row["cell_type_assignment"])] += int(row["count"])
    total = sum(counts.values())
    if total == 0:
        raise ValueError("Cannot plot a zero-total PBMC distribution")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.titlesize": 16,
            "axes.labelsize": 13,
        }
    )
    fig, ax = plt.subplots(figsize=(6.3, 6.6))
    bottom = 0.0
    for assignment in ("Lymphocyte group", "Monocyte"):
        count = counts[assignment]
        percent = 100 * count / total
        ax.bar(
            [0],
            [percent],
            bottom=[bottom],
            width=0.5,
            color=COLORS[assignment],
            edgecolor="white",
            linewidth=1.5,
        )
        ax.text(
            0,
            bottom + percent / 2,
            f"{percent:.0f}%",
            ha="center",
            va="center",
            color="white",
            weight="bold",
            fontsize=13,
        )
        label = "Lymphocytes (T+NK+B)" if assignment == "Lymphocyte group" else assignment
        ax.text(
            0.29,
            bottom + percent / 2,
            f"{label} (n={count})",
            ha="left",
            va="center",
            color=COLORS[assignment],
            weight="bold",
            fontsize=11,
        )
        bottom += percent

    ax.set_xticks([0], [f"PBMC\n(n={total})"])
    ax.set_ylim(0, 100)
    ax.set_xlim(-0.65, 1.55)
    ax.set_ylabel("% of cell-type-associated variant calls")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E8E8E8", linewidth=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(
            RESULTS_DIR / f"figure_2f_pbmc_component.{suffix}",
            dpi=350,
            bbox_inches="tight",
            facecolor="white",
        )
    plt.close(fig)


def main() -> None:
    germline = load_germline()
    celltype = load_celltype()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    write_tsv(
        RESULTS_DIR / "supplementary_table_1_pbmc_germline_counts.tsv",
        germline,
        ["sample_id", "tissue", "small_variant_count", "sv_count", "tr_loci_different_from_reference"],
    )
    write_tsv(
        RESULTS_DIR / "supplementary_table_2_pbmc_celltype_counts.tsv",
        build_celltype_table(celltype),
        ["sample_id", "variant_class", "cell_type_assignment", "count"],
    )
    write_tsv(
        RESULTS_DIR / "manuscript_count_reconciliation.tsv",
        build_reconciliation(germline, celltype),
        [
            "scope",
            "sample_id",
            "metric",
            "manuscript_reported_count",
            "validated_source_count",
            "difference",
            "status",
        ],
    )
    plot_pbmc_component(celltype)

    total_calls = sum(int(row["count"]) for row in celltype)
    total_tr = sum(int(row["count"]) for row in celltype if row["variant_class"] == "TR")
    total_sv = sum(int(row["count"]) for row in celltype if row["variant_class"] == "SV")
    print(f"PBMC cell-type-associated calls: total={total_calls}; TR={total_tr}; SV={total_sv}")
    print(f"Wrote tables and Figure 2f PBMC component under {RESULTS_DIR}")


if __name__ == "__main__":
    main()
