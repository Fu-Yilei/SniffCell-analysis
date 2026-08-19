#!/usr/bin/env python3
"""Validate Zenodo scDNA coverage summaries and reproduce Extended Data Fig. 2a-b."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


SAMPLE_ORDER = ["bcontrol1", "bcontrol2", "bparkinsons", "bmsa"]
CELL_TYPE_ORDER = ["Oligo", "Neuron"]
DISPLAY_CELL_TYPES = ["Oligodendrocyte", "Neuron"]
EXPECTED_COUNTS = {
    ("bcontrol1", "Oligo"): 14,
    ("bcontrol1", "Neuron"): 90,
    ("bcontrol2", "Oligo"): 64,
    ("bcontrol2", "Neuron"): 93,
    ("bparkinsons", "Oligo"): 128,
    ("bparkinsons", "Neuron"): 93,
    ("bmsa", "Oligo"): 46,
}
COLORS = {"Oligodendrocyte": "#1B9E77", "Neuron": "#E67E22"}
SUMMARY_SUFFIX = ".fig13_all_122.mosdepth.summary.txt"
REQUIRED_COLUMNS = {
    "sample",
    "bam_id",
    "scdna_cell_type",
    "unique_variant_loci",
    "inside_variant_bp",
    "inside_variant_mean_depth",
    "genome_coverage_bases",
    "mosdepth_reported_genome_mean",
    "outside_variant_mean_depth",
}


def parse_args() -> argparse.Namespace:
    section_dir = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="Extracted scdna_coverage_summaries archive root.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=section_dir / "results",
        help="Output directory (default: this section's results directory).",
    )
    parser.add_argument(
        "--skip-mosdepth-check",
        action="store_true",
        help="Render from the per-run TSV without checking all raw mosdepth summaries.",
    )
    return parser.parse_args()


def load_per_run(data_dir: Path) -> pd.DataFrame:
    path = data_dir / "data" / "scdna_coverage_per_run.tsv"
    if not path.is_file():
        raise FileNotFoundError(path)
    data = pd.read_csv(path, sep="\t")
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    forbidden = [column for column in data.columns if "path" in column.lower()]
    if forbidden:
        raise ValueError(f"Path-bearing columns are not permitted: {forbidden}")
    if len(data) != sum(EXPECTED_COUNTS.values()):
        raise ValueError(f"Expected 528 runs, found {len(data)}")
    if data.duplicated(["sample", "bam_id"]).any():
        raise ValueError("Duplicate sample/run identifiers found")

    observed = data.groupby(["sample", "scdna_cell_type"]).size().to_dict()
    if observed != EXPECTED_COUNTS:
        raise ValueError(f"Unexpected sample/cell-type run counts: {observed}")
    if not (data["unique_variant_loci"] == 105).all():
        raise ValueError("Expected 105 deduplicated variant intervals per run")
    if not (data["inside_variant_bp"] == 71_154).all():
        raise ValueError("Expected a 71,154-bp inside-interval denominator per run")

    numeric = data.select_dtypes(include=[np.number])
    if (numeric < 0).any().any() or not np.isfinite(numeric.to_numpy()).all():
        raise ValueError("Coverage metrics must be finite and non-negative")
    data["sample"] = pd.Categorical(data["sample"], SAMPLE_ORDER, ordered=True)
    data["scdna_cell_type"] = pd.Categorical(
        data["scdna_cell_type"], CELL_TYPE_ORDER, ordered=True
    )
    return data.sort_values(["sample", "scdna_cell_type", "bam_id"]).reset_index(
        drop=True
    )


def validate_mosdepth_summaries(data: pd.DataFrame, data_dir: Path) -> None:
    seen: set[Path] = set()
    for row in data.itertuples(index=False):
        path = (
            data_dir
            / "mosdepth"
            / str(row.sample)
            / str(row.scdna_cell_type)
            / f"{row.bam_id}{SUMMARY_SUFFIX}"
        )
        if not path.is_file():
            raise FileNotFoundError(path)
        summary = pd.read_csv(path, sep="\t")
        total = summary.loc[summary["chrom"] == "total"]
        if len(total) != 1:
            raise ValueError(f"{path}: expected exactly one total row")
        total_row = total.iloc[0]
        if int(total_row["bases"]) != int(row.genome_coverage_bases):
            raise ValueError(f"{path}: genome coverage-bases mismatch")
        if not np.isclose(
            float(total_row["mean"]),
            float(row.mosdepth_reported_genome_mean),
            atol=0.001,
        ):
            raise ValueError(f"{path}: genome mean-depth mismatch")
        seen.add(path)
    archived = set((data_dir / "mosdepth").rglob(f"*{SUMMARY_SUFFIX}"))
    if archived != seen:
        raise ValueError(
            f"Expected exactly {len(seen)} selected mosdepth files; found {len(archived)}"
        )


def make_long(data: pd.DataFrame) -> pd.DataFrame:
    shared = ["sample", "bam_id", "scdna_cell_type"]
    outside = data[shared].copy()
    outside["coverage_context"] = "Outside SniffCell variant"
    outside["average_depth"] = data["outside_variant_mean_depth"]
    inside = data[shared].copy()
    inside["coverage_context"] = "Inside SniffCell variant"
    inside["average_depth"] = data["inside_variant_mean_depth"]
    return pd.concat([outside, inside], ignore_index=True)


def build_totals(data: pd.DataFrame) -> pd.DataFrame:
    return (
        data.groupby(["sample", "scdna_cell_type"], observed=True, as_index=False)
        .agg(
            n_scdna_runs=("bam_id", "size"),
            cumulative_outside_variant_depth=("outside_variant_mean_depth", "sum"),
            cumulative_inside_variant_depth=("inside_variant_mean_depth", "sum"),
            mean_outside_variant_depth=("outside_variant_mean_depth", "mean"),
            median_outside_variant_depth=("outside_variant_mean_depth", "median"),
            mean_inside_variant_depth=("inside_variant_mean_depth", "mean"),
            median_inside_variant_depth=("inside_variant_mean_depth", "median"),
        )
        .sort_values(["sample", "scdna_cell_type"])
    )


def sample_labels() -> list[str]:
    sample_counts = {
        sample: sum(count for (name, _), count in EXPECTED_COUNTS.items() if name == sample)
        for sample in SAMPLE_ORDER
    }
    return [f"{sample}\nn={sample_counts[sample]}" for sample in SAMPLE_ORDER]


def render_figure(long_data: pd.DataFrame, totals: pd.DataFrame, output_dir: Path) -> None:
    sns.set_theme(style="whitegrid", context="talk", font_scale=1.05)
    plot_data = long_data.copy()
    plot_data["Cell type"] = plot_data["scdna_cell_type"].astype(str).map(
        {"Oligo": "Oligodendrocyte", "Neuron": "Neuron"}
    )

    total_frames: list[pd.DataFrame] = []
    for context, metric in {
        "Outside SniffCell variant": "cumulative_outside_variant_depth",
        "Inside SniffCell variant": "cumulative_inside_variant_depth",
    }.items():
        frame = totals[["sample", "scdna_cell_type", metric]].rename(
            columns={metric: "total_depth"}
        )
        frame["coverage_context"] = context
        total_frames.append(frame)
    total_data = pd.concat(total_frames, ignore_index=True)
    total_data["Cell type"] = total_data["scdna_cell_type"].astype(str).map(
        {"Oligo": "Oligodendrocyte", "Neuron": "Neuron"}
    )

    fig, axes = plt.subplots(2, 2, figsize=(17, 12))
    contexts = ["Outside SniffCell variant", "Inside SniffCell variant"]
    for column, context in enumerate(contexts):
        average_ax = axes[0, column]
        sns.boxplot(
            data=plot_data.loc[plot_data["coverage_context"] == context],
            x="sample",
            y="average_depth",
            hue="Cell type",
            order=SAMPLE_ORDER,
            hue_order=DISPLAY_CELL_TYPES,
            palette=COLORS,
            saturation=1,
            fliersize=2.5,
            ax=average_ax,
        )
        average_ax.set_title(context)
        average_ax.set_xlabel("")
        average_ax.set_ylabel("Average coverage per scDNA run (×)" if column == 0 else "")
        average_ax.set_xticks(range(len(SAMPLE_ORDER)), sample_labels())
        average_ax.set_ylim(bottom=0)
        average_ax.legend_.remove()

        total_ax = axes[1, column]
        sns.barplot(
            data=total_data.loc[total_data["coverage_context"] == context],
            x="sample",
            y="total_depth",
            hue="Cell type",
            order=SAMPLE_ORDER,
            hue_order=DISPLAY_CELL_TYPES,
            palette=COLORS,
            saturation=1,
            errorbar=None,
            ax=total_ax,
        )
        total_ax.set_title(context)
        total_ax.set_xlabel("Sample")
        total_ax.set_ylabel("Total coverage (×)" if column == 0 else "")
        total_ax.set_xticks(range(len(SAMPLE_ORDER)), sample_labels())
        total_ax.set_ylim(bottom=0)
        total_ax.legend_.remove()

    axes[0, 0].text(
        -0.16, 1.10, "a", transform=axes[0, 0].transAxes,
        fontsize=28, fontweight="bold", va="top"
    )
    axes[1, 0].text(
        -0.16, 1.10, "b", transform=axes[1, 0].transAxes,
        fontsize=28, fontweight="bold", va="top"
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        title="Cell type",
        loc="upper center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=2,
        frameon=False,
    )
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.94), h_pad=2.3, w_pad=2.0)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / "extended_data_figure_2ab_reproduced"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    data = load_per_run(args.data_dir)
    if not args.skip_mosdepth_check:
        validate_mosdepth_summaries(data, args.data_dir)
    totals = build_totals(data)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    totals.to_csv(
        args.output_dir / "scdna_coverage_by_sample_cell_type.tsv",
        sep="\t",
        index=False,
    )
    render_figure(make_long(data), totals, args.output_dir)
    print(f"validated_runs={len(data)}")
    print(f"validated_mosdepth_files={0 if args.skip_mosdepth_check else len(data)}")
    print(totals.to_string(index=False))


if __name__ == "__main__":
    main()
