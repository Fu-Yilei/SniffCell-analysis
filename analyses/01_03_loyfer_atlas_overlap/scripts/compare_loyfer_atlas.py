#!/usr/bin/env python3
"""Compare the public Loyfer ctDMR atlas with eight SniffCell tissue atlases."""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


TISSUE_ORDER = [
    "Brain",
    "PBMC",
    "Breast",
    "Colon",
    "Kidney",
    "Liver",
    "Lung",
    "Pancreas",
]
TISSUE_COLORS = {
    "Brain": "#1f77b4",
    "PBMC": "#ff7f0e",
    "Breast": "#e377c2",
    "Colon": "#8c564b",
    "Kidney": "#17becf",
    "Liver": "#bcbd22",
    "Lung": "#9467bd",
    "Pancreas": "#2ca02c",
}


def parse_args() -> argparse.Namespace:
    analysis_dir = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--loyfer-atlas",
        type=Path,
        default=analysis_dir / "data/loyfer_atlas/Atlas.U250.l4.hg38.full.tsv",
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        default=analysis_dir / "data/lofyer_celltype_mapping.tsv",
    )
    parser.add_argument(
        "--sniffcell-atlas-dir",
        type=Path,
        default=(
            analysis_dir.parent
            / "01_02_eight_tissue_ctdmr_atlases/data/sniffcell_find_atlases"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=analysis_dir / "results",
    )
    parser.add_argument("--bedtools", default="bedtools")
    return parser.parse_args()


def contains_leaf(value: object, group: str) -> bool:
    return group in str(value).split("|")


def atlas_path(root: Path, key: str) -> Path:
    filename = f"{key}.combined_dt04.ctdmr.tsv"
    nested = root / key / filename
    flat = root / filename
    if nested.is_file():
        return nested
    if flat.is_file():
        return flat
    raise FileNotFoundError(f"missing SniffCell atlas: {nested}")


def validate_inputs(loyfer: pd.DataFrame, mapping: pd.DataFrame) -> None:
    required = {"chr", "start", "end", "target", "direction"}
    missing = sorted(required - set(loyfer.columns))
    if missing:
        raise ValueError(f"Loyfer atlas is missing columns: {missing}")
    if len(loyfer) != 9520:
        raise ValueError(f"expected 9,520 Loyfer markers, found {len(loyfer):,}")
    directions = set(loyfer["direction"].dropna().astype(str))
    if directions != {"U"}:
        raise ValueError(f"expected only U-direction Loyfer markers, found {directions}")
    if len(mapping) != 27:
        raise ValueError(f"expected 27 cell-type mappings, found {len(mapping)}")
    missing_targets = sorted(set(mapping["atlas_target"]) - set(loyfer["target"]))
    if missing_targets:
        raise ValueError(f"mapping targets absent from Loyfer atlas: {missing_targets}")


def run_bedtools(
    bedtools: str,
    loyfer_bed: Path,
    catalog_no_header: Path,
) -> str:
    completed = subprocess.run(
        [
            bedtools,
            "intersect",
            "-a",
            str(loyfer_bed),
            "-b",
            str(catalog_no_header),
            "-wa",
            "-wb",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    return completed.stdout


def compare(
    loyfer: pd.DataFrame,
    mapping: pd.DataFrame,
    atlas_dir: Path,
    bedtools: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    loyfer = loyfer.copy()
    # The published full atlas uses a 1-based start, whereas BEDTools expects
    # zero-based, half-open BED coordinates. This is the conversion used by the
    # manuscript analysis and by the companion hg38 BED in the upstream repo.
    loyfer["bed_start"] = pd.to_numeric(loyfer["start"], errors="raise") - 1
    loyfer["loyfer_id"] = (
        loyfer["chr"].astype(str)
        + ":"
        + loyfer["bed_start"].astype(str)
        + "-"
        + loyfer["end"].astype(str)
        + "|"
        + loyfer["target"].astype(str)
    )

    coverage_rows: list[dict[str, object]] = []
    direction_rows: list[dict[str, object]] = []
    detail_rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="sniffcell_loyfer_") as temporary:
        temporary_dir = Path(temporary)
        catalog_cache: dict[str, tuple[list[str], Path]] = {}

        for rec in mapping.itertuples(index=False):
            tissue = str(rec.tissue)
            key = str(rec.json_key)
            group = str(rec.sniffcell_group)
            target = str(rec.atlas_target)

            if tissue not in catalog_cache:
                catalog = pd.read_csv(atlas_path(atlas_dir, key), sep="\t")
                required = {
                    "chr",
                    "start",
                    "end",
                    "name",
                    "best_group",
                    "best_dir",
                    "hyper_group_leaves",
                    "hypo_group_leaves",
                }
                missing = sorted(required - set(catalog.columns))
                if missing:
                    raise ValueError(f"{tissue} atlas is missing columns: {missing}")
                no_header = temporary_dir / f"{key}.no_header.tsv"
                catalog.to_csv(no_header, sep="\t", header=False, index=False)
                catalog_cache[tissue] = (list(catalog.columns), no_header)
            catalog_columns, catalog_no_header = catalog_cache[tissue]

            subset = loyfer[loyfer["target"] == target].copy()
            loyfer_bed = temporary_dir / f"{key}.{group}.loyfer.bed"
            subset[
                ["chr", "bed_start", "end", "target", "direction", "loyfer_id"]
            ].to_csv(loyfer_bed, sep="\t", header=False, index=False)
            stdout = run_bedtools(bedtools, loyfer_bed, catalog_no_header)

            interval_status: dict[str, set[str]] = {}
            overlap_count = 0
            for line in stdout.splitlines():
                overlap_count += 1
                parts = line.split("\t")
                loyfer_id = parts[5]
                sniffcell = dict(zip(catalog_columns, parts[6:]))
                if contains_leaf(sniffcell.get("hypo_group_leaves"), group):
                    status = "concordant"
                elif contains_leaf(sniffcell.get("hyper_group_leaves"), group):
                    status = "discordant"
                else:
                    status = "group_not_in_sniffcell_contrast"
                interval_status.setdefault(loyfer_id, set()).add(status)
                detail_rows.append(
                    {
                        "tissue": tissue,
                        "sniffcell_group": group,
                        "atlas_target": target,
                        "loyfer_direction": "U",
                        "expected_sniffcell_group_direction": "hypo",
                        "loyfer_id": loyfer_id,
                        "sniffcell_region": (
                            f"{sniffcell['chr']}:{sniffcell['start']}-{sniffcell['end']}"
                        ),
                        "sniffcell_name": sniffcell["name"],
                        "sniffcell_best_group": sniffcell["best_group"],
                        "sniffcell_best_dir": sniffcell["best_dir"],
                        "sniffcell_hyper_group_leaves": sniffcell[
                            "hyper_group_leaves"
                        ],
                        "sniffcell_hypo_group_leaves": sniffcell[
                            "hypo_group_leaves"
                        ],
                        "direction_status": status,
                    }
                )

            n_total = len(subset)
            n_overlapped = len(interval_status)
            n_concordant = sum(
                "concordant" in statuses for statuses in interval_status.values()
            )
            n_discordant = sum(
                "discordant" in statuses and "concordant" not in statuses
                for statuses in interval_status.values()
            )
            n_both = sum(
                "concordant" in statuses and "discordant" in statuses
                for statuses in interval_status.values()
            )
            n_other = sum(
                statuses == {"group_not_in_sniffcell_contrast"}
                for statuses in interval_status.values()
            )
            coverage_rows.append(
                {
                    "tissue": tissue,
                    "json_key": key,
                    "sniffcell_group": group,
                    "atlas_target": target,
                    "n_loyfer_regions": n_total,
                    "n_overlapped_by_sniffcell": n_overlapped,
                    "coverage_fraction": n_overlapped / n_total,
                    "coverage_percent": 100 * n_overlapped / n_total,
                }
            )
            direction_rows.append(
                {
                    "tissue": tissue,
                    "sniffcell_group": group,
                    "atlas_target": target,
                    "loyfer_direction": "U",
                    "expected_sniffcell_group_direction": "hypo",
                    "n_loyfer_regions": n_total,
                    "n_loyfer_regions_with_any_sniffcell_overlap": n_overlapped,
                    "n_concordant_intervals": n_concordant,
                    "n_discordant_only_intervals": n_discordant,
                    "n_intervals_with_both_concordant_and_discordant_overlaps": n_both,
                    "n_group_not_in_sniffcell_contrast_intervals": n_other,
                    "n_pairwise_sniffcell_overlaps": overlap_count,
                    "concordant_fraction_of_overlapped": (
                        n_concordant / n_overlapped if n_overlapped else 0.0
                    ),
                    "discordant_only_fraction_of_overlapped": (
                        n_discordant / n_overlapped if n_overlapped else 0.0
                    ),
                }
            )

    coverage = pd.DataFrame(coverage_rows)
    direction = pd.DataFrame(direction_rows)
    detail = pd.DataFrame(detail_rows)
    by_tissue = (
        coverage.drop_duplicates(["tissue", "atlas_target"])
        .groupby("tissue", as_index=False, sort=False)
        .agg(
            n_loyfer_regions=("n_loyfer_regions", "sum"),
            n_overlapped_by_sniffcell=("n_overlapped_by_sniffcell", "sum"),
        )
    )
    by_tissue["coverage_fraction"] = (
        by_tissue["n_overlapped_by_sniffcell"] / by_tissue["n_loyfer_regions"]
    )
    by_tissue["coverage_percent"] = 100 * by_tissue["coverage_fraction"]
    return coverage, by_tissue, direction, detail


def setup_style() -> None:
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8,
        }
    )
    sns.set_theme(
        style="ticks",
        context="paper",
        rc={
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#444444",
            "grid.color": "#e0e0e0",
            "grid.linewidth": 0.45,
        },
    )


def plot_panel_d(ax, coverage: pd.DataFrame, direction: pd.DataFrame) -> None:
    panel = coverage.merge(
        direction[
            [
                "tissue",
                "sniffcell_group",
                "n_concordant_intervals",
                "n_loyfer_regions_with_any_sniffcell_overlap",
            ]
        ],
        on=["tissue", "sniffcell_group"],
        validate="one_to_one",
    )
    y = np.arange(len(panel))
    colors = [TISSUE_COLORS[tissue] for tissue in panel["tissue"]]
    ax.hlines(
        y,
        60,
        panel["coverage_percent"],
        color="#d7d7d7",
        linewidth=0.8,
        zorder=1,
    )
    ax.scatter(
        panel["coverage_percent"],
        y,
        s=30,
        c=colors,
        edgecolor="#333333",
        linewidth=0.35,
        zorder=3,
    )
    ax.axvline(100, color="#777777", linewidth=0.6, linestyle=":")
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{row.tissue}: {row.sniffcell_group}" for row in panel.itertuples()],
        fontsize=7.2,
    )
    ax.invert_yaxis()
    ax.set_xlim(60, 102)
    ax.set_xticks([60, 70, 80, 90, 100])
    ax.set_xlabel("Loyfer ctDMRs overlapped by SniffCell (%)")
    ax.xaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for yi, percent in zip(y, panel["coverage_percent"]):
        ax.text(
            min(percent + 0.8, 101.0),
            yi,
            f"{percent:.1f}",
            va="center",
            ha="left",
            fontsize=6.5,
        )
    n_concordant = int(panel["n_concordant_intervals"].sum())
    n_overlap = int(panel["n_loyfer_regions_with_any_sniffcell_overlap"].sum())
    if (n_concordant, n_overlap) != (6050, 6050):
        raise ValueError(
            f"expected 6,050/6,050 directional concordance; found "
            f"{n_concordant:,}/{n_overlap:,}"
        )
    ax.text(
        0.0,
        1.03,
        f"Direction concordance: {n_concordant}/{n_overlap} overlapped intervals "
        "(U atlas -> SniffCell hypo)",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
    )
    ax.text(
        -0.06,
        1.10,
        "D",
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="top",
        ha="left",
    )


def main() -> int:
    args = parse_args()
    loyfer = pd.read_csv(args.loyfer_atlas, sep="\t")
    mapping = pd.read_csv(args.mapping, sep="\t")
    validate_inputs(loyfer, mapping)
    coverage, by_tissue, direction, detail = compare(
        loyfer,
        mapping,
        args.sniffcell_atlas_dir,
        args.bedtools,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(
        args.output_dir / "lofyer_overlap_coverage_by_celltype.tsv",
        sep="\t",
        index=False,
    )
    by_tissue.to_csv(
        args.output_dir / "lofyer_overlap_coverage_by_tissue.tsv",
        sep="\t",
        index=False,
    )
    direction.to_csv(
        args.output_dir / "lofyer_directional_concordance_by_celltype.tsv",
        sep="\t",
        index=False,
    )
    detail.to_csv(
        args.output_dir / "lofyer_directional_concordance_overlaps.tsv",
        sep="\t",
        index=False,
    )

    setup_style()
    figure, axis = plt.subplots(figsize=(6.8, 5.1))
    plot_panel_d(axis, coverage, direction)
    figure.subplots_adjust(left=0.34, right=0.98, top=0.91, bottom=0.11)
    output = args.output_dir / "supplementary_figure_2d.png"
    figure.savefig(output, dpi=300, facecolor="white")
    figure.savefig(output.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)

    print(f"wrote {len(coverage)} mapped cell-type summaries")
    print("directional concordance: 6050/6050")
    print(f"saved {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
