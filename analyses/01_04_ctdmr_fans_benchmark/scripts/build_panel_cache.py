#!/usr/bin/env python3
"""Build Supplementary Figure 3 panel arrays from public atlases and pileups."""

from __future__ import annotations

import argparse
import shlex
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


HG38_SIZES = [
    ("chr1", 248956422), ("chr2", 242193529), ("chr3", 198295559),
    ("chr4", 190214555), ("chr5", 181538259), ("chr6", 170805979),
    ("chr7", 159345973), ("chr8", 145138636), ("chr9", 138394717),
    ("chr10", 133797422), ("chr11", 135086622), ("chr12", 133275309),
    ("chr13", 114364328), ("chr14", 107043718), ("chr15", 101991189),
    ("chr16", 90338345), ("chr17", 83257441), ("chr18", 80373285),
    ("chr19", 58617616), ("chr20", 64444167), ("chr21", 46709983),
    ("chr22", 50818468), ("chrX", 156040895), ("chrY", 57227415),
    ("chrM", 16569),
]
CHROMS = [chrom for chrom, _ in HG38_SIZES]
CHROM_ORDER = {chrom: index for index, chrom in enumerate(CHROMS)}
BRAIN_SAMPLES = ["b1", "b4", "b2", "b3"]
PBMC_SAMPLES = ["pbmc1", "pbmc2"]
BRAIN_LABELS = ["Neuron", "Oligodendrocyte"]
PBMC_LABELS = ["T-cell", "NK-cell", "Monocyte", "B-cell"]


def parse_args() -> argparse.Namespace:
    analysis_dir = Path(__file__).resolve().parents[1]
    analyses_dir = analysis_dir.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pileup-root",
        type=Path,
        default=analyses_dir / "shared_data/modkit_cpg_pileups_01_04_to_01_06/data/sniffcell_modkit_cpg_pileups",
    )
    parser.add_argument(
        "--brain-atlas",
        type=Path,
        default=analysis_dir / "data/brain_cereb.maxgap3000_dt04.ctdmr.tsv",
        help="Exact brain benchmark atlas mirrored in GitHub and Zenodo.",
    )
    parser.add_argument(
        "--pbmc-atlas",
        type=Path,
        default=analyses_dir / "01_02_eight_tissue_ctdmr_atlases/data/sniffcell_find_atlases/pbmc/pbmc.combined_dt04.ctdmr.tsv",
    )
    parser.add_argument("--output", type=Path, default=analysis_dir / "data/supplementary_figure_3_panel_data.npz")
    return parser.parse_args()


def chromosome_key(value: object) -> tuple[int, str]:
    chrom = str(value)
    return CHROM_ORDER.get(chrom, len(CHROM_ORDER)), chrom


def write_regions(catalog: pd.DataFrame, path: Path) -> None:
    missing = {"chr", "start", "end"} - set(catalog.columns)
    if missing:
        raise ValueError(f"Atlas is missing columns: {sorted(missing)}")
    regions = catalog[["chr", "start", "end"]].copy()
    regions["region_index"] = np.arange(len(regions), dtype=np.int64)
    regions["_order"] = regions["chr"].map(lambda value: chromosome_key(value)[0])
    regions.sort_values(["_order", "start", "end"], kind="stable").drop(columns="_order").to_csv(
        path, sep="\t", header=False, index=False
    )


def aggregate(pileup: Path, regions: Path, genome: Path, expected_rows: int, temporary: Path) -> tuple[np.ndarray, np.ndarray]:
    mapped = temporary / f"mapped_{abs(hash(str(pileup)))}.tsv"
    command = (
        f"gzip -dc {shlex.quote(str(pileup))} | "
        "awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$2+1,$3,$4}' | "
        f"bedtools map -g {shlex.quote(str(genome))} -a {shlex.quote(str(regions))} "
        f"-b - -c 4,5 -o sum,sum -null 0 > {shlex.quote(str(mapped))}"
    )
    subprocess.run(["bash", "-c", command], check=True)
    frame = pd.read_csv(mapped, sep="\t", header=None, names=["chr", "start", "end", "region_index", "coverage", "modified"])
    frame.sort_values("region_index", kind="stable", inplace=True)
    if len(frame) != expected_rows:
        raise RuntimeError(f"{pileup} mapped to {len(frame)} regions; expected {expected_rows}")
    coverage = pd.to_numeric(frame["coverage"], errors="coerce").fillna(0).to_numpy(dtype=float)
    modified = pd.to_numeric(frame["modified"], errors="coerce").fillna(0).to_numpy(dtype=float)
    beta = np.divide(modified, coverage, out=np.full(expected_rows, np.nan), where=coverage > 0)
    return coverage, beta


def leaves(value: object) -> list[str]:
    return [] if pd.isna(value) else [label for label in str(value).split("|") if label]


def grouped_effect(catalog: pd.DataFrame, values: dict[str, np.ndarray]) -> np.ndarray:
    output = np.full(len(catalog), np.nan)
    best = catalog["best_group_leaves"].fillna("").astype(str)
    other = catalog["other_group_leaves"].fillna("").astype(str)
    for pair in pd.DataFrame({"best": best, "other": other}).drop_duplicates().itertuples(index=False):
        best_labels = [label for label in leaves(pair.best) if label in values]
        other_labels = [label for label in leaves(pair.other) if label in values]
        if not best_labels or not other_labels:
            continue
        mask = ((best == pair.best) & (other == pair.other)).to_numpy()
        output[mask] = np.nanmean(np.column_stack([values[label][mask] for label in best_labels]), axis=1) - np.nanmean(
            np.column_stack([values[label][mask] for label in other_labels]), axis=1
        )
    return output


def clean_pair(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    keep = np.isfinite(x) & np.isfinite(y)
    return x[keep], y[keep]


def main() -> int:
    args = parse_args()
    for required in (args.brain_atlas, args.pbmc_atlas, args.pileup_root):
        if not required.exists():
            raise FileNotFoundError(required)
    brain = pd.read_csv(args.brain_atlas, sep="\t")
    pbmc = pd.read_csv(args.pbmc_atlas, sep="\t")
    output: dict[str, np.ndarray] = {}

    with tempfile.TemporaryDirectory(prefix="sniffcell_ctdmr_benchmark_") as name:
        temporary = Path(name)
        genome = temporary / "hg38.genome"
        genome.write_text("".join(f"{chrom}\t{size}\n" for chrom, size in HG38_SIZES), encoding="utf-8")
        brain_regions, pbmc_regions = temporary / "brain.regions.bed", temporary / "pbmc.regions.bed"
        write_regions(brain, brain_regions)
        write_regions(pbmc, pbmc_regions)

        brain_catalog_values = {label: pd.to_numeric(brain[f"mean_{label}"], errors="coerce").to_numpy() for label in BRAIN_LABELS}
        brain_x = brain_catalog_values["Neuron"] - brain_catalog_values["Oligodendrocyte"]
        for sample in BRAIN_SAMPLES:
            sample_values = {}
            for label in BRAIN_LABELS:
                _, sample_values[label] = aggregate(
                    args.pileup_root / f"brain/{sample}/fans/{label}.cpg_simple.tsv.gz",
                    brain_regions,
                    genome,
                    len(brain),
                    temporary,
                )
            x, y = clean_pair(brain_x, sample_values["Neuron"] - sample_values["Oligodendrocyte"])
            output[f"{sample}_x"], output[f"{sample}_y"], output[f"{sample}_n_total"] = x, y, np.asarray(len(brain))

        pbmc_catalog_values = {label: pd.to_numeric(pbmc[f"mean_{label}"], errors="coerce").to_numpy() for label in PBMC_LABELS}
        pbmc_x = grouped_effect(pbmc, pbmc_catalog_values)
        for sample in PBMC_SAMPLES:
            sample_values = {}
            for label in PBMC_LABELS:
                _, sample_values[label] = aggregate(
                    args.pileup_root / f"pbmc/{sample}/fans/{label}.cpg_simple.tsv.gz",
                    pbmc_regions,
                    genome,
                    len(pbmc),
                    temporary,
                )
            x, y = clean_pair(pbmc_x, grouped_effect(pbmc, sample_values))
            output[f"{sample}_x"], output[f"{sample}_y"], output[f"{sample}_n_total"] = x, y, np.asarray(len(pbmc))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **output)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
