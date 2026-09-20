#!/usr/bin/env python3
"""Table 5 motifs (or independent catalog matches) -> unique-locus GAA tests.

Applies a uniform motif definition to both cell types and calculates
Fisher exact tests from per-locus counts.
"""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import sys

SHARED = Path(__file__).resolve().parents[2] / "03_01_brain_variant_counts_fans_validation/scripts"
sys.path.insert(0, str(SHARED))
from paper_source import read_candidates, records_sha256
from rebuild_candidate_counts import fingerprint, write_table

CELLS = ("Neuron", "Oligodendrocyte")


def is_gaa_element(motif):
    motif = motif.upper()
    if not motif or set(motif) - set("ACGT"):
        raise ValueError("Catalog motif must contain only A/C/G/T; resolve multiple motifs explicitly")
    for sequence in (motif, motif.translate(str.maketrans("ACGT", "TGCA"))[::-1]):
        wrapped = (sequence * 3)[:len(sequence) + 2]
        if any(wrapped[i:i + 3] == "GAA" for i in range(len(sequence))):
            return True
    return False


def fisher(a, b, c, d):
    n1, n2, successes = a + b, c + d, a + c
    if not n1 or not n2:
        raise ValueError("Both cell types need eligible loci")
    def probability(x):
        return math.comb(n1, x) * math.comb(n2, successes - x) / math.comb(n1 + n2, successes)
    observed = probability(a)
    two = sum(probability(x) for x in range(max(0, successes - n2), min(successes, n1) + 1)
              if probability(x) <= observed + 1e-15)
    greater = sum(probability(x) for x in range(a, min(successes, n1) + 1))
    odds = a * d / (b * c) if b * c else (math.inf if a * d else math.nan)
    return odds, min(1., greater), min(1., two)


def match_catalog(rows, path):
    # Stream the full catalog once; retain only overlaps with eligible candidates.
    by_chrom, hits = defaultdict(list), defaultdict(list)
    for r in rows:
        by_chrom[r["chrom"]].append(r)
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end, info = line.rstrip().split("\t")[:4]
            if chrom not in by_chrom:
                continue
            start, end = int(start), int(end)
            for row in by_chrom[chrom]:
                overlap = min(end, int(row["end"])) - max(start, int(row["start"]))
                if overlap > 0:
                    attrs = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
                    hits[row["variant_id"]].append((overlap, start, end, attrs["MOTIFS"].upper()))
    out = []
    for r in rows:
        matches = sorted(hits[r["variant_id"]], reverse=True)
        if not matches or (len(matches) > 1 and matches[0][0] == matches[1][0]):
            raise ValueError(f"Missing or tied catalog match for source row {r['source_row']}")
        overlap, start, end, motif = matches[0]
        out.append(dict(variant_id=r["variant_id"], source_row=r["source_row"], donor=r["donor"],
            cell_type=r["cell_type"], chrom=r["chrom"], call_start=r["start"], call_end=r["end"],
            catalog_start=start, catalog_end=end, catalog_motif=motif, overlap_bp=overlap,
            is_gaa_element=is_gaa_element(motif), motif_source="Adotto catalog"))
    return out


def analyze(candidates, catalog, outdir):
    all_rows = read_candidates(candidates)
    rows = [r for r in all_rows if r["variant_type"] == "TR" and r["tissue"] == "brain"
            and r["observed_in_paired_FANS"] == "Yes"]
    if catalog is not None:
        annotations = match_catalog(rows, catalog)
    else:
        annotations = []
        for r in rows:
            motif = r.get("motif", "").upper()
            qualifies = is_gaa_element(motif)
            annotations.append(dict(variant_id=r["variant_id"], source_row=r["source_row"], donor=r["donor"],
                cell_type=r["cell_type"], chrom=r["chrom"], call_start=r["start"], call_end=r["end"],
                catalog_start=int(r["start"]), catalog_end=int(r["end"]), catalog_motif=motif,
                overlap_bp=int(r["end"]) - int(r["start"]), is_gaa_element=qualifies,
                motif_source="Table 5 supplied motif and interval"))
    primary_unit = "unique_catalog_locus" if catalog is not None else "unique_table_locus"
    results = []
    for unit in (primary_unit, "donor_call"):
        counts = {}
        for cell in CELLS:
            group = [r for r in annotations if r["cell_type"] == cell]
            if unit == primary_unit:
                group = list({(r["chrom"], r["catalog_start"], r["catalog_end"]): r for r in group}.values())
            n = sum(r["is_gaa_element"] for r in group)
            counts[cell] = (n, len(group) - n)
        a, b = counts[CELLS[0]]
        c, d = counts[CELLS[1]]
        odds, greater, two = fisher(a, b, c, d)
        results.append(dict(analysis_unit=unit, neuron_gaa=a, neuron_non_gaa=b,
            oligodendrocyte_gaa=c, oligodendrocyte_non_gaa=d,
            neuron_gaa_percent=100 * a / (a + b), oligodendrocyte_gaa_percent=100 * c / (c + d),
            odds_ratio=odds, fisher_exact_p_greater=greater, fisher_exact_p_two_sided=two))
    outdir = Path(outdir)
    private = outdir / "private_inputs"
    private.mkdir(parents=True, exist_ok=True)
    write_table(private / "gaa_per_locus_annotations.tsv", list(annotations[0]), annotations)
    write_table(outdir / "fans_validated_brain_tr_gaa_tests.tsv", list(results[0]), results)
    provenance = dict(candidate_sha256=fingerprint(candidates), record_sha256=records_sha256(all_rows),
        catalog_sha256=fingerprint(catalog) if catalog is not None else None,
        input_rows=len(all_rows), eligible_donor_rows=len(rows),
        selection="Brain TR; manual paired FANS Yes; same definition for both cell types",
        motif_match="Largest interval overlap with Adotto v2; ties or absent matches fail" if catalog is not None
                    else "Motif and genomic interval provided in Supplementary Table 5",
        gaa_definition="Any-length motif containing cyclic GAA on either strand",
        primary_unit=primary_unit + " within each cell type",
        test="Two-sided Fisher exact; donor-row and one-sided results are sensitivity outputs",
        analysis_output="GAA-element enrichment test statistics and per-locus motif annotations")
    (outdir / "gaa_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidates", type=Path, required=True)
    p.add_argument("--catalog", type=Path, help="Optional independent Adotto matching; default uses Table 5 motif/interval")
    p.add_argument("--outdir", type=Path, required=True)
    a = p.parse_args()
    analyze(a.candidates, a.catalog, a.outdir)
