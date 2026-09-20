#!/usr/bin/env python3
"""Derive Table-5-based Figure 2e/f/g from detailed candidate records.

Does not call variants or infer manual FANS validation from summary counts.
Candidate details remain in private_inputs; only aggregate outputs are public.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from paper_source import read_candidates, depth_counts, summarize, records_sha256
from rebuild_candidate_counts import write_table, fingerprint
from plot_paper_figure2 import render

HERE = Path(__file__).resolve().parent


def rebuild(candidates, outdir, catalog=None, motif_outdir=None):
    rows = read_candidates(candidates)
    depth = depth_counts(rows)
    validation, composition, unique = summarize(rows)
    outdir = Path(outdir)
    private = outdir / "private_inputs"
    private.mkdir(parents=True, exist_ok=True)
    # Full record-to-panel trace is deliberately not a public output.
    write_table(private / "normalized_candidates.tsv", list(rows[0]), rows)
    ledger = [dict(variant_id=r["variant_id"], source_row=r["source_row"], donor=r["donor"],
                   tissue=r["tissue"], variant_class=r["variant_type"], cell_type=r["cell_type"],
                   coverage=coverage, fans_observed=r["observed_in_paired_FANS"])
              for r in rows for coverage in r["coverages_detected"].split(";")]
    write_table(private / "panel_membership.tsv", list(ledger[0]), ledger)
    for name, records in (("figure_2e_counts.tsv", depth),
                          ("figure_2f_validation.tsv", validation),
                          ("figure_2g_composition.tsv", composition),
                          ("brain_unique_locus_validation.tsv", [unique])):
        write_table(outdir / name, list(records[0]), records)
    classes = []
    for tissue, cells in (("pbmc", ("T-cell+NK-cell+B-cell", "Monocyte")),
                          ("brain", ("Neuron", "Oligodendrocyte"))):
        for kind in ("TR", "SV"):
            for cell in cells:
                classes.append(dict(tissue=tissue, variant_class=kind, cell_type=cell,
                    count=sum(r["tissue"] == tissue and r["variant_type"] == kind
                              and r["cell_type"] == cell for r in rows)))
    write_table(outdir / "candidate_classes_by_celltype.tsv", list(classes[0]), classes)
    provenance = dict(source="Supplementary Table 5 SniffCell", sheet_gid=1888934510,
        source_url="https://docs.google.com/spreadsheets/d/1o31UEyJkIofEj7v_beWsKueVgT5R-KwZPNxyhfV1ubw/edit#gid=1888934510",
        input_sha256=fingerprint(candidates), record_sha256=records_sha256(rows), candidate_rows=len(rows),
        figure_2e="Listed candidate membership at each coverage; not independent discovery totals",
        figure_2f="Candidate/donor rows; Yes and No are manual matched-FANS outcomes",
        figure_2g="All listed TR and SV candidate rows by tissue/cell type",
        brain_validation="Unique genomic intervals across donors; observed if any donor is Yes")
    (outdir / "paper_panel_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    render(outdir, outdir)
    if motif_outdir is not None:
        motif_script = HERE.parents[1] / "03_02_brain_tr_motif_composition/scripts/analyze_gaa_enrichment.py"
        command = [sys.executable, str(motif_script), "--candidates", str(candidates),
                   "--outdir", str(motif_outdir)]
        if catalog is not None:
            command.extend(["--catalog", str(catalog)])
        subprocess.run(command, check=True)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidates", type=Path, required=True, help="Detailed Table 5 TSV/CSV or connector JSON")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--catalog", type=Path)
    p.add_argument("--motif-outdir", type=Path)
    a = p.parse_args()
    if a.catalog is not None and a.motif_outdir is None:
        p.error("--catalog requires --motif-outdir")
    rebuild(a.candidates, a.outdir, a.catalog, a.motif_outdir)


if __name__ == "__main__":
    main()
