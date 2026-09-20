# 3.2 TR motif composition from candidate-level input

`scripts/derive_tr_motif_composition.py` implements the original
`build_figure13_tr_motif_composition.py` analysis: select TR rows, classify each
motif, and calculate counts and percentages. The checked-in category counts
are an output snapshot, not the only supported starting point.

```bash
python scripts/derive_tr_motif_composition.py \
  --candidates private_inputs/full_depth_annotated_candidates.tsv \
  --outdir results/rebuilt
python scripts/plot_tr_motif_composition.py \
  --input results/rebuilt/tr_motif_composition.tsv \
  --outdir results/rebuilt
```

The section 3.1 `scripts/run_analysis.py` runs both sections. Install
`requirements.txt` for plotting; derivation uses the Python standard library.

## Definitions

- Input TSV columns: `variant_type, tissue, donor, chrom, start, end, motif`.
  Tissue is `brain` or `pbmc`; variant type is `TR` or `SV`.
- Only TRs are counted. The unit is a **donor-locus candidate**, not a unique
  genomic locus. Repeated loci across donors count once per donor. Duplicate
  rows within a donor/locus fail explicitly.
- Motif length is annotated string length, not minimal repeat period or a
  motif inferred from expanded read sequence. Bins: 2, 3, 4, 5, 6–10, >10 bp.
- GC fraction is `(G + C) / length`. AT-rich: GC <= 0.20; GC-rich: GC >= 0.80;
  Mixed: between these bounds. Case is normalized. Empty, ambiguous, and
  one-base motifs fail because silently dropping them changes denominators.
- Counts and denominators are calculated from rows, not hardcoded to 10/105.
  Percentages are within tissue and dimension, rounded to one decimal place.
  Zero categories are retained.

## Verified archived input

`Figure13_all_full_depth_variant_detailed_summary.tsv` contains 122 candidates,
including 10 PBMC and 105 brain TR rows. Derivation from these source rows
exactly reproduces all 18 rows of `data/tr_motif_composition.tsv`. The source
SHA256 and input/TR row counts are recorded in the generated provenance table.

This full-depth set differs from the older FANS table, which includes subset
candidates. They must not be treated as identical cohorts. Motif annotations
are an explicit input; this step does not infer motifs or rerun discovery.

Only category statistics and a checksum are exported. Source calls, coordinates,
and motif-by-locus records remain private. The current public deposition alone
supports the bundled plot, but not private candidate-level derivation.

The archived manuscript PNG remains `results/tr_motif_composition_source.png`.
To plot the bundled snapshot alone, run `python scripts/plot_tr_motif_composition.py`.
