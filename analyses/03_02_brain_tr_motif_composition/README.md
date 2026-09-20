# 3.2 GAA-element enrichment from per-locus Table 5 records

This section analyzes GAA-element enrichment in FANS-validated brain TRs.

The input is the same [Supplementary Table 5](https://docs.google.com/spreadsheets/d/1o31UEyJkIofEj7v_beWsKueVgT5R-KwZPNxyhfV1ubw/edit#gid=1888934510)
snapshot as section 3.1. The calculation starts with individual candidate
intervals, donor IDs, cell-type assignments, and manual FANS annotations,
then classifies their supplied motifs. Independent matching to the Adotto v2
catalog is also supported.

## Run

```bash
python scripts/analyze_gaa_enrichment.py \
  --candidates ../03_01_brain_variant_counts_fans_validation/private_inputs/supplementary_table5.csv \
  --outdir results/table5
```

Or run both sections with section 3.1's `scripts/run_analysis.py`.
No precomputed motif category counts or extra reference download is required.
Add `--catalog private_inputs/adotto.v2.trgt.bed` to independently match catalog
motifs rather than use the Table 5 annotations.

## Analysis steps

1. Validate the candidate table using the same parser as section 3.1.
2. Select brain TR rows annotated Yes in paired FANS. The source snapshot has
   92 eligible donor records: 44 neuron and 48 oligodendrocyte records.
3. Use Table 5's motif and interval. In optional catalog-matching mode,
   match to the Adotto v2 interval with greatest overlap. Missing or tied best
   matches fail.
4. Define a GAA element uniformly in both cell types: a catalog motif of any
   length containing GAA as a cyclic three-base substring on either strand.
   GAA, AAG, AGA and reverse complements qualify; long compound motifs can also
   qualify. Ambiguous or comma-separated motifs require explicit
   resolution and are rejected.
5. Deduplicate chromosome/start/end within each cell type for
   the primary analysis. A recurrent locus counts once per cell type, not once
   per donor. Do not sum these denominators to estimate the union of brain loci.
6. Compare GAA versus non-GAA counts between cell types with two-sided Fisher's
   exact test. Also export the one-sided test and donor-row sensitivity analysis.

## Results

| Unit | Neuron GAA / total | Oligodendrocyte GAA / total | Two-sided Fisher P |
| --- | --- | --- | --- |
| Unique Table 5 locus | 18/41 (43.9%) | 13/43 (30.2%) | 0.258800690627 |
| Donor record | 19/44 (43.2%) | 13/48 (27.1%) | 0.127803765100 |

Results are calculated from the Table 5 records using the motif rule above.
The analysis exports test statistics and per-locus motif annotations.

## Outputs and publication boundaries

- [Test results](results/table5/fans_validated_brain_tr_gaa_tests.tsv)
- [Source hashes and definitions](results/table5/gaa_provenance.json)
- Local-only `results/table5/private_inputs/gaa_per_locus_annotations.tsv`: source row,
  candidate ID, donor, cell type, call/catalog intervals, catalog motif, overlap,
  and GAA classification. This is git-ignored and must not be uploaded as calls.

The default reproduction uses Table 5 alone. When a catalog is supplied, its
SHA256 is recorded in the source metadata. No private server path, callset,
or age metadata is included in public outputs.

The supplementary `derive_tr_motif_composition.py`,
`plot_tr_motif_composition.py`, `data/`, and results outside
`results/table5/` contain archived motif-length and GC-composition analyses.
