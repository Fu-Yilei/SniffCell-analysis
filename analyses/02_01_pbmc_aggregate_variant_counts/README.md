# 2.1 PBMC aggregate variant counts

This analysis generates the PBMC count-level inputs for Supplementary Tables 1
and 2 and the PBMC component of Figure 2f. It summarizes:

- bulk germline small-variant, structural-variant, and tandem-repeat totals for
  `pbmc1` and `pbmc2`; and
- full-depth cell-type-associated TR and SV counts by donor and assigned cell
  group.

Only aggregate counts are included. The public files contain no variant loci,
coordinates, alleles, VCF records, read identifiers, or server paths.

Run:

```bash
python scripts/summarize_pbmc_variant_counts.py
```

The script writes:

- `results/supplementary_table_1_pbmc_germline_counts.tsv`;
- `results/supplementary_table_2_pbmc_celltype_counts.tsv`;
- `results/figure_2f_pbmc_component.png` and `.pdf`; and
- `results/manuscript_count_reconciliation.tsv`.

## Count reconciliation

The committed count inputs reflect the current validated source summaries. Two
sentences in the manuscript draft use older or mismatched totals:

- the PBMC TR totals are 616,479 (`pbmc1`) and 616,405 (`pbmc2`), whereas the
  draft lists 613,802 and 616,454; and
- the current full-depth PBMC candidate table contains 10 TRs (8 lymphocyte
  group and 2 monocyte) plus 3 lymphocyte-group SVs. Thus the Figure 2f PBMC
  total is 13 calls (11 lymphocyte group and 2 monocyte), rather than 12 TRs
  split 10 and 2.

The script recomputes these comparisons in
`results/manuscript_count_reconciliation.tsv` so the discrepancy remains
visible and testable.
