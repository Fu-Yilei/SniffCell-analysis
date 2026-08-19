# 1.6 read-level disagreement at ctDMRs

This analysis reproduces Supplementary Figure 5: the percentage of classified
reads whose methylation differs from the matched FANS reference by more than
the ctDMR assignment threshold (`|Delta beta| > 0.4`).

The brain and monocyte reference betas were calculated from the shared
modkit-derived FANS CpG pileups in
`../shared_data/modkit_cpg_pileups_01_04_to_01_06/`. PBMC lymphocyte reads were
evaluated against their assigned FANS leaf (T-cell, NK-cell, or B-cell); the
multi/unresolved bucket used the mean of the applicable FANS leaves.

The read-level calculation requires the SniffCell split BAMs and is therefore
not a public-data-only step. The path-free aggregate result needed for the
paper figure is committed as `data/read_level_fans_disagreement.tsv`. It
contains counts only, with no read identifiers, variant calls, or server paths.

The original calculation used seed 7, sampled at most 2,000 ctDMRs per
brain/Monocyte comparison and 500 ctDMRs per PBMC leaf-resolved comparison,
and counted disagreement strictly greater than 0.4. `n_obs` is the number of
assessed read-by-ctDMR observations; `total_reads` is the full classified-read
denominator displayed above each bar.

Run:

```bash
python scripts/plot_supplementary_figure_5.py
```

Outputs are written to `results/supplementary_figure_5.png` and `.pdf`.
