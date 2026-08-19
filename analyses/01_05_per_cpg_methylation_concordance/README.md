# 1.5 per-CpG methylation concordance

This analysis reproduces Supplementary Figure 4, comparing methylation beta
values at shared CpGs with coverage greater than 5 in both the SniffCell-
assigned bulk reads and the matched FANS reference.

The 12 benchmark comparisons pool four brain samples (Neuron and
Oligodendrocyte) and two PBMC samples (Monocyte and Lymphocyte). For each PBMC
Lymphocyte comparison, the T-cell, NK-cell, and B-cell FANS pileups are summed
at each CpG before comparison.

The committed `data/*.npz` files are compact, path-free histogram/statistic
caches used to render the exact manuscript panels. To recompute them from the
28 shared Zenodo pileups, first download the common data under
`../shared_data/modkit_cpg_pileups_01_04_to_01_06/`, then run:

```bash
python scripts/plot_supplementary_figure_4.py --rebuild
```

Without `--rebuild`, the script renders directly from the committed caches.
Outputs are written to `results/supplementary_figure_4_brain.*` and
`results/supplementary_figure_4_pbmc.*`.
