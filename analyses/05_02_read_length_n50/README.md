# Read-length N50 cohort comparison

This analysis reproduces Supplementary Figure 17 from the per-sample N50
source table. It compares four FANS-validation brain samples, 44 multi-region
brain samples, and 155 large-scale CARD/HBCC brain samples.

The public table contains anonymous within-cohort row identifiers and
sequencing summary metrics. Original sample identifiers, server paths, and
alignment filenames were removed. The script recomputes the cohort summaries,
one-sided Mann-Whitney tests, and the figure.

Run:

```bash
python scripts/plot_supplementary_figure_17.py
```

`results/supplementary_figure_17_source.png` and `.pdf` are the exact retained
source figure. No age-at-death metadata are included.
