# 1.8 FACS/FANS versus SniffCell ctDMR PCA

This analysis reproduces Supplementary Figure 7 from the committed PCA scores.
The PCA compares ctDMR methylation profiles from sorted FACS/FANS reads with
the corresponding SniffCell-assigned reads. Grey segments join each matched
source pair.

The public data contain PCA coordinates, explained variance, and aggregate
pair distances. They contain no BAM paths, read identifiers, or variant calls.
The expensive per-read methylation extraction and SniffCell run do not need to
be repeated to verify the figure.
The corresponding original per-read assignment tables are available through
`../shared_data/read_assignments_01_07_to_01_08/`.

Run:

```bash
python scripts/plot_supplementary_figure_7.py
```

Outputs are written to `results/supplementary_figure_7.png` and `.pdf`.
