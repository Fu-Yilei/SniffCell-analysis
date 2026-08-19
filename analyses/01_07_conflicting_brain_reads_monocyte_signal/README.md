# 1.7 Monocyte-like signal among conflicting brain reads

This analysis reproduces Supplementary Figure 6. It asks whether brain reads
with conflicting Neuron-like and Oligodendrocyte-like ctDMR assignments carry
at least one ctDMR whose methylation is closer to the Monocyte reference.

The committed table contains aggregate counts only. It excludes read names,
server paths, BAMs, and variant calls. Because the public figure can be
reproduced from these aggregates, rerunning SniffCell is not required.
The six bulk and sixteen purified FACS/FANS original per-read assignment
tables are available through
`../shared_data/read_assignments_01_07_to_01_08/`.

Run:

```bash
python scripts/plot_supplementary_figure_6.py
```

Outputs are written to `results/supplementary_figure_6.png` and `.pdf`.
