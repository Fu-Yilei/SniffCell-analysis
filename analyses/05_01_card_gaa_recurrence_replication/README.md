# CARD GAA-family recurrence and cross-cohort replication

This analysis reproduces Figure 5 from aggregate recurrence and replication
tables. Panel A compares Neuron and Oligodendrocyte recurrence across GAA
motif families, panel B shows directional cross-cohort replication, and panel
C shows six identical-coordinate loci observed in all three cohorts.

The public three-cohort table retains locus coordinates, motifs, genes, and
aggregate donor counts only. Donor/sample lists were removed because they are
not required to reproduce the figure.

Run:

```bash
python scripts/plot_figure_5.py
```

`results/figure_5_source.png` and `.pdf` are the exact retained source figure.
No age-at-death metadata or age-association inputs are included.
