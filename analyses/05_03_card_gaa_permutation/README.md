# CARD GAA-family within-donor permutation tests

This analysis reproduces the age-independent permutation component of
Supplementary Figure 18. The included 20-row result table contains observed
Neuron/Oligodendrocyte recurrence counts and permutation statistics for five
GAA motif families at four recurrence thresholds.

Run:

```bash
python scripts/plot_supplementary_figure_18a.py
```

`results/supplementary_figure_18a_source.png` and `.pdf` are the exact retained
source panel. The source result table is precomputed from within-donor label
swaps; donor-level calls are not included. The age-association panel, its
sample metadata, and all age-analysis scripts are intentionally excluded.
