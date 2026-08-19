# Multi-region brain regional prevalence

This analysis covers Figure 4b-c and the donor-stratified regional permutation
tests. The public inputs are aggregate locus counts, donor-normalized regional
prevalence, and the reported one-million-permutation results (seed 20260729).
No locus-level calls are included.

Run:

```bash
python scripts/plot_regional_prevalence.py
```

The sampled donor counts are 12 cerebellum, 9 frontal lobe, 8 temporal lobe,
and 8 hippocampus. Global heterogeneity was not significant (one-sided
permutation P=0.565 for neuron and P=0.355 for oligodendrocyte).
`results/figure_4bc_source.png` is the exact retained source panel.
