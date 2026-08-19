# 3.1 Brain aggregate counts and FANS validation

This section rebuilds the brain component of Figure 2f and a count-level FANS
validation panel from aggregate tables only. No locus-level callset or VCF is
included.

The exact manuscript-source count panel is retained as
`results/figure_2f_cell_type_counts_source.png`. Run:

```bash
python scripts/plot_brain_counts_and_validation.py
```

The current Figure 13-matched source contains 56 neuron and 53
oligodendrocyte candidates (109 total). The current FANS validation table
contains 73 TR rows, of which 47 are observed in paired FANS. These do not
match the draft prose's class-separated 53/52 TR and 9/2 SV counts or its
44/65 FANS-locus statement. The script preserves the current source totals and
writes a machine-readable reconciliation table instead of hiding the mismatch.
