# 2.2 V(D)J rearrangement analysis

This section reproduces Figure 3a-b from the sanitized read-level event tables
in the Zenodo archive `vdj_rearrangement_events.tar.zst` at
doi:10.5281/zenodo.22003085. It contains no BAMs or variant-call files.

After extracting the archive, run:

```bash
python scripts/plot_vdj_rearrangements.py \
  --events-dir /path/to/vdj_rearrangement_events
```

The exact manuscript-source PNG is retained as
`results/figure_3ab_vdj_rearrangements_source.png`; the script rebuilds the
same panel from the public archive. The event tables retain read names and alignment evidence, but the original
server-path column was removed before upload. The plot applies the manuscript
filters (event size >=1,000 bp; MAPQ >=20; read length >=1,000 bp; strict split
read geometry) and writes PNG, PDF, and a plotted-event count table to
`results/`.
