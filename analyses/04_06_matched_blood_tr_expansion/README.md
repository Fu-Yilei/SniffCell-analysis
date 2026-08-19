# Matched-blood TR expansion comparison

This analysis covers Figure 4g-h. The sanitized source table is in the shared
Zenodo archive `multi_region_brain_tr_figure_data.tar.zst`, under
`multi_region_brain_tr_figure_data/matched_blood_expansion/`. It contains only
the values used by the figure; VCF paths, genotypes, support fields, and server
paths are excluded.

After the Zenodo record is published, download the input and run:

```bash
python scripts/download_zenodo_data.py
python scripts/plot_matched_blood_expansion.py
```

For the unpublished draft only, the downloader also accepts
`--token-file PATH`; the token is read at runtime and is never stored here.

The current source has 83 donor-locus rows across 35 loci and 12 donors (35
neuron-expanded rows and 48 oligodendrocyte-expanded rows). The 88-row count
in an earlier caption is stale. `results/figure_4gh_source.{png,pdf}` is the
exact retained source figure.
