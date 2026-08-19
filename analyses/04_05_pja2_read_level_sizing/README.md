# PJA2 read-level repeat sizing

This analysis covers Figure 4f at `chr5:109340592-109340704` (motif CCCTCC),
near PJA2. The exact panel is retained in `results/figure_4f_source.{png,pdf}`.

The sanitized point-level input is in the Zenodo archive
`multi_region_brain_tr_figure_data.tar.zst`, under
`multi_region_brain_tr_figure_data/pja2_read_lengths/`. It contains donor,
tissue, platform/source, cell type, and repeat length. Read identifiers and
server paths were removed.

After the Zenodo record is published, download the input and run:

```bash
python scripts/download_zenodo_data.py
python scripts/plot_pja2_read_lengths.py
```

For the unpublished draft only, the downloader also accepts
`--token-file PATH`; the token is read at runtime and is never stored here.

The six displayed donors are SMHT018, SMHT020, SMHT022, SMHT023, SMHT024,
and SMHT031.
