# 3.3 ONT low-quality SNV-like artifact

Supplementary Figure 12 is the manuscript's IGV view of the FBXL7
`chr5:15569886 C>T` SNV-like signal in `bcontrol1`. The signal is enriched in
SniffCell- and FANS-neuron reads, but the uniformly low base quality makes it
unreliable and illustrates why no high-confidence cell-type-specific SNV was
reported.

The exact image embedded in the manuscript is retained here. This is a manual
review result, so the executable component verifies the asset checksum rather
than attempting to regenerate an IGV session without its BAMs. No VCF, call
table, BAM, read identifier, or server path is included.

Run:

```bash
python scripts/verify_figure_asset.py
```
