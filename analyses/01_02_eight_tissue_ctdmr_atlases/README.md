# 1.2 eight-tissue ctDMR atlases

This analysis reproduces Supplementary Figure 2 of the SniffCell paper and
provides the exact `sniffcell find` workflow used to generate ctDMR atlases for
eight tissue groups.

## Data organization

The large inputs and generated atlases are stored in the Zenodo record
[10.5281/zenodo.22003085](https://doi.org/10.5281/zenodo.22003085). The download
script organizes the record's flat files into two local folders:

```text
data/
├── sniffcell_find_inputs/
│   ├── all_celltypes_blocks.npy
│   ├── all_celltypes_blocks.index.gz
│   ├── all_celltypes.txt
│   └── index_to_major_celltypes.json
└── sniffcell_find_atlases/
    ├── brain_cereb/
    ├── pbmc/
    ├── breast/
    ├── colon/
    ├── kidney/
    ├── liver/
    ├── lung/
    └── pancreas/
```

The repository contains no private server paths. It also contains no variant
calls or variant-level call tables.

## Download the Zenodo data

After the Zenodo draft has been published, run from this analysis directory:

```bash
python scripts/download_zenodo_data.py
```

Every file is checked against `data/zenodo_manifest.tsv` using its byte count
and SHA-256 digest.

## Regenerate the eight atlases

SniffCell v0.9.4 was used with the following parameters:

- `diff_threshold=0.40`
- `min_rows=1`
- `min_cpgs=3`
- `max_gap_bp=2000`

Run:

```bash
bash scripts/run_sniffcell_find.sh
```

The script writes one final atlas to each tissue subfolder under
`data/sniffcell_find_atlases/`. The historical SniffCell source revision used
for these results was `68ec30efc6d54a7989be253b6527b2d7275c2cb2`.

## Reproduce Supplementary Figure 2

Once the atlases have been downloaded or regenerated, run:

```bash
python scripts/plot_supplementary_figure_2.py
```

The command writes:

- `results/supplementary_figure_2.png`
- `results/supplementary_figure_2.pdf`

Panels A and B are calculated directly from the eight public atlas TSVs.
Panels C and D use the path-free aggregate tables committed under
`data/figure_ready/`. Those tables contain region counts and overlap summaries,
not variant calls.

## Provenance

The NPY matrix and companion files are the reprocessed SniffCell `find` inputs.
The eight atlas TSVs are the manuscript outputs generated from that matrix. The
figure script is a path-independent cleanup of the script used to render the
manuscript figure, with the manuscript's panel order, colors, labels, dimensions,
and plotting parameters retained.
