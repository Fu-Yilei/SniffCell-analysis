# 1.3 Loyfer atlas overlap

This analysis reproduces the Loyfer-atlas comparison shown in Supplementary
Figure 2d. It measures how often the eight SniffCell tissue atlases overlap the
corresponding cell-type markers from the normal human DNA methylation atlas and
checks whether the methylation direction agrees.

## Inputs

- `data/loyfer_atlas/Atlas.U250.l4.hg38.full.tsv`: the public hg38 Loyfer UXM
  atlas, committed directly to this repository.
- `data/lofyer_celltype_mapping.tsv`: the explicit mapping between SniffCell
  tissue groups and Loyfer atlas targets used in the manuscript analysis.
- The eight SniffCell atlas TSVs from
  `../01_02_eight_tissue_ctdmr_atlases/data/sniffcell_find_atlases/`.

The Loyfer file contains 9,520 ctDMR markers and no variant calls. Its values
represent unmethylated-fragment signal. Every marker in this file has direction
`U`, so agreement requires the mapped SniffCell cell type to occur in the
SniffCell region's hypomethylated group.

## Run

Install the Python requirements and BEDTools, then run from this analysis
directory:

```bash
python scripts/compare_loyfer_atlas.py
```

The command writes path-free overlap and directionality tables plus:

- `results/supplementary_figure_2d.png`
- `results/supplementary_figure_2d.pdf`

Use `--sniffcell-atlas-dir` if the eight SniffCell atlas folders are stored
elsewhere.

## Interpretation

Coverage is the percentage of Loyfer ctDMR intervals for a mapped cell type
that overlap at least one ctDMR in the corresponding SniffCell tissue atlas.
Directional concordance is evaluated only among overlapped Loyfer intervals.
The manuscript result is 66.0%-100.0% overlap across mapped cell types and
6,050/6,050 directionally concordant overlapped intervals.

## Public-atlas provenance and license

The atlas was downloaded from the public
[`nloyfer/UXM_deconv`](https://github.com/nloyfer/UXM_deconv) repository at
revision `8d0bb456742ece21802bc8e651e080869c67c944`. Its SHA-256 digest is
`5cd05666dc1fbf9686b30a022d6178979dd996a91f4d84cf0ffd7c46dba6a18b`.

Use and distribution of this Technology is subject to the Software Research
License included herein. See `LICENSE.LOYFER.md`; the upstream license restricts
use and redistribution to Research Use as defined in that license.

Please cite Loyfer, Magenheim, Peretz et al. (2023), *A DNA methylation atlas of
normal human cell types*, Nature.
