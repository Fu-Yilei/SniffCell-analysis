#!/usr/bin/env bash
set -euo pipefail

analysis_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
input_dir=${SNIFFCELL_FIND_INPUT_DIR:-"${analysis_dir}/data/sniffcell_find_inputs"}
atlas_dir=${SNIFFCELL_FIND_ATLAS_DIR:-"${analysis_dir}/data/sniffcell_find_atlases"}
sniffcell_bin=${SNIFFCELL_BIN:-sniffcell}

npy="${input_dir}/all_celltypes_blocks.npy"
index="${input_dir}/all_celltypes_blocks.index.gz"
meta="${input_dir}/all_celltypes.txt"
celltypes="${input_dir}/index_to_major_celltypes.json"

for required in "$npy" "$index" "$meta" "$celltypes"; do
    if [[ ! -f "$required" ]]; then
        echo "missing required input: $required" >&2
        exit 1
    fi
done

tissues=(brain_cereb pbmc breast colon kidney liver lung pancreas)
for tissue in "${tissues[@]}"; do
    tissue_dir="${atlas_dir}/${tissue}"
    output="${tissue_dir}/${tissue}.combined_dt04.ctdmr.tsv"
    mkdir -p "$tissue_dir"
    "$sniffcell_bin" find \
        -n "$npy" \
        -i "$index" \
        -m "$meta" \
        -cf "$celltypes" \
        -ck "$tissue" \
        -o "$output" \
        --diff_threshold 0.40 \
        --min_rows 1 \
        --min_cpgs 3 \
        --max_gap_bp 2000
done
