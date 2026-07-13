#!/bin/bash
set -euo pipefail

base="HRR031_3_lam17D_ts_002.mrc_nova"
tiltFile="${base}.tlt"

n=$(ls ${base}.defocus_* | wc -l)

for ((i=0; i<n; i++)); do
    corrected="${base}_corrected.mrc_${i}"
    flipped="${base}_corrected_flipped.mrc_${i}"
    filtered="${base}_filtered.mrc_${i}"

    clip flipyz "$corrected" "$flipped"

    novaCTF -Algorithm filterProjections \
        -InputProjections "$flipped" \
        -OutputFile "$filtered" \
        -TILTFILE "$tiltFile" \
        -StackOrientation xz
done
