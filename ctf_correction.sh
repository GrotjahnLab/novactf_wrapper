#!/bin/bash
set -euo pipefail

base="HRR031_3_lam17D_ts_002.mrc_nova"
stack="${base}.mrc"
tiltFile="${base}.tlt"

# microscope parameters -- fill in your actual values
amplitudeContrast=0.07
cs=2.7
volt=300

n=$(ls ${base}.defocus_* | wc -l)

for ((i=0; i<n; i++)); do
    defocus="${base}.defocus_${i}"
    corrected="${base}_corrected.mrc_${i}"

    novaCTF -Algorithm ctfCorrection \
        -InputProjections "$stack" \
        -OutputFile "$corrected" \
        -DefocusFile "$defocus" \
        -TILTFILE "$tiltFile" \
        -CorrectionType phaseflip \
        -DefocusFileFormat imod \
        -CorrectAstigmatism 1 \
        -PixelSize 0.332 \
        -AmplitudeContrast "$amplitudeContrast" \
        -Cs "$cs" \
        -Volt "$volt"
done
