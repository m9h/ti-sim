#!/bin/bash
# Check a MIDA model download (IT'IS Foundation / US FDA; Iacono et al. 2015) for what
# ti_sim.mida needs: a NIfTI label volume and the label-name text file. MIDA is not
# fetched by this script and must never be committed, packaged or redistributed: its
# licence (LicenseAgreementMIDA, IT'IS) allows use but forbids distributing the data or
# derivatives, requires disguised faces on published images and a citation of the FDA,
# IT'IS and Iacono et al. in publications.
#
#   scripts/check-mida.sh [MIDA_DIR]      (default: $TI_SIM_MIDA or ./data/mida)
#
# To obtain it: https://itis.swiss/virtual-population/regional-human-models/mida-model/
# (DOI 10.13099/ViP-MIDA-V1.0; register with the form on the "MIDA V1.0" page), unpack
# the voxel model (MAT/RAW/NII/TXT) under MIDA_DIR, then `export TI_SIM_MIDA=MIDA_DIR`.
set -euo pipefail
DIR="${1:-${TI_SIM_MIDA:-$PWD/data/mida}}"
if [ ! -d "$DIR" ]; then
    echo "no MIDA directory at $DIR (see the comment in this script for how to obtain it)" >&2
    exit 2
fi
nii=$(find "$DIR" -iname '*.nii' -o -iname '*.nii.gz' | head -1)
txt=$(find "$DIR" -iname '*.txt' | head -1)
[ -n "$nii" ] || { echo "no .nii/.nii.gz label volume under $DIR" >&2; exit 1; }
[ -n "$txt" ] || { echo "no .txt label list under $DIR" >&2; exit 1; }
echo "label volume: $nii"
echo "label names:  $txt"
python3 -I - "$DIR" <<'PY'
import sys
sys.path.insert(0, "src")
from ti_sim import mida
nii, txt = mida.find_files(sys.argv[1])
names = mida.parse_label_file(txt)
mapping = mida.classify(names)
print(f"{len(names)} structures parsed, all mapped to SimNIBS tags; hippocampus label(s): "
      f"{[k for k, v in mapping.items() if v == mida.HIPPOCAMPUS]}")
PY
echo "MIDA ok: export TI_SIM_MIDA=$DIR"
