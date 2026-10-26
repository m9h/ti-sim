#!/bin/bash
# Fetch the SimNIBS "ernie" low-resolution example head (charm m2m folder, ~413 MB zip)
# into a data directory, verify its checksum, and unpack it. Idempotent.
#
#   scripts/fetch-ernie.sh [DATA_DIR]      -> DATA_DIR/m2m_ernie  (default: ./data)
#
# Third-party data (CC BY-NC 4.0): downloaded into its own directory and only unpacked
# after the sha256 matches. Never packaged in the free COPR.
set -euo pipefail

DATA="${1:-$PWD/data}"
URL="https://github.com/simnibs/example-dataset/releases/download/v4.0-lowres/ernie_lowres_V2.zip"
SHA256="399f184f5121818ca7c2c7f06dbb2967f69dc3aff6ce076e78d1e834debfa6ab"
ZIP="$DATA/download/ernie_lowres_V2.zip"

if [ -f "$DATA/m2m_ernie/ernie.msh" ]; then
    echo "ernie already present: $DATA/m2m_ernie"
    exit 0
fi
mkdir -p "$DATA/download"
if [ ! -f "$ZIP" ] || ! echo "$SHA256  $ZIP" | sha256sum --quiet -c - 2>/dev/null; then
    echo "downloading $URL"
    curl -fL --retry 3 -o "$ZIP.part" "$URL"
    mv "$ZIP.part" "$ZIP"
fi
echo "$SHA256  $ZIP" | sha256sum -c -
UNPACK="$DATA/unpack.$$"
mkdir -p "$UNPACK"
unzip -q "$ZIP" -d "$UNPACK"
if [ ! -d "$UNPACK/m2m_ernie" ]; then
    echo "unexpected archive layout under $UNPACK" >&2
    exit 1
fi
mv "$UNPACK/m2m_ernie" "$DATA/m2m_ernie"
[ -f "$UNPACK/readme.txt" ] && mv "$UNPACK/readme.txt" "$DATA/ernie-readme.txt"
rmdir "$UNPACK" 2>/dev/null || true
echo "ernie ready: $DATA/m2m_ernie"
