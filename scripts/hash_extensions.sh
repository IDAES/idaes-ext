#!/usr/bin/env bash

# Usage: ./hash_all_tarballs.sh <release> <directory>

set -euo pipefail

RELEASE="$1"
DIR="$2"

OUTFILE="${DIR}/sha256sum_${RELEASE}.txt"

echo "Writing hashes to: ${OUTFILE}"
rm -f "$OUTFILE"
touch "$OUTFILE"

shopt -s nullglob
for f in "$DIR"/*.tar.gz; do
    hash=$(sha256sum "$f" | awk '{print $1}')
    filename=$(basename "$f")
    echo "${hash}  ${filename}" >> "$OUTFILE"
done

echo "Done."

