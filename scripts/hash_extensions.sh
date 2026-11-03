#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") <release> <directory>

Description:
  Computes SHA-256 hashes for all .tar.gz files in <directory>
  and writes them to:

    <directory>/sha256sum_<release>.txt

Arguments:
  <release>    Release label used in the output filename
  <directory>  Directory containing .tar.gz files

Options:
  -h, --help   Show this help message and exit

Example:
  $(basename "$0") 4.0.0 /path/to/artifacts

  Copy the resulting file into the releases/ directory of this repository.
  idaes-pse downloads it from there to verify the binaries.
EOF
}

if [[ $# -eq 0 ]]; then
    usage
    exit 1
fi

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ $# -ne 2 ]]; then
    echo "Error: expected 2 arguments." >&2
    usage
    exit 1
fi

RELEASE="$1"
DIR="$2"

OUTFILE="${DIR}/sha256sum_${RELEASE}.txt"

# macOS ships shasum rather than sha256sum
if command -v sha256sum >/dev/null 2>&1; then
    SHA256="sha256sum"
else
    SHA256="shasum -a 256"
fi

echo "Writing hashes to: ${OUTFILE}"
rm -f "$OUTFILE"
touch "$OUTFILE"

shopt -s nullglob
for f in "$DIR"/*.tar.gz; do
    hash=$($SHA256 "$f" | awk '{print $1}')
    filename=$(basename "$f")
    echo "${hash}  ${filename}" >> "$OUTFILE"
done

echo "Done."

