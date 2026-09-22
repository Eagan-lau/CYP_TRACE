#!/usr/bin/env bash
# Run with the absolute path of a prepared analysis directory.
# Install third-party binaries/environments first; numerical arguments are retained.
set -euo pipefail
TASK="${1:?Usage: bash this_script.sh /absolute/path/to/workspace/analysis}"
cd "$TASK"
TASK="$PWD"
PROJECT="$(dirname "$TASK")"
PYTHON="${CYPTRACE_PYTHON:-python}"
MMSEQS="${CYPTRACE_MMSEQS:-mmseqs}"
FOLDSEEK="${CYPTRACE_FOLDSEEK:-foldseek}"
FS="$FOLDSEEK"
mkdir -p logs

set -euo pipefail
test -s reverse_alignment_01/mmseqs_raw.tsv
test -s reverse_alignment_01/blast_raw.tsv
test ! -e reverse_alignment_01/audit.json
test ! -e reverse_alignment_validation_01
"$PYTHON" normalize_reverse_alignments.py
"$PYTHON" verify_reverse_alignments.py
