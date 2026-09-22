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
"$PYTHON" prepare_ordered_heme_anchor.py
"$FS" filterdb foldseek_raw_01/alignments ordered_anchor_01/to_reference_alignments --filter-file ordered_anchor_01/reference_foldseek_key.txt --filter-column 1 --positive-filter 1 --threads 4
"$FS" convertalis foldseek_raw_01/structures foldseek_raw_01/structures ordered_anchor_01/to_reference_alignments ordered_anchor_01/structure_to_reference.tsv --threads 4 --format-output query,target,qstart,qend,tstart,tend,qlen,tlen,qcov,tcov,evalue,bits,qaln,taln
sha256sum ordered_anchor_01/structure_to_reference.tsv > ordered_anchor_01/alignment.sha256
