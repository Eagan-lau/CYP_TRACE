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
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
"$PYTHON" build_dataset.py --raw-run run_03 --output dataset_01
mkdir -p similarity_01
"$MMSEQS" version > similarity_01/mmseqs_version.txt
"$MMSEQS" easy-search dataset_01/core_sequences.fasta dataset_01/core_sequences.fasta \
  similarity_01/all_vs_all.tsv similarity_01/tmp \
  --threads 4 -s 7.5 -e 100 --max-seqs 100000 --alignment-mode 3 \
  --format-output query,target,fident,qcov,tcov,alnlen,evalue,bits
