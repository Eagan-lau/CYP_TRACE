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
# Search only. Grouping and models require a separate inspected analysis stage.
"$PYTHON" -c 'import json,hashlib,pathlib; p=pathlib.Path("dataset_02"); a=json.loads((p/"dataset_audit.json").read_text()); assert a["status"]=="PASS" and all(a["checks"].values()); c=json.loads((p/"output_checksums.json").read_text()); assert all(hashlib.sha256((p/f).read_bytes()).hexdigest()==h for f,h in c.items())'
if [ -e similarity_02 ]; then
  echo 'Refusing to overwrite similarity_02; inspect existing status before choosing a new output.' >&2
  exit 1
fi
mkdir similarity_02
sha256sum dataset_02/core_sequences.fasta > similarity_02/input_fasta.sha256
"$MMSEQS" version > similarity_02/mmseqs_version.txt
"$MMSEQS" easy-search dataset_02/core_sequences.fasta dataset_02/core_sequences.fasta \
  similarity_02/all_vs_all.tsv similarity_02/tmp \
  --threads 4 -s 7.5 -e 100 --max-seqs 100000 --alignment-mode 3 \
  --format-output query,target,fident,qcov,tcov,alnlen,evalue,bits
sha256sum similarity_02/all_vs_all.tsv > similarity_02/output_alignment.sha256
touch similarity_02/SEARCH_COMMAND_COMPLETED_NOT_GROUPING_OR_MODELS
