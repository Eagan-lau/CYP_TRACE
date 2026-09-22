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
if [ -e blast_main_01 ]; then
  echo 'Refusing to overwrite blast_main_01' >&2
  exit 1
fi
mkdir blast_main_01
blastp -version > blast_main_01/blast_version.txt
sha256sum dataset_02/core_sequences.fasta > blast_main_01/input_fasta.sha256
makeblastdb -in dataset_02/core_sequences.fasta -dbtype prot -out blast_main_01/main_reference
# Keep aligned sequences and spans for later coordinate mapping; do not infer SRS boundaries here.
blastp -query dataset_02/core_sequences.fasta -db blast_main_01/main_reference \
  -out blast_main_01/all_vs_all.tsv -num_threads 4 -evalue 100 \
  -max_target_seqs 100000 -max_hsps 1 -seg yes -comp_based_stats 2 \
  -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qseq sseq'
sha256sum blast_main_01/all_vs_all.tsv > blast_main_01/output_alignment.sha256
touch blast_main_01/SEARCH_COMPLETE_NOT_MODEL_EVALUATION
