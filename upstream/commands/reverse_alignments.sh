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
if [[ -e reverse_alignment_01 ]]; then
  echo 'Refusing to overwrite reverse_alignment_01' >&2
  exit 1
fi
mkdir reverse_alignment_01
"$MMSEQS" version > reverse_alignment_01/mmseqs_version.txt
"$MMSEQS" easy-search reverse_candidate_census_01/candidate_sequences.fasta dataset_02/core_sequences.fasta \
  reverse_alignment_01/mmseqs_raw.tsv reverse_alignment_01/mmseqs_tmp \
  --threads 8 -s 7.5 -e 100 --max-seqs 100000 --alignment-mode 3 \
  --format-output query,target,fident,qcov,tcov,alnlen,evalue,bits
blastp -version > reverse_alignment_01/blast_version.txt
makeblastdb -in dataset_02/core_sequences.fasta -dbtype prot -out reverse_alignment_01/core_reference
blastp -query reverse_candidate_census_01/candidate_sequences.fasta -db reverse_alignment_01/core_reference \
  -out reverse_alignment_01/blast_raw.tsv -num_threads 8 -evalue 100 \
  -max_target_seqs 100000 -max_hsps 1 -seg yes -comp_based_stats 2 \
  -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen'
"$PYTHON" normalize_reverse_alignments.py
"$PYTHON" verify_reverse_alignments.py
