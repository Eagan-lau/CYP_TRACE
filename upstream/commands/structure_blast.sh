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
cd $TASK
mkdir structure_blast_01
makeblastdb -in structure_mapping_01/resolved_chain_sequences.fasta -dbtype prot -out structure_blast_01/reference
blastp -query dataset_02/core_sequences.fasta -db structure_blast_01/reference \
  -out structure_blast_01/core_to_resolved_chains.tsv -num_threads 4 -evalue 100 \
  -max_target_seqs 100000 -max_hsps 1 -seg yes -comp_based_stats 2 \
  -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qseq sseq'
sha256sum dataset_02/core_sequences.fasta structure_mapping_01/resolved_chain_sequences.fasta > structure_blast_01/input.sha256
sha256sum structure_blast_01/core_to_resolved_chains.tsv > structure_blast_01/output.sha256
blastp -version > structure_blast_01/version.txt
