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
mkdir foldseek_raw_01
"$PYTHON" verify_structure_source.py
"$FOLDSEEK" version > foldseek_raw_01/version.txt
"$FOLDSEEK" createdb "$PROJECT/data/raw/rcsb_cyp_ligand_templates_v1/mmcif" foldseek_raw_01/structures --threads 8
"$FOLDSEEK" convert2fasta foldseek_raw_01/structures foldseek_raw_01/structure_sequences.fasta
"$FOLDSEEK" search foldseek_raw_01/structures foldseek_raw_01/structures \
  foldseek_raw_01/alignments foldseek_raw_01/tmp --threads 8 -s 9.5 -e 100 \
  --max-seqs 10000 --alignment-type 2 --exhaustive-search 1 -a 1
"$FOLDSEEK" convertalis foldseek_raw_01/structures foldseek_raw_01/structures \
  foldseek_raw_01/alignments foldseek_raw_01/all_vs_all.tsv --threads 8 \
  --format-output 'query,target,fident,qcov,tcov,alnlen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,alntmscore,lddt'
sha256sum foldseek_raw_01/all_vs_all.tsv foldseek_raw_01/structure_sequences.fasta > foldseek_raw_01/output.sha256
