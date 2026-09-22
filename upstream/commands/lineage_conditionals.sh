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
# Execute once for every task index from 0 through 1.
CYPTRACE_TASK_INDEX="${2:?Supply the task index as the second argument (0-1)}"
[[ "$CYPTRACE_TASK_INDEX" =~ ^[0-9]+$ ]] || exit 2
(( CYPTRACE_TASK_INDEX >= 0 && CYPTRACE_TASK_INDEX <= 1 )) || exit 2
set -euo pipefail
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
CHANNELS=(sequence structure)
"$PYTHON" run_lineage_conditionals.py --channel "${CHANNELS[$CYPTRACE_TASK_INDEX]}"
