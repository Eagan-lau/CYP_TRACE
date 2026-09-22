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
# Execute once for every task index from 0 through 3.
CYPTRACE_TASK_INDEX="${2:?Supply the task index as the second argument (0-3)}"
[[ "$CYPTRACE_TASK_INDEX" =~ ^[0-9]+$ ]] || exit 2
(( CYPTRACE_TASK_INDEX >= 0 && CYPTRACE_TASK_INDEX <= 3 )) || exit 2

set -euo pipefail
case "$CYPTRACE_TASK_INDEX" in
  0) CONTROL=interaction_control_protein_position_permutation_01; CHANNEL=sequence ;;
  1) CONTROL=interaction_control_protein_position_permutation_01; CHANNEL=structure ;;
  2) CONTROL=interaction_control_reaction_center_row_permutation_01; CHANNEL=sequence ;;
  3) CONTROL=interaction_control_reaction_center_row_permutation_01; CHANNEL=structure ;;
esac
cd "$TASK/$CONTROL"
"$PYTHON" run_interaction_conditionals.py --channel "$CHANNEL"
"$PYTHON" verify_interaction_conditionals.py --channel "$CHANNEL"
"$PYTHON" interaction_paired_comparisons.py --channel "$CHANNEL"
"$PYTHON" verify_interaction_paired.py --channel "$CHANNEL"
