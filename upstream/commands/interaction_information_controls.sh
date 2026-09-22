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
if [ "$CYPTRACE_TASK_INDEX" -eq 0 ]; then CHANNEL=sequence; else CHANNEL=structure; fi
"$PYTHON" interaction_information_controls.py --channel "$CHANNEL"
"$PYTHON" verify_interaction_information_controls.py --channel "$CHANNEL"
