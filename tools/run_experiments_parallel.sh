#!/usr/bin/env bash
# Helper wrapper that uses GNU parallel to execute BuffCut experiments with per-run
# memory limits. Usage:
#   tools/run_experiments_parallel.sh <config> <binary> <jobs> <memory_mb> [task_dir]

set -euo pipefail

if [ $# -lt 4 ]; then
    echo "Usage: $0 <config> <binary> <jobs> <memory_mb> [task_dir]" >&2
    exit 1
fi

if ! command -v parallel >/dev/null 2>&1; then
    echo "GNU parallel is required but not found in PATH." >&2
    exit 1
fi

CONFIG=$1
BINARY=$2
JOBS=$3
MEM_MB=$4
TASK_DIR=${5:-"$(mktemp -d experiments/tasks_XXXXXX)"}

echo "Planning tasks into ${TASK_DIR}"
python3 tools/run_experiments.py \
    --config "${CONFIG}" \
    --binary "${BINARY}" \
    --plan-output-dir "${TASK_DIR}" \
    --plan-only

TASK_COUNT=$(ls "${TASK_DIR}"/*.json 2>/dev/null | wc -l | tr -d ' ')
if [ "${TASK_COUNT}" -eq 0 ]; then
    echo "No tasks were generated; nothing to run."
    exit 0
fi

echo "Launching ${TASK_COUNT} tasks with GNU parallel (${JOBS} jobs, limit ${MEM_MB} MB)..."
parallel --jobs "${JOBS}" \
    python3 tools/run_experiments.py --run-task {} --memory-limit-mb "${MEM_MB}" \
    ::: "${TASK_DIR}"/*.json

echo "All parallel tasks finished."
