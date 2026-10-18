#!/bin/bash
#SBATCH --job-name=cisc_12_filter
#SBATCH --output=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%A_%a.out
#SBATCH --error=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%A_%a.err
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=0-3

set -euo pipefail


# ── PATHS ─────────────────────────────────────────────────────────────────

ROOT=${CISC_ROOT:-/proj/sebastianlab/users/hquintal/cisc}

PYTHON_ENV=${CISC_PYTHON_ENV:-/proj/sebastianlab/users/hquintal/venvs/eca-hi}

PY_SCRIPT=${ROOT}/scripts/12_filter_clustered_rp.py

INPUT_DIR=${ROOT}/data/out/11_clusters/clustered_rp5


# ── SETTINGS ──────────────────────────────────────────────────────────────

WORKER_COUNT=${CISC_FILTER_WORKER_COUNT:-4}


# ── VALIDATION ────────────────────────────────────────────────────────────

[[ -f "${PY_SCRIPT}" ]] || {
    echo "ERROR: Missing Python script: ${PY_SCRIPT}" >&2
    exit 1
}

[[ -d "${INPUT_DIR}" ]] || {
    echo "ERROR: Missing clustered RP5 directory: ${INPUT_DIR}" >&2
    exit 1
}

[[ -f "${PYTHON_ENV}/bin/activate" ]] || {
    echo "ERROR: Missing Python environment: ${PYTHON_ENV}" >&2
    exit 1
}

[[ "${WORKER_COUNT}" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: Invalid CISC_FILTER_WORKER_COUNT: ${WORKER_COUNT}" >&2
    exit 1
}

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not defined." >&2
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= WORKER_COUNT )); then
    echo \
        "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} is outside " \
        "0-$((WORKER_COUNT - 1))." \
        >&2
    exit 1
fi


# ── ENVIRONMENT ───────────────────────────────────────────────────────────

module purge
module load python

# shellcheck disable=SC1090
source "${PYTHON_ENV}/bin/activate"

export CISC_ROOT="${ROOT}"
export CISC_FILTER_WORKER_COUNT="${WORKER_COUNT}"


# ── PRE-FLIGHT ────────────────────────────────────────────────────────────

python -c \
'import pandas, pyarrow; print("Python dependencies verified"); print("  pandas:", pandas.__version__); print("  pyarrow:", pyarrow.__version__)'

input_count=$(
    find "${INPUT_DIR}" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name "cluster_rp5_*.parquet" \
        | wc -l
)

if [[ "${input_count}" -eq 0 ]]; then
    echo "ERROR: No clustered RP5 Parquet files found under ${INPUT_DIR}" >&2
    exit 1
fi


# ── RUN ───────────────────────────────────────────────────────────────────

echo "=================================================="
echo "12_filter_clustered_rp.sh"
echo "=================================================="
echo "ROOT=${ROOT}"
echo "PY_SCRIPT=${PY_SCRIPT}"
echo "INPUT_DIR=${INPUT_DIR}"
echo "INPUT_FILES=${input_count}"
echo "WORKER_COUNT=${WORKER_COUNT}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID:-NA}"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"
echo "HOSTNAME=$(hostname)"
echo "START_TIME=$(date --iso-8601=seconds)"
echo "=================================================="

python "${PY_SCRIPT}"

echo "=================================================="
echo "END_TIME=$(date --iso-8601=seconds)"
echo "Completed worker ${SLURM_ARRAY_TASK_ID}"
echo "=================================================="