#!/bin/bash
#SBATCH --job-name=cisc_14_model_eval
#SBATCH --output=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%j.out
#SBATCH --error=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%j.err
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -euo pipefail

ROOT=${CISC_ROOT:-/proj/sebastianlab/users/hquintal/cisc}
DATA_OUT=${CISC_DATA_OUT:-${ROOT}/data/out}
PERF_DIR=${CISC_PERF_DIR:-${DATA_OUT}/13_performance}

PYTHON_ENV=${CISC_PYTHON_ENV:-/proj/sebastianlab/users/hquintal/venvs/eca-hi}
PY_SCRIPT=${ROOT}/scripts/14_model_evaluation_postprocess.py

[[ -f "${PY_SCRIPT}" ]] || {
    echo "ERROR: Missing Python script: ${PY_SCRIPT}" >&2
    exit 1
}

[[ -f "${PYTHON_ENV}/bin/activate" ]] || {
    echo "ERROR: Missing Python environment: ${PYTHON_ENV}" >&2
    exit 1
}

module purge
module load python

# shellcheck disable=SC1090
source "${PYTHON_ENV}/bin/activate"

export CISC_ROOT="${ROOT}"
export CISC_DATA_OUT="${DATA_OUT}"
export CISC_PERF_DIR="${PERF_DIR}"

python -c '
import numpy
import openpyxl
import pandas
import pyarrow
import scipy

print("Python dependencies verified")
print("  numpy:", numpy.__version__)
print("  pandas:", pandas.__version__)
print("  pyarrow:", pyarrow.__version__)
print("  scipy:", scipy.__version__)
print("  openpyxl:", openpyxl.__version__)
'

HEAT_FILE=${PERF_DIR}/metrics_county/metrics_heatindex_all_heat_multiwindow.parquet
PRECIP_FILE=${PERF_DIR}/metrics_county/metrics_precipitation_all_precip_multiwindow.parquet

[[ -f "${HEAT_FILE}" ]] || {
    echo "ERROR: Missing heat-index metrics: ${HEAT_FILE}" >&2
    exit 1
}

[[ -f "${PRECIP_FILE}" ]] || {
    echo "ERROR: Missing precipitation metrics: ${PRECIP_FILE}" >&2
    exit 1
}

echo "=================================================="
echo "14_model_evaluation_postprocess.sh"
echo "=================================================="
echo "ROOT=${ROOT}"
echo "PERF_DIR=${PERF_DIR}"
echo "WINDOW_TAG=20020101_20231231"
echo "RETURN_PERIODS=5,10,25,50"
echo "HOSTNAME=$(hostname)"
echo "START_TIME=$(date --iso-8601=seconds)"
echo "=================================================="

python "${PY_SCRIPT}"

echo "=================================================="
echo "END_TIME=$(date --iso-8601=seconds)"
echo "Completed Stage 14 model evaluation"
echo "=================================================="