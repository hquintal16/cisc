#!/bin/bash

# =============================================================================
# Mesonet daily union and Ar correlation workflow
#
# Launch from the CISC project root with:
#
#   bash slurm/10_mesonet_warning_timeseries.sh
#
# The driver submits:
#
#   1. One array task per year, 2002-2023.
#   2. A dependent combine, correlation, and plotting job.
#
# Internal stages:
#
#   CISC_MESONET_STAGE=annual
#   CISC_MESONET_STAGE=combine
# =============================================================================

set -euo pipefail

ROOT=${CISC_ROOT:-/proj/sebastianlab/users/hquintal/cisc}

SELF="${ROOT}/slurm/10_mesonet_warning_timeseries.sh"
PY_SCRIPT="${ROOT}/scripts/10_mesonet_warning_timeseries.py"

LOG_DIR="${ROOT}/logs"
OUT_DIR="${ROOT}/data/out/10_mesonet"
ANNUAL_DIR="${OUT_DIR}/annual_combined"
CORRELATION_DIR="${OUT_DIR}/ar_correlations"
FIGURE_DIR="${ROOT}/figures/fig_s09_mesonet_ar"

STAGE=${CISC_MESONET_STAGE:-driver}

mkdir -p \
  "${LOG_DIR}" \
  "${ANNUAL_DIR}" \
  "${CORRELATION_DIR}" \
  "${FIGURE_DIR}"

if [[ ! -f "${SELF}" ]]; then
    echo "ERROR: SLURM script not found: ${SELF}" >&2
    exit 1
fi

if [[ ! -f "${PY_SCRIPT}" ]]; then
    echo "ERROR: Python script not found: ${PY_SCRIPT}" >&2
    exit 1
fi


# =============================================================================
# Driver stage
# =============================================================================

if [[ "${STAGE}" == "driver" ]]; then
    echo "Submitting Mesonet annual daily-union array..."

    ANNUAL_JOB=$(
        sbatch \
          --parsable \
          --job-name=cisc_10_mesonet_union \
          --output="${LOG_DIR}/%x_%A_%a.out" \
          --error="${LOG_DIR}/%x_%A_%a.err" \
          --partition=general \
          --nodes=1 \
          --ntasks=1 \
          --cpus-per-task=1 \
          --mem=64G \
          --time=1-00:00:00 \
          --array=0-21%6 \
          --export=ALL,CISC_ROOT="${ROOT}",CISC_MESONET_STAGE=annual \
          "${SELF}"
    )

    echo "Annual array job: ${ANNUAL_JOB}"

    COMBINE_JOB=$(
        sbatch \
          --parsable \
          --job-name=cisc_10_mesonet_ar \
          --output="${LOG_DIR}/%x_%j.out" \
          --error="${LOG_DIR}/%x_%j.err" \
          --partition=general \
          --nodes=1 \
          --ntasks=1 \
          --cpus-per-task=1 \
          --mem=32G \
          --time=04:00:00 \
          --dependency="afterok:${ANNUAL_JOB}" \
          --export=ALL,CISC_ROOT="${ROOT}",CISC_MESONET_STAGE=combine \
          "${SELF}"
    )

    echo "Combine job: ${COMBINE_JOB}"
    echo
    echo "Monitor:"
    echo "  squeue -j ${ANNUAL_JOB},${COMBINE_JOB}"
    echo
    echo "Accounting:"
    echo "  sacct -j ${ANNUAL_JOB} \\"
    echo "    --format=JobID,State,ExitCode,Elapsed,MaxRSS,ReqMem"

    exit 0
fi


# =============================================================================
# Runtime environment
# =============================================================================

module purge
module load python/3.12.4

source \
  /proj/sebastianlab/users/hquintal/venvs/eca-hi/bin/activate

export CISC_ROOT="${ROOT}"

echo "============================================================"
echo "Stage: ${STAGE}"
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-not_set}"
echo "Python: $(which python)"
echo "Started: $(date)"
echo "============================================================"

python - <<'PY'
import geopandas
import matplotlib
import numpy
import pandas
import pyarrow
import pyogrio
import pyproj
import scipy
import shapely

print("Python geospatial environment loaded.")
print("geopandas:", geopandas.__version__)
print("shapely:", shapely.__version__)
print("pyogrio:", pyogrio.__version__)
print("pyarrow:", pyarrow.__version__)
print("scipy:", scipy.__version__)
PY


# =============================================================================
# Annual stage
# =============================================================================

if [[ "${STAGE}" == "annual" ]]; then
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        echo \
          "ERROR: SLURM_ARRAY_TASK_ID is not defined." \
          >&2
        exit 1
    fi

    YEARS=($(seq 2002 2023))
    YEAR=${YEARS[$SLURM_ARRAY_TASK_ID]}

    echo "Array task: ${SLURM_ARRAY_TASK_ID}"
    echo "Processing year: ${YEAR}"

    python -u "${PY_SCRIPT}" \
      --stage annual \
      --year "${YEAR}" \
      --root "${ROOT}" \
      --overwrite

    echo "Completed year: ${YEAR}"
    echo "Finished: $(date)"

    exit 0
fi


# =============================================================================
# Combine, correlate, and plot stage
# =============================================================================

if [[ "${STAGE}" == "combine" ]]; then
    echo "Combining annual files and correlating warnings with Ar."

    python -u "${PY_SCRIPT}" \
      --stage combine \
      --root "${ROOT}"

    echo "Finished: $(date)"

    exit 0
fi


echo "ERROR: Unknown CISC_MESONET_STAGE: ${STAGE}" >&2
exit 1