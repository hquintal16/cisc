#!/bin/bash
#SBATCH --job-name=cisc_cluster_tifs_pngs
#SBATCH --output=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%j.out
#SBATCH --error=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%j.err
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --time=04:00:00

set -euo pipefail

ROOT=${CISC_ROOT:-/proj/sebastianlab/users/hquintal/cisc}
PY_SCRIPT=${ROOT}/scripts/f_cluster_tifs_pngs_2023.py
YEAR=${CISC_TIF_YEAR:-2016}
REFERENCE_RASTER=${CISC_TIF_REFERENCE_RASTER:-${ROOT}/data/out/04_rp/heatindex_idf_hi1_return_periods_bbox.nc}
WORKERS=${CISC_TIF_WORKERS:-${SLURM_CPUS_PER_TASK:-10}}
OVERWRITE=${CISC_TIF_OVERWRITE:-0}

mkdir -p "${ROOT}/logs"

module purge
module load python

source /proj/sebastianlab/users/hquintal/venvs/eca-hi/bin/activate

export CISC_ROOT="${ROOT}"
export MPLBACKEND=Agg
export GDAL_CACHEMAX=4096
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1

if [[ ! -f "${PY_SCRIPT}" ]]; then
  echo "ERROR: Python script not found: ${PY_SCRIPT}" >&2
  exit 1
fi

if [[ ! -f "${REFERENCE_RASTER}" ]]; then
  echo "ERROR: Reference raster not found: ${REFERENCE_RASTER}" >&2
  exit 1
fi

if [[ ! -d "${ROOT}/data/out/06_clusters/filtered_rp25" ]]; then
  echo "ERROR: filtered_rp25 directory not found:" >&2
  echo "  ${ROOT}/data/out/06_clusters/filtered_rp25" >&2
  exit 1
fi

ARGS=(
  --root "${ROOT}"
  --year "${YEAR}"
  --reference-raster "${REFERENCE_RASTER}"
  --workers "${WORKERS}"
)

# Optional environment variables:
#   CISC_STATE_SHP=/path/to/state_boundaries.shp
#   CISC_TIF_OVERWRITE=1
#   CISC_TIF_DATASETS="heatindex_0.25 precipitation_0.25"
if [[ -n "${CISC_STATE_SHP:-}" ]]; then
  ARGS+=(--state-shapefile "${CISC_STATE_SHP}")
fi

if [[ "${OVERWRITE}" == "1" ]]; then
  ARGS+=(--overwrite)
fi

if [[ -n "${CISC_TIF_DATASETS:-}" ]]; then
  read -r -a DATASET_ARRAY <<< "${CISC_TIF_DATASETS}"
  ARGS+=(--datasets "${DATASET_ARRAY[@]}")
fi

echo "Python: $(which python)"
python --version
echo "Script: ${PY_SCRIPT}"
echo "Reference: ${REFERENCE_RASTER}"
echo "Input root: ${ROOT}/data/out/06_clusters/filtered_rp25"
echo "Year: ${YEAR}"
echo "Workers: ${WORKERS}"
echo "Overwrite: ${OVERWRITE}"

python -u "${PY_SCRIPT}" "${ARGS[@]}"