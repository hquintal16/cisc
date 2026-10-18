#!/bin/bash
#SBATCH --job-name=cisc_11_mesonet_corr
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

PY_SCRIPT="${ROOT}/scripts/11_mesonet_warning_correlations.py"

mkdir -p \
  "${ROOT}/logs" \
  "${ROOT}/data/out/10_mesonet/covariance_correlations" \
  "${ROOT}/figures/fig_s09_mesonet_covariance"

module purge
module load python/3.12.4

source \
  /proj/sebastianlab/users/hquintal/venvs/eca-hi/bin/activate

export CISC_ROOT="${ROOT}"

python -m py_compile \
  "${PY_SCRIPT}"

python -u "${PY_SCRIPT}" \
  --root "${ROOT}"