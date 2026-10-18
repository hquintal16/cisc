#!/bin/bash
#SBATCH --job-name=cisc_13_driver
#SBATCH --output=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%A_%a.out
#SBATCH --error=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%A_%a.err
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

set -euo pipefail

ROOT=${CISC_ROOT:-/proj/sebastianlab/users/hquintal/cisc}
DATA_OUT=${CISC_DATA_OUT:-${ROOT}/data/out}
WWA_DIR=${CISC_WWA_DIR:-${ROOT}/figures/fig_02_wwa}
SELF=${ROOT}/slurm/13_wwa_county_evaluation.sh
PY_SCRIPT=${ROOT}/scripts/13_wwa_county_evaluation.py
PYTHON_ENV=${CISC_PYTHON_ENV:-/proj/sebastianlab/users/hquintal/venvs/eca-hi}
PERF_DIR=${CISC_PERF_DIR:-${DATA_OUT}/13_performance}
STAGE=${CISC_WWA_EVAL_STAGE:-driver}
LOOKUP_WORKERS=${CISC_WWA_LOOKUP_WORKERS:-4}
MODELED_WORKERS=${CISC_WWA_MODELED_WORKERS:-8}
EVAL_WORKERS=${CISC_WWA_EVAL_WORKERS:-8}
OVERWRITE_LOOKUPS=${CISC_OVERWRITE_LOOKUPS:-false}
OVERWRITE_OBSERVED=${CISC_OVERWRITE_OBSERVED:-false}
OVERWRITE_MODELED=${CISC_OVERWRITE_MODELED:-false}
OVERWRITE_METRICS=${CISC_OVERWRITE_METRICS:-true}
EVALUATE_INDIVIDUAL=${CISC_EVALUATE_INDIVIDUAL_WWA_TAGS:-false}
EVAL_SCENARIOS=filtered_rp5,filtered_rp10,filtered_rp25,filtered_rp50
R_LIBS_USER=${R_LIBS_USER:-/nas/longleaf/home/hquintal/R/x86_64-pc-linux-gnu-library/4.5}

die() { echo "ERROR: $*" >&2; exit 1; }
valid_int() { [[ "$1" =~ ^[1-9][0-9]*$ ]] || die "$2 must be a positive integer: $1"; }

load_environment() {
    module purge
    module load python
    module load r

    export R_LIBS_USER=/nas/longleaf/home/hquintal/R/x86_64-pc-linux-gnu-library/4.5

    [[ -f "${PYTHON_ENV}/bin/activate" ]] || {
        die "Missing Python environment: ${PYTHON_ENV}"
    }

    source "${PYTHON_ENV}/bin/activate"

    command -v python >/dev/null 2>&1 || {
        die "python is unavailable"
    }

    command -v Rscript >/dev/null 2>&1 || {
        die "Rscript is unavailable"
    }
}

check_deps() {
    python -c 'import geopandas, numpy, pandas, pyarrow, rasterio, requests, shapely; print("Python dependencies verified"); print("  geopandas:", geopandas.__version__); print("  pandas:", pandas.__version__); print("  pyarrow:", pyarrow.__version__); print("  requests:", requests.__version__)'
    Rscript -e 'pkgs <- c("qs", "arrow", "data.table"); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop("Missing R packages: ", paste(missing, collapse=", ")); cat("R dependencies verified:", paste(pkgs, collapse=", "), "\n")'
}

export_common() {
    export CISC_ROOT="${ROOT}"
    export CISC_DATA_OUT="${DATA_OUT}"
    export CISC_PERF_DIR="${PERF_DIR}"
    export CISC_EVAL_SCENARIOS="${EVAL_SCENARIOS}"
    export CISC_OVERWRITE_LOOKUPS="${OVERWRITE_LOOKUPS}"
    export CISC_OVERWRITE_OBSERVED="${OVERWRITE_OBSERVED}"
    export CISC_OVERWRITE_MODELED="${OVERWRITE_MODELED}"
    export CISC_OVERWRITE_METRICS="${OVERWRITE_METRICS}"
    export CISC_EVALUATE_INDIVIDUAL_WWA_TAGS="${EVALUATE_INDIVIDUAL}"
    export CISC_WWA_DIR="${WWA_DIR}"
}

set_mode() {
    case "${STAGE}" in
        setup) export CISC_WWA_EVAL_MODE=setup CISC_WWA_WORKER_COUNT=1 ;;
        lookup_worker) export CISC_WWA_EVAL_MODE=lookup_worker CISC_WWA_WORKER_COUNT="${LOOKUP_WORKERS}" ;;
        modeled_worker) export CISC_WWA_EVAL_MODE=modeled_worker CISC_WWA_WORKER_COUNT="${MODELED_WORKERS}" ;;
        eval_worker) export CISC_WWA_EVAL_MODE=eval_worker CISC_WWA_WORKER_COUNT="${EVAL_WORKERS}" ;;
        summarize) export CISC_WWA_EVAL_MODE=summarize CISC_WWA_WORKER_COUNT=1 ;;
        *) die "Unsupported stage: ${STAGE}" ;;
    esac
}

run_stage() {
    [[ -f "${PY_SCRIPT}" ]] || die "Missing Python script: ${PY_SCRIPT}"
    load_environment
    check_deps
    export_common
    set_mode
    mkdir -p "${ROOT}/logs" "${PERF_DIR}/_control"
    cd "${ROOT}"
    echo "STAGE=${STAGE} MODE=${CISC_WWA_EVAL_MODE} WORKER=${SLURM_ARRAY_TASK_ID:-NA}/${CISC_WWA_WORKER_COUNT}"
    echo "SCENARIOS=${CISC_EVAL_SCENARIOS}"
    python "${PY_SCRIPT}"
}

submit_pipeline() {
    mkdir -p "${ROOT}/logs" "${PERF_DIR}/_control"
    [[ -f "${SELF}" ]] || die "Missing shell script: ${SELF}"
    [[ -f "${PY_SCRIPT}" ]] || die "Missing Python script: ${PY_SCRIPT}"
    valid_int "${LOOKUP_WORKERS}" CISC_WWA_LOOKUP_WORKERS
    valid_int "${MODELED_WORKERS}" CISC_WWA_MODELED_WORKERS
    valid_int "${EVAL_WORKERS}" CISC_WWA_EVAL_WORKERS
    lmax=$((LOOKUP_WORKERS-1)); mmax=$((MODELED_WORKERS-1)); emax=$((EVAL_WORKERS-1))
    EXP="ALL,CISC_ROOT=${ROOT},CISC_DATA_OUT=${DATA_OUT},CISC_WWA_DIR=${WWA_DIR},CISC_PERF_DIR=${PERF_DIR},CISC_PYTHON_ENV=${PYTHON_ENV},CISC_OVERWRITE_LOOKUPS=${OVERWRITE_LOOKUPS},CISC_OVERWRITE_OBSERVED=${OVERWRITE_OBSERVED},CISC_OVERWRITE_MODELED=${OVERWRITE_MODELED},CISC_OVERWRITE_METRICS=${OVERWRITE_METRICS},CISC_EVALUATE_INDIVIDUAL_WWA_TAGS=${EVALUATE_INDIVIDUAL},R_LIBS_USER=${R_LIBS_USER}"
    
    setup=$(sbatch --parsable --job-name=cisc_13_setup --partition=general --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=32G --time=04:00:00 --output="${ROOT}/logs/cisc_13_setup_%j.out" --error="${ROOT}/logs/cisc_13_setup_%j.err" --export="${EXP},CISC_WWA_EVAL_STAGE=setup" "${SELF}")
    lookup=$(sbatch --parsable --job-name=cisc_13_lookup --partition=general --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=64G --time=12:00:00 --array="0-${lmax}" --dependency="afterok:${setup}" --output="${ROOT}/logs/cisc_13_lookup_%A_%a.out" --error="${ROOT}/logs/cisc_13_lookup_%A_%a.err" --export="${EXP},CISC_WWA_EVAL_STAGE=lookup_worker,CISC_WWA_LOOKUP_WORKERS=${LOOKUP_WORKERS}" "${SELF}")
    modeled=$(sbatch --parsable --job-name=cisc_13_modeled --partition=general --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=64G --time=1-00:00:00 --array="0-${mmax}" --dependency="afterok:${lookup}" --output="${ROOT}/logs/cisc_13_modeled_%A_%a.out" --error="${ROOT}/logs/cisc_13_modeled_%A_%a.err" --export="${EXP},CISC_WWA_EVAL_STAGE=modeled_worker,CISC_WWA_MODELED_WORKERS=${MODELED_WORKERS}" "${SELF}")
    evalj=$(sbatch --parsable --job-name=cisc_13_eval --partition=general --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=32G --time=12:00:00 --array="0-${emax}" --dependency="afterok:${modeled}" --output="${ROOT}/logs/cisc_13_eval_%A_%a.out" --error="${ROOT}/logs/cisc_13_eval_%A_%a.err" --export="${EXP},CISC_WWA_EVAL_STAGE=eval_worker,CISC_WWA_EVAL_WORKERS=${EVAL_WORKERS}" "${SELF}")
    summary=$(sbatch --parsable --job-name=cisc_13_summary --partition=general --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=32G --time=04:00:00 --dependency="afterok:${evalj}" --output="${ROOT}/logs/cisc_13_summary_%j.out" --error="${ROOT}/logs/cisc_13_summary_%j.err" --export="${EXP},CISC_WWA_EVAL_STAGE=summarize" "${SELF}")

    echo "Submitted: setup=${setup} lookup=${lookup} modeled=${modeled} eval=${evalj} summary=${summary}"
    echo "Monitor: squeue -j ${setup},${lookup},${modeled},${evalj},${summary}"
}

if [[ "${STAGE}" == driver ]]; then submit_pipeline; else run_stage; fi