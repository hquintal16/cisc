#!/bin/bash
#SBATCH --job-name=cisc_11_driver
#SBATCH --output=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%A_%a.out
#SBATCH --error=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%A_%a.err
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=03:00:00

set -euo pipefail


# =============================================================================
# CISC 11: Cluster RP >= 5 exceedances, then retain clusters reaching
#           RP >= 5, 10, 25, and 50 years
#
# Stages:
#   driver
#   build_manifest
#   cluster_worker
#   summarize
# =============================================================================


# ── PATHS ───────────────────────────────────────────────────────────────────

ROOT=${CISC_ROOT:-/proj/sebastianlab/users/hquintal/cisc}
DATA_OUT=${CISC_DATA_OUT:-${ROOT}/data/out}

SELF=${ROOT}/slurm/11_cluster_rp5_filter.sh
PY_SCRIPT=${ROOT}/scripts/11_cluster_rp5_filter.py

EXCEEDANCE_DIR=${CISC_EXCEEDANCE_DIR:-${DATA_OUT}/05_exceedances}
CLUSTER_DIR=${CISC_CLUSTER_DIR:-${DATA_OUT}/11_clusters}

MANIFEST=${CISC_CLUSTER_MANIFEST:-${CLUSTER_DIR}/_control/cluster_manifest.csv}


# ── PIPELINE SETTINGS ───────────────────────────────────────────────────────

STAGE=${CISC_CLUSTER_STAGE:-driver}

EPS=${CISC_DBSCAN_EPS:-4.0}
MIN_SAMPLES=${CISC_DBSCAN_MIN_SAMPLES:-4}
TIME_STEP_HOURS=${CISC_CLUSTER_TIME_STEP_HOURS:-24}

CLUSTER_FLOOR=${CISC_CLUSTER_FLOOR_YEARS:-5}
FILTER_RPS=${CISC_CLUSTER_FILTER_RPS:-5,10,25,50}

OVERWRITE=${CISC_OVERWRITE:-false}

WORKER_COUNT=${CISC_CLUSTER_WORKER_COUNT:-4}
WORKER_MEM=${CISC_CLUSTER_WORKER_MEM:-150G}
WORKER_TIME=${CISC_CLUSTER_WORKER_TIME:-2-00:00:00}

PYTHON_ENV=${CISC_PYTHON_ENV:-/proj/sebastianlab/users/hquintal/venvs/eca-hi}


# ── HELPERS ─────────────────────────────────────────────────────────────────

die() {
    echo "ERROR: $*" >&2
    exit 1
}


validate_integer() {
    local value=$1
    local name=$2

    [[ "${value}" =~ ^[1-9][0-9]*$ ]] || {
        die "${name} must be a positive integer; received: ${value}"
    }
}


load_environment() {
    module purge
    module load python
    module load r
    
    export R_LIBS_USER=${CISC_R_LIBS_USER:-/nas/longleaf/home/hquintal/R/x86_64-pc-linux-gnu-library/4.5}

    [[ -f "${PYTHON_ENV}/bin/activate" ]] || {
        die "Missing Python environment: ${PYTHON_ENV}"
    }

    # shellcheck disable=SC1090
    source "${PYTHON_ENV}/bin/activate"

    command -v python >/dev/null 2>&1 || {
        die "python is unavailable after activating ${PYTHON_ENV}"
    }

    command -v Rscript >/dev/null 2>&1 || {
        die "Rscript is unavailable after loading the R module"
    }
}


validate_dependencies() {
    python -c \
        'import numpy, pandas, pyarrow, sklearn; print("Python dependencies verified"); print("  numpy:", numpy.__version__); print("  pandas:", pandas.__version__); print("  pyarrow:", pyarrow.__version__); print("  sklearn:", sklearn.__version__)'

    Rscript -e \
        'pkgs <- c("qs", "arrow", "data.table"); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) stop("Missing required R packages: ", paste(missing, collapse = ", ")); cat("R dependencies verified:", paste(pkgs, collapse = ", "), "\n")'
}


export_environment() {
    export CISC_ROOT="${ROOT}"
    export CISC_DATA_OUT="${DATA_OUT}"

    export CISC_EXCEEDANCE_DIR="${EXCEEDANCE_DIR}"
    export CISC_CLUSTER_DIR="${CLUSTER_DIR}"
    export CISC_CLUSTER_MANIFEST="${MANIFEST}"

    export CISC_DBSCAN_EPS="${EPS}"
    export CISC_DBSCAN_MIN_SAMPLES="${MIN_SAMPLES}"
    export CISC_CLUSTER_TIME_STEP_HOURS="${TIME_STEP_HOURS}"

    export CISC_CLUSTER_FLOOR_YEARS="${CLUSTER_FLOOR}"
    export CISC_CLUSTER_FILTER_RPS="${FILTER_RPS}"

    export CISC_CLUSTER_WORKER_COUNT="${WORKER_COUNT}"
    export CISC_OVERWRITE="${OVERWRITE}"
}


set_mode() {
    case "${STAGE}" in
        build_manifest)
            export CISC_CLUSTER_MODE=build_manifest
            ;;

        cluster_worker)
            export CISC_CLUSTER_MODE=cluster_worker
            ;;

        summarize)
            export CISC_CLUSTER_MODE=summarize
            ;;

        *)
            die "Unsupported CISC_CLUSTER_STAGE: ${STAGE}"
            ;;
    esac
}


print_stage_settings() {
    echo "=================================================="
    echo "11_cluster_rp5_filter.sh"
    echo "=================================================="
    echo "CISC_CLUSTER_STAGE=${STAGE}"
    echo "CISC_CLUSTER_MODE=${CISC_CLUSTER_MODE}"
    echo
    echo "CISC_ROOT=${CISC_ROOT}"
    echo "CISC_DATA_OUT=${CISC_DATA_OUT}"
    echo "CISC_EXCEEDANCE_DIR=${CISC_EXCEEDANCE_DIR}"
    echo "CISC_CLUSTER_DIR=${CISC_CLUSTER_DIR}"
    echo "CISC_CLUSTER_MANIFEST=${CISC_CLUSTER_MANIFEST}"
    echo
    echo "CISC_DBSCAN_EPS=${CISC_DBSCAN_EPS}"
    echo "CISC_DBSCAN_MIN_SAMPLES=${CISC_DBSCAN_MIN_SAMPLES}"
    echo "CISC_CLUSTER_TIME_STEP_HOURS=${CISC_CLUSTER_TIME_STEP_HOURS}"
    echo "CISC_CLUSTER_FLOOR_YEARS=${CISC_CLUSTER_FLOOR_YEARS}"
    echo "CISC_CLUSTER_FILTER_RPS=${CISC_CLUSTER_FILTER_RPS}"
    echo "CISC_CLUSTER_WORKER_COUNT=${CISC_CLUSTER_WORKER_COUNT}"
    echo "CISC_OVERWRITE=${CISC_OVERWRITE}"
    echo
    echo "SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"
    echo "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID:-NA}"
    echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-NA}"
    echo "HOSTNAME=$(hostname)"
    echo "START_TIME=$(date --iso-8601=seconds)"
    echo "=================================================="
}


# ── RUN ONE PIPELINE STAGE ─────────────────────────────────────────────────

run_stage() {
    [[ -f "${PY_SCRIPT}" ]] || {
        die "Missing Python script: ${PY_SCRIPT}"
    }

    load_environment
    validate_dependencies
    export_environment
    set_mode

    mkdir -p "${ROOT}/logs"
    mkdir -p "${CLUSTER_DIR}/_control"

    cd "${ROOT}"

    print_stage_settings

    python "${PY_SCRIPT}"

    echo "=================================================="
    echo "END_TIME=$(date --iso-8601=seconds)"
    echo "Completed stage: ${STAGE}"
    echo "=================================================="
}


# ── SUBMIT COMPLETE PIPELINE ───────────────────────────────────────────────

submit_pipeline() {
    mkdir -p "${ROOT}/logs"
    mkdir -p "${CLUSTER_DIR}/_control"

    [[ -f "${SELF}" ]] || {
        die "Missing shell script: ${SELF}"
    }

    [[ -f "${PY_SCRIPT}" ]] || {
        die "Missing Python script: ${PY_SCRIPT}"
    }

    [[ -d "${EXCEEDANCE_DIR}" ]] || {
        die "Missing exceedance directory: ${EXCEEDANCE_DIR}"
    }

    validate_integer "${WORKER_COUNT}" "CISC_CLUSTER_WORKER_COUNT"
    validate_integer "${MIN_SAMPLES}" "CISC_DBSCAN_MIN_SAMPLES"
    validate_integer "${TIME_STEP_HOURS}" "CISC_CLUSTER_TIME_STEP_HOURS"
    validate_integer "${CLUSTER_FLOOR}" "CISC_CLUSTER_FLOOR_YEARS"

    worker_max=$((WORKER_COUNT - 1))

    EXPORTS="ALL"
    EXPORTS+=",CISC_ROOT=${ROOT}"
    EXPORTS+=",CISC_DATA_OUT=${DATA_OUT}"
    EXPORTS+=",CISC_EXCEEDANCE_DIR=${EXCEEDANCE_DIR}"
    EXPORTS+=",CISC_CLUSTER_DIR=${CLUSTER_DIR}"
    EXPORTS+=",CISC_CLUSTER_MANIFEST=${MANIFEST}"
    EXPORTS+=",CISC_DBSCAN_EPS=${EPS}"
    EXPORTS+=",CISC_DBSCAN_MIN_SAMPLES=${MIN_SAMPLES}"
    EXPORTS+=",CISC_CLUSTER_TIME_STEP_HOURS=${TIME_STEP_HOURS}"
    EXPORTS+=",CISC_CLUSTER_FLOOR_YEARS=${CLUSTER_FLOOR}"
    EXPORTS+=",CISC_CLUSTER_FILTER_RPS=${FILTER_RPS}"
    EXPORTS+=",CISC_CLUSTER_WORKER_COUNT=${WORKER_COUNT}"
    EXPORTS+=",CISC_CLUSTER_WORKER_MEM=${WORKER_MEM}"
    EXPORTS+=",CISC_CLUSTER_WORKER_TIME=${WORKER_TIME}"
    EXPORTS+=",CISC_PYTHON_ENV=${PYTHON_ENV}"
    EXPORTS+=",CISC_OVERWRITE=${OVERWRITE}"

    manifest_job=$(
        sbatch \
            --parsable \
            --job-name=cisc_1_manifest \
            --partition=general \
            --nodes=1 \
            --ntasks=1 \
            --cpus-per-task=1 \
            --mem=8G \
            --time=00:30:00 \
            --output="${ROOT}/logs/cisc_11_manifest_%j.out" \
            --error="${ROOT}/logs/cisc_11_manifest_%j.err" \
            --export="${EXPORTS},CISC_CLUSTER_STAGE=build_manifest" \
            "${SELF}"
    )

    cluster_job=$(
        sbatch \
            --parsable \
            --job-name=cisc_11_cluster \
            --partition=general \
            --nodes=1 \
            --ntasks=1 \
            --cpus-per-task=1 \
            --mem="${WORKER_MEM}" \
            --time="${WORKER_TIME}" \
            --array="0-${worker_max}" \
            --dependency="afterok:${manifest_job}" \
            --output="${ROOT}/logs/cisc_11_cluster_%A_%a.out" \
            --error="${ROOT}/logs/cisc_11_cluster_%A_%a.err" \
            --export="${EXPORTS},CISC_CLUSTER_STAGE=cluster_worker" \
            "${SELF}"
    )

    summary_job=$(
        sbatch \
            --parsable \
            --job-name=cisc_11_summary \
            --partition=general \
            --nodes=1 \
            --ntasks=1 \
            --cpus-per-task=1 \
            --mem=16G \
            --time=01:00:00 \
            --dependency="afterok:${cluster_job}" \
            --output="${ROOT}/logs/cisc_11_summary_%j.out" \
            --error="${ROOT}/logs/cisc_11_summary_%j.err" \
            --export="${EXPORTS},CISC_CLUSTER_STAGE=summarize" \
            "${SELF}"
    )

    echo
    echo "Submitted CISC RP5 clustering pipeline:"
    echo
    echo "  Manifest job:"
    echo "    ${manifest_job}"
    echo
    echo "  Cluster array:"
    echo "    ${cluster_job}"
    echo "    Array indices: 0-${worker_max}"
    echo "    Workers: ${WORKER_COUNT}"
    echo "    Memory per worker: ${WORKER_MEM}"
    echo "    Time per worker: ${WORKER_TIME}"
    echo "    Dependency: afterok:${manifest_job}"
    echo
    echo "  Summary job:"
    echo "    ${summary_job}"
    echo "    Dependency: afterok:${cluster_job}"
    echo
    echo "  Clustering floor:"
    echo "    RP >= ${CLUSTER_FLOOR} years"
    echo
    echo "  Retained cluster products:"
    echo "    RP >= ${FILTER_RPS} years"
    echo
    echo "Monitor with:"
    echo "  squeue -j ${manifest_job},${cluster_job},${summary_job}"
}


# ── ENTRY POINT ─────────────────────────────────────────────────────────────

if [[ "${STAGE}" == "driver" ]]; then
    submit_pipeline
else
    run_stage
fi