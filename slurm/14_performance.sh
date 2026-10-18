#!/usr/bin/env bash
#SBATCH --job-name=cisc_warn_perf
#SBATCH --output=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%j.out
#SBATCH --error=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%j.err
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=08:00:00

set -euo pipefail

ROOT="${CISC_ROOT:-/proj/sebastianlab/users/hquintal/cisc}"
SCRIPT="${CISC_MODEL_PERFORMANCE_PYTHON:-${ROOT}/scripts/14_performance.py}"
VENV="${CISC_VENV:-/proj/sebastianlab/users/hquintal/venvs/eca-hi}"
ACTION="${1:-submit_all}"
OUT_DIR="${CISC_PERF_DIR:-${ROOT}/data/out/14_performance}"
R_USER_LIB="${CISC_R_LIBS_USER:-/nas/longleaf/home/hquintal/R/x86_64-pc-linux-gnu-library/4.5}"

mkdir -p "${ROOT}/logs" "${OUT_DIR}" "${R_USER_LIB}"

load_environment() {
    module purge
    module load python/3.12.4

    if module load r/4.5.0 2>/dev/null; then
        :
    elif module load r/4.5.1 2>/dev/null; then
        :
    else
        module load r
    fi

    if [[ ! -f "${VENV}/bin/activate" ]]; then
        echo "ERROR: missing Python environment: ${VENV}" >&2
        exit 1
    fi

    # shellcheck disable=SC1090
    source "${VENV}/bin/activate"

    export R_LIBS_USER="${R_USER_LIB}"
    export CISC_ROOT="${ROOT}"
    export CISC_PERF_DIR="${OUT_DIR}"
    export CISC_RSCRIPT="$(command -v Rscript)"

    python - <<'PY'
required = {
    "geopandas": "geopandas",
    "numpy": "numpy",
    "pandas": "pandas",
    "pyarrow": "pyarrow",
    "matplotlib": "matplotlib",
    "openpyxl": "openpyxl",
    "requests": "requests",
    "scipy": "scipy",
    "shapely": "shapely",
}
missing = []
for package, module in required.items():
    try:
        __import__(module)
    except Exception:
        missing.append(package)
if missing:
    raise SystemExit("Missing Python packages: " + ", ".join(missing))
PY

    Rscript --vanilla - "${R_USER_LIB}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
user_lib <- args[[1]]
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

repos <- "https://cloud.r-project.org"
archive <- c(
  stringfish = "https://cran.r-project.org/src/contrib/Archive/stringfish/stringfish_0.16.0.tar.gz",
  qs = "https://cran.r-project.org/src/contrib/Archive/qs/qs_0.27.3.tar.gz"
)

if (!requireNamespace("RApiSerialize", quietly = TRUE)) {
  install.packages("RApiSerialize", repos = repos, type = "source", lib = user_lib)
}
if (!requireNamespace("stringfish", quietly = TRUE)) {
  install.packages(archive[["stringfish"]], repos = NULL, type = "source", lib = user_lib)
}
if (!requireNamespace("qs", quietly = TRUE)) {
  install.packages(archive[["qs"]], repos = NULL, type = "source", lib = user_lib)
}

required <- c("RApiSerialize", "stringfish", "qs", "data.table", "arrow")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "))
}
RSCRIPT

    if [[ ! -f "${SCRIPT}" ]]; then
        echo "ERROR: missing Python script: ${SCRIPT}" >&2
        exit 1
    fi
}

run_mode() {
    local mode="$1"
    export CISC_WWA_COMBO_MODE="${mode}"
    python "${SCRIPT}"
}

if [[ "${ACTION}" == "submit_all" ]]; then
    if [[ ! -f "${SCRIPT}" ]]; then
        echo "ERROR: missing Python script: ${SCRIPT}" >&2
        exit 1
    fi
    PREP_JOB=$(sbatch --parsable "${BASH_SOURCE[0]}" prepare_shared)
    EVAL_JOB=$(sbatch --parsable --dependency="afterok:${PREP_JOB}" "${BASH_SOURCE[0]}" evaluate_all)
    SUMMARY_JOB=$(sbatch --parsable --dependency="afterok:${EVAL_JOB}" "${BASH_SOURCE[0]}" summarize_all)
    PLOT_JOB=$(sbatch --parsable --dependency="afterok:${SUMMARY_JOB}" "${BASH_SOURCE[0]}" plot_spatial)

    echo "Submitted preparation: ${PREP_JOB}"
    echo "Submitted evaluation: ${EVAL_JOB}"
    echo "Submitted summary: ${SUMMARY_JOB}"
    echo "Submitted spatial plotting: ${PLOT_JOB}"
    exit 0
fi

load_environment

case "${ACTION}" in
    check_environment)
        echo "Python: $(command -v python)"
        echo "Rscript: $(command -v Rscript)"
        echo "Script: ${SCRIPT}"
        echo "Output: ${OUT_DIR}"
        ;;

    prepare_shared)
        run_mode prepare_shared
        ;;

    evaluate_all)
        run_mode evaluate_all
        ;;

    summarize_all)
        run_mode summarize_all
        ;;

    plot_spatial)
        run_mode plot_spatial
        ;;

    run_all)
        run_mode prepare_shared
        run_mode evaluate_all
        run_mode summarize_all
        run_mode plot_spatial
        ;;

    *)
        echo "Usage: $0 {check_environment|submit_all|run_all|prepare_shared|evaluate_all|summarize_all|plot_spatial}" >&2
        exit 2
        ;;
esac