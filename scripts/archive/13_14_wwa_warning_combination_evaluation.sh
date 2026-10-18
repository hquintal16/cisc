#!/usr/bin/env bash
#SBATCH --job-name=cisc_wwa_combo
#SBATCH --output=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%A_%a.out
#SBATCH --error=/proj/sebastianlab/users/hquintal/cisc/logs/%x_%A_%a.err
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

ROOT="${CISC_ROOT:-/proj/sebastianlab/users/hquintal/cisc}"
PYTHON_SCRIPT="${CISC_COMBO_PYTHON:-${ROOT}/scripts/13_14_wwa_warning_combination_evaluation.py}"
VENV="${CISC_VENV:-/proj/sebastianlab/users/hquintal/venvs/eca-hi}"
ACTION="${1:-submit_all}"

COMBINATIONS_PER_TASK="${CISC_COMBINATIONS_PER_TASK:-25}"
MAX_CONCURRENT="${CISC_MAX_CONCURRENT:-50}"

R_USER_LIB="${CISC_R_LIBS_USER:-/nas/longleaf/home/hquintal/R/x86_64-pc-linux-gnu-library/4.5}"

RAPISERIALIZE_REPO="https://cloud.r-project.org"
STRINGFISH_URL="https://cran.r-project.org/src/contrib/Archive/stringfish/stringfish_0.16.0.tar.gz"
QS_URL="https://cran.r-project.org/src/contrib/Archive/qs/qs_0.27.3.tar.gz"

mkdir -p "${ROOT}/logs" "${R_USER_LIB}"

load_modules() {
    module purge

    if module load python/3.12.4 2>/dev/null; then
        :
    elif module load python 2>/dev/null; then
        :
    else
        echo "ERROR: no usable Python module was found." >&2
        exit 1
    fi

    if module load r/4.5.0 2>/dev/null; then
        :
    elif module load r/4.5.1 2>/dev/null; then
        :
    elif module load r 2>/dev/null; then
        :
    else
        echo "ERROR: no usable R module was found." >&2
        exit 1
    fi

    if [[ ! -f "${VENV}/bin/activate" ]]; then
        echo "ERROR: missing Python environment: ${VENV}" >&2
        exit 1
    fi

    # shellcheck disable=SC1090
    source "${VENV}/bin/activate"
}

remove_stale_r_locks() {
    find "${R_USER_LIB}" \
        -maxdepth 1 \
        -type d \
        -name '00LOCK*' \
        -exec rm -rf {} +
}

ensure_r_packages() {
    export R_LIBS_USER="${R_USER_LIB}"

    local rscript_path
    rscript_path="$(command -v Rscript || true)"

    if [[ -z "${rscript_path}" ]]; then
        echo "ERROR: Rscript was not found after loading the R module." >&2
        exit 1
    fi

    export CISC_RSCRIPT="${rscript_path}"

    remove_stale_r_locks

    "${CISC_RSCRIPT}" --vanilla - \
        "${R_USER_LIB}" \
        "${RAPISERIALIZE_REPO}" \
        "${STRINGFISH_URL}" \
        "${QS_URL}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)

user_lib <- args[[1]]
cran_repo <- args[[2]]
stringfish_url <- args[[3]]
qs_url <- args[[4]]

dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

message("[R] User library: ", user_lib)
message("[R] Library paths: ", paste(.libPaths(), collapse = " | "))

install_cran_if_missing <- function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    message(
      "[R] ", package, " already installed: ",
      as.character(packageVersion(package))
    )
    return(invisible(TRUE))
  }

  message("[R] Installing CRAN dependency: ", package)

  install.packages(
    package,
    repos = cran_repo,
    type = "source",
    lib = user_lib,
    dependencies = TRUE
  )

  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package ", package, " was not available after installation.")
  }

  message(
    "[R] Installed ", package, ": ",
    as.character(packageVersion(package))
  )
}

install_archive_if_missing <- function(package, required_version, url) {
  if (requireNamespace(package, quietly = TRUE)) {
    installed_version <- as.character(packageVersion(package))

    if (identical(installed_version, required_version)) {
      message("[R] ", package, " already installed: ", installed_version)
      return(invisible(TRUE))
    }

    message(
      "[R] Reinstalling ", package,
      "; found ", installed_version,
      ", require ", required_version
    )
  } else {
    message("[R] Installing ", package, " from ", url)
  }

  install.packages(
    url,
    repos = NULL,
    type = "source",
    lib = user_lib
  )

  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package ", package, " was not available after installation.")
  }

  installed_version <- as.character(packageVersion(package))

  if (!identical(installed_version, required_version)) {
    stop(
      "Installed ", package, " version ", installed_version,
      "; expected ", required_version, "."
    )
  }

  message("[R] Installed ", package, ": ", installed_version)
}

install_cran_if_missing("RApiSerialize")

install_archive_if_missing(
  "stringfish",
  "0.16.0",
  stringfish_url
)

install_archive_if_missing(
  "qs",
  "0.27.3",
  qs_url
)

required <- c(
  "RApiSerialize",
  "stringfish",
  "qs",
  "data.table",
  "arrow"
)

missing <- required[
  !vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing)) {
  stop(
    "Required R packages remain unavailable: ",
    paste(missing, collapse = ", "),
    ". Install them into ",
    user_lib
  )
}

message(
  "[R] Required packages available: ",
  paste(required, collapse = ", ")
)
RSCRIPT

    echo "Python: $(command -v python)"
    echo "Rscript: ${CISC_RSCRIPT}"
    echo "R_LIBS_USER: ${R_LIBS_USER}"
}

validate_inputs() {
    if [[ ! -f "${PYTHON_SCRIPT}" ]]; then
        echo "ERROR: missing Python workflow: ${PYTHON_SCRIPT}" >&2
        exit 1
    fi

    if ! [[ "${COMBINATIONS_PER_TASK}" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: CISC_COMBINATIONS_PER_TASK must be a positive integer." >&2
        exit 1
    fi

    if ! [[ "${MAX_CONCURRENT}" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: CISC_MAX_CONCURRENT must be a positive integer." >&2
        exit 1
    fi
}

run_python() {
    python "${PYTHON_SCRIPT}"
}

read_combination_counts() {
    read -r HEAT_COMBINATIONS PRECIP_COMBINATIONS < <(
        python - "${PYTHON_SCRIPT}" <<'PYTHON'
import runpy
import sys

namespace = runpy.run_path(sys.argv[1])

heat = tuple(namespace["HEAT_WARNINGS"])
precipitation = tuple(namespace["PRECIPITATION_WARNINGS"])

if not heat:
    raise RuntimeError("HEAT_WARNINGS is empty.")

if not precipitation:
    raise RuntimeError("PRECIPITATION_WARNINGS is empty.")

print(
    (2 ** len(heat)) - 1,
    (2 ** len(precipitation)) - 1,
)
PYTHON
    )

    if [[ -z "${HEAT_COMBINATIONS:-}" || -z "${PRECIP_COMBINATIONS:-}" ]]; then
        echo "ERROR: failed to derive warning-combination counts from Python." >&2
        exit 1
    fi
}

load_modules
ensure_r_packages
validate_inputs

export CISC_ROOT="${ROOT}"
export CISC_COMBINATIONS_PER_TASK="${COMBINATIONS_PER_TASK}"

case "${ACTION}" in
    check_environment)
        read_combination_counts
        echo "Heat combinations: ${HEAT_COMBINATIONS}"
        echo "Precipitation combinations: ${PRECIP_COMBINATIONS}"
        echo "Environment check complete."
        ;;

    prepare_shared)
        export CISC_WWA_COMBO_MODE="prepare_shared"
        run_python
        ;;

    evaluate_heat)
        export CISC_WWA_COMBO_MODE="evaluate_worker"
        export CISC_COMBINATION_VARIABLE="heat"
        run_python
        ;;

    evaluate_precipitation)
        export CISC_WWA_COMBO_MODE="evaluate_worker"
        export CISC_COMBINATION_VARIABLE="precipitation"
        run_python
        ;;

    summarize_all)
        export CISC_WWA_COMBO_MODE="summarize_all"
        run_python
        ;;

    submit_all)
        read_combination_counts

        HEAT_TASKS=$(( (HEAT_COMBINATIONS + COMBINATIONS_PER_TASK - 1) / COMBINATIONS_PER_TASK ))
        PRECIP_TASKS=$(( (PRECIP_COMBINATIONS + COMBINATIONS_PER_TASK - 1) / COMBINATIONS_PER_TASK ))

        echo "Heat combinations: ${HEAT_COMBINATIONS}; array tasks: ${HEAT_TASKS}"
        echo "Precipitation combinations: ${PRECIP_COMBINATIONS}; array tasks: ${PRECIP_TASKS}"

        PREP_JOB=$(sbatch --parsable "${BASH_SOURCE[0]}" prepare_shared)
        echo "Submitted shared preparation: ${PREP_JOB}"

        HEAT_JOB=$(sbatch --parsable \
            --dependency="afterok:${PREP_JOB}" \
            --array="0-$((HEAT_TASKS - 1))%${MAX_CONCURRENT}" \
            --export="ALL,CISC_COMBINATIONS_PER_TASK=${COMBINATIONS_PER_TASK}" \
            "${BASH_SOURCE[0]}" evaluate_heat)
        echo "Submitted heat combinations: ${HEAT_JOB}"

        PRECIP_JOB=$(sbatch --parsable \
            --dependency="afterok:${PREP_JOB}" \
            --array="0-$((PRECIP_TASKS - 1))%${MAX_CONCURRENT}" \
            --export="ALL,CISC_COMBINATIONS_PER_TASK=${COMBINATIONS_PER_TASK}" \
            "${BASH_SOURCE[0]}" evaluate_precipitation)
        echo "Submitted precipitation combinations: ${PRECIP_JOB}"

        SUMMARY_JOB=$(sbatch --parsable \
            --dependency="afterok:${HEAT_JOB}:${PRECIP_JOB}" \
            "${BASH_SOURCE[0]}" summarize_all)
        echo "Submitted final summary: ${SUMMARY_JOB}"
        ;;

    *)
        echo "Usage: $0 {check_environment|submit_all|prepare_shared|evaluate_heat|evaluate_precipitation|summarize_all}" >&2
        exit 2
        ;;
esac