# scripts/summarise_confusion_matrices_by_variable_rp.R
# =============================================================================
# Read county-level performance CSVs for heat index and precipitation,
# then sum TP, FP, FN, TN across rows by return period and variable.
#
# Input:
#   V:/users/hquintal/phd1_cluster_southeast_v2/data/out/09_performance/metrics_county/
#     - metrics_heatindex_all_heat_multiwindow.csv
#     - metrics_precipitation_all_precip_multiwindow.csv
#
# Output:
#   V:/users/hquintal/phd1_cluster_southeast_v2/data/out/09_performance/confusion_matrices.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
in_dir <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/09_performance/metrics_county"
out_file <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/09_performance/confusion_matrices.csv"

files <- c(
  file.path(in_dir, "metrics_heatindex_all_heat_multiwindow.csv"),
  file.path(in_dir, "metrics_precipitation_all_precip_multiwindow.csv")
)

# -----------------------------------------------------------------------------
# Helper: extract return period from scenario string like hi10_pr10
# -----------------------------------------------------------------------------
extract_rp <- function(x, variable_name) {
  if (tolower(variable_name) == "heatindex") {
    return(as.integer(str_extract(x, "(?<=hi)\\d+")))
  }
  if (tolower(variable_name) == "precipitation") {
    return(as.integer(str_extract(x, "(?<=pr)\\d+")))
  }
  return(NA_integer_)
}

# -----------------------------------------------------------------------------
# Read + combine
# -----------------------------------------------------------------------------
dt_list <- lapply(files, fread)
dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# -----------------------------------------------------------------------------
# Filter to full-period window only
# -----------------------------------------------------------------------------
dt <- dt[window_tag == "20020101_20231231"]

# -----------------------------------------------------------------------------
# Standardize variable names
# -----------------------------------------------------------------------------
dt[, var_group := tolower(trimws(as.character(var_group)))]
dt[var_group %in% c("heat index", "heat_index"), var_group := "heatindex"]
dt[var_group %in% c("precip", "precipitation"), var_group := "precipitation"]

dt[, scenario := as.character(scenario)]
dt[, return_period := mapply(extract_rp, scenario, var_group)]

# Keep only requested return periods and variables
dt <- dt[
  var_group %in% c("heatindex", "precipitation") &
    return_period %in% c(2, 5, 10, 25)
]

# -----------------------------------------------------------------------------
# Sum confusion matrix elements across rows
# -----------------------------------------------------------------------------
conf_mat <- dt[
  ,
  .(
    TP = sum(as.numeric(TP), na.rm = TRUE),
    FP = sum(as.numeric(FP), na.rm = TRUE),
    FN = sum(as.numeric(FN), na.rm = TRUE),
    TN = sum(as.numeric(TN), na.rm = TRUE),
    N_rows = .N
  ),
  by = .(variable = var_group, dataset_tag, return_period)
]

# Optional: add totals/derived rates from the summed confusion matrix
conf_mat[
  ,
  `:=`(
    total = TP + FP + FN + TN,
    precision = fifelse((TP + FP) > 0, TP / (TP + FP), NA_real_),
    recall = fifelse((TP + FN) > 0, TP / (TP + FN), NA_real_),
    specificity = fifelse((TN + FP) > 0, TN / (TN + FP), NA_real_),
    accuracy = fifelse((TP + TN + FP + FN) > 0, (TP + TN) / (TP + TN + FP + FN), NA_real_),
    f1 = fifelse((2 * TP + FP + FN) > 0, (2 * TP) / (2 * TP + FP + FN), NA_real_),
    f2 = fifelse((5 * TP + 4 * FN + FP) > 0, (5 * TP) / (5 * TP + 4 * FN + FP), NA_real_)
  )
]

setorder(conf_mat, variable, dataset_tag, return_period)

# -----------------------------------------------------------------------------
# Write output
# -----------------------------------------------------------------------------
fwrite(conf_mat, out_file)

cat("✅ Wrote confusion matrix summary to:\n", out_file, "\n\n")
print(conf_mat)