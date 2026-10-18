# # =============================================================================
# # Summarize county-level F1 for a targeted subset only
# #
# # FILTERS:
# #   - window_tag must equal "20190101_20231231"
# #   - heat index:    scenario == "hi2_pr2" and dataset_tag starts with "heatindex"
# #   - precipitation: scenario == "hi5_pr5" and dataset_tag starts with "precipitation"
# #
# # OUTPUT:
# #   - one row per dataset_tag
# #   - mean and median F1 across counties
# #
# # Writes:
# #   - summary_f1_heatindex_hi2_pr2_20190101_20231231.csv
# #   - summary_f1_precipitation_hi5_pr5_20190101_20231231.csv
# # =============================================================================
# 
# suppressPackageStartupMessages({
#   library(data.table)
#   library(fs)
#   library(stringr)
# })
# 
# # -----------------------------------------------------------------------------
# # Paths
# # -----------------------------------------------------------------------------
# ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
# IN_DIR <- file.path(ROOT, "data/out/09_performance/metrics_county")
# 
# HEAT_FILE <- file.path(IN_DIR, "metrics_heatindex_all_heat_multiwindow.csv")
# PREC_FILE <- file.path(IN_DIR, "metrics_precipitation_all_precip_multiwindow.csv")
# 
# TARGET_WINDOW <- "20190101_20231231"
# HEAT_SCENARIO <- "hi25_pr25"
# PREC_SCENARIO <- "hi25_pr25"
# 
# OUT_HEAT <- file.path(
#   IN_DIR,
#   paste0("summary_f1_heatindex_", HEAT_SCENARIO, "_", TARGET_WINDOW, ".csv")
# )
# 
# OUT_PREC <- file.path(
#   IN_DIR,
#   paste0("summary_f1_precipitation_", PREC_SCENARIO, "_", TARGET_WINDOW, ".csv")
# )
# 
# # -----------------------------------------------------------------------------
# # Helper: read + validate
# # -----------------------------------------------------------------------------
# .read_metrics <- function(f) {
#   if (!file_exists(f)) stop("Missing file: ", f)
#   
#   dt <- fread(f)
#   
#   required_cols <- c("scenario", "dataset_tag", "window_tag", "f1")
#   missing_cols <- setdiff(required_cols, names(dt))
#   if (length(missing_cols) > 0) {
#     stop(
#       "File is missing required columns: ", basename(f), " -> ",
#       paste(missing_cols, collapse = ", ")
#     )
#   }
#   
#   dt[, f1 := as.numeric(f1)]
#   dt[]
# }
# 
# # -----------------------------------------------------------------------------
# # Helper: summarize
# # -----------------------------------------------------------------------------
# .summarize_f1 <- function(dt) {
#   out <- dt[
#     ,
#     .(
#       n_counties = sum(!is.na(f1)),
#       mean_f1    = mean(f1, na.rm = TRUE),
#       median_f1  = median(f1, na.rm = TRUE),
#       sd_f1      = sd(f1, na.rm = TRUE),
#       min_f1     = suppressWarnings(min(f1, na.rm = TRUE)),
#       max_f1     = suppressWarnings(max(f1, na.rm = TRUE))
#     ),
#     by = .(scenario, dataset_tag, window_tag)
#   ]
#   
#   out[!is.finite(mean_f1), `:=`(
#     mean_f1   = NA_real_,
#     median_f1 = NA_real_,
#     sd_f1     = NA_real_,
#     min_f1    = NA_real_,
#     max_f1    = NA_real_
#   )]
#   
#   setorder(out, dataset_tag)
#   out[]
# }
# 
# # -----------------------------------------------------------------------------
# # Read
# # -----------------------------------------------------------------------------
# heat_dt <- .read_metrics(HEAT_FILE)
# prec_dt <- .read_metrics(PREC_FILE)
# 
# # -----------------------------------------------------------------------------
# # Filter
# # -----------------------------------------------------------------------------
# heat_sub <- heat_dt[
#   window_tag == TARGET_WINDOW &
#     scenario == HEAT_SCENARIO &
#     stringr::str_starts(dataset_tag, "heatindex")
# ]
# 
# prec_sub <- prec_dt[
#   window_tag == TARGET_WINDOW &
#     scenario == PREC_SCENARIO &
#     stringr::str_starts(dataset_tag, "precipitation")
# ]
# 
# # Optional safety checks
# if (nrow(heat_sub) == 0) {
#   warning("No heat index rows matched the requested filters.")
# }
# if (nrow(prec_sub) == 0) {
#   warning("No precipitation rows matched the requested filters.")
# }
# 
# # -----------------------------------------------------------------------------
# # Summarize
# # -----------------------------------------------------------------------------
# heat_summary <- .summarize_f1(heat_sub)
# prec_summary <- .summarize_f1(prec_sub)
# 
# # -----------------------------------------------------------------------------
# # Write
# # -----------------------------------------------------------------------------
# fwrite(heat_summary, OUT_HEAT)
# fwrite(prec_summary, OUT_PREC)
# 
# # -----------------------------------------------------------------------------
# # Report
# # -----------------------------------------------------------------------------
# cat("\nDone.\n")
# cat("Target window: ", TARGET_WINDOW, "\n", sep = "")
# cat("Heat scenario: ", HEAT_SCENARIO, "\n", sep = "")
# cat("Precip scenario: ", PREC_SCENARIO, "\n\n", sep = "")
# 
# cat("Heat index summary rows: ", nrow(heat_summary), "\n", sep = "")
# cat("Precipitation summary rows: ", nrow(prec_summary), "\n", sep = "")
# cat("Wrote:\n")
# cat("  ", OUT_HEAT, "\n", sep = "")
# cat("  ", OUT_PREC, "\n\n", sep = "")
# 
# cat("Heat index summary:\n")
# print(heat_summary)
# 
# cat("\nPrecipitation summary:\n")
# print(prec_summary)
# 
# # ANOVA ----
# # =============================================================================
# # Compare F1 distributions against native resolution and plot differences
# #
# # Filters:
# #   - window_tag == "20190101_20231231"
# #   - heat index:    scenario == "hi2_pr2"
# #   - precipitation: scenario == "hi5_pr5"
# #
# # Output tables:
# #   - f1_comparison_vs_native_heatindex.csv
# #   - f1_comparison_vs_native_precipitation.csv
# #
# # Output figure:
# #   - V:/users/hquintal/phd1_cluster_southeast_v2/figures/fig_07_evaluation_stats/
# #       f1_difference_vs_native_heatindex_precipitation.png
# # =============================================================================
# 
# suppressPackageStartupMessages({
#   library(data.table)
#   library(stringr)
#   library(ggplot2)
#   library(ragg)
#   library(patchwork)
# })
# 
# # -----------------------------------------------------------------------------
# # Paths
# # -----------------------------------------------------------------------------
# ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
# 
# IN_DIR <- file.path(ROOT, "data/out/09_performance/metrics_county")
# OUT_DIR <- file.path(ROOT, "figures/fig_07_evaluation_stats")
# dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
# 
# HEAT_FILE <- file.path(IN_DIR, "metrics_heatindex_all_heat_multiwindow.csv")
# PREC_FILE <- file.path(IN_DIR, "metrics_precipitation_all_precip_multiwindow.csv")
# 
# OUT_HEAT_CSV <- file.path(IN_DIR, "f1_comparison_vs_native_heatindex.csv")
# OUT_PREC_CSV <- file.path(IN_DIR, "f1_comparison_vs_native_precipitation.csv")
# OUT_FIG <- file.path(
#   OUT_DIR,
#   "f1_difference_vs_native_heatindex_precipitation.png"
# )
# 
# # -----------------------------------------------------------------------------
# # Filters
# # -----------------------------------------------------------------------------
# TARGET_WINDOW <- "20190101_20231231"
# HEAT_SCENARIO <- "hi2_pr2"
# PREC_SCENARIO <- "hi5_pr5"
# 
# # -----------------------------------------------------------------------------
# # Label mapping
# # -----------------------------------------------------------------------------
# .label_heat <- function(x) {
#   out <- x
#   out[x == "heatindex_0.30_bilinear"] <- "STM1 bilinear"
#   out[x == "heatindex_0.30_max"]      <- "STM1 max"
#   out[x == "heatindex_0.36_bilinear"] <- "STM2 bilinear"
#   out[x == "heatindex_0.36_max"]      <- "STM2 max"
#   out
# }
# 
# .label_prec <- function(x) {
#   out <- x
#   out[x == "precipitation_0.11_bilinear"] <- "STM1 bilinear"
#   out[x == "precipitation_0.11_max"]      <- "STM1 max"
#   out[x == "precipitation_0.14_bilinear"] <- "STM2 bilinear"
#   out[x == "precipitation_0.14_max"]      <- "STM2 max"
#   out
# }
# 
# # -----------------------------------------------------------------------------
# # Helper: format p-values for plotting
# # -----------------------------------------------------------------------------
# .format_p <- function(p) {
#   out <- ifelse(
#     is.na(p), "p = NA",
#     ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, format = "f", digits = 3)))
#   )
#   out
# }
# 
# # -----------------------------------------------------------------------------
# # Function: compare datasets to native
# # -----------------------------------------------------------------------------
# compare_to_native <- function(dt, native_tag, label_fun) {
#   native_vals <- dt[dataset_tag == native_tag, f1]
#   
#   others <- setdiff(sort(unique(dt$dataset_tag)), native_tag)
#   
#   res <- rbindlist(lapply(others, function(tag) {
#     vals <- dt[dataset_tag == tag, f1]
#     
#     t_res <- t.test(vals, native_vals)
#     w_res <- wilcox.test(vals, native_vals)
#     
#     data.table(
#       dataset_tag    = tag,
#       label          = label_fun(tag),
#       mean_dataset   = mean(vals, na.rm = TRUE),
#       median_dataset = median(vals, na.rm = TRUE),
#       mean_native    = mean(native_vals, na.rm = TRUE),
#       median_native  = median(native_vals, na.rm = TRUE),
#       diff_mean      = mean(vals, na.rm = TRUE) - mean(native_vals, na.rm = TRUE),
#       ci_low         = t_res$conf.int[1],
#       ci_high        = t_res$conf.int[2],
#       t_p_value      = t_res$p.value,
#       wilcox_p_value = w_res$p.value,
#       n_dataset      = length(vals),
#       n_native       = length(native_vals)
#     )
#   }))
#   
#   desired_order <- c("STM1 bilinear", "STM1 max", "STM2 bilinear", "STM2 max")
#   
#   res[, label := factor(label, levels = rev(desired_order))]
#   res[, order_id := match(as.character(label), rev(desired_order))]
#   setorder(res, order_id)
#   res[, order_id := NULL]
#   
#   res[, p_label := .format_p(t_p_value)]
#   res[]
# }
# 
# # -----------------------------------------------------------------------------
# # Function: panel plot
# # -----------------------------------------------------------------------------
# plot_diff_panel <- function(df, show_y = TRUE, x_limits = NULL) {
#   
#   x_span <- diff(x_limits)
#   label_x_nudge <- 0.03 * x_span
#   
#   ggplot(df, aes(x = diff_mean, y = label)) +
#     geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
#     geom_errorbarh(
#       aes(xmin = ci_low, xmax = ci_high),
#       height = 0.18,
#       linewidth = 0.6,
#       color = "grey35"
#     ) +
#     geom_point(size = 2.8) +
#     geom_text(
#       aes(x = ci_high + label_x_nudge, label = p_label),
#       hjust = 0,
#       size = 3.5
#     ) +
#     scale_x_continuous(limits = x_limits) +
#     labs(
#       x = "Difference in F1 relative to native grid",
#       y = NULL
#     ) +
#     theme_minimal(base_size = 13) +
#     theme(
#       plot.title = element_blank(),
#       panel.grid.minor = element_blank(),
#       axis.text.y = if (show_y) element_text(color = "black") else element_blank(),
#       axis.ticks.y = if (show_y) element_line(color = "black") else element_blank()
#     )
# }
# 
# # -----------------------------------------------------------------------------
# # Read data
# # -----------------------------------------------------------------------------
# heat_dt <- fread(HEAT_FILE)
# prec_dt <- fread(PREC_FILE)
# 
# heat_dt[, f1 := as.numeric(f1)]
# prec_dt[, f1 := as.numeric(f1)]
# 
# # -----------------------------------------------------------------------------
# # Apply filters
# # -----------------------------------------------------------------------------
# heat_sub <- heat_dt[
#   window_tag == TARGET_WINDOW &
#     scenario == HEAT_SCENARIO &
#     str_starts(dataset_tag, "heatindex") &
#     !is.na(f1)
# ]
# 
# prec_sub <- prec_dt[
#   window_tag == TARGET_WINDOW &
#     scenario == PREC_SCENARIO &
#     str_starts(dataset_tag, "precipitation") &
#     !is.na(f1)
# ]
# 
# # -----------------------------------------------------------------------------
# # Run comparisons
# # -----------------------------------------------------------------------------
# heat_results <- compare_to_native(
#   dt = heat_sub,
#   native_tag = "heatindex_0.25",
#   label_fun = .label_heat
# )
# 
# prec_results <- compare_to_native(
#   dt = prec_sub,
#   native_tag = "precipitation_0.25",
#   label_fun = .label_prec
# )
# 
# # -----------------------------------------------------------------------------
# # Write tables
# # -----------------------------------------------------------------------------
# fwrite(heat_results, OUT_HEAT_CSV)
# fwrite(prec_results, OUT_PREC_CSV)
# 
# # -----------------------------------------------------------------------------
# # Shared x-axis limits, expanded to make room for p-value text
# # -----------------------------------------------------------------------------
# x_min <- min(c(heat_results$ci_low, prec_results$ci_low), na.rm = TRUE)
# x_max <- max(c(heat_results$ci_high, prec_results$ci_high), na.rm = TRUE)
# x_pad <- 0.05 * max(abs(c(x_min, x_max)))
# x_span_core <- (x_max + x_pad) - (x_min - x_pad)
# 
# # extra room on right for p-value labels
# x_limits <- c(
#   x_min - x_pad,
#   x_max + x_pad + 0.30 * x_span_core
# )
# 
# # -----------------------------------------------------------------------------
# # Build panels
# # -----------------------------------------------------------------------------
# p_heat <- plot_diff_panel(
#   heat_results,
#   show_y = TRUE,
#   x_limits = x_limits
# )
# 
# p_prec <- plot_diff_panel(
#   prec_results,
#   show_y = FALSE,
#   x_limits = x_limits
# )
# 
# combined_plot <- p_heat + p_prec + plot_layout(ncol = 2)
# 
# # -----------------------------------------------------------------------------
# # Save figure
# # -----------------------------------------------------------------------------
# agg_png(
#   OUT_FIG,
#   width = 2800,
#   height = 1400,
#   res = 300
# )
# print(combined_plot)
# dev.off()
# 
# # -----------------------------------------------------------------------------
# # Report
# # -----------------------------------------------------------------------------
# cat("\nDone.\n")
# cat("Wrote tables:\n")
# cat("  ", OUT_HEAT_CSV, "\n", sep = "")
# cat("  ", OUT_PREC_CSV, "\n", sep = "")
# cat("Wrote figure:\n")
# cat("  ", OUT_FIG, "\n\n", sep = "")
# 
# cat("Heat index results:\n")
# print(heat_results)
# 
# cat("\nPrecipitation results:\n")
# print(prec_results)

# =============================================================================
# Summarize county-level F1 and compare resampled vs native across:
#   - year ranges: 2002-2023, 2014-2023, 2019-2023
#   - return periods: 2, 5, 10, 25
# for both heat index and precipitation
#
# OUTPUT TABLES
#   1) summary_f1_heatindex_all_windows_all_rps.csv
#   2) summary_f1_precipitation_all_windows_all_rps.csv
#   3) f1_comparison_vs_native_heatindex_all_windows_all_rps.csv
#   4) f1_comparison_vs_native_precipitation_all_windows_all_rps.csv
#
# NOTES
#   - Scenario is assumed to match the RP as: hi<rp>_pr<rp>
#     e.g., rp=2  -> hi2_pr2
#           rp=5  -> hi5_pr5
#           rp=10 -> hi10_pr10
#           rp=25 -> hi25_pr25
#   - Native-vs-resampled tests are paired by GEOID
#   - p_label is based on paired t-test p-value
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(stringr)
})

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
IN_DIR <- file.path(ROOT, "data/out/09_performance/metrics_county")

HEAT_FILE <- file.path(IN_DIR, "metrics_heatindex_all_heat_multiwindow.csv")
PREC_FILE <- file.path(IN_DIR, "metrics_precipitation_all_precip_multiwindow.csv")

OUT_HEAT_SUMMARY <- file.path(
  IN_DIR, "summary_f1_heatindex_all_windows_all_rps.csv"
)
OUT_PREC_SUMMARY <- file.path(
  IN_DIR, "summary_f1_precipitation_all_windows_all_rps.csv"
)

OUT_HEAT_COMPARE <- file.path(
  IN_DIR, "f1_comparison_vs_native_heatindex_all_windows_all_rps.csv"
)
OUT_PREC_COMPARE <- file.path(
  IN_DIR, "f1_comparison_vs_native_precipitation_all_windows_all_rps.csv"
)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
WINDOWS <- c(
  "20020101_20231231",
  "20140101_20231231",
  "20190101_20231231"
)

RETURN_PERIODS <- c(2, 5, 10, 25)

scenario_from_rp <- function(rp) {
  paste0("hi", rp, "_pr", rp)
}

# -----------------------------------------------------------------------------
# Label mapping
# -----------------------------------------------------------------------------
.label_heat <- function(x) {
  out <- x
  out[x == "heatindex_0.30_bilinear"] <- "STM1 bilinear"
  out[x == "heatindex_0.30_max"]      <- "STM1 max"
  out[x == "heatindex_0.36_bilinear"] <- "STM2 bilinear"
  out[x == "heatindex_0.36_max"]      <- "STM2 max"
  out
}

.label_prec <- function(x) {
  out <- x
  out[x == "precipitation_0.11_bilinear"] <- "STM1 bilinear"
  out[x == "precipitation_0.11_max"]      <- "STM1 max"
  out[x == "precipitation_0.14_bilinear"] <- "STM2 bilinear"
  out[x == "precipitation_0.14_max"]      <- "STM2 max"
  out
}

# -----------------------------------------------------------------------------
# Helper: format p-values
# -----------------------------------------------------------------------------
.format_p <- function(p) {
  fifelse(
    is.na(p), "p = NA",
    fifelse(p < 0.001, "p < 0.001",
            paste0("p = ", formatC(p, format = "f", digits = 3)))
  )
}

# -----------------------------------------------------------------------------
# Helper: read + validate
# -----------------------------------------------------------------------------
.read_metrics <- function(f) {
  if (!file_exists(f)) stop("Missing file: ", f)
  
  dt <- fread(f)
  
  required_cols <- c("scenario", "dataset_tag", "window_tag", "GEOID", "f1")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop(
      "File is missing required columns: ", basename(f), " -> ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  dt[, `:=`(
    scenario   = as.character(scenario),
    dataset_tag = as.character(dataset_tag),
    window_tag = as.character(window_tag),
    GEOID      = as.character(GEOID),
    f1         = as.numeric(f1)
  )]
  
  dt[]
}

# -----------------------------------------------------------------------------
# Helper: summarize
# -----------------------------------------------------------------------------
.summarize_f1 <- function(dt, variable_name, rp) {
  out <- dt[
    ,
    .(
      variable      = variable_name,
      return_period = rp,
      n_counties    = sum(!is.na(f1)),
      mean_f1       = mean(f1, na.rm = TRUE),
      median_f1     = median(f1, na.rm = TRUE),
      sd_f1         = sd(f1, na.rm = TRUE),
      min_f1        = if (all(is.na(f1))) NA_real_ else min(f1, na.rm = TRUE),
      max_f1        = if (all(is.na(f1))) NA_real_ else max(f1, na.rm = TRUE)
    ),
    by = .(scenario, dataset_tag, window_tag)
  ]
  
  setcolorder(
    out,
    c("variable", "return_period", "scenario", "dataset_tag", "window_tag",
      "n_counties", "mean_f1", "median_f1", "sd_f1", "min_f1", "max_f1")
  )
  
  setorder(out, return_period, window_tag, dataset_tag)
  out[]
}

# -----------------------------------------------------------------------------
# Helper: paired comparison vs native
# -----------------------------------------------------------------------------
compare_to_native <- function(dt, native_tag, label_fun, variable_name, rp) {
  
  native_dt <- dt[
    dataset_tag == native_tag & !is.na(f1),
    .(GEOID, native_f1 = f1)
  ]
  
  if (nrow(native_dt) == 0) {
    return(data.table(
      variable       = character(),
      return_period  = integer(),
      scenario       = character(),
      window_tag     = character(),
      dataset_tag    = character(),
      label          = character(),
      mean_dataset   = numeric(),
      median_dataset = numeric(),
      mean_native    = numeric(),
      median_native  = numeric(),
      diff_mean      = numeric(),
      ci_low         = numeric(),
      ci_high        = numeric(),
      t_p_value      = numeric(),
      wilcox_p_value = numeric(),
      n_dataset      = integer(),
      n_native       = integer(),
      n_paired       = integer(),
      p_label        = character()
    ))
  }
  
  others <- setdiff(sort(unique(dt$dataset_tag)), native_tag)
  
  res_list <- lapply(others, function(tag) {
    
    ds_dt <- dt[
      dataset_tag == tag & !is.na(f1),
      .(GEOID, dataset_f1 = f1)
    ]
    
    paired <- merge(native_dt, ds_dt, by = "GEOID", all = FALSE)
    
    if (nrow(paired) == 0) {
      return(data.table(
        variable       = variable_name,
        return_period  = rp,
        scenario       = unique(dt$scenario)[1],
        window_tag     = unique(dt$window_tag)[1],
        dataset_tag    = tag,
        label          = label_fun(tag),
        mean_dataset   = NA_real_,
        median_dataset = NA_real_,
        mean_native    = NA_real_,
        median_native  = NA_real_,
        diff_mean      = NA_real_,
        ci_low         = NA_real_,
        ci_high        = NA_real_,
        t_p_value      = NA_real_,
        wilcox_p_value = NA_real_,
        n_dataset      = nrow(ds_dt),
        n_native       = nrow(native_dt),
        n_paired       = 0L,
        p_label        = "p = NA"
      ))
    }
    
    t_res <- tryCatch(
      t.test(paired$dataset_f1, paired$native_f1, paired = TRUE),
      error = function(e) NULL
    )
    
    w_res <- tryCatch(
      wilcox.test(paired$dataset_f1, paired$native_f1, paired = TRUE, exact = FALSE),
      error = function(e) NULL
    )
    
    out <- data.table(
      variable       = variable_name,
      return_period  = rp,
      scenario       = unique(dt$scenario)[1],
      window_tag     = unique(dt$window_tag)[1],
      dataset_tag    = tag,
      label          = label_fun(tag),
      mean_dataset   = mean(paired$dataset_f1, na.rm = TRUE),
      median_dataset = median(paired$dataset_f1, na.rm = TRUE),
      mean_native    = mean(paired$native_f1, na.rm = TRUE),
      median_native  = median(paired$native_f1, na.rm = TRUE),
      diff_mean      = mean(paired$dataset_f1 - paired$native_f1, na.rm = TRUE),
      ci_low         = if (!is.null(t_res)) t_res$conf.int[1] else NA_real_,
      ci_high        = if (!is.null(t_res)) t_res$conf.int[2] else NA_real_,
      t_p_value      = if (!is.null(t_res)) t_res$p.value else NA_real_,
      wilcox_p_value = if (!is.null(w_res)) w_res$p.value else NA_real_,
      n_dataset      = nrow(ds_dt),
      n_native       = nrow(native_dt),
      n_paired       = nrow(paired)
    )
    
    out[, p_label := .format_p(t_p_value)]
    out[]
  })
  
  res <- rbindlist(res_list, fill = TRUE)
  
  desired_order <- c("STM1 bilinear", "STM1 max", "STM2 bilinear", "STM2 max")
  res[, label := factor(label, levels = desired_order)]
  setorder(res, window_tag, return_period, label)
  
  res[]
}

# -----------------------------------------------------------------------------
# Read source data
# -----------------------------------------------------------------------------
heat_dt <- .read_metrics(HEAT_FILE)
prec_dt <- .read_metrics(PREC_FILE)

# -----------------------------------------------------------------------------
# Loop over windows and return periods
# -----------------------------------------------------------------------------
heat_summary_list <- list()
prec_summary_list <- list()
heat_compare_list <- list()
prec_compare_list <- list()

idx_hs <- 1L
idx_ps <- 1L
idx_hc <- 1L
idx_pc <- 1L

for (w in WINDOWS) {
  for (rp in RETURN_PERIODS) {
    
    scenario_i <- scenario_from_rp(rp)
    
    cat("\n--------------------------------------------------\n")
    cat("Processing window: ", w, "\n", sep = "")
    cat("Processing RP:     ", rp, "\n", sep = "")
    cat("Scenario:          ", scenario_i, "\n", sep = "")
    
    # ---- Heat subset ----
    heat_sub <- heat_dt[
      window_tag == w &
        scenario == scenario_i &
        str_starts(dataset_tag, "heatindex")
    ]
    
    if (nrow(heat_sub) == 0) {
      warning("No heat index rows matched window=", w, ", scenario=", scenario_i)
    } else {
      heat_summary_list[[idx_hs]] <- .summarize_f1(
        dt = heat_sub,
        variable_name = "heatindex",
        rp = rp
      )
      idx_hs <- idx_hs + 1L
      
      heat_compare_list[[idx_hc]] <- compare_to_native(
        dt = heat_sub,
        native_tag = "heatindex_0.25",
        label_fun = .label_heat,
        variable_name = "heatindex",
        rp = rp
      )
      idx_hc <- idx_hc + 1L
    }
    
    # ---- Precip subset ----
    prec_sub <- prec_dt[
      window_tag == w &
        scenario == scenario_i &
        str_starts(dataset_tag, "precipitation")
    ]
    
    if (nrow(prec_sub) == 0) {
      warning("No precipitation rows matched window=", w, ", scenario=", scenario_i)
    } else {
      prec_summary_list[[idx_ps]] <- .summarize_f1(
        dt = prec_sub,
        variable_name = "precipitation",
        rp = rp
      )
      idx_ps <- idx_ps + 1L
      
      prec_compare_list[[idx_pc]] <- compare_to_native(
        dt = prec_sub,
        native_tag = "precipitation_0.25",
        label_fun = .label_prec,
        variable_name = "precipitation",
        rp = rp
      )
      idx_pc <- idx_pc + 1L
    }
  }
}

# -----------------------------------------------------------------------------
# Combine outputs
# -----------------------------------------------------------------------------
heat_summary <- rbindlist(heat_summary_list, fill = TRUE)
prec_summary <- rbindlist(prec_summary_list, fill = TRUE)
heat_compare <- rbindlist(heat_compare_list, fill = TRUE)
prec_compare <- rbindlist(prec_compare_list, fill = TRUE)

setorder(heat_summary, return_period, window_tag, dataset_tag)
setorder(prec_summary, return_period, window_tag, dataset_tag)

setorder(heat_compare, return_period, window_tag, label)
setorder(prec_compare, return_period, window_tag, label)

# -----------------------------------------------------------------------------
# Write outputs
# -----------------------------------------------------------------------------
fwrite(heat_summary, OUT_HEAT_SUMMARY)
fwrite(prec_summary, OUT_PREC_SUMMARY)
fwrite(heat_compare, OUT_HEAT_COMPARE)
fwrite(prec_compare, OUT_PREC_COMPARE)

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------
cat("\n==================================================\n")
cat("Done.\n")
cat("Wrote summary tables:\n")
cat("  ", OUT_HEAT_SUMMARY, "\n", sep = "")
cat("  ", OUT_PREC_SUMMARY, "\n", sep = "")
cat("Wrote comparison tables:\n")
cat("  ", OUT_HEAT_COMPARE, "\n", sep = "")
cat("  ", OUT_PREC_COMPARE, "\n", sep = "")

cat("\nHeat summary rows: ", nrow(heat_summary), "\n", sep = "")
cat("Precip summary rows: ", nrow(prec_summary), "\n", sep = "")
cat("Heat comparison rows: ", nrow(heat_compare), "\n", sep = "")
cat("Precip comparison rows: ", nrow(prec_compare), "\n", sep = "")

cat("\nHeat comparison preview:\n")
print(heat_compare)

cat("\nPrecipitation comparison preview:\n")
print(prec_compare)
