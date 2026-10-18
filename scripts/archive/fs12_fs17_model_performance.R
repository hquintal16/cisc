#!/usr/bin/env Rscript

# scripts/fig_05_plot_f1_native_vs_resampled_16panel_bothvars.R
# =============================================================================
# County-level native-versus-resampled model-performance comparisons using
# outputs from scripts/13_14_wwa_warning_combination_evaluation.py.
#
# Six 16-panel figures are written:
#   1. Heat-index native versus resampled F1
#   2. Precipitation native versus resampled F1
#   3. Heat-index recall difference: resampled - native
#   4. Precipitation recall difference: resampled - native
#   5. Heat-index precision difference: resampled - native
#   6. Precipitation precision difference: resampled - native
#
# Facets:
#   rows    : return period (5, 10, 25, 50 year)
#   columns : STM1/STM2 x bilinear/max resampling
#
# Fill:
#   county count in each two-dimensional bin
#   heat-index figures use RColorBrewer Reds
#   precipitation figures use RColorBrewer Blues
#
# Default warning combinations:
#   heat          : EH_W
#   precipitation : FF_W
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
  library(scales)
})

# =============================================================================
# Configuration
# =============================================================================

ROOT <- Sys.getenv(
  "CISC_ROOT",
  "/proj/sebastianlab/users/hquintal/cisc"
)

PERF_BASE <- Sys.getenv(
  "CISC_PERF_DIR",
  file.path(ROOT, "data", "out", "20260802_model_performance")
)

HEAT_COMBINATION <- Sys.getenv(
  "CISC_HEAT_COMBINATION",
  "EH_W"
)

PRECIPITATION_COMBINATION <- Sys.getenv(
  "CISC_PRECIPITATION_COMBINATION",
  "FF_W"
)

WINDOW_TAG <- "20020101_20231231"
RETURN_PERIODS <- c(5L, 10L, 25L, 50L)

OUT_DIR <- file.path(ROOT, "figures/fs12_fs17_performance_1to1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_BINS_XY <- 45L
N_FILL_BINS <- 5L

# =============================================================================
# Helpers
# =============================================================================

.log <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf(...),
    "\n",
    sep = ""
  )
}

.safe_file_tag <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", x)
  out <- gsub("^_+|_+$", "", out)
  tolower(out)
}

.theme_panel <- function(base_size = 20) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial", size = 20),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(
        fill = "white",
        color = "grey35",
        linewidth = 0.7
      ),
      strip.text.x = element_text(
        face = "bold",
        size = 20,
        margin = margin(3, 3, 3, 3)
      ),
      strip.text.y = element_text(
        face = "bold",
        size = 20,
        angle = -90,
        margin = margin(3, 3, 3, 3)
      ),
      axis.title = element_text(size = 20, family = "Arial"),
      axis.text = element_text(size = 12, family = "Arial"),
      axis.ticks.length = unit(1.5, "mm"),
      legend.title = element_text(size = 20, family = "Arial"),
      legend.text = element_text(size = 13, family = "Arial"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(6, 6, 6, 6),
      panel.spacing.x = unit(0, "mm"),
      panel.spacing.y = unit(0, "mm")
    )
}

.read_metrics <- function(variable_dir, combination_tag) {
  input_dir <- file.path(PERF_BASE, variable_dir, combination_tag)
  csv_file <- file.path(input_dir, "metrics_county.csv")
  parquet_file <- file.path(input_dir, "metrics_county.parquet")
  
  if (file.exists(csv_file)) {
    .log("Reading %s", csv_file)
    return(fread(csv_file))
  }
  
  if (file.exists(parquet_file)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop(
        "Found only the Parquet input but package 'arrow' is unavailable: ",
        parquet_file
      )
    }
    .log("Reading %s", parquet_file)
    return(as.data.table(arrow::read_parquet(parquet_file)))
  }
  
  stop(
    "Missing metrics_county.csv and metrics_county.parquet under: ",
    input_dir
  )
}

.validate_metrics <- function(dt, expected_var_group, expected_combination) {
  required <- c(
    "var_group", "combination_tag", "dataset_tag", "return_period",
    "window_tag", "GEOID", "f1", "recall", "precision"
  )
  
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0L) {
    stop(
      "Metrics table is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
  
  character_columns <- intersect(
    c(
      "variable", "var_group", "combination_tag", "warnings",
      "dataset_tag", "window_tag", "GEOID"
    ),
    names(dt)
  )
  
  for (column_name in character_columns) {
    set(dt, j = column_name, value = as.character(dt[[column_name]]))
  }
  
  set(dt, j = "var_group", value = tolower(dt[["var_group"]]))
  set(dt, j = "return_period", value = as.integer(dt[["return_period"]]))
  
  for (metric_name in c("f1", "recall", "precision")) {
    set(dt, j = metric_name, value = as.numeric(dt[[metric_name]]))
  }
  
  # Deliberately use ordinary vectors rather than data.table NSE here. This
  # prevents a column named combination_tag from masking a function argument
  # or triggering stale parent-scope lookup errors.
  keep <-
    dt[["var_group"]] == tolower(as.character(expected_var_group)) &
    dt[["combination_tag"]] == as.character(expected_combination) &
    dt[["window_tag"]] == as.character(WINDOW_TAG) &
    dt[["return_period"]] %in% RETURN_PERIODS
  
  keep[is.na(keep)] <- FALSE
  out <- dt[which(keep)]
  
  if (nrow(out) == 0L) {
    stop(
      "No rows remain for var_group=", expected_var_group,
      ", combination_tag=", expected_combination,
      ", window_tag=", WINDOW_TAG,
      ", and return periods ",
      paste(RETURN_PERIODS, collapse = ", ")
    )
  }
  
  observed_combinations <- unique(out[["combination_tag"]])
  if (length(observed_combinations) != 1L ||
      !identical(observed_combinations, as.character(expected_combination))) {
    stop(
      "Unexpected combination_tag values after filtering: ",
      paste(observed_combinations, collapse = ", ")
    )
  }
  
  out[]
}

.add_labels <- function(dt, var_group_value) {
  out <- copy(dt)
  
  if (var_group_value == "heatindex") {
    native_tag <- "heatindex_0.25"
    brewer_palette <- "Reds"
    strategy_lookup <- c(
      "heatindex_0.30_bilinear" = "STM1 + Bilinear",
      "heatindex_0.30_max" = "STM1 + Max",
      "heatindex_0.36_bilinear" = "STM2 + Bilinear",
      "heatindex_0.36_max" = "STM2 + Max"
    )
  } else if (var_group_value == "precipitation") {
    native_tag <- "precipitation_0.25"
    brewer_palette <- "Blues"
    strategy_lookup <- c(
      "precipitation_0.11_bilinear" = "STM1 + Bilinear",
      "precipitation_0.11_max" = "STM1 + Max",
      "precipitation_0.14_bilinear" = "STM2 + Bilinear",
      "precipitation_0.14_max" = "STM2 + Max"
    )
  } else {
    stop("Unsupported var_group: ", var_group_value)
  }
  
  out[, strategy_lab := unname(strategy_lookup[dataset_tag])]
  out[, rp_lab := sprintf("%d-year", return_period)]
  
  list(
    data = out,
    native_tag = native_tag,
    brewer_palette = brewer_palette,
    strategy_levels = unname(strategy_lookup)
  )
}

.make_fill_colors <- function(
    palette_name,
    n_colors = 5L,
    start_fraction = 0.50
) {
  full_ramp <- colorRampPalette(
    brewer.pal(9L, palette_name),
    space = "Lab"
  )(256L)
  
  start_index <- max(
    1L,
    floor(start_fraction * length(full_ramp))
  )
  
  selected_indices <- round(
    seq(
      from = start_index,
      to = length(full_ramp),
      length.out = n_colors
    )
  )
  
  full_ramp[selected_indices]
}

.nice_legend_step <- function(max_count, n_bins = 10L) {
  raw_step <- max_count / n_bins
  
  if (!is.finite(raw_step) || raw_step <= 0) {
    return(1)
  }
  
  magnitude <- 10 ^ floor(log10(raw_step))
  normalized <- raw_step / magnitude
  candidates <- c(1, 2, 2.5, 5, 7.5, 10)
  selected <- candidates[which(candidates >= normalized)[1L]]
  
  if (is.na(selected)) {
    selected <- 10
  }
  
  selected * magnitude
}

.make_legend_spec <- function(max_count, n_bins = 10L) {
  if (!is.finite(max_count) || max_count < 1) {
    stop("max_count must be at least 1")
  }
  
  step <- .nice_legend_step(max_count, n_bins)
  upper <- step * n_bins
  boundaries <- seq(step, upper, by = step)
  
  list(
    step = step,
    upper = upper,
    breaks = boundaries,
    labels = comma(boundaries, accuracy = 1)
  )
}

.build_bin2d_counts <- function(
    plot_dt,
    x_column,
    y_column,
    x_limits,
    y_limits,
    n_bins_xy = 45L
) {
  dt <- copy(plot_dt)
  
  x_width <- diff(x_limits) / n_bins_xy
  y_width <- diff(y_limits) / n_bins_xy
  
  dt[, x_value := pmin(pmax(get(x_column), x_limits[1L]), x_limits[2L])]
  dt[, y_value := pmin(pmax(get(y_column), y_limits[1L]), y_limits[2L])]
  
  dt[, xbin := pmin(
    floor((x_value - x_limits[1L]) / x_width),
    n_bins_xy - 1L
  )]
  
  dt[, ybin := pmin(
    floor((y_value - y_limits[1L]) / y_width),
    n_bins_xy - 1L
  )]
  
  output <- dt[
    ,
    .(count = .N),
    by = .(rp_lab, strategy_lab, xbin, ybin)
  ]
  
  output[, x := x_limits[1L] + (xbin + 0.5) * x_width]
  output[, y := y_limits[1L] + (ybin + 0.5) * y_width]
  output[, `:=`(tile_width = x_width, tile_height = y_width)]
  output[]
}

.prepare_paired_data <- function(
    variable_dir,
    var_group_value,
    combination_tag,
    metric_name
) {
  dt <- .read_metrics(variable_dir, combination_tag)
  dt <- .validate_metrics(dt, var_group_value, combination_tag)
  labels <- .add_labels(dt, var_group_value)
  
  dt <- labels$data
  native_tag <- labels$native_tag
  strategy_levels <- labels$strategy_levels
  rp_levels <- sprintf("%d-year", RETURN_PERIODS)
  
  dt[, rp_lab := factor(rp_lab, levels = rp_levels)]
  dt[, strategy_lab := factor(strategy_lab, levels = strategy_levels)]
  
  native_dt <- dt[
    dataset_tag == native_tag,
    .(
      return_period,
      window_tag,
      GEOID,
      rp_lab,
      native_value = get(metric_name),
      native_f1 = f1
    )
  ]
  
  resampled_dt <- dt[
    dataset_tag != native_tag & !is.na(strategy_lab),
    .(
      return_period,
      window_tag,
      GEOID,
      dataset_tag,
      strategy_lab,
      rp_lab,
      resampled_value = get(metric_name)
    )
  ]
  
  duplicate_native <- native_dt[
    ,
    .N,
    by = .(return_period, window_tag, GEOID)
  ][N > 1L]
  
  if (nrow(duplicate_native) > 0L) {
    stop("Native metrics contain duplicate county/return-period rows")
  }
  
  duplicate_resampled <- resampled_dt[
    ,
    .N,
    by = .(return_period, window_tag, GEOID, dataset_tag)
  ][N > 1L]
  
  if (nrow(duplicate_resampled) > 0L) {
    stop("Resampled metrics contain duplicate county/dataset/return-period rows")
  }
  
  plot_dt <- merge(
    x = resampled_dt,
    y = native_dt,
    by = c("return_period", "window_tag", "GEOID", "rp_lab"),
    all = FALSE,
    sort = FALSE
  )
  
  plot_dt <- plot_dt[
    is.finite(native_value) &
      is.finite(resampled_value) &
      between(native_value, 0, 1) &
      between(resampled_value, 0, 1)
  ]
  
  if (nrow(plot_dt) == 0L) {
    stop(
      "No finite paired native/resampled county values remain for ",
      var_group_value,
      " / ",
      metric_name
    )
  }
  
  plot_dt[, difference := resampled_value - native_value]
  
  expected_panels <- CJ(
    rp_lab = factor(rp_levels, levels = rp_levels),
    strategy_lab = factor(strategy_levels, levels = strategy_levels),
    unique = TRUE
  )
  
  missing_panels <- fsetdiff(
    expected_panels,
    unique(plot_dt[, .(rp_lab, strategy_lab)])
  )
  
  if (nrow(missing_panels) > 0L) {
    stop(
      "Missing paired data for these panels:\n",
      paste(capture.output(print(missing_panels)), collapse = "\n")
    )
  }
  
  list(
    data = plot_dt,
    brewer_palette = labels$brewer_palette,
    strategy_levels = strategy_levels
  )
}

# =============================================================================
# Figure builder
# =============================================================================

.build_metric_figure <- function(
    variable_dir,
    var_group_value,
    combination_tag,
    metric_name,
    output_stem,
    plot_type = c("comparison", "difference")
) {
  plot_type <- match.arg(plot_type)
  
  prepared <- .prepare_paired_data(
    variable_dir = variable_dir,
    var_group_value = var_group_value,
    combination_tag = combination_tag,
    metric_name = metric_name
  )
  
  plot_dt <- prepared$data
  metric_label <- switch(
    metric_name,
    f1 = "F1",
    recall = "Recall",
    precision = "Precision",
    stop("Unsupported metric: ", metric_name)
  )
  
  if (plot_type == "comparison") {
    x_column <- "native_value"
    y_column <- "resampled_value"
    x_limits <- c(0, 1)
    y_limits <- c(0, 1)
    x_breaks <- seq(0, 1, by = 0.2)
    y_breaks <- seq(0, 1, by = 0.2)
    x_labels <- number(x_breaks, accuracy = 0.1)
    y_labels <- number(y_breaks, accuracy = 0.1)
    x_title <- sprintf("Native %s", metric_label)
    y_title <- sprintf("Resampled %s", metric_label)
  } else {
    x_column <- "native_value"
    y_column <- "difference"
    x_limits <- c(0, 1)
    y_limits <- c(-1, 1)
    x_breaks <- seq(0, 1, by = 0.2)
    y_breaks <- seq(-1, 1, by = 0.2)
    x_labels <- number(x_breaks, accuracy = 0.1)
    y_labels <- number(y_breaks, accuracy = 0.1)
    x_title <- sprintf("Native %s", metric_label)
    y_title <- sprintf("Difference in %s (resampled - native)", metric_label)
  }
  
  tile_dt <- .build_bin2d_counts(
    plot_dt = plot_dt,
    x_column = x_column,
    y_column = y_column,
    x_limits = x_limits,
    y_limits = y_limits,
    n_bins_xy = N_BINS_XY
  )
  
  max_count <- max(tile_dt$count, na.rm = TRUE)
  legend_spec <- .make_legend_spec(max_count, N_FILL_BINS)
  fill_colors <- .make_fill_colors(
    prepared$brewer_palette,
    N_FILL_BINS
  )
  
  panel_stats <- plot_dt[
    ,
    .(
      mean_native_value = mean(native_value),
      mean_resampled_value = mean(resampled_value),
      mean_difference = mean(difference),
      mean_native_f1 = mean(native_f1, na.rm = TRUE),
      n_counties = .N
    ),
    by = .(rp_lab, strategy_lab)
  ]
  
  if (plot_type == "comparison") {
    panel_stats[
      ,
      `:=`(
        top_left_label = sprintf(
          "Resampled %s = %.3f\nDifference = %+.3f\nn = %s",
          metric_label,
          mean_resampled_value,
          mean_difference,
          comma(n_counties)
        ),
        bottom_right_label = sprintf(
          "Native F1 = %.3f",
          mean_native_f1
        )
      )
    ]
  } else {
    panel_stats[
      ,
      `:=`(
        top_left_label = sprintf(
          "Mean difference = %+.3f\nn = %s",
          mean_difference,
          comma(n_counties)
        ),
        bottom_right_label = sprintf(
          "Native %s = %.3f",
          metric_label,
          mean_native_value
        )
      )
    ]
  }
  
  panel_stats[, `:=`(
    x_left = x_limits[1L] + 0.03 * diff(x_limits),
    y_top = y_limits[2L] - 0.03 * diff(y_limits),
    x_right = x_limits[2L] - 0.03 * diff(x_limits),
    y_bottom = y_limits[1L] + 0.03 * diff(y_limits)
  )]
  
  p <- ggplot(tile_dt, aes(x = x, y = y, fill = count)) +
    geom_tile(
      aes(width = tile_width, height = tile_height)
    )
  
  if (plot_type == "comparison") {
    p <- p + geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    )
  } else {
    p <- p + geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    )
  }
  
  p <- p +
    geom_text(
      data = panel_stats,
      aes(x = x_left, y = y_top, label = top_left_label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 4.5,
      lineheight = 0.95,
      family = "Arial"
    ) +
    geom_text(
      data = panel_stats,
      aes(x = x_right, y = y_bottom, label = bottom_right_label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0,
      size = 4.5,
      family = "Arial"
    ) +
    facet_grid(
      rows = vars(rp_lab),
      cols = vars(strategy_lab),
      drop = FALSE
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = x_breaks,
      labels = x_labels,
      expand = c(0, 0),
      sec.axis = dup_axis(
        breaks = x_breaks,
        labels = x_labels,
        name = NULL
      )
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      labels = y_labels,
      expand = c(0, 0),
      sec.axis = dup_axis(
        breaks = y_breaks,
        labels = y_labels,
        name = NULL
      )
    ) +
    scale_fill_stepsn(
      colours = fill_colors,
      limits = c(1, legend_spec$upper),
      breaks = legend_spec$breaks,
      labels = legend_spec$labels,
      na.value = fill_colors[1L],
      oob = scales::squish,
      name = "County count",
      guide = guide_colorsteps(
        title.position = "top",
        title.hjust = 0.5,
        direction = "horizontal",
        barwidth = unit(8.5, "in"),
        barheight = unit(0.18, "in"),
        even.steps = TRUE,
        show.limits = TRUE
      )
    ) +
    labs(x = x_title, y = y_title) +
    coord_cartesian(clip = "off") +
    .theme_panel() +
    theme(
      axis.text.x = element_text(size = 12, family = "Arial"),
      axis.text.y = element_text(size = 12, family = "Arial"),
      axis.text.x.top = element_text(size = 12, family = "Arial"),
      axis.text.y.right = element_text(size = 12, family = "Arial"),
      axis.ticks.x.top = element_blank(),
      axis.ticks.y.right = element_blank()
    )
  
  output_file <- file.path(
    OUT_DIR,
    paste0(output_stem, ".png")
  )
  
  ggsave(
    filename = output_file,
    plot = p,
    width = 12.8,
    height = 13.4,
    units = "in",
    dpi = 200,
    bg = "white",
    limitsize = FALSE
  )
  
  .log(
    "%s %s %s: paired %s county-dataset rows; max bin count=%s; legend step=%s",
    var_group_value,
    metric_name,
    plot_type,
    comma(nrow(plot_dt)),
    comma(max_count),
    comma(legend_spec$step)
  )
  .log("Wrote %s", output_file)
  
  invisible(output_file)
}

# =============================================================================
# Build all figures
# =============================================================================

outputs <- c(
  .build_metric_figure(
    variable_dir = "heat",
    var_group_value = "heatindex",
    combination_tag = HEAT_COMBINATION,
    metric_name = "f1",
    output_stem = "fs12_1to1_heat_f1",
    plot_type = "comparison"
  ),
  .build_metric_figure(
    variable_dir = "precipitation",
    var_group_value = "precipitation",
    combination_tag = PRECIPITATION_COMBINATION,
    metric_name = "f1",
    output_stem = "fs13_1to1_precip_f1",
    plot_type = "comparison"
  ),
  .build_metric_figure(
    variable_dir = "heat",
    var_group_value = "heatindex",
    combination_tag = HEAT_COMBINATION,
    metric_name = "recall",
    output_stem = "fs14_1to1_heat_recall_difference",
    plot_type = "difference"
  ),
  .build_metric_figure(
    variable_dir = "precipitation",
    var_group_value = "precipitation",
    combination_tag = PRECIPITATION_COMBINATION,
    metric_name = "recall",
    output_stem = "fs15_1to1_precip_recall_difference",
    plot_type = "difference"
  ),
  .build_metric_figure(
    variable_dir = "heat",
    var_group_value = "heatindex",
    combination_tag = HEAT_COMBINATION,
    metric_name = "precision",
    output_stem = "fs16_1to1_heat_precision_difference",
    plot_type = "difference"
  ),
  .build_metric_figure(
    variable_dir = "precipitation",
    var_group_value = "precipitation",
    combination_tag = PRECIPITATION_COMBINATION,
    metric_name = "precision",
    output_stem = "fs17_1to1_precip_precision_difference",
    plot_type = "difference"
  )
)

cat("\nSaved figures:\n")
cat(outputs, sep = "\n")
cat("\n")