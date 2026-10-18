#!/usr/bin/env Rscript

# scripts/fs08_heat_cluster_time_series.R
# =============================================================================
# Supplemental Figure FS08
#
# Annual heat-index cluster-characteristic time series for:
#   5-, 10-, 25-, and 50-year return periods
#
# Columns:
#   Native 0.25 degrees
#   STM1 + Bilinear 0.30 degrees
#   STM1 + Max 0.30 degrees
#   STM2 + Bilinear 0.36 degrees
#   STM2 + Max 0.36 degrees
#
# Rows:
#   Peak intensity
#   Duration
#   Frequency
#   Extent
# =============================================================================


suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(stringr)
  library(lubridate)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(ggh4x)
})

.log <- function(fmt, ...) {
  message(
    sprintf(
      "[%s] %s",
      format(Sys.time(), "%H:%M:%S"),
      sprintf(fmt, ...)
    )
  )
}

.get_year_col <- function(dt) {
  candidates <- c(
    "start_time",
    "start_datetime",
    "start_date",
    "end_time",
    "end_datetime",
    "end_date"
  )
  
  hit <- candidates[candidates %in% names(dt)]
  
  if (!length(hit)) {
    stop(
      "Could not identify a cluster time column. Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }
  
  hit[[1L]]
}

.parse_event_year <- function(x) {
  x <- as.character(x)
  
  parsed <- suppressWarnings(
    parse_date_time(
      x,
      orders = c(
        "Ymd HMS",
        "Ymd HM",
        "Ymd",
        "Y-m-d H:M:S",
        "Y-m-d H:M",
        "Y-m-d"
      ),
      tz = "UTC",
      quiet = TRUE
    )
  )
  
  year(parsed)
}

.get_duration_days <- function(dt) {
  if ("duration_days" %in% names(dt)) {
    return(as.numeric(dt$duration_days))
  }
  
  if ("duration_hr" %in% names(dt)) {
    return(as.numeric(dt$duration_hr) / 24)
  }
  
  if ("duration_hours" %in% names(dt)) {
    return(as.numeric(dt$duration_hours) / 24)
  }
  
  stop(
    "Missing duration field. Expected duration_days, duration_hr, ",
    "or duration_hours."
  )
}

.get_cluster_id_col <- function(dt) {
  candidates <- c(
    "cluster_id",
    "cluster",
    "event_id"
  )
  
  hit <- candidates[candidates %in% names(dt)]
  
  if (!length(hit)) {
    stop(
      "Could not identify cluster identifier column. Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }
  
  hit[[1L]]
}

.get_peak_intensity_col <- function(dt) {
  candidates <- c(
    "peak_intensity",
    "peak_intensity_C",
    "peak_intensity_mm"
  )
  
  hit <- candidates[candidates %in% names(dt)]
  
  if (!length(hit)) {
    stop(
      "Could not identify peak-intensity column. Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }
  
  hit[[1L]]
}

.get_extent_col <- function(dt, preferred = "max_t_extent_km2") {
  candidates <- unique(
    c(
      preferred,
      "max_hourly_extent_km2",
      "max_t_extent_km2",
      "total_extent_km2"
    )
  )
  
  hit <- candidates[candidates %in% names(dt)]
  
  if (!length(hit)) {
    stop(
      "Could not identify an extent column. Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }
  
  hit[[1L]]
}

.safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  
  if (!length(x)) {
    return(NA_real_)
  }
  
  mean(x)
}

.read_characteristic_sets <- function(
    root,
    return_periods,
    cluster_sets,
    variable_name
) {
  pieces <- list()
  piece_index <- 0L
  
  for (rp in return_periods) {
    for (cluster_set in cluster_sets) {
      input_csv <- file.path(
        root,
        "data",
        "out",
        "07_characteristics",
        paste0("filtered_rp", rp),
        cluster_set,
        "cluster_characteristics.csv"
      )
      
      if (!file.exists(input_csv)) {
        warning(
          sprintf(
            "Missing input; skipping RP%d %s: %s",
            rp,
            cluster_set,
            input_csv
          ),
          call. = FALSE
        )
        next
      }
      
      .log(
        "Reading RP%d %s: %s",
        rp,
        cluster_set,
        input_csv
      )
      
      dt <- fread(input_csv)
      
      if (!nrow(dt)) {
        warning(
          sprintf(
            "Input contains no rows; skipping: %s",
            input_csv
          ),
          call. = FALSE
        )
        next
      }
      
      dt[, var_group := variable_name]
      dt[, rp_years := as.integer(rp)]
      dt[, cluster_set := cluster_set]
      dt[, source_csv := input_csv]
      
      piece_index <- piece_index + 1L
      pieces[[piece_index]] <- dt
    }
  }
  
  if (!length(pieces)) {
    stop(
      "No readable cluster-characteristics files were found for ",
      variable_name,
      "."
    )
  }
  
  rbindlist(
    pieces,
    use.names = TRUE,
    fill = TRUE
  )
}

.make_annual_summary <- function(
    dt,
    strategy_map,
    extent_metric = "max_t_extent_km2"
) {
  dt <- copy(dt)
  
  time_col <- .get_year_col(dt)
  id_col <- .get_cluster_id_col(dt)
  intensity_col <- .get_peak_intensity_col(dt)
  extent_col <- .get_extent_col(dt, extent_metric)
  
  dt[, event_year := .parse_event_year(get(time_col))]
  dt[, duration_days_use := .get_duration_days(.SD)]
  
  dt <- dt[
    !is.na(event_year) &
      is.finite(event_year) &
      rp_years %in% RETURN_PERIODS
  ]
  
  if (!nrow(dt)) {
    stop(
      "No valid rows remained after parsing event years and filtering ",
      "return periods."
    )
  }
  
  dt <- merge(
    dt,
    strategy_map,
    by = "cluster_set",
    all.x = TRUE,
    sort = FALSE
  )
  
  if (anyNA(dt$strategy)) {
    missing_sets <- unique(dt[is.na(strategy), cluster_set])
    
    stop(
      "Missing strategy-map entries for: ",
      paste(missing_sets, collapse = ", ")
    )
  }
  
  annual <- dt[
    ,
    .(
      annual_peak_intensity = .safe_mean(get(intensity_col)),
      annual_duration_days = .safe_mean(duration_days_use),
      annual_frequency = uniqueN(get(id_col)),
      annual_extent_km2 = .safe_mean(get(extent_col))
    ),
    by = .(
      var_group,
      rp_years,
      event_year,
      cluster_set,
      strategy,
      strategy_label
    )
  ]
  
  setorder(
    annual,
    strategy,
    rp_years,
    event_year
  )
  
  annual[]
}

.to_long_characteristics <- function(
    annual,
    intensity_label
) {
  long <- melt(
    annual,
    id.vars = c(
      "var_group",
      "rp_years",
      "event_year",
      "cluster_set",
      "strategy",
      "strategy_label"
    ),
    measure.vars = c(
      "annual_peak_intensity",
      "annual_duration_days",
      "annual_frequency",
      "annual_extent_km2"
    ),
    variable.name = "characteristic",
    value.name = "value"
  )
  
  long[
    ,
    characteristic_label := factor(
      characteristic,
      levels = c(
        "annual_peak_intensity",
        "annual_duration_days",
        "annual_frequency",
        "annual_extent_km2"
      ),
      labels = c(
        intensity_label,
        "Duration (days)",
        "Frequency",
        "Extent (km²)"
      )
    )
  ]
  
  long[
    ,
    rp_label := factor(
      rp_years,
      levels = RETURN_PERIODS,
      labels = paste0(
        RETURN_PERIODS,
        "-year"
      )
    )
  ]
  
  long[]
}

.make_time_series_plot <- function(
    long,
    rp_palette,
    log_duration = FALSE
) {
  dt_plot <- copy(long)
  
  dt_plot <- dt_plot[
    !is.na(value) &
      is.finite(value)
  ]
  
  if (!nrow(dt_plot)) {
    stop(
      "No finite annual summary values are available to plot."
    )
  }
  
  dt_plot[
    ,
    strategy := factor(
      as.character(strategy),
      levels = STRATEGY_LEVELS
    )
  ]
  
  dt_plot[
    ,
    strategy_label := factor(
      as.character(strategy_label),
      levels = STRATEGY_LABELS
    )
  ]
  
  dt_plot[
    ,
    rp_label := factor(
      as.character(rp_label),
      levels = names(rp_palette)
    )
  ]
  
  duration_label <- "Duration (days)"
  extent_label <- "Extent (km²)"
  
  # Logarithmic scales cannot include zero or negative values.
  dt_plot[
    characteristic_label == extent_label &
      value <= 0,
    value := NA_real_
  ]
  
  if (isTRUE(log_duration)) {
    dt_plot[
      characteristic_label == duration_label &
        value <= 0,
      value := NA_real_
    ]
  }
  
  dt_plot <- dt_plot[
    !is.na(value) &
      is.finite(value)
  ]
  
  if (!nrow(dt_plot)) {
    stop(
      "No finite positive values remain after applying the ",
      "row-specific logarithmic-scale requirements."
    )
  }
  
  p <- ggplot(
    dt_plot,
    aes(
      x = event_year,
      y = value,
      color = rp_label,
      fill = rp_label,
      group = rp_label
    )
  ) +
    geom_point(
      size = 1.6,
      alpha = 0.85
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      level = 0.95,
      linewidth = 0.95,
      alpha = 0.20
    ) +
    facet_grid(
      rows = vars(characteristic_label),
      cols = vars(strategy_label),
      scales = "free_y",
      drop = FALSE
    ) +
    scale_color_manual(
      values = rp_palette,
      breaks = names(rp_palette),
      limits = names(rp_palette),
      drop = FALSE,
      name = "Return period"
    ) +
    scale_fill_manual(
      values = rp_palette,
      breaks = names(rp_palette),
      limits = names(rp_palette),
      drop = FALSE,
      name = "Return period"
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(
        n = 6
      )
    ) +
    labs(
      x = "Year",
      y = NULL
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          linewidth = 1.0,
          alpha = 1
        )
      ),
      fill = "none"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(
        size = 11
      ),
      legend.text = element_text(
        size = 10
      ),
      strip.background = element_rect(
        fill = "white"
      ),
      strip.text.x = element_text(
        size = 10.5,
        face = "bold"
      ),
      strip.text.y = element_text(
        size = 10.5,
        face = "bold"
      ),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_blank()
    )
  
  row_labels <- levels(
    dt_plot$characteristic_label
  )
  
  y_scales <- lapply(
    row_labels,
    function(label) {
      if (identical(label, extent_label)) {
        return(
          scale_y_log10(
            labels = label_number()
          )
        )
      }
      
      if (
        isTRUE(log_duration) &&
        identical(label, duration_label)
      ) {
        return(
          scale_y_log10(
            labels = label_number()
          )
        )
      }
      
      scale_y_continuous(
        labels = label_number()
      )
    }
  )
  
  p +
    ggh4x::facetted_pos_scales(
      y = y_scales
    )
}


# =============================================================================
# SETTINGS
# =============================================================================

ROOT <- Sys.getenv(
  "CISC_ROOT",
  "/proj/sebastianlab/users/hquintal/cisc"
)

RETURN_PERIODS <- c(
  5L,
  10L,
  25L,
  50L
)

CLUSTER_SETS <- c(
  "heatindex_0.25",
  "heatindex_0.30_bilinear",
  "heatindex_0.30_max",
  "heatindex_0.36_bilinear",
  "heatindex_0.36_max"
)

STRATEGY_LEVELS <- c(
  "native",
  "stm1_bilinear",
  "stm1_max",
  "stm2_bilinear",
  "stm2_max"
)

STRATEGY_LABELS <- c(
  "Native\n(0.25°)",
  "STM1 + Bilinear\n(0.30°)",
  "STM1 + Max\n(0.30°)",
  "STM2 + Bilinear\n(0.36°)",
  "STM2 + Max\n(0.36°)"
)

STRATEGY_MAP <- data.table(
  cluster_set = CLUSTER_SETS,
  strategy = factor(
    STRATEGY_LEVELS,
    levels = STRATEGY_LEVELS
  ),
  strategy_label = factor(
    STRATEGY_LABELS,
    levels = STRATEGY_LABELS
  )
)

# Use four ordered colors sampled from the full nine-color Brewer Reds scale.
RP_PALETTE <- brewer.pal(
  n = 9,
  name = "Reds"
)[c(3, 5, 7, 9)]

names(RP_PALETTE) <- paste0(
  RETURN_PERIODS,
  "-year"
)

EXTENT_METRIC <- Sys.getenv(
  "CISC_EXTENT_METRIC",
  "max_t_extent_km2"
)

OUT_DATA <- file.path(
  ROOT,
  "data",
  "out",
  "07_characteristics",
  "annual_summaries",
  "fs08_heat_cluster_time_series.csv"
)

OUT_FIG <- file.path(
  ROOT,
  "figures",
  "fs08_heat_cluster_time_series",
  "fs08_heat_cluster_time_series.png"
)

dir_create(
  path_dir(OUT_DATA),
  recurse = TRUE
)

dir_create(
  path_dir(OUT_FIG),
  recurse = TRUE
)

# =============================================================================
# READ, SUMMARIZE, AND PLOT
# =============================================================================

heat <- .read_characteristic_sets(
  root = ROOT,
  return_periods = RETURN_PERIODS,
  cluster_sets = CLUSTER_SETS,
  variable_name = "heatindex"
)

.log(
  "Loaded %s heat-index cluster rows.",
  comma(nrow(heat))
)

annual <- .make_annual_summary(
  dt = heat,
  strategy_map = STRATEGY_MAP,
  extent_metric = EXTENT_METRIC
)

fwrite(
  annual,
  OUT_DATA
)

.log(
  "Wrote annual heat-index summaries: %s",
  OUT_DATA
)

long <- .to_long_characteristics(
  annual = annual,
  intensity_label = "Peak intensity (°C)"
)

plot_fs08 <- .make_time_series_plot(
  long = long,
  rp_palette = RP_PALETTE,
  log_duration = TRUE
)

ggsave(
  filename = OUT_FIG,
  plot = plot_fs08,
  width = 17,
  height = 12,
  dpi = 300,
  bg = "white"
)

.log(
  "Wrote FS08 heat-index figure: %s",
  OUT_FIG
)

.log("Done.")