# scripts/fig_s02_cluster_characteristics_16panel_byStrategy_rpcolor.R
# ------------------------------------------------------------------------------
# Annual cluster-characteristic summaries + 16-panel figures by variable
#
# Figure layout:
#   rows    = characteristics
#   columns = resampling strategies
#
# Within each panel:
#   - one annual time series per return period
#   - one linear trend per return period
#
# Rows have consistent y-axis scales across columns.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(stringr)
  library(lubridate)
  library(ggplot2)
  library(viridis)
  library(scales)
  library(ggh4x)
})

# -----------------------------
# HELPERS
# -----------------------------
.log <- function(fmt, ...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))
}

.parse_resolution <- function(x) {
  m <- str_match(as.character(x), "_([0-9]+\\.[0-9]+)(?:_|$)")
  suppressWarnings(as.numeric(m[, 2]))
}

.parse_method <- function(x) {
  x <- as.character(x)
  out <- rep("native", length(x))
  out[grepl("_bilinear$", x)] <- "bilinear"
  out[grepl("_max$", x)] <- "max"
  out
}

.detect_variable_from_file <- function(path) {
  b <- basename(path)
  if (grepl("heatindex", b, ignore.case = TRUE)) return("heatindex")
  if (grepl("precipitation", b, ignore.case = TRUE)) return("precipitation")
  NA_character_
}

.get_rp_col <- function(dt, var_group) {
  if ("rp" %in% names(dt) && is.numeric(dt$rp)) return("rp")
  if (tolower(var_group) == "heatindex" && "hi_floor_yr" %in% names(dt)) return("hi_floor_yr")
  if (tolower(var_group) == "precipitation" && "pr_floor_yr" %in% names(dt)) return("pr_floor_yr")
  stop("Could not identify return-period column for var_group = ", var_group)
}

.get_year_col <- function(dt) {
  if ("start_time" %in% names(dt)) return("start_time")
  if ("end_time" %in% names(dt)) return("end_time")
  stop("Could not identify a time column (start_time/end_time).")
}

.get_duration_days <- function(dt) {
  if ("duration_days" %in% names(dt)) return(dt$duration_days)
  if ("duration_hr" %in% names(dt)) return(dt$duration_hr / 24)
  stop("Missing duration_days and duration_hr.")
}

.get_cluster_id_for_count <- function(dt) {
  if ("cluster_id" %in% names(dt)) return("cluster_id")
  if ("cluster" %in% names(dt)) return("cluster")
  stop("Could not identify cluster identifier column.")
}

.make_strategy_map <- function(dt_var) {
  tmp <- unique(dt_var[, .(cluster_set)])
  tmp[, resolution_deg := .parse_resolution(cluster_set)]
  tmp[, method := .parse_method(cluster_set)]
  tmp[, is_native := method == "native"]
  
  res_levels <- sort(unique(tmp[!is_native & is.finite(resolution_deg), resolution_deg]))
  res_rank_dt <- data.table(
    resolution_deg = res_levels,
    stm_rank = paste0("stm", seq_along(res_levels))
  )
  
  tmp <- merge(tmp, res_rank_dt, by = "resolution_deg", all.x = TRUE)
  
  tmp[, strategy := fifelse(
    is_native,
    "native",
    paste(stm_rank, method, sep = "_")
  )]
  
  tmp[, strategy_label := fifelse(
    strategy == "native",
    paste0("Native\n(", number(resolution_deg, accuracy = 0.01), "°)"),
    fifelse(
      grepl("^stm\\d+_bilinear$", strategy),
      paste0(
        toupper(sub("_.*$", "", strategy)),
        " + Bilinear\n(",
        number(resolution_deg, accuracy = 0.01), "°)"
      ),
      fifelse(
        grepl("^stm\\d+_max$", strategy),
        paste0(
          toupper(sub("_.*$", "", strategy)),
          " + Max\n(",
          number(resolution_deg, accuracy = 0.01), "°)"
        ),
        strategy
      )
    )
  )]
  
  desired <- c("native", "stm1_bilinear", "stm1_max", "stm2_bilinear", "stm2_max")
  tmp[, strategy := factor(strategy, levels = desired)]
  setorder(tmp, strategy)
  
  tmp[]
}

.make_annual_summary <- function(dt, var_group, extent_metric = "max_t_extent_km2") {
  rp_col   <- .get_rp_col(dt, var_group)
  time_col <- .get_year_col(dt)
  id_col   <- .get_cluster_id_for_count(dt)
  
  if (!(extent_metric %in% names(dt))) stop("Extent metric not found: ", extent_metric)
  if (!("peak_intensity" %in% names(dt))) stop("Missing peak_intensity")
  
  dt <- copy(dt)
  dt[, var_group := tolower(as.character(var_group))]
  dt[, cluster_set := as.character(cluster_set)]
  dt[, rp_years := suppressWarnings(as.integer(get(rp_col)))]
  dt[, event_year := year(ymd_hms(get(time_col), quiet = TRUE))]
  
  if (all(is.na(dt$event_year))) {
    dt[, event_year := year(as.POSIXct(get(time_col), tz = "UTC"))]
  }
  
  dt[, duration_days_use := .get_duration_days(.SD)]
  
  dt <- dt[
    !is.na(rp_years) &
      !is.na(event_year) &
      is.finite(event_year)
  ]
  
  strat_map <- .make_strategy_map(dt)
  dt <- merge(
    dt,
    strat_map[, .(cluster_set, strategy, strategy_label)],
    by = "cluster_set",
    all.x = TRUE
  )
  
  out <- dt[, .(
    annual_peak_intensity = mean(peak_intensity, na.rm = TRUE),
    annual_duration_days  = mean(duration_days_use, na.rm = TRUE),
    annual_frequency      = uniqueN(get(id_col)),
    annual_extent_km2     = mean(get(extent_metric), na.rm = TRUE)
  ), by = .(var_group, rp_years, event_year, cluster_set, strategy, strategy_label)]
  
  out[]
}

.to_long_characteristics <- function(dt_annual, var_group_name) {
  intensity_lab <- if (tolower(var_group_name) == "heatindex") {
    "Peak intensity (°C)"
  } else if (tolower(var_group_name) == "precipitation") {
    "Peak intensity (mm)"
  } else {
    "Peak intensity"
  }
  
  melt(
    dt_annual,
    id.vars = c("var_group", "rp_years", "event_year", "cluster_set", "strategy", "strategy_label"),
    measure.vars = c(
      "annual_peak_intensity",
      "annual_duration_days",
      "annual_frequency",
      "annual_extent_km2"
    ),
    variable.name = "characteristic",
    value.name = "value"
  )[
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
        intensity_lab,
        "Duration (days)",
        "Frequency",
        "Extent (km²)"
      )
    )
  ][
    ,
    rp_label := factor(
      rp_years,
      levels = c(2, 5, 10, 25),
      labels = c("2-year", "5-year", "10-year", "25-year")
    )
  ][]
}

.make_rp_palette <- function() {
  cols <- viridis::viridis(4, option = "D", direction = -1)
  names(cols) <- c("2-year", "5-year", "10-year", "25-year")
  cols
}

.make_16panel_plot <- function(dt_long, var_group_name) {
  dt_plot <- copy(dt_long)
  dt_plot <- dt_plot[!is.na(value) & is.finite(value)]
  if (!nrow(dt_plot)) return(NULL)
  
  strategy_levels <- c("native", "stm1_bilinear", "stm1_max", "stm2_bilinear", "stm2_max")
  dt_plot[, strategy := factor(as.character(strategy), levels = strategy_levels)]
  
  strategy_lab_map <- unique(dt_plot[, .(strategy, strategy_label)])
  strategy_lab_map <- strategy_lab_map[match(strategy_levels, as.character(strategy))]
  strategy_label_levels <- strategy_lab_map$strategy_label
  strategy_label_levels <- strategy_label_levels[!is.na(strategy_label_levels)]
  
  dt_plot[, strategy_label := factor(strategy_label, levels = strategy_label_levels)]
  dt_plot[, rp_label := factor(rp_label, levels = c("2-year", "5-year", "10-year", "25-year"))]
  
  rp_pal <- .make_rp_palette()
  
  # --------------------------------
  # log-axis row controls
  # --------------------------------
  intensity_label <- levels(dt_plot$characteristic_label)[grepl("^Peak intensity", levels(dt_plot$characteristic_label))]
  duration_label  <- "Duration (days)"
  extent_label    <- "Extent (km²)"
  
  dt_plot[, use_log_y := FALSE]
  dt_plot[characteristic_label == extent_label, use_log_y := TRUE]
  
  if (tolower(var_group_name) == "heatindex") {
    dt_plot[characteristic_label == duration_label, use_log_y := TRUE]
  }
  
  # log scales cannot include <= 0
  dt_plot[use_log_y == TRUE & value <= 0, value := NA_real_]
  dt_plot <- dt_plot[!is.na(value)]
  
  p <- ggplot(
    dt_plot,
    aes(
      x = event_year,
      y = value,
      color = rp_label,
      group = rp_label
    )
  ) +
    geom_line(linewidth = 0.45, alpha = 0.45) +
    geom_point(size = 0.9, alpha = 0.75) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.95
    ) +
    facet_grid(
      rows = vars(characteristic_label),
      cols = vars(strategy_label),
      scales = "free_y",
      drop = FALSE
    ) +
    scale_color_manual(
      values = rp_pal,
      breaks = c("2-year", "5-year", "10-year", "25-year"),
      drop = FALSE,
      name = "Return period"
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(n = 6)
    ) +
    labs(
      x = "Year",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      strip.background = element_rect(fill = "white"),
      strip.text.x = element_text(size = 10.5, face = "bold"),
      strip.text.y = element_text(size = 10.5, face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
  
  # apply log y for the chosen rows only using faceted scale transform
  # requires ggh4x
  if (!requireNamespace("ggh4x", quietly = TRUE)) {
    stop("Please install ggh4x to use row-specific y-axis transforms: install.packages('ggh4x')")
  }
  
  y_scales_list <- list()
  
  for (lab in strategy_label_levels) {
    # intensity row
    y_scales_list[[paste(as.character(intensity_label), lab, sep = "::")]] <- scale_y_continuous()
    
    # duration row
    if (tolower(var_group_name) == "heatindex") {
      y_scales_list[[paste(duration_label, lab, sep = "::")]] <- scale_y_log10(labels = label_number())
    } else {
      y_scales_list[[paste(duration_label, lab, sep = "::")]] <- scale_y_continuous()
    }
    
    # frequency row
    y_scales_list[[paste("Frequency", lab, sep = "::")]] <- scale_y_continuous()
    
    # extent row
    y_scales_list[[paste(extent_label, lab, sep = "::")]] <- scale_y_log10(labels = label_number())
  }
  
  p +
    ggh4x::facetted_pos_scales(
      y = unname(y_scales_list)
    )
}

# -----------------------------
# USER SETTINGS
# -----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

IN_DIR <- file.path(ROOT, "data/out/07_characteristics")
SCENARIO_DIRS <- c("hi2_pr2", "hi5_pr5", "hi10_pr10", "hi25_pr25")

OUT_DATA_DIR <- file.path(IN_DIR, "annual_summaries")
OUT_FIG_DIR  <- file.path(ROOT, "figures/fig_s02")

dir_create(OUT_DATA_DIR, recurse = TRUE)
dir_create(OUT_FIG_DIR, recurse = TRUE)

RPS_TO_KEEP <- c(2, 5, 10, 25)

# Change to "total_extent_km2" if preferred
EXTENT_METRIC <- "max_t_extent_km2"

# -----------------------------
# READ ALL CHARACTERISTIC FILES
# -----------------------------
all_csvs <- unlist(lapply(
  file.path(IN_DIR, SCENARIO_DIRS),
  function(d) dir_ls(d, glob = "*.csv", type = "file")
), use.names = FALSE)

char_files <- all_csvs[
  grepl("^cluster_characteristics_(heatindex|precipitation)_.*\\.csv$",
        basename(all_csvs), ignore.case = TRUE)
]

if (!length(char_files)) {
  stop("No cluster characteristics CSVs found under: ", IN_DIR)
}

.log("Found %d cluster characteristics file(s).", length(char_files))

all_dt <- rbindlist(lapply(char_files, function(f) {
  .log("Reading %s", f)
  dt <- fread(f)
  dt[, source_csv := f]
  dt[, scenario_dir := basename(dirname(f))]
  dt[, var_group := .detect_variable_from_file(f)]
  dt
}), use.names = TRUE, fill = TRUE)

all_dt <- all_dt[tolower(var_group) %in% c("heatindex", "precipitation")]

# -----------------------------
# BUILD ANNUAL SUMMARIES
# -----------------------------
annual_list <- list()

for (vg in c("heatindex", "precipitation")) {
  dsub <- all_dt[tolower(var_group) == vg]
  if (!nrow(dsub)) next
  
  .log("Summarizing %s", vg)
  
  annual_dt <- .make_annual_summary(
    dt = dsub,
    var_group = vg,
    extent_metric = EXTENT_METRIC
  )
  
  annual_dt <- annual_dt[!is.na(rp_years) & rp_years %in% RPS_TO_KEEP]
  
  out_csv <- file.path(
    OUT_DATA_DIR,
    paste0("annual_cluster_characteristics_", vg, ".csv")
  )
  fwrite(annual_dt, out_csv)
  .log("Wrote %s", out_csv)
  
  annual_list[[vg]] <- annual_dt
}

# -----------------------------
# PLOT FIGURES
# -----------------------------
for (vg in names(annual_list)) {
  annual_dt <- annual_list[[vg]]
  if (!nrow(annual_dt)) next
  
  long_dt <- .to_long_characteristics(annual_dt, vg)
  
  p <- .make_16panel_plot(long_dt, vg)
  if (is.null(p)) next
  
  out_png <- file.path(
    OUT_FIG_DIR,
    paste0("fig_s02_cluster_characteristics_16panel_byStrategy_", vg, ".png")
  )
  
  ggsave(
    filename = out_png,
    plot = p,
    width = 17,
    height = 12,
    dpi = 300
  )
  .log("Wrote %s", out_png)
}

.log("Done.")