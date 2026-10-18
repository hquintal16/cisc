# scripts/fig_06_plot_f1diff_spatial_county_specificscenario_twovars.R
# =============================================================================
# Fig 06: 2-panel county maps of F1 difference relative to native resolution
#
# Left  panel: Heat Index    ΔF1 = F1(selected resampled dataset) - F1(native)
# Right panel: Precipitation ΔF1 = F1(selected resampled dataset) - F1(native)
#
# Updated to:
# - use clean lon/lat axes in EPSG:4326
# - avoid repeated graticule labels
# - keep grid lines from showing through NA land areas
# - use 20 pt Arial for all text
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tigris)
  library(ragg)
  library(patchwork)
  library(terra)
  library(scales)
  library(grid)
})

options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)

# ----------------------------
# USER PATHS
# ----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

MET_HEAT <- file.path(
  ROOT,
  "data/out/09_performance/metrics_county/metrics_heatindex_all_heat_multiwindow.csv"
)

MET_PREC <- file.path(
  ROOT,
  "data/out/09_performance/metrics_county/metrics_precipitation_all_precip_multiwindow.csv"
)

WFO_SHP <- "V:/users/hquintal/phd1_cluster_southeast/data/input/regional_aggregation/nws_wfo_cws/w_05mr24/w_05mr24.shp"

OUT_DIR <- file.path(ROOT, "figures/fig_00_evaluation_spatial")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# USER SETTINGS
# ----------------------------
CRS_PLOT <- 4326
BBOX_LONLAT <- c(xmin = -95, xmax = -75, ymin = 25, ymax = 40)

# Choose scenario independently for each variable
HEAT_SCENARIO_TO_PLOT <- "hi5_pr5"
PREC_SCENARIO_TO_PLOT <- "hi2_pr2"

# Choose one window independently if needed
HEAT_WINDOW_TYPE_TO_PLOT <- "fixed"
HEAT_WINDOW_TAG_TO_PLOT  <- "20020101_20231231"

PREC_WINDOW_TYPE_TO_PLOT <- "fixed"
PREC_WINDOW_TAG_TO_PLOT  <- "20020101_20231231"

# Choose one dataset independently for each variable
HEAT_DATASET_TO_PLOT <- "heatindex_0.36_max"
PREC_DATASET_TO_PLOT <- "precipitation_0.14_max"

# Matching native datasets
HEAT_NATIVE_DATASET <- "heatindex_0.25"
PREC_NATIVE_DATASET <- "precipitation_0.25"

# Figure settings
RES_MAIN      <- 200
WIDTH_PANELS  <- 3600
HEIGHT_PANELS <- 1900

# ----------------------------
# Logging
# ----------------------------
.log <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
}

# ----------------------------
# Helpers
# ----------------------------
.norm_geoid5 <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_extract(x, "\\d+")
  x <- ifelse(is.na(x) | !nzchar(x), NA_character_, x)
  str_pad(x, 5, pad = "0")
}

.read_metrics <- function(path) {
  if (!file.exists(path)) stop("Missing metrics file: ", path)
  
  dt <- fread(path, showProgress = FALSE)
  
  req <- c(
    "scenario", "dataset_tag", "var_group", "event",
    "window_type", "window_tag", "GEOID", "f1"
  )
  miss <- setdiff(req, names(dt))
  if (length(miss)) {
    stop("Metrics missing columns: ", paste(miss, collapse = ", "), "\nFile: ", path)
  }
  
  dt[, GEOID := .norm_geoid5(GEOID)]
  dt[, scenario    := as.character(scenario)]
  dt[, dataset_tag := as.character(dataset_tag)]
  dt[, var_group   := tolower(as.character(var_group))]
  dt[, event       := as.character(event)]
  dt[, window_type := as.character(window_type)]
  dt[, window_tag  := as.character(window_tag)]
  dt[, f1          := as.numeric(f1)]
  
  dt[!is.na(GEOID) & nzchar(GEOID)]
}

.clip_to_mask <- function(x, mask_poly) {
  x <- st_make_valid(x)
  suppressWarnings(st_intersection(x, mask_poly))
}

.make_bbox_poly <- function(bbox_lonlat, crs_out = 4326) {
  bb_ll <- st_as_sfc(st_bbox(bbox_lonlat, crs = st_crs(4326)))
  if (!is.null(crs_out) && crs_out != 4326) {
    bb_ll <- st_transform(bb_ll, crs_out)
  }
  st_make_valid(bb_ll)
}

.extract_heat_rp <- function(scenario) {
  x <- str_match(scenario, "(?:^|_)hi(\\d+)(?:_|$)")[, 2]
  ifelse(is.na(x), NA_character_, paste0(x, "-yr"))
}

.extract_prec_rp <- function(scenario) {
  x <- str_match(scenario, "(?:^|_)pr(\\d+)(?:_|$)")[, 2]
  ifelse(is.na(x), NA_character_, paste0(x, "-yr"))
}

.parse_dataset_info <- function(dataset_tag, var_group) {
  out <- list(
    resolution = NA_character_,
    stm        = "Native",
    method     = "Native"
  )
  
  if (tolower(var_group) == "heatindex") {
    if (dataset_tag == "heatindex_0.25") {
      out$resolution <- "0.25°"
      out$stm        <- "Native"
      out$method     <- "Native"
      return(out)
    }
    if (dataset_tag %in% c("heatindex_0.30_bilinear", "heatindex_0.30_max")) {
      out$resolution <- "0.30°"
      out$stm        <- "STM1"
    }
    if (dataset_tag %in% c("heatindex_0.36_bilinear", "heatindex_0.36_max")) {
      out$resolution <- "0.36°"
      out$stm        <- "STM2"
    }
    if (dataset_tag %in% c("heatindex_0.34_bilinear", "heatindex_0.34_max")) {
      out$resolution <- "0.34°"
      out$stm        <- "STM1"
    }
    if (dataset_tag %in% c("heatindex_0.44_bilinear", "heatindex_0.44_max")) {
      out$resolution <- "0.44°"
      out$stm        <- "STM2"
    }
  }
  
  if (tolower(var_group) == "precipitation") {
    if (dataset_tag == "precipitation_0.25") {
      out$resolution <- "0.25°"
      out$stm        <- "Native"
      out$method     <- "Native"
      return(out)
    }
    if (dataset_tag %in% c("precipitation_0.11_bilinear", "precipitation_0.11_max")) {
      out$resolution <- "0.11°"
      out$stm        <- "STM1"
    }
    if (dataset_tag %in% c("precipitation_0.14_bilinear", "precipitation_0.14_max")) {
      out$resolution <- "0.14°"
      out$stm        <- "STM2"
    }
    if (dataset_tag %in% c("precipitation_0.12_bilinear", "precipitation_0.12_max")) {
      out$resolution <- "0.12°"
      out$stm        <- "STM1"
    }
    if (dataset_tag %in% c("precipitation_0.17_bilinear", "precipitation_0.17_max")) {
      out$resolution <- "0.17°"
      out$stm        <- "STM2"
    }
  }
  
  if (grepl("_bilinear$", dataset_tag)) out$method <- "Bilinear"
  if (grepl("_max$", dataset_tag))      out$method <- "Max"
  
  out
}

.make_panel_title_diff <- function(var_group, scenario, dataset_tag) {
  info <- .parse_dataset_info(dataset_tag, var_group)
  
  rp <- if (tolower(var_group) == "heatindex") {
    .extract_heat_rp(scenario)
  } else {
    .extract_prec_rp(scenario)
  }
  
  vg_lab <- if (tolower(var_group) == "heatindex") "Heat Index" else "Precipitation"
  
  paste0(
    vg_lab,
    " | RP ", rp,
    " | ", info$stm,
    " | ", info$method,
    " | F1 Improvement (%) vs Native"
  )
}

.build_f1_dt <- function(dt, var_group_value, event_name, scenario_value,
                         window_type_value, window_tag_value, dataset_tag_value) {
  
  d <- dt[
    var_group   == tolower(var_group_value) &
      event       == event_name &
      scenario    == scenario_value &
      window_type == window_type_value &
      window_tag  == window_tag_value &
      dataset_tag == dataset_tag_value
  ]
  
  if (!nrow(d)) {
    return(data.table(GEOID = character(), f1 = numeric()))
  }
  
  d <- d[
    is.finite(f1) &
      f1 >= 0 &
      f1 <= 1
  ]
  
  d[, .(f1 = mean(f1, na.rm = TRUE)), by = GEOID]
}

.build_f1_diff_dt <- function(dt, var_group_value, event_name, scenario_value,
                              window_type_value, window_tag_value,
                              selected_dataset_tag, native_dataset_tag) {
  
  sel <- .build_f1_dt(
    dt                = dt,
    var_group_value   = var_group_value,
    event_name        = event_name,
    scenario_value    = scenario_value,
    window_type_value = window_type_value,
    window_tag_value  = window_tag_value,
    dataset_tag_value = selected_dataset_tag
  )
  
  nat <- .build_f1_dt(
    dt                = dt,
    var_group_value   = var_group_value,
    event_name        = event_name,
    scenario_value    = scenario_value,
    window_type_value = window_type_value,
    window_tag_value  = window_tag_value,
    dataset_tag_value = native_dataset_tag
  )
  
  if (!nrow(sel)) {
    stop(
      "No selected rows found for:\n",
      "var_group=", var_group_value,
      " | scenario=", scenario_value,
      " | window_type=", window_type_value,
      " | window_tag=", window_tag_value,
      " | dataset_tag=", selected_dataset_tag
    )
  }
  
  if (!nrow(nat)) {
    stop(
      "No native rows found for:\n",
      "var_group=", var_group_value,
      " | scenario=", scenario_value,
      " | window_type=", window_type_value,
      " | window_tag=", window_tag_value,
      " | dataset_tag=", native_dataset_tag
    )
  }
  
  out <- merge(
    sel, nat,
    by = "GEOID",
    all = TRUE,
    suffixes = c("_sel", "_nat")
  )
  
  out[, f1_diff := f1_sel - f1_nat]
  
  out[, f1_improve_pct := fifelse(
    is.finite(f1_nat) & f1_nat > 0,
    100 * (f1_sel - f1_nat) / f1_nat,
    NA_real_
  )]
  
  out[, .(GEOID, f1_sel, f1_nat, f1_diff, f1_improve_pct)]
}

.theme_map <- function() {
  theme_minimal(base_size = 20, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial", size = 20),
      panel.grid.major = element_line(color = "grey70", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.title = element_blank(),
      axis.text = element_text(size = 20, family = "Arial", color = "grey25"),
      axis.ticks = element_line(color = "grey25", linewidth = 0.3),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 20, family = "Arial"),
      legend.text  = element_text(size = 20, family = "Arial"),
      legend.key.width  = grid::unit(1.2, "in"),
      legend.key.height = grid::unit(0.18, "in"),
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.justification = "center",
      plot.title = element_text(size = 20, family = "Arial"),
      plot.tag = element_text(size = 20, family = "Arial", face = "bold"),
      plot.tag.position = c(0.02, 0.95),
      plot.margin = margin(4, 4, 4, 4)
    )
}

.make_plot_diff <- function(counties_sf, states_sf, wfo_sf, mask_poly,
                            fill_limits, title_text, panel_tag, show_legend = TRUE) {
  
  bb <- st_bbox(mask_poly)
  se_land <- st_union(states_sf)
  
  p <- ggplot() +
    geom_sf(data = se_land, fill = "white", color = NA) +
    geom_sf(data = counties_sf, aes(fill = f1_improve_pct), color = NA) +
    geom_sf(data = wfo_sf, fill = NA, color = "grey55", linewidth = 0.25) +
    geom_sf(data = states_sf, fill = NA, color = "black", linewidth = 0.35) +
    scale_fill_steps2(
      low = "red",
      mid = "grey80",
      high = "blue",
      midpoint = 0,
      limits = fill_limits,
      breaks = seq(fill_limits[1], fill_limits[2], length.out = 8),
      labels = function(x) paste0(round(x, 0), "%"),
      oob = scales::squish,
      na.value = NA,
      name = "F1 Improvement (%)",
      guide = guide_colorsteps(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(8.0, "in"),
        barheight = grid::unit(0.22, "in"),
        even.steps = TRUE,
        show.limits = TRUE
      )
    ) +
    coord_sf(
      crs   = st_crs(CRS_PLOT),
      datum = st_crs(4326),
      xlim  = c(bb["xmin"], bb["xmax"]),
      ylim  = c(bb["ymin"], bb["ymax"]),
      expand = FALSE,
      label_axes = list(bottom = "E", left = "N")
    ) +
    labs(tag = panel_tag) +
    .theme_map()
  
  if (!isTRUE(show_legend)) {
    p <- p + theme(legend.position = "none")
  }
  
  p
}
# ----------------------------
# TRUE bbox polygon in lon/lat
# ----------------------------
bbox_poly <- .make_bbox_poly(BBOX_LONLAT, CRS_PLOT)
plot_mask_poly <- st_make_valid(bbox_poly)

# ----------------------------
# Geometries
# ----------------------------
COUNTIES <- counties(cb = TRUE, year = 2022, class = "sf")
STATES   <- states(cb = TRUE, year = 2022, class = "sf")

drop_statefp <- c("02", "15", "72", "60", "66", "69", "78")
COUNTIES <- COUNTIES[!(COUNTIES$STATEFP %in% drop_statefp), ]
STATES   <- STATES[!(STATES$STATEFP %in% drop_statefp), ]

COUNTIES <- st_transform(COUNTIES, CRS_PLOT) |> st_make_valid()
STATES   <- st_transform(STATES, CRS_PLOT)   |> st_make_valid()

if (!file.exists(WFO_SHP)) stop("Missing WFO shapefile: ", WFO_SHP)
WFO <- st_read(WFO_SHP, quiet = TRUE) |> st_make_valid()
WFO <- st_transform(WFO, CRS_PLOT)

COUNTIES_SE <- .clip_to_mask(COUNTIES, plot_mask_poly)
STATES_SE   <- .clip_to_mask(STATES,   plot_mask_poly)
WFO_SE      <- .clip_to_mask(WFO,      plot_mask_poly)

# ----------------------------
# Read metrics
# ----------------------------
heat_dt <- .read_metrics(MET_HEAT)
prec_dt <- .read_metrics(MET_PREC)

# ----------------------------
# Build one F1-difference dataset per variable
# ----------------------------
heat_diff_dt <- .build_f1_diff_dt(
  dt                   = heat_dt,
  var_group_value      = "heatindex",
  event_name           = "all_heat",
  scenario_value       = HEAT_SCENARIO_TO_PLOT,
  window_type_value    = HEAT_WINDOW_TYPE_TO_PLOT,
  window_tag_value     = HEAT_WINDOW_TAG_TO_PLOT,
  selected_dataset_tag = HEAT_DATASET_TO_PLOT,
  native_dataset_tag   = HEAT_NATIVE_DATASET
)

prec_diff_dt <- .build_f1_diff_dt(
  dt                   = prec_dt,
  var_group_value      = "precipitation",
  event_name           = "all_precip",
  scenario_value       = PREC_SCENARIO_TO_PLOT,
  window_type_value    = PREC_WINDOW_TYPE_TO_PLOT,
  window_tag_value     = PREC_WINDOW_TAG_TO_PLOT,
  selected_dataset_tag = PREC_DATASET_TO_PLOT,
  native_dataset_tag   = PREC_NATIVE_DATASET
)

# ----------------------------
# Join to counties
# ----------------------------
counties_heat <- COUNTIES_SE %>%
  left_join(as.data.frame(heat_diff_dt), by = "GEOID")

counties_prec <- COUNTIES_SE %>%
  left_join(as.data.frame(prec_diff_dt), by = "GEOID")

# ----------------------------
# Shared symmetric legend range across both panels
# ----------------------------
vals_heat <- counties_heat$f1_improve_pct[is.finite(counties_heat$f1_improve_pct)]
vals_prec <- counties_prec$f1_improve_pct[is.finite(counties_prec$f1_improve_pct)]
all_vals  <- c(vals_heat, vals_prec)

if (!length(all_vals)) stop("No finite F1 improvement (%) values available to plot.")

max_abs <- max(abs(all_vals), na.rm = TRUE)
if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
fill_limits <- c(-max_abs, max_abs)

.log("Heat counties with F1 improvement (%%): %d", sum(is.finite(counties_heat$f1_improve_pct)))
.log("Precip counties with F1 improvement (%%): %d", sum(is.finite(counties_prec$f1_improve_pct)))
.log("Shared F1 improvement (%%) scale: [%.2f, %.2f]", fill_limits[1], fill_limits[2])

# ----------------------------
# Plot objects
# ----------------------------
heat_title <- .make_panel_title_diff("heatindex", HEAT_SCENARIO_TO_PLOT, HEAT_DATASET_TO_PLOT)
prec_title <- .make_panel_title_diff("precipitation", PREC_SCENARIO_TO_PLOT, PREC_DATASET_TO_PLOT)

p_heat <- .make_plot_diff(
  counties_sf = counties_heat,
  states_sf   = STATES_SE,
  wfo_sf      = WFO_SE,
  mask_poly   = plot_mask_poly,
  fill_limits = fill_limits,
  title_text  = heat_title,
  panel_tag   = "(a)",
  show_legend = TRUE
)

p_prec <- .make_plot_diff(
  counties_sf = counties_prec,
  states_sf   = STATES_SE,
  wfo_sf      = WFO_SE,
  mask_poly   = plot_mask_poly,
  fill_limits = fill_limits,
  title_text  = prec_title,
  panel_tag   = "(b)",
  show_legend = FALSE
)

p_panels <- (p_heat + p_prec) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    theme = theme(
      text = element_text(family = "Arial", size = 20),
      legend.position = "bottom",
      legend.box.just = "center",
      legend.justification = "center",
      panel.spacing = unit(0.15, "lines")
    )
  )

# ----------------------------
# Output
# ----------------------------
safe_heat_scen <- gsub("[^A-Za-z0-9]+", "_", HEAT_SCENARIO_TO_PLOT)
safe_prec_scen <- gsub("[^A-Za-z0-9]+", "_", PREC_SCENARIO_TO_PLOT)
safe_heat_tag  <- gsub("[^A-Za-z0-9]+", "_", HEAT_DATASET_TO_PLOT)
safe_prec_tag  <- gsub("[^A-Za-z0-9]+", "_", PREC_DATASET_TO_PLOT)

OUT_PANELS <- file.path(
  OUT_DIR,
  paste0(
    "fig_06_f1diff_spatial_",
    "heat_", safe_heat_scen, "_", safe_heat_tag,
    "__prec_", safe_prec_scen, "_", safe_prec_tag,
    ".png"
  )
)

ragg::agg_png(OUT_PANELS, width = WIDTH_PANELS, height = HEIGHT_PANELS, res = RES_MAIN)
print(p_panels)
dev.off()

.log("Wrote: %s", OUT_PANELS)
.log("Done.")