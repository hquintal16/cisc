#!/usr/bin/env Rscript

# HCQ Leftoff generating individual events

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(scales)
  library(tigris)
  library(grid)
})

options(
  tigris_use_cache = TRUE,
  tigris_class = "sf"
)

# =============================================================================
# Create heat-index and precipitation comparison figures from selected
# cluster GeoTIFFs.
#
# Layout:
#
#                 Native       Bilinear       Maximum
#   gamma2 row    native        gamma2         gamma2
#   gamma3 row    native        gamma3         gamma3
#
# The native panel spans both map rows in the first column.
#
# A single shared discrete legend is placed in a dedicated bottom guide row.
#
# Each panel includes:
#   1. Return-period raster cells
#   2. Black dissolved outline around all finite raster cells > 0
#   3. State boundaries plotted last
#
# Inputs:
#   data/out/08_tifs_pngs/geotiffs/<dataset_tag>/<cluster_id>.tif
#
# Outputs:
#   figures/08_tifs_pngs_panels/
#     fig_heatindex_cluster_resampling_comparison.png
#     fig_precipitation_cluster_resampling_comparison.png
# =============================================================================

# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------

ROOT <- Sys.getenv(
  "CISC_ROOT",
  "/proj/sebastianlab/users/hquintal/cisc"
)

GEOTIFF_ROOT <- file.path(
  ROOT,
  "data/out/08_tifs_pngs/geotiffs"
)

OUT_DIR <- file.path(
  ROOT,
  "figures/08_tifs_pngs_panels"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

HEAT_CLUSTER_IDS <- c(
  "heatindex_0.25"          = "heatindex-heatindex_0.25-2023-06-15-001",
  "heatindex_0.30_bilinear" = "heatindex-heatindex_0.30_bilinear-2023-06-15-002",
  "heatindex_0.30_max"      = "heatindex-heatindex_0.30_max-2023-06-15-001",
  "heatindex_0.36_bilinear" = "heatindex-heatindex_0.36_bilinear-2023-06-15-001",
  "heatindex_0.36_max"      = "heatindex-heatindex_0.36_max-2023-06-15-001"
)

PRECIP_CLUSTER_IDS <- c(
  "precipitation_0.25"          = "precipitation-precipitation_0.25-2023-08-30-001",
  "precipitation_0.11_bilinear" = "precipitation-precipitation_0.11_bilinear-2023-08-30-002",
  "precipitation_0.11_max"      = "precipitation-precipitation_0.11_max-2023-08-30-002",
  "precipitation_0.14_bilinear" = "precipitation-precipitation_0.14_bilinear-2023-08-30-002",
  "precipitation_0.14_max"      = "precipitation-precipitation_0.14_max-2023-08-30-001"
)

# Optional CSV override.
#
# Required columns:
#   variable,dataset_tag,cluster_id
SELECTION_CSV <- Sys.getenv(
  "CISC_PANEL_SELECTION_CSV",
  ""
)

FIG_WIDTH_IN <- 8.5
FIG_HEIGHT_IN <- 6.15
DPI <- 300

BASE_SIZE <- 8

# -----------------------------------------------------------------------------
# FIXED RETURN-PERIOD CLASSES
# -----------------------------------------------------------------------------

RP_LEVELS <- c(
  "2",
  "5",
  "10",
  "25",
  "50",
  "100",
  "200",
  "500",
  "1000"
)

RP_CUT_BREAKS <- c(
  2,
  5,
  10,
  25,
  50,
  100,
  200,
  500,
  1000,
  Inf
)

RP_COLORS <- c(
  "2"    = "#d9d9d9",
  "5"    = "#fff5f0",
  "10"   = "#fee0d2",
  "25"   = "#fcbba1",
  "50"   = "#fc9272",
  "100"  = "#fb6a4a",
  "200"  = "#ef3b2c",
  "500"  = "#cb181d",
  "1000" = "#99000d"
)

# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------

.log <- function(fmt, ...) {
  message(
    sprintf(
      "[%s] %s",
      format(Sys.time(), "%H:%M:%S"),
      sprintf(fmt, ...)
    )
  )
}

.validate_selection <- function(
    selection,
    expected_tags,
    variable_name
) {
  missing_tags <- setdiff(
    expected_tags,
    names(selection)
  )
  
  if (length(missing_tags)) {
    stop(
      variable_name,
      " selection is missing dataset tags: ",
      paste(missing_tags, collapse = ", "),
      call. = FALSE
    )
  }
  
  bad <- grepl(
    "^REPLACE_WITH_",
    unname(selection[expected_tags])
  )
  
  if (any(bad)) {
    stop(
      "Replace placeholder cluster IDs for ",
      variable_name,
      ": ",
      paste(expected_tags[bad], collapse = ", "),
      call. = FALSE
    )
  }
}

.read_selection_csv <- function(path) {
  tab <- readr::read_csv(
    path,
    show_col_types = FALSE
  )
  
  required <- c(
    "variable",
    "dataset_tag",
    "cluster_id"
  )
  
  missing_cols <- setdiff(
    required,
    names(tab)
  )
  
  if (length(missing_cols)) {
    stop(
      "Selection CSV is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  tab <- tab |>
    mutate(
      variable = trimws(as.character(variable)),
      dataset_tag = trimws(as.character(dataset_tag)),
      cluster_id = trimws(as.character(cluster_id))
    )
  
  heat <- tab |>
    filter(variable == "heatindex")
  
  precipitation <- tab |>
    filter(variable == "precipitation")
  
  list(
    heat = setNames(
      heat$cluster_id,
      heat$dataset_tag
    ),
    precipitation = setNames(
      precipitation$cluster_id,
      precipitation$dataset_tag
    )
  )
}

.cluster_tif_path <- function(
    geotiff_root,
    dataset_tag,
    cluster_id
) {
  file.path(
    geotiff_root,
    dataset_tag,
    paste0(cluster_id, ".tif")
  )
}

.load_states <- function(target_crs) {
  states_sf <- tigris::states(
    cb = TRUE,
    year = 2023,
    class = "sf"
  )
  
  states_sf <- states_sf |>
    dplyr::filter(
      !STUSPS %in% c(
        "AK",
        "HI",
        "PR",
        "VI",
        "GU",
        "MP",
        "AS"
      )
    )
  
  sf::st_transform(
    states_sf,
    target_crs
  )
}

.read_cluster_raster <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Cluster GeoTIFF not found: ",
      path,
      call. = FALSE
    )
  }
  
  r <- terra::rast(path)[[1]]
  
  names(r) <- "rp"
  
  list(
    raster = r,
    extent = terra::ext(r),
    crs = terra::crs(r)
  )
}

.raster_to_df <- function(r) {
  df <- as.data.frame(
    r,
    xy = TRUE,
    na.rm = TRUE
  )
  
  names(df) <- c(
    "x",
    "y",
    "rp"
  )
  
  df <- df |>
    mutate(
      rp_class = cut(
        rp,
        breaks = RP_CUT_BREAKS,
        labels = RP_LEVELS,
        right = FALSE,
        include.lowest = TRUE
      ),
      rp_class = factor(
        rp_class,
        levels = RP_LEVELS,
        ordered = TRUE
      )
    )
  
  df
}

.nonzero_outline_sf <- function(r) {
  footprint_mask <- terra::ifel(
    is.finite(r) & r > 0,
    1,
    NA
  )
  
  if (!any(!is.na(terra::values(footprint_mask)))) {
    return(NULL)
  }
  
  footprint <- terra::as.polygons(
    footprint_mask,
    dissolve = TRUE,
    na.rm = TRUE,
    values = FALSE
  )
  
  if (nrow(footprint) == 0L) {
    return(NULL)
  }
  
  footprint_sf <- sf::st_as_sf(footprint)
  
  footprint_sf <- sf::st_make_valid(
    footprint_sf
  )
  
  footprint_sf
}

.return_period_scale <- function() {
  scale_fill_manual(
    values = RP_COLORS,
    breaks = RP_LEVELS,
    limits = RP_LEVELS,
    drop = FALSE,
    na.value = "white",
    name = "Maximum return period (yr)",
    guide = guide_legend(
      direction = "horizontal",
      nrow = 1,
      byrow = TRUE,
      title.position = "bottom",
      title.hjust = 0.5,
      label.position = "left",
      keywidth = unit(1.10, "cm"),
      keyheight = unit(0.28, "cm"),
      override.aes = list(
        colour = "black",
        linewidth = 0.25
      )
    )
  )
}

.panel_plot <- function(
    raster_path,
    states,
    panel_title,
    show_x = TRUE,
    show_y = TRUE
) {
  obj <- .read_cluster_raster(
    raster_path
  )
  
  r <- obj$raster
  ex <- obj$extent
  
  df <- .raster_to_df(r)
  
  nonzero_outline <- .nonzero_outline_sf(r)
  
  dx <- terra::res(r)[1]
  dy <- terra::res(r)[2]
  
  p <- ggplot() +
    geom_tile(
      data = df,
      mapping = aes(
        x = x,
        y = y,
        fill = rp_class
      ),
      width = dx,
      height = dy,
      inherit.aes = FALSE,
      show.legend = TRUE
    )
  
  # Draw the dissolved raster footprint above the raster cells.
  if (!is.null(nonzero_outline)) {
    p <- p +
      geom_sf(
        data = nonzero_outline,
        mapping = aes(),
        fill = NA,
        colour = "black",
        linewidth = 0.45,
        linejoin = "mitre",
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
  
  # Draw state boundaries last.
  if (!is.null(states)) {
    p <- p +
      geom_sf(
        data = states,
        mapping = aes(),
        fill = NA,
        colour = "black",
        linewidth = 0.22,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
  
  p +
    .return_period_scale() +
    coord_sf(
      xlim = c(
        terra::xmin(ex),
        terra::xmax(ex)
      ),
      ylim = c(
        terra::ymin(ex),
        terra::ymax(ex)
      ),
      expand = FALSE,
      datum = NA
    ) +
    labs(
      title = panel_title,
      x = if (show_x) "Longitude" else NULL,
      y = if (show_y) "Latitude" else NULL
    ) +
    theme_bw(
      base_size = BASE_SIZE
    ) +
    theme(
      plot.title = element_text(
        size = BASE_SIZE,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 3
        )
      ),
      axis.title = element_text(
        size = BASE_SIZE
      ),
      axis.text = element_text(
        size = BASE_SIZE - 1
      ),
      axis.ticks = element_line(
        linewidth = 0.3
      ),
      panel.grid = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.4
      ),
      legend.position = "bottom",
      legend.title = element_text(
        size = BASE_SIZE
      ),
      legend.text = element_text(
        size = BASE_SIZE - 1
      ),
      legend.spacing.x = unit(0, "pt"),
      plot.margin = margin(
        t = 2,
        r = 2,
        b = 2,
        l = 2
      )
    )
}

.build_variable_figure <- function(
    variable,
    selection,
    geotiff_root,
    out_file
) {
  if (variable == "heatindex") {
    native_tag <- "heatindex_0.25"
    gamma2_bilinear_tag <- "heatindex_0.30_bilinear"
    gamma2_max_tag <- "heatindex_0.30_max"
    gamma3_bilinear_tag <- "heatindex_0.36_bilinear"
    gamma3_max_tag <- "heatindex_0.36_max"
    
    gamma2_label <- "γ\u2082 (0.30°)"
    gamma3_label <- "γ\u2083 (0.36°)"
  } else if (variable == "precipitation") {
    native_tag <- "precipitation_0.25"
    gamma2_bilinear_tag <- "precipitation_0.11_bilinear"
    gamma2_max_tag <- "precipitation_0.11_max"
    gamma3_bilinear_tag <- "precipitation_0.14_bilinear"
    gamma3_max_tag <- "precipitation_0.14_max"
    
    gamma2_label <- "γ\u2082 (0.11°)"
    gamma3_label <- "γ\u2083 (0.14°)"
  } else {
    stop(
      "Unsupported variable: ",
      variable,
      call. = FALSE
    )
  }
  
  expected_tags <- c(
    native_tag,
    gamma2_bilinear_tag,
    gamma2_max_tag,
    gamma3_bilinear_tag,
    gamma3_max_tag
  )
  
  .validate_selection(
    selection,
    expected_tags,
    variable
  )
  
  paths <- setNames(
    vapply(
      expected_tags,
      function(tag) {
        .cluster_tif_path(
          geotiff_root,
          tag,
          selection[[tag]]
        )
      },
      character(1)
    ),
    expected_tags
  )
  
  missing_files <- paths[
    !file.exists(paths)
  ]
  
  if (length(missing_files)) {
    stop(
      "Missing cluster GeoTIFFs:\n",
      paste(
        missing_files,
        collapse = "\n"
      ),
      call. = FALSE
    )
  }
  
  native_raster <- terra::rast(
    paths[[native_tag]]
  )[[1]]
  
  states <- .load_states(
    terra::crs(native_raster)
  )
  
  p_native <- .panel_plot(
    raster_path = paths[[native_tag]],
    states = states,
    panel_title = "Native (0.25°)",
    show_x = TRUE,
    show_y = TRUE
  )
  
  p_g2_bilinear <- .panel_plot(
    raster_path = paths[[gamma2_bilinear_tag]],
    states = states,
    panel_title = paste0(
      "Bilinear\n",
      gamma2_label
    ),
    show_x = FALSE,
    show_y = TRUE
  )
  
  p_g2_max <- .panel_plot(
    raster_path = paths[[gamma2_max_tag]],
    states = states,
    panel_title = paste0(
      "Maximum\n",
      gamma2_label
    ),
    show_x = FALSE,
    show_y = FALSE
  )
  
  p_g3_bilinear <- .panel_plot(
    raster_path = paths[[gamma3_bilinear_tag]],
    states = states,
    panel_title = paste0(
      "Bilinear\n",
      gamma3_label
    ),
    show_x = TRUE,
    show_y = TRUE
  )
  
  p_g3_max <- .panel_plot(
    raster_path = paths[[gamma3_max_tag]],
    states = states,
    panel_title = paste0(
      "Maximum\n",
      gamma3_label
    ),
    show_x = TRUE,
    show_y = FALSE
  )
  
  # A = native panel spanning both map rows
  # B-E = four resampling panels
  # F = dedicated shared legend area
  design <- "
    ABC
    ADE
    FFF
  "
  
  combined <- (
    p_native +
      p_g2_bilinear +
      p_g2_max +
      p_g3_bilinear +
      p_g3_max +
      guide_area() +
      plot_layout(
        design = design,
        guides = "collect",
        widths = c(
          1,
          1,
          1
        ),
        heights = c(
          1,
          1,
          0.12
        )
      )
  ) &
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.box.just = "center",
      legend.direction = "horizontal",
      legend.spacing.x = unit(0, "pt"),
      legend.margin = margin(
        t = 2,
        r = 0,
        b = 0,
        l = 0
      ),
      legend.box.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      )
    )
  
  ggsave(
    filename = out_file,
    plot = combined,
    width = FIG_WIDTH_IN,
    height = FIG_HEIGHT_IN,
    units = "in",
    dpi = DPI,
    bg = "white"
  )
  
  .log(
    "Wrote: %s",
    out_file
  )
}

# -----------------------------------------------------------------------------
# RUN
# -----------------------------------------------------------------------------

heat_selection <- HEAT_CLUSTER_IDS
precip_selection <- PRECIP_CLUSTER_IDS

if (nzchar(SELECTION_CSV)) {
  if (!file.exists(SELECTION_CSV)) {
    stop(
      "Selection CSV does not exist: ",
      SELECTION_CSV,
      call. = FALSE
    )
  }
  
  selection <- .read_selection_csv(
    SELECTION_CSV
  )
  
  heat_selection <- selection$heat
  precip_selection <- selection$precipitation
}

.build_variable_figure(
  variable = "heatindex",
  selection = heat_selection,
  geotiff_root = GEOTIFF_ROOT,
  out_file = file.path(
    OUT_DIR,
    "fig_heatindex_cluster_resampling_comparison.png"
  )
)

.build_variable_figure(
  variable = "precipitation",
  selection = precip_selection,
  geotiff_root = GEOTIFF_ROOT,
  out_file = file.path(
    OUT_DIR,
    "fig_precipitation_cluster_resampling_comparison.png"
  )
)

.log("Complete.")