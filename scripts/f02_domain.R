# #!/usr/bin/env Rscript
# 
# # scripts/f02_domain.R
# # =============================================================================
# # Study-domain figure
# #
# # Main map:
# #   Southeast study domain filled by aggregated ecoregion, with state outlines
# #   and state abbreviations.
# #
# # Inset:
# #   CONUS with the Southeast study-domain bounding box outlined in red.
# #
# # Figure-only workflow. No intermediate CSV, Parquet, QS, or spatial files are
# # written.
# # =============================================================================
# 
# suppressPackageStartupMessages({
#   library(sf)
#   library(ggplot2)
#   library(dplyr)
#   library(cowplot)
#   library(scales)
#   library(ragg)
#   library(RColorBrewer)
#   library(tigris)
#   library(fs)
# })
# 
# options(tigris_use_cache = TRUE)
# sf::sf_use_s2(FALSE)
# 
# # =============================================================================
# # User settings
# # =============================================================================
# 
# ROOT <- Sys.getenv(
#   "CISC_ROOT",
#   "/proj/sebastianlab/users/hquintal/cisc"
# )
# 
# OUT_DIR <- file.path(ROOT, "figures", "f02_domain")
# OUT_PNG <- file.path(OUT_DIR, "f02_domain.png")
# OUT_PDF <- file.path(OUT_DIR, "f02_domain.pdf")
# 
# WRITE_PDF <- FALSE
# CRS_PLOT <- 4326
# 
# BBOX_LONLAT <- c(
#   xmin = -95,
#   xmax = -75,
#   ymin = 25,
#   ymax = 40
# )
# 
# RES_MAIN <- 300
# WIDTH_IN <- 12
# HEIGHT_IN <- 8.1
# 
# ECO_SHP <- file.path(
#   ROOT,
#   "data/in/shapefile/Aggr_Ecoregions_2015/Aggr_Ecoregions_2015.shp"
# )
# 
# # =============================================================================
# # Helpers
# # =============================================================================
# 
# .log <- function(fmt, ...) {
#   message(
#     sprintf(
#       "[%s] %s",
#       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
#       sprintf(fmt, ...)
#     )
#   )
# }
# 
# .make_bbox <- function(bbox_values) {
#   sf::st_as_sfc(
#     sf::st_bbox(
#       c(
#         xmin = as.numeric(bbox_values[["xmin"]]),
#         ymin = as.numeric(bbox_values[["ymin"]]),
#         xmax = as.numeric(bbox_values[["xmax"]]),
#         ymax = as.numeric(bbox_values[["ymax"]])
#       ),
#       crs = sf::st_crs(4326)
#     )
#   )
# }
# 
# .prepare_states <- function() {
#   states <- tigris::states(cb = TRUE, year = 2022, class = "sf")
#   drop_statefp <- c("02", "15", "60", "66", "69", "72", "78")
#   states <- states[!(as.character(states$STATEFP) %in% drop_statefp), ]
#   states <- sf::st_make_valid(states)
#   states <- sf::st_transform(states, CRS_PLOT)
#   states$abb <- as.character(states$STUSPS)
#   states
# }
# 
# .prepare_state_labels <- function(states_clip) {
#   label_points <- states_clip |>
#     sf::st_transform(5070) |>
#     sf::st_point_on_surface() |>
#     sf::st_transform(4326)
#   
#   label_points[
#     !is.na(label_points$abb) &
#       label_points$abb != "",
#   ]
# }
# 
# .build_ecoregion_palette <- function(levels) {
#   levels <- as.character(levels)
#   levels <- levels[!is.na(levels) & levels != ""]
#   levels <- sort(unique(levels))
#   
#   if (!length(levels)) {
#     stop("No ecoregion levels were found inside the study domain.")
#   }
#   
#   base_palette <- RColorBrewer::brewer.pal(8, "Accent")
#   colours <- if (length(levels) <= length(base_palette)) {
#     base_palette[seq_along(levels)]
#   } else {
#     grDevices::colorRampPalette(base_palette)(length(levels))
#   }
#   
#   stats::setNames(colours, levels)
# }
# 
# .theme_main_map <- function() {
#   ggplot2::theme_bw(base_size = 18, base_family = "sans") +
#     ggplot2::theme(
#       text = ggplot2::element_text(family = "sans", size = 18),
#       panel.grid.major = ggplot2::element_line(color = "grey82", linewidth = 0.35),
#       panel.grid.minor = ggplot2::element_blank(),
#       panel.background = ggplot2::element_rect(fill = "white", color = NA),
#       plot.background = ggplot2::element_rect(fill = "white", color = NA),
#       axis.title = ggplot2::element_blank(),
#       axis.text = ggplot2::element_text(size = 14, color = "grey25"),
#       axis.ticks = ggplot2::element_line(color = "grey25", linewidth = 0.3),
#       legend.position = "bottom",
#       legend.direction = "horizontal",
#       legend.title = ggplot2::element_text(size = 14),
#       legend.text = ggplot2::element_text(size = 11),
#       # adjust
#       legend.key.width = grid::unit(0.5, "cm"),
#       legend.key.height = grid::unit(0.5, "cm"),
#       legend.spacing.x = grid::unit(0.8, "cm"),
#       # adjust
#       legend.margin = ggplot2::margin(2, 0, 0, 0),
#       legend.background = ggplot2::element_blank(),
#       legend.box.background = ggplot2::element_blank(),
#       legend.box.just = "center",
#       legend.justification = "center",
#       plot.margin = ggplot2::margin(8, 8, 2, 8)
#     )
# }
# 
# # =============================================================================
# # Resolve inputs
# # =============================================================================
# 
# fs::dir_create(OUT_DIR, recurse = TRUE)
# 
# if (!file.exists(ECO_SHP)) {
#   stop("Missing ecoregion shapefile: ", ECO_SHP)
# }
# 
# .log("Ecoregion shapefile: %s", ECO_SHP)
# .log("Output directory: %s", OUT_DIR)
# 
# bbox_poly <- .make_bbox(BBOX_LONLAT)
# 
# # =============================================================================
# # Ecoregions
# # =============================================================================
# 
# ecoregions_all <- sf::st_read(ECO_SHP, quiet = TRUE)
# 
# if (!("WSA9_NAME" %in% names(ecoregions_all))) {
#   stop(
#     "Ecoregion shapefile is missing WSA9_NAME. Available fields: ",
#     paste(names(ecoregions_all), collapse = ", ")
#   )
# }
# 
# ecoregions_all <- ecoregions_all |>
#   sf::st_make_valid() |>
#   sf::st_transform(CRS_PLOT)
# 
# ecoregions_clip <- suppressWarnings(
#   sf::st_intersection(ecoregions_all, bbox_poly)
# )
# 
# ecoregions_clip <- ecoregions_clip[!sf::st_is_empty(ecoregions_clip), ]
# ecoregions_clip$WSA9_NAME <- trimws(as.character(ecoregions_clip$WSA9_NAME))
# ecoregions_clip <- ecoregions_clip[
#   !is.na(ecoregions_clip$WSA9_NAME) & ecoregions_clip$WSA9_NAME != "",
# ]
# 
# if (!nrow(ecoregions_clip)) {
#   stop("No ecoregions intersect the requested study-domain bounding box.")
# }
# 
# # =============================================================================
# # States and labels
# # =============================================================================
# 
# states_all <- .prepare_states()
# states_clip <- suppressWarnings(sf::st_intersection(states_all, bbox_poly))
# states_clip <- states_clip[!sf::st_is_empty(states_clip), ]
# state_labels <- .prepare_state_labels(states_clip)
# 
# # =============================================================================
# # Palette
# # =============================================================================
# 
# eco_levels <- sort(unique(ecoregions_clip$WSA9_NAME))
# eco_palette <- .build_ecoregion_palette(eco_levels)
# 
# # =============================================================================
# # Inset map
# # =============================================================================
# 
# p_inset <- ggplot2::ggplot() +
#   ggplot2::geom_sf(
#     data = states_all,
#     fill = "white",
#     color = "grey55",
#     linewidth = 0.25
#   ) +
#   ggplot2::geom_sf(
#     data = bbox_poly,
#     fill = NA,
#     color = "red",
#     linewidth = 0.9
#   ) +
#   ggplot2::coord_sf(
#     crs = sf::st_crs(4326),
#     datum = NA,
#     xlim = c(-125, -66),
#     ylim = c(24, 50),
#     expand = FALSE
#   ) +
#   ggplot2::theme_void() +
#   ggplot2::theme(
#     panel.background = ggplot2::element_rect(
#       fill = "white",
#       color = "black",
#       linewidth = 0.8
#     ),
#     plot.background = ggplot2::element_rect(fill = "white", color = NA),
#     plot.margin = ggplot2::margin(0, 0, 0, 0)
#   )
# 
# # =============================================================================
# # Main map
# # =============================================================================
# 
# p_main <- ggplot2::ggplot() +
#   ggplot2::geom_sf(
#     data = ecoregions_clip,
#     ggplot2::aes(fill = WSA9_NAME),
#     color = "grey35",
#     linewidth = 0.25
#   ) +
#   ggplot2::geom_sf(
#     data = states_clip,
#     fill = NA,
#     color = "black",
#     linewidth = 0.45
#   ) +
#   ggplot2::geom_sf_text(
#     data = state_labels,
#     ggplot2::aes(label = abb),
#     size = 4.2,
#     color = "grey20",
#     family = "sans"
#   ) +
#   ggplot2::coord_sf(
#     crs = sf::st_crs(4326),
#     datum = sf::st_crs(4326),
#     xlim = c(BBOX_LONLAT[["xmin"]], BBOX_LONLAT[["xmax"]]),
#     ylim = c(BBOX_LONLAT[["ymin"]], BBOX_LONLAT[["ymax"]]),
#     expand = FALSE,
#     label_axes = list(bottom = "E", left = "N")
#   ) +
#   ggplot2::scale_fill_manual(
#     values = eco_palette,
#     limits = eco_levels,
#     name = "Ecoregion",
#     drop = FALSE
#   ) +
#   ggplot2::guides(
#     fill = ggplot2::guide_legend(
#       nrow = 1,
#       byrow = TRUE,
#       title.position = "top",
#       title.hjust = 0.5,
#       keywidth = grid::unit(0.6, "cm"),
#       keyheight = grid::unit(0.6, "cm")
#     )
#   ) +
#   .theme_main_map()
# 
# # =============================================================================
# # Combine and save
# # =============================================================================
# 
# combined <- cowplot::ggdraw() +
#   cowplot::draw_plot(p_main, x = 0, y = 0, width = 1, height = 1) +
#   cowplot::draw_plot(
#     p_inset,
#     x = 0.205,
#     y = 0.125,
#     width = 0.27,
#     height = 0.25
#   )
# 
# ragg::agg_png(
#   filename = OUT_PNG,
#   width = WIDTH_IN,
#   height = HEIGHT_IN,
#   units = "in",
#   res = RES_MAIN,
#   background = "white"
# )
# print(combined)
# grDevices::dev.off()
# 
# .log("Wrote: %s", OUT_PNG)
# 
# if (WRITE_PDF) {
#   grDevices::cairo_pdf(
#     filename = OUT_PDF,
#     width = WIDTH_IN,
#     height = HEIGHT_IN,
#     family = "sans",
#     onefile = TRUE
#   )
#   print(combined)
#   grDevices::dev.off()
#   .log("Wrote: %s", OUT_PDF)
# }
# 
# .log("Done.")

#!/usr/bin/env Rscript

# scripts/f02_domain.R
# =============================================================================
# Study-domain figure
#
# Main map:
#   Southeast study domain filled by aggregated ecoregion, with state outlines
#   and state abbreviations.
#
# Inset:
#   CONUS with the Southeast study-domain bounding box outlined in red.
#
# PNG-only figure workflow. No intermediate CSV, Parquet, QS, PDF, or spatial
# files are written.
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(scales)
  library(ragg)
  library(RColorBrewer)
  library(tigris)
  library(fs)
})

options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)

# =============================================================================
# User settings
# =============================================================================

ROOT <- Sys.getenv(
  "CISC_ROOT",
  "/proj/sebastianlab/users/hquintal/cisc"
)

OUT_DIR <- file.path(ROOT, "figures", "f02_domain")
OUT_PNG <- file.path(OUT_DIR, "f02_domain.png")
OUT_PDF <- file.path(OUT_DIR, "f02_domain.pdf")

WRITE_PDF <- TRUE
CRS_PLOT <- 4326

BBOX_LONLAT <- c(
  xmin = -95,
  xmax = -75,
  ymin = 25,
  ymax = 40
)

RES_MAIN <- 300
WIDTH_IN <- 12
HEIGHT_IN <- 8.1

ECO_SHP <- file.path(
  ROOT,
  "data/in/shapefile/Aggr_Ecoregions_2015/Aggr_Ecoregions_2015.shp"
)

# =============================================================================
# Helpers
# =============================================================================

.log <- function(fmt, ...) {
  message(
    sprintf(
      "[%s] %s",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      sprintf(fmt, ...)
    )
  )
}

.make_bbox <- function(bbox_values) {
  sf::st_as_sfc(
    sf::st_bbox(
      c(
        xmin = as.numeric(bbox_values[["xmin"]]),
        ymin = as.numeric(bbox_values[["ymin"]]),
        xmax = as.numeric(bbox_values[["xmax"]]),
        ymax = as.numeric(bbox_values[["ymax"]])
      ),
      crs = sf::st_crs(4326)
    )
  )
}

.prepare_states <- function() {
  states <- tigris::states(cb = TRUE, year = 2022, class = "sf")
  drop_statefp <- c("02", "15", "60", "66", "69", "72", "78")
  states <- states[!(as.character(states$STATEFP) %in% drop_statefp), ]
  states <- sf::st_make_valid(states)
  states <- sf::st_transform(states, CRS_PLOT)
  states$abb <- as.character(states$STUSPS)
  states
}

.prepare_state_labels <- function(states_clip) {
  label_points <- states_clip |>
    sf::st_transform(5070) |>
    sf::st_point_on_surface() |>
    sf::st_transform(4326)
  
  label_points[
    !is.na(label_points$abb) &
      label_points$abb != "",
  ]
}

.build_ecoregion_palette <- function(levels) {
  levels <- as.character(levels)
  levels <- levels[!is.na(levels) & levels != ""]
  levels <- sort(unique(levels))
  
  if (!length(levels)) {
    stop("No ecoregion levels were found inside the study domain.")
  }
  
  base_palette <- RColorBrewer::brewer.pal(8, "Accent")
  colours <- if (length(levels) <= length(base_palette)) {
    base_palette[seq_along(levels)]
  } else {
    grDevices::colorRampPalette(base_palette)(length(levels))
  }
  
  stats::setNames(colours, levels)
}

.theme_main_map <- function() {
  ggplot2::theme_bw(base_size = 18, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans", size = 18),
      panel.grid.major = ggplot2::element_line(color = "grey82", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 14, color = "grey25"),
      axis.ticks = ggplot2::element_line(color = "grey25", linewidth = 0.3),
      # adjust
      legend.position = c(0.97, 0.40),
      # adjust
      legend.justification = c(1, 0.5),
      legend.direction = "vertical",
      legend.title = ggplot2::element_text(
        size = 13,
        hjust = 0.5
      ),
      legend.text = ggplot2::element_text(
        size = 10
      ),
      legend.key.width = grid::unit(0.42, "cm"),
      legend.key.height = grid::unit(0.42, "cm"),
      legend.spacing.y = grid::unit(0.25, "cm"),
      legend.margin = ggplot2::margin(4, 4, 4, 4),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.82),
        color = NA
      ),
      legend.box.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

# =============================================================================
# Resolve inputs
# =============================================================================

fs::dir_create(OUT_DIR, recurse = TRUE)

if (!file.exists(ECO_SHP)) {
  stop("Missing ecoregion shapefile: ", ECO_SHP)
}

.log("Ecoregion shapefile: %s", ECO_SHP)
.log("Output directory: %s", OUT_DIR)

bbox_poly <- .make_bbox(BBOX_LONLAT)

# =============================================================================
# Ecoregions
# =============================================================================

ecoregions_all <- sf::st_read(ECO_SHP, quiet = TRUE)

if (!("WSA9_NAME" %in% names(ecoregions_all))) {
  stop(
    "Ecoregion shapefile is missing WSA9_NAME. Available fields: ",
    paste(names(ecoregions_all), collapse = ", ")
  )
}

ecoregions_all <- ecoregions_all |>
  sf::st_make_valid() |>
  sf::st_transform(CRS_PLOT)

ecoregions_clip <- suppressWarnings(
  sf::st_intersection(ecoregions_all, bbox_poly)
)

ecoregions_clip <- ecoregions_clip[!sf::st_is_empty(ecoregions_clip), ]
ecoregions_clip$WSA9_NAME <- trimws(as.character(ecoregions_clip$WSA9_NAME))
ecoregions_clip <- ecoregions_clip[
  !is.na(ecoregions_clip$WSA9_NAME) & ecoregions_clip$WSA9_NAME != "",
]

if (!nrow(ecoregions_clip)) {
  stop("No ecoregions intersect the requested study-domain bounding box.")
}

# =============================================================================
# States and labels
# =============================================================================

states_all <- .prepare_states()
states_clip <- suppressWarnings(sf::st_intersection(states_all, bbox_poly))
states_clip <- states_clip[!sf::st_is_empty(states_clip), ]
state_labels <- .prepare_state_labels(states_clip)

# =============================================================================
# Palette
# =============================================================================

eco_levels <- sort(unique(ecoregions_clip$WSA9_NAME))
eco_palette <- .build_ecoregion_palette(eco_levels)

# =============================================================================
# Inset map
# =============================================================================

p_inset <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = states_all,
    fill = "white",
    color = "grey55",
    linewidth = 0.25
  ) +
  ggplot2::geom_sf(
    data = bbox_poly,
    fill = NA,
    color = "red",
    linewidth = 0.9
  ) +
  ggplot2::coord_sf(
    crs = sf::st_crs(4326),
    datum = NA,
    xlim = c(-125, -66),
    ylim = c(24, 50),
    expand = FALSE
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.8
    ),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

# =============================================================================
# Main map
# =============================================================================

p_main <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = ecoregions_clip,
    ggplot2::aes(fill = WSA9_NAME),
    color = "grey35",
    linewidth = 0.25
  ) +
  ggplot2::geom_sf(
    data = states_clip,
    fill = NA,
    color = "black",
    linewidth = 0.45
  ) +
  ggplot2::geom_sf_text(
    data = state_labels,
    ggplot2::aes(label = abb),
    size = 4.2,
    color = "grey20",
    family = "sans"
  ) +
  ggplot2::coord_sf(
    crs = sf::st_crs(4326),
    datum = sf::st_crs(4326),
    xlim = c(BBOX_LONLAT[["xmin"]], BBOX_LONLAT[["xmax"]]),
    ylim = c(BBOX_LONLAT[["ymin"]], BBOX_LONLAT[["ymax"]]),
    expand = FALSE,
    label_axes = list(bottom = "E", left = "N")
  ) +
  ggplot2::scale_fill_manual(
    values = eco_palette,
    limits = eco_levels,
    name = "Ecoregion",
    drop = FALSE
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      ncol = 1,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.5,
      keywidth = grid::unit(0.42, "cm"),
      keyheight = grid::unit(0.42, "cm")
    )
  ) +
  .theme_main_map()

# =============================================================================
# Combine and save
# =============================================================================

combined <- cowplot::ggdraw() +
  cowplot::draw_plot(p_main, x = 0, y = 0, width = 1, height = 1) +
  cowplot::draw_plot(
    p_inset,
    x = 0.18,
    y = 0.05,
    width = 0.27,
    height = 0.25
  )

ragg::agg_png(
  filename = OUT_PNG,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  res = RES_MAIN,
  background = "white"
)
print(combined)
grDevices::dev.off()

.log("Wrote: %s", OUT_PNG)


.log("Done.")