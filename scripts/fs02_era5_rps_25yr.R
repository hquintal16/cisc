# Southeast ----
# scripts/fig_02_era5_return_periods_layer04_southeast.R
# =============================================================================
# ERA5 Return Period Intensity Maps (Layer 4) -- SOUTHEAST VERSION
#
# Updates:
# - plot in lon/lat (EPSG:4326), no Albers projection
# - keep raster as raster (no polygon conversion)
# - crop to Southeast lon/lat bbox
# - apply land mask in native lon/lat grid so ocean is NA
# - render NA as empty
# - 20 pt Arial font throughout
# - panel labels (a), (b)
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(tigris)
  library(ragg)
  library(magick)
  library(patchwork)
  library(grid)
  library(scales)
  library(RColorBrewer)
})

options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)

# ----------------------------
# INPUT FILES
# ----------------------------
HI_NC <- "/proj/sebastianlab/users/hquintal/cisc/data/out/04_rp/heatindex_idf_hi1_return_periods_bbox.nc"
PR_NC <- "/proj/sebastianlab/users/hquintal/cisc/data/out/04_rp/precip_idf_pr24_return_periods_bbox.nc"
# layer:return period
# 1 : 1
# 2 : 2
# 3 : 5
# 4 : 10
# 5 : 25
# 6 : 50
# 7 : 100
# 8 : 200
# 9 : 500
# 10: 1000

WFO_SHP <- "/proj/sebastianlab/users/hquintal/cisc/data/in/shapefile/nws_wfo_cws/w_05mr24/w_05mr24.shp"

LAYER_I <- 5
CRS_USE <- 4326
BBOX_LONLAT <- c(xmin = -95, xmax = -75, ymin = 24.75, ymax = 40)

# ----------------------------
# OUTPUT DIRECTORY + FILES
# ----------------------------
OUT_BASE <- "/proj/sebastianlab/users/hquintal/cisc/figures/fs02_era5_rps"
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

OUT_HI          <- file.path(OUT_BASE, "fs02_rp_25yr_heatindex_southeast.png")
OUT_PR          <- file.path(OUT_BASE, "fs02_rp_25yr_precipitation_southeast.png")
OUT_BOTH_STACK  <- file.path(OUT_BASE, "fs02_rp_25yr_both_southeast.png")
OUT_BOTH_PANELS <- file.path(OUT_BASE, "fs02_rp_25yr_panels_southeast.png")

# ----------------------------
# FIGURE SETTINGS
# ----------------------------
RES_MAIN <- 200

WIDTH_ALL  <- 1700
HEIGHT_ALL <- 1700

WIDTH_PANELS  <- 3800
HEIGHT_PANELS <- 1900

# ----------------------------
# HELPERS
# ----------------------------
.log <- function(fmt, ...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))
}

.clip_to_mask <- function(x, mask_poly) {
  x <- st_make_valid(x)
  mask_poly <- st_make_valid(mask_poly)
  suppressWarnings(st_intersection(x, mask_poly))
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
      legend.key.width  = grid::unit(0.9, "in"),
      legend.key.height = grid::unit(0.15, "in"),
      
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.justification = "center",
      
      plot.margin = margin(1, 1, 1, 1)
    )
}

# ----------------------------
# SOUTHEAST bbox polygon in lon/lat
# ----------------------------
bbox_poly <- st_as_sfc(st_bbox(BBOX_LONLAT, crs = st_crs(CRS_USE)))
bbox_poly <- st_make_valid(bbox_poly)

bb <- st_bbox(bbox_poly)

# ----------------------------
# STATES + WFO in lon/lat, clipped to Southeast bbox
# ----------------------------
STATES <- tigris::states(cb = TRUE, year = 2022, class = "sf")
drop_statefp <- c("02", "15", "72", "60", "66", "69", "78")
STATES <- STATES[!(STATES$STATEFP %in% drop_statefp), ] |>
  st_make_valid() |>
  st_transform(st_crs(bbox_poly))

STATES_crop <- .clip_to_mask(STATES, bbox_poly)

if (!file.exists(WFO_SHP)) stop("Missing WFO_SHP: ", WFO_SHP)
WFO <- st_read(WFO_SHP, quiet = TRUE) |>
  st_make_valid() |>
  st_transform(st_crs(bbox_poly))

WFO_crop <- .clip_to_mask(WFO, bbox_poly)

# ----------------------------
# LAND MASK
#   build from clipped states so ocean becomes NA
# ----------------------------
land_mask_vect <- STATES_crop |>
  st_union() |>
  st_make_valid() |>
  terra::vect()

# ----------------------------
# Raster -> data frame
#   crop -> land mask -> data frame
# ----------------------------
raster_layer_to_df_landmasked <- function(r, layer_i, mask_bb, land_mask_vect) {
  lyr <- r[[layer_i]]
  
  lyr <- terra::crop(
    lyr,
    terra::ext(mask_bb["xmin"], mask_bb["xmax"], mask_bb["ymin"], mask_bb["ymax"])
  )
  
  # make sure mask vector uses same CRS as raster
  if (!terra::same.crs(lyr, land_mask_vect)) {
    land_mask_vect <- terra::project(land_mask_vect, terra::crs(lyr))
  }
  
  land_mask_r <- terra::rasterize(
    land_mask_vect,
    lyr,
    field = 1,
    background = NA,
    touches = FALSE
  )
  
  lyr <- terra::mask(
    lyr,
    land_mask_r,
    maskvalues = NA,
    updatevalue = NA
  )
  
  df <- as.data.frame(lyr, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", "value")
  df
}

# ----------------------------
# Plot builder
# ----------------------------
make_rp_plot <- function(df, states_sf, wfo_sf, bbox_poly,
                         legend_title, fill_colours, panel_tag = NULL,
                         fill_limits = NULL) {
  
  bb <- st_bbox(bbox_poly)
  
  se_land <- st_union(states_sf)
  
  p <- ggplot() +
    geom_sf(data = se_land, fill = "white", color = NA) +
    geom_raster(
      data = df,
      aes(x = x, y = y, fill = value),
      interpolate = FALSE
    ) +
    geom_sf(data = wfo_sf, fill = NA, color = "grey55", linewidth = 0.25) +
    geom_sf(data = states_sf, fill = NA, color = "black", linewidth = 0.35) +
    coord_sf(
      crs = st_crs(4326),
      datum = st_crs(4326),
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE,
      label_axes = list(bottom = "E", left = "N")
    ) +
    labs(tag = panel_tag) +
    .theme_map() +
    theme(
      plot.tag = element_text(family = "Arial", size = 20, face = "bold"),
      plot.tag.position = c(0.02, 0.95),
      
      legend.position = "bottom",
      legend.box.just = "center",
      legend.justification = "center"
    )
  
  if (is.null(fill_limits)) {
    fill_limits <- range(df$value, na.rm = TRUE, finite = TRUE)
  }
  
  p <- p +
    scale_fill_stepsn(
      colours = fill_colours,
      limits = fill_limits,
      breaks = seq(
        fill_limits[1],
        fill_limits[2],
        length.out = 8
      ),
      labels = scales::label_number(
        accuracy = 1
      ),
      oob = scales::squish,
      na.value = NA,
      name = legend_title,
      guide = guide_colorsteps(
        barheight = grid::unit(0.18, "in"),
        barwidth = grid::unit(6.0, "in"),
        even.steps = TRUE,
        show.limits = TRUE
      )
    )
  
  p
}

# ----------------------------
# READ RASTERS
# ----------------------------
if (!file.exists(HI_NC)) stop("Missing HI netCDF: ", HI_NC)
if (!file.exists(PR_NC)) stop("Missing PR netCDF: ", PR_NC)

hi_r <- terra::rast(HI_NC)
pr_r <- terra::rast(PR_NC)

hi_df <- raster_layer_to_df_landmasked(
  r = hi_r,
  layer_i = LAYER_I,
  mask_bb = bb,
  land_mask_vect = land_mask_vect
)

pr_df <- raster_layer_to_df_landmasked(
  r = pr_r,
  layer_i = LAYER_I,
  mask_bb = bb,
  land_mask_vect = land_mask_vect
)

# ----------------------------
# BUILD PLOTS
# ----------------------------
p_hi <- make_rp_plot(
  df = hi_df,
  states_sf = STATES_crop,
  wfo_sf = WFO_crop,
  bbox_poly = bbox_poly,
  legend_title = "°C",
  fill_colours = RColorBrewer::brewer.pal(
    n = 9,
    name = "Reds"
  ),
  panel_tag = "(a)"
)

p_pr <- make_rp_plot(
  df = pr_df,
  states_sf = STATES_crop,
  wfo_sf = WFO_crop,
  bbox_poly = bbox_poly,
  legend_title = "mm",
  fill_colours = RColorBrewer::brewer.pal(
    n = 9,
    name = "Blues"
  ),
  panel_tag = "(b)"
)

# ----------------------------
# SAVE INDIVIDUAL FIGURES
# ----------------------------
ragg::agg_png(OUT_HI, width = WIDTH_ALL, height = HEIGHT_ALL, res = RES_MAIN)
print(p_hi)
dev.off()
.log("Wrote: %s", OUT_HI)

ragg::agg_png(OUT_PR, width = WIDTH_ALL, height = HEIGHT_ALL, res = RES_MAIN)
print(p_pr)
dev.off()
.log("Wrote: %s", OUT_PR)

# ----------------------------
# STACK WITHOUT RESIZING
# ----------------------------
img_hi <- magick::image_read(OUT_HI)
img_pr <- magick::image_read(OUT_PR)

info_hi <- magick::image_info(img_hi)
info_pr <- magick::image_info(img_pr)

if (info_hi$width != info_pr$width || info_hi$height != info_pr$height) {
  stop(
    "HI and PR PNG sizes differ: ",
    "hi=", info_hi$width, "x", info_hi$height,
    " pr=", info_pr$width, "x", info_pr$height
  )
}

img_both <- magick::image_append(c(img_hi, img_pr), stack = TRUE)
magick::image_write(img_both, OUT_BOTH_STACK)
.log("Wrote stacked: %s", OUT_BOTH_STACK)

# ----------------------------
# TRUE SIDE-BY-SIDE PANELS
#   one shared centered legend
# ----------------------------
p_panels <- (p_hi + p_pr) +
  plot_layout(ncol = 2) +
  plot_annotation(
    theme = theme(
      text = element_text(family = "Arial", size = 20),
      panel.spacing = unit(0.15, "lines")
    )
  )

ragg::agg_png(OUT_BOTH_PANELS, width = WIDTH_PANELS, height = HEIGHT_PANELS, res = RES_MAIN)
print(p_panels)
dev.off()
.log("Wrote panels: %s", OUT_BOTH_PANELS)

.log("Done.")
