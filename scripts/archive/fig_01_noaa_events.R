# scripts/fig_01_noaa_events_from_definitive_countydays.R
# =============================================================================
# Plot NOAA county-day frequency maps from definitive county-day cache files
#
# Reads:
#   data/out/09_performance/noaa_county_days_definitive/noaa_county_days__all_heat.{qs,csv}
#   data/out/09_performance/noaa_county_days_definitive/noaa_county_days__all_precip.{qs,csv}
#
# Writes:
#   figures/fig_01_noaa_events/noaa_freq_allheat.png
#   figures/fig_01_noaa_events/noaa_freq_allprecip.png
#   figures/fig_01_noaa_events/noaa_freq_combo_heat_vs_precip.png
#   figures/fig_01_noaa_events/noaa_freq_panels_heat_vs_precip.png
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tigris)
  library(ragg)
  library(magick)
  library(patchwork)
  library(terra)
  library(qs)
})

options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)

# ----------------------------
# SHARED HELPERS (identical across scripts 1–5)
# ----------------------------
.norm_digits <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_extract(x, "\\d+")
  x[!nzchar(x)] <- NA_character_
  x
}
.norm_geoid5 <- function(x) {
  x <- .norm_digits(x)
  ok <- !is.na(x) & nchar(x) <= 5
  out <- rep(NA_character_, length(x))
  out[ok] <- stringr::str_pad(x[ok], 5, pad = "0")
  out
}
.norm_statefp2 <- function(x) {
  x <- .norm_digits(x)
  ok <- !is.na(x) & nchar(x) <= 2
  out <- rep(NA_character_, length(x))
  out[ok] <- stringr::str_pad(x[ok], 2, pad = "0")
  out
}
.norm_wfo <- function(x) {
  x <- toupper(as.character(x))
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "[^A-Z0-9]+", "")
  x[!nzchar(x)] <- NA_character_
  x
}
event_tag_norm <- function(x) {
  x <- tolower(as.character(x))
  x <- stringr::str_replace_all(x, "[^a-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  x
}
.log <- function(fmt, ...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))

# ----------------------------
# USER SETTINGS
# ----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
DEF_DIR <- file.path(ROOT, "data/out/09_performance/noaa_county_days_definitive")

WFO_SHP <- "V:/users/hquintal/phd1_cluster_southeast/data/input/regional_aggregation/nws_wfo_cws/w_05mr24/w_05mr24.shp"

OUT_BASE <- file.path(ROOT, "figures/fig_01_noaa_events")
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

CRS_ALBERS <- 5070

EXTENT_ANCHOR_NC    <- file.path(ROOT, "data/out/04_rp/heatindex_idf_hi3_return_periods.nc")
EXTENT_ANCHOR_LAYER <- 4

RES_MAIN <- 200
WIDTH_ALL  <- 2600
HEIGHT_ALL <- 1700
WIDTH_PANELS  <- 5200
HEIGHT_PANELS <- 1700

OUT_ALLHEAT      <- file.path(OUT_BASE, "noaa_freq_allheat.png")
OUT_ALLPRECIP    <- file.path(OUT_BASE, "noaa_freq_allprecip.png")
OUT_COMBO_STACK  <- file.path(OUT_BASE, "noaa_freq_combo_heat_vs_precip.png")
OUT_COMBO_PANELS <- file.path(OUT_BASE, "noaa_freq_panels_heat_vs_precip.png")

# ----------------------------
# Read definitive days and compute counts
# ----------------------------
.read_def_group <- function(name) {
  nm <- event_tag_norm(name)
  q <- file.path(DEF_DIR, sprintf("noaa_county_days__%s.qs", nm))
  c <- file.path(DEF_DIR, sprintf("noaa_county_days__%s.csv", nm))
  f <- if (file.exists(q)) q else c
  if (!file.exists(f)) stop("Missing definitive group file: ", f)
  
  x <- if (grepl("\\.qs$", f, ignore.case = TRUE)) {
    as.data.table(qread(f))
  } else {
    fread(f, showProgress = FALSE)
  }
  
  stopifnot(all(c("GEOID","date") %in% names(x)))
  x[, GEOID := .norm_geoid5(GEOID)]
  x[, date := as.Date(date)]
  x <- x[!is.na(GEOID) & !is.na(date)]
  x
}

all_heat   <- .read_def_group("all_heat")
all_precip <- .read_def_group("all_precip")

counts_allheat   <- all_heat[, .(freq = .N), by = GEOID]
counts_allprecip <- all_precip[, .(freq = .N), by = GEOID]

# ----------------------------
# GEOMETRIES
# ----------------------------
COUNTIES <- counties(cb = TRUE, year = 2022, class = "sf")
STATES   <- states(cb = TRUE, year = 2022, class = "sf")
drop_statefp <- c("02","15","72","60","66","69","78")
COUNTIES <- COUNTIES[!(COUNTIES$STATEFP %in% drop_statefp), ]
STATES   <- STATES[!(STATES$STATEFP %in% drop_statefp), ]
COUNTIES <- st_transform(COUNTIES, CRS_ALBERS)
STATES   <- st_transform(STATES, CRS_ALBERS)

if (!file.exists(WFO_SHP)) stop("Missing WFO shapefile: ", WFO_SHP)
WFO <- st_read(WFO_SHP, quiet = TRUE)
WFO <- st_make_valid(WFO)
WFO <- st_transform(WFO, CRS_ALBERS)
conus_union <- st_union(st_geometry(STATES))
WFO <- WFO[st_intersects(WFO, conus_union, sparse = FALSE), ]

# Join diagnostics (surface GEOID mismatches immediately)
.check_join <- function(count_dt, label) {
  miss <- setdiff(unique(count_dt$GEOID), COUNTIES$GEOID)
  message(sprintf("[JOIN CHECK] %s: %d unmatched GEOIDs", label, length(miss)))
  if (length(miss)) message("  Examples: ", paste(head(miss, 20), collapse=", "))
}
.check_join(counts_allheat, "all_heat")
.check_join(counts_allprecip, "all_precip")

# ----------------------------
# EXTENT ALIGNMENT
# ----------------------------
anchor_r <- terra::rast(EXTENT_ANCHOR_NC)[[EXTENT_ANCHOR_LAYER]]
anchor_r <- terra::project(anchor_r, paste0("EPSG:", CRS_ALBERS), method="near")
anchor_r <- terra::crop(anchor_r, terra::ext(st_bbox(STATES)))
ex <- terra::ext(anchor_r)
XLIM <- c(ex[1], ex[2]); YLIM <- c(ex[3], ex[4])

# ----------------------------
# Theme + plot
# ----------------------------
theme_map_bottom <- theme_void(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 20),
    legend.text  = element_text(size = 20),
    legend.key.width  = grid::unit(1.6, "in"),
    legend.key.height = grid::unit(0.18, "in")
  )

make_plot_from_counts <- function(count_dt, high_col, panel_tag=NULL) {
  geo <- COUNTIES %>% left_join(as.data.frame(count_dt), by="GEOID")
  ggplot() +
    geom_sf(data = geo, aes(fill = freq), color = NA) +
    geom_sf(data = WFO, fill = NA, color = "grey55", linewidth = 0.25) +
    geom_sf(data = STATES, fill = NA, color = "black", linewidth = 0.35) +
    scale_fill_gradient(low="grey85", high=high_col, na.value="white", name="County-days") +
    coord_sf(crs=st_crs(CRS_ALBERS), xlim=XLIM, ylim=YLIM, expand=FALSE) +
    theme_map_bottom +
    labs(tag=panel_tag) +
    theme(plot.tag = element_text(size=20, face="bold"),
          plot.tag.position = c(0.02, 0.98))
}

p_allheat   <- make_plot_from_counts(counts_allheat, "red",  "(a)")
p_allprecip <- make_plot_from_counts(counts_allprecip, "blue", "(b)")

ragg::agg_png(OUT_ALLHEAT, width=WIDTH_ALL, height=HEIGHT_ALL, res=RES_MAIN); print(p_allheat); dev.off()
ragg::agg_png(OUT_ALLPRECIP, width=WIDTH_ALL, height=HEIGHT_ALL, res=RES_MAIN); print(p_allprecip); dev.off()

img_combo <- magick::image_append(
  c(magick::image_read(OUT_ALLHEAT), magick::image_read(OUT_ALLPRECIP)),
  stack=TRUE
)
magick::image_write(img_combo, OUT_COMBO_STACK)

p_panels <- p_allheat + p_allprecip + plot_layout(ncol=2)
ragg::agg_png(OUT_COMBO_PANELS, width=WIDTH_PANELS, height=HEIGHT_PANELS, res=RES_MAIN); print(p_panels); dev.off()

message("Done.")
