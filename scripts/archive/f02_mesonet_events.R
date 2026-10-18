# # scripts/fig_01_wwa_two_panel_heat_vs_precip_union_restartsafe.R
# # =============================================================================
# # Build 3-panel WWA frequency figure:
# #   (a) heat union (deduped per zone-day)
# #   (b) precip union (deduped per zone-day)
# #   (c) CONUS inset showing Southeast bbox
# #
# # Updated to:
# # - use 20 pt Arial for all text
# # - avoid Albers projection and plot in lon/lat (EPSG:4326)
# # - add panel tags (a), (b), (c)
# # =============================================================================
# 
# suppressPackageStartupMessages({
#   library(data.table)
#   library(sf)
#   library(dplyr)
#   library(ggplot2)
#   library(ragg)
#   library(tigris)
#   library(fs)
#   library(qs)
#   library(magick)
#   library(patchwork)
#   library(scales)
#   library(grid)
#   library(viridis)
# })
# 
# options(tigris_use_cache = TRUE)
# sf::sf_use_s2(FALSE)
# 
# # ----------------------------
# # USER SETTINGS
# # ----------------------------
# ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
# 
# OUT_DIR <- file.path(ROOT, "figures/fig_02_wwa")
# 
# ZONE_SHP <- file.path(ROOT, "data/in/shapefile/nws_public_forecast_zone/z_18mr25.shp")
# WFO_SHP  <- "V:/users/hquintal/phd1_cluster_southeast/data/input/regional_aggregation/nws_wfo_cws/w_05mr24/w_05mr24.shp"
# 
# BBOX_LONLAT <- c(xmin = -95, xmax = -75, ymin = 25, ymax = 40)
# CRS_PLOT <- 4326
# 
# RES_MAIN <- 200
# WIDTH_ALL  <- 1700
# HEIGHT_ALL <- 1700
# WIDTH_PANELS  <- 3800
# HEIGHT_PANELS <- 2600
# 
# REBUILD_COMBOS <- TRUE
# 
# # Heat (left)
# HEAT_PHENOM <- c("EH")
# HEAT_SIG    <- c("W")
# 
# # Precip (middle)
# PREC_PHENOM <- c("SV", "FF", "HI", "HU", "TI", "TR")
# PREC_SIG    <- c("W")
# 
# # Outputs
# OUT_HEAT_PNG <- file.path(OUT_DIR, "wwa_freq_combo_heat_AW.png")
# OUT_HEAT_QS  <- file.path(OUT_DIR, "wwa_freq_combo_heat_AW.qs")
# 
# OUT_PREC_PNG <- file.path(OUT_DIR, "wwa_freq_combo_precip_AW.png")
# OUT_PREC_QS  <- file.path(OUT_DIR, "wwa_freq_combo_precip_AW.qs")
# 
# OUT_PANELS   <- file.path(OUT_DIR, "wwa_freq_panels_heat_vs_precip.png")
# OUT_INSET_PNG <- file.path(OUT_DIR, "wwa_freq_inset_conus.png")
# OUT_FINAL_PNG <- file.path(OUT_DIR, "wwa_freq_panels_heat_vs_precip_with_inset.png")
# 
# # ----------------------------
# # Helpers
# # ----------------------------
# .log <- function(fmt, ...) {
#   message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))
# }
# 
# .auto_find_out_dir <- function(root, out_dir) {
#   if (dir_exists(out_dir)) {
#     f0 <- dir_ls(out_dir, recurse = TRUE, type = "file", regexp = "wwa_freq_.*\\.qs$")
#     if (length(f0)) return(out_dir)
#   }
#   
#   candidates <- c(
#     file.path(root, "figures"),
#     file.path(root, "figures/fig_01_wwa"),
#     file.path(root, "figures/fig_01_wwa_events"),
#     file.path(root, "figures/fig_01_noaa_events")
#   )
#   candidates <- candidates[dir_exists(candidates)]
#   
#   for (d in candidates) {
#     f <- dir_ls(d, recurse = TRUE, type = "file", regexp = "wwa_freq_.*\\.qs$")
#     if (length(f)) {
#       .log("[FIX] Found wwa_freq_*.qs under: %s", d)
#       return(d)
#     }
#   }
#   
#   stop("Could not find any wwa_freq_*.qs under ROOT/figures. Check drive mapping or output location.")
# }
# 
# idx_qs <- function(out_dir) {
#   qs_files <- dir_ls(out_dir, recurse = TRUE, type = "file", regexp = "wwa_freq_.*\\.qs$")
#   if (!length(qs_files)) {
#     return(data.table(file = character(), fname = character(), PHENOM = character(), SIG = character()))
#   }
#   
#   dt <- data.table(file = qs_files)
#   dt[, fname := path_file(file)]
#   
#   toks <- strsplit(dt$fname, "_", fixed = TRUE)
#   n_tok <- vapply(toks, length, integer(1))
#   dt <- dt[n_tok >= 4]
#   toks <- toks[n_tok >= 4]
#   
#   dt[, PHENOM := toupper(vapply(toks, `[`, character(1), 3))]
#   dt[, SIG    := toupper(vapply(toks, `[`, character(1), 4))]
#   
#   dt <- dt[nchar(PHENOM) == 2 & nchar(SIG) == 1]
#   dt[, .(file, fname, PHENOM, SIG)]
# }
# 
# .read_union_ts <- function(index_dt, phenoms, sigs) {
#   phenoms <- toupper(phenoms)
#   sigs <- toupper(sigs)
#   wanted <- CJ(PHENOM = phenoms, SIG = sigs, unique = TRUE)
#   
#   comp <- index_dt[PHENOM %in% phenoms & SIG %in% sigs]
#   have <- unique(comp[, .(PHENOM, SIG)])
#   missing <- wanted[!have, on = .(PHENOM, SIG)]
#   
#   if (nrow(missing)) {
#     .log(
#       "[WARN] Missing %d component TS files: %s",
#       nrow(missing),
#       paste(sprintf("%s.%s", missing$PHENOM, missing$SIG), collapse = ", ")
#     )
#   }
#   
#   if (!nrow(comp)) {
#     .log("[WARN] No component TS files found for requested union; returning EMPTY ts.")
#     return(data.table(STATE_ZONE = character(), date = as.Date(character()), present = integer()))
#   }
#   
#   ts_list <- lapply(comp$file, function(f) {
#     x <- tryCatch(as.data.table(qread(f)), error = function(e) NULL)
#     if (is.null(x)) return(NULL)
#     if (!all(c("STATE_ZONE", "date") %in% names(x))) return(NULL)
#     x <- x[, .(STATE_ZONE, date)]
#     x[, STATE_ZONE := toupper(trimws(as.character(STATE_ZONE)))]
#     x[, date := as.Date(date)]
#     x
#   })
#   
#   ts_list <- Filter(Negate(is.null), ts_list)
#   if (!length(ts_list)) {
#     return(data.table(STATE_ZONE = character(), date = as.Date(character()), present = integer()))
#   }
#   
#   ts <- rbindlist(ts_list, use.names = TRUE, fill = TRUE)
#   if (!nrow(ts)) {
#     return(data.table(STATE_ZONE = character(), date = as.Date(character()), present = integer()))
#   }
#   
#   setkey(ts, STATE_ZONE, date)
#   ts <- unique(ts)
#   ts[, present := 1L]
#   ts
# }
# 
# .freq_from_ts <- function(ts) {
#   if (is.null(ts) || !nrow(ts)) {
#     return(data.table(STATE_ZONE = character(), freq_days = integer()))
#   }
#   ts[, .(freq_days = .N), by = STATE_ZONE]
# }
# 
# .clip_to_mask <- function(x, mask_poly) {
#   x <- st_make_valid(x)
#   suppressWarnings(st_intersection(x, mask_poly))
# }
# 
# .make_bbox_poly <- function(bbox_lonlat, crs_out = 4326) {
#   bb_ll <- st_as_sfc(st_bbox(bbox_lonlat, crs = st_crs(4326)))
#   if (!is.null(crs_out) && crs_out != 4326) {
#     bb_ll <- st_transform(bb_ll, crs_out)
#   }
#   st_make_valid(bb_ll)
# }
# 
# .theme_map <- function() {
#   theme_minimal(base_size = 20, base_family = "Arial") +
#     theme(
#       text = element_text(family = "Arial", size = 20),
#       panel.grid.major = element_line(color = "grey70", linewidth = 0.35),
#       panel.grid.minor = element_blank(),
#       panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
#       axis.title = element_blank(),
#       axis.text = element_text(size = 20, family = "Arial", color = "grey25"),
#       axis.ticks = element_line(color = "grey25", linewidth = 0.3),
#       legend.position = "bottom",
#       legend.direction = "horizontal",
#       legend.title = element_text(size = 20, family = "Arial"),
#       legend.text  = element_text(size = 20, family = "Arial"),
#       legend.key.width  = grid::unit(0.9, "in"),
#       legend.key.height = grid::unit(0.15, "in"),
#       legend.box = "horizontal",
#       legend.box.just = "center",
#       legend.justification = "center",
#       plot.tag = element_text(size = 20, family = "Arial", face = "bold"),
#       plot.tag.position = c(0.02, 0.95),
#       plot.margin = margin(1, 1, 1, 1)
#     )
# }
# 
# .make_plot <- function(zones_sf, states_sf, wfo_sf, bbox_poly, fill_limits,
#                        legend_title = "Count", panel_tag = NULL) {
#   bb <- st_bbox(bbox_poly)
#   se_land <- st_union(states_sf)
#   ggplot() +
#     geom_sf(data = se_land, fill = "white", color = NA) +
#     geom_sf(data = zones_sf, aes(fill = freq_days), color = NA) +
#     geom_sf(data = wfo_sf, fill = NA, color = "grey55", linewidth = 0.25) +
#     geom_sf(data = states_sf, fill = NA, color = "black", linewidth = 0.35) +
#     scale_fill_viridis_b(
#       option = "viridis",
#       limits = fill_limits,
#       breaks = seq(fill_limits[1], fill_limits[2], length.out = 8),
#       labels = scales::label_number(accuracy = 1),
#       oob = scales::squish,
#       na.value = NA,
#       name = legend_title,
#       guide = guide_colorsteps(
#         direction = "horizontal",
#         title.position = "top",
#         title.hjust = 0.5,
#         barwidth = grid::unit(8.0, "in"),
#         barheight = grid::unit(0.22, "in"),
#         even.steps = TRUE,
#         show.limits = TRUE
#       )
#     ) +
#     coord_sf(
#       crs   = st_crs(CRS_PLOT),
#       datum = st_crs(4326),
#       xlim  = c(bb["xmin"], bb["xmax"]),
#       ylim  = c(bb["ymin"], bb["ymax"]),
#       expand = FALSE,
#       label_axes = list(bottom = "E", left = "N")
#     ) +
#     labs(tag = panel_tag) +
#     .theme_map()
# }
# 
# .make_conus_inset <- function(states_conus, se_bbox, panel_tag = NULL) {
#   ggplot() +
#     geom_sf(data = states_conus, fill = "grey95", color = "grey70", linewidth = 0.50) +
#     geom_sf(data = se_bbox, fill = NA, color = "red", linewidth = 0.85) +
#     coord_sf(
#       crs = st_crs(CRS_PLOT),
#       datum = NA,
#       expand = FALSE
#     ) +
#     labs(tag = panel_tag) +
#     theme_void(base_family = "Arial") +
#     theme(
#       text = element_text(family = "Arial", size = 20),
#       plot.tag = element_text(size = 20, family = "Arial", face = "bold"),
#       plot.tag.position = c(0.02, 0.98),
#       legend.position = "none",
#       panel.background = element_rect(fill = "white", color = NA),
#       plot.background = element_rect(fill = "white", color = NA),
#       plot.margin = margin(5, 5, 5, 5)
#     )
# }
# 
# # ----------------------------
# # Resolve OUT_DIR to where the per-code .qs actually live
# # ----------------------------
# OUT_DIR <- .auto_find_out_dir(ROOT, OUT_DIR)
# dir_create(OUT_DIR, recurse = TRUE)
# 
# OUT_HEAT_PNG <- file.path(OUT_DIR, "wwa_freq_combo_heat_AW.png")
# OUT_HEAT_QS  <- file.path(OUT_DIR, "wwa_freq_combo_heat_AW.qs")
# OUT_PREC_PNG <- file.path(OUT_DIR, "wwa_freq_combo_precip_AW.png")
# OUT_PREC_QS  <- file.path(OUT_DIR, "wwa_freq_combo_precip_AW.qs")
# OUT_PANELS   <- file.path(OUT_DIR, "wwa_freq_panels_heat_vs_precip.png")
# 
# # ----------------------------
# # BBOX mask polygon in lon/lat
# # ----------------------------
# bbox_poly <- .make_bbox_poly(BBOX_LONLAT, CRS_PLOT)
# 
# # ----------------------------
# # Geometries transformed + TRUE polygon clip in lon/lat
# # ----------------------------
# if (!file.exists(ZONE_SHP)) stop("Missing ZONE_SHP: ", ZONE_SHP)
# zones <- st_read(ZONE_SHP, quiet = TRUE) |> st_make_valid()
# zones$STATE_ZONE <- toupper(trimws(as.character(zones$STATE_ZONE)))
# zones <- st_transform(zones, CRS_PLOT)
# zones_crop <- .clip_to_mask(zones, bbox_poly)
# 
# STATES <- tigris::states(cb = TRUE, year = 2022, class = "sf")
# drop_statefp <- c("02", "15", "72", "60", "66", "69", "78")
# STATES <- STATES[!(STATES$STATEFP %in% drop_statefp), ] |> st_make_valid()
# STATES <- st_transform(STATES, CRS_PLOT)
# STATES_crop <- .clip_to_mask(STATES, bbox_poly)
# 
# if (!file.exists(WFO_SHP)) stop("Missing WFO_SHP: ", WFO_SHP)
# WFO <- st_read(WFO_SHP, quiet = TRUE) |> st_make_valid()
# WFO <- st_transform(WFO, CRS_PLOT)
# WFO_crop <- .clip_to_mask(WFO, bbox_poly)
# 
# # ----------------------------
# # CONUS states + Southeast bbox for inset
# # ----------------------------
# STATES_conus <- tigris::states(cb = TRUE, year = 2022, class = "sf")
# STATES_conus <- STATES_conus[!(STATES_conus$STATEFP %in% drop_statefp), ] |>
#   st_make_valid()
# STATES_conus <- st_transform(STATES_conus, CRS_PLOT)
# 
# SE_BBOX <- st_as_sf(
#   data.frame(id = 1),
#   geometry = st_geometry(bbox_poly)
# )
# 
# # ----------------------------
# # Index per-code TS .qs files
# # ----------------------------
# index_dt <- idx_qs(OUT_DIR)
# .log("OUT_DIR: %s", OUT_DIR)
# .log("Indexed per-code TS files: %d", nrow(index_dt))
# if (nrow(index_dt)) {
#   .log("Available SIG: %s", paste(sort(unique(index_dt$SIG)), collapse = ", "))
#   .log("Available PHENOM (sample): %s", paste(head(sort(unique(index_dt$PHENOM)), 40), collapse = ", "))
# }
# 
# req_pairs <- rbind(
#   CJ(PHENOM = toupper(HEAT_PHENOM), SIG = toupper(HEAT_SIG), group = "heat", unique = TRUE),
#   CJ(PHENOM = toupper(PREC_PHENOM), SIG = toupper(PREC_SIG), group = "precip", unique = TRUE)
# )
# req_pairs[, exists := (paste(PHENOM, SIG) %chin% paste(index_dt$PHENOM, index_dt$SIG))]
# .log("Requested pair availability:")
# print(req_pairs[order(group, PHENOM, SIG)])
# 
# # ----------------------------
# # Heat union: build or read
# # ----------------------------
# if (REBUILD_COMBOS || !file.exists(OUT_HEAT_QS)) {
#   ts_heat <- .read_union_ts(index_dt, HEAT_PHENOM, HEAT_SIG)
#   qsave(ts_heat, OUT_HEAT_QS)
#   .log("Wrote: %s  (rows=%d)", OUT_HEAT_QS, nrow(ts_heat))
# } else {
#   ts_heat <- as.data.table(qread(OUT_HEAT_QS))
#   .log("Loaded: %s  (rows=%d)", OUT_HEAT_QS, nrow(ts_heat))
# }
# 
# freq_heat <- .freq_from_ts(ts_heat)
# .log("Heat union: zones=%d | total zone-days=%d", length(unique(freq_heat$STATE_ZONE)), nrow(ts_heat))
# 
# zones_heat <- zones_crop %>%
#   left_join(as.data.frame(freq_heat), by = c("STATE_ZONE" = "STATE_ZONE"))
# 
# # ----------------------------
# # Precip union: build or read
# # ----------------------------
# if (REBUILD_COMBOS || !file.exists(OUT_PREC_QS)) {
#   ts_prec <- .read_union_ts(index_dt, PREC_PHENOM, PREC_SIG)
#   qsave(ts_prec, OUT_PREC_QS)
#   .log("Wrote: %s  (rows=%d)", OUT_PREC_QS, nrow(ts_prec))
# } else {
#   ts_prec <- as.data.table(qread(OUT_PREC_QS))
#   .log("Loaded: %s  (rows=%d)", OUT_PREC_QS, nrow(ts_prec))
# }
# 
# freq_prec <- .freq_from_ts(ts_prec)
# .log("Precip union: zones=%d | total zone-days=%d", length(unique(freq_prec$STATE_ZONE)), nrow(ts_prec))
# 
# zones_prec <- zones_crop %>%
#   left_join(as.data.frame(freq_prec), by = c("STATE_ZONE" = "STATE_ZONE"))
# 
# # ----------------------------
# # Shared fill limits
# # ----------------------------
# max_heat <- if ("freq_days" %in% names(zones_heat)) {
#   suppressWarnings(max(zones_heat$freq_days, na.rm = TRUE))
# } else {
#   0
# }
# max_prec <- if ("freq_days" %in% names(zones_prec)) {
#   suppressWarnings(max(zones_prec$freq_days, na.rm = TRUE))
# } else {
#   0
# }
# fill_limits <- c(0, max(max_heat, max_prec, 1, na.rm = TRUE))
# 
# # ----------------------------
# # Build three panel plots
# # ----------------------------
# p_heat <- .make_plot(
#   zones_heat,
#   STATES_crop,
#   WFO_crop,
#   bbox_poly,
#   fill_limits,
#   legend_title = "Count",
#   panel_tag = "(b)"
# )
# 
# p_prec <- .make_plot(
#   zones_prec,
#   STATES_crop,
#   WFO_crop,
#   bbox_poly,
#   fill_limits,
#   legend_title = "Count",
#   panel_tag = "(c)"
# )
# 
# p_inset <- .make_conus_inset(
#   STATES_conus,
#   SE_BBOX,
#   panel_tag = "(a)"
# )
# 
# # ----------------------------
# # Combine figure
# #   top row: centered inset
# #   bottom row: heat + precip
# # ----------------------------
# lay <- c(
#   patchwork::area(t = 1, l = 2, b = 1, r = 3),  # inset centered on top
#   patchwork::area(t = 2, l = 1, b = 2, r = 2),  # heat bottom left
#   patchwork::area(t = 2, l = 3, b = 2, r = 4)   # precip bottom right
# )
# 
# p_panels <- p_inset + p_heat + p_prec +
#   plot_layout(
#     design = lay,
#     widths = c(1, 1, 1, 1),
#     heights = c(0.38, 1),
#     guides = "collect"
#   ) +
#   plot_annotation(
#     theme = theme(
#       text = element_text(family = "Arial", size = 20),
#       legend.position = "bottom",
#       legend.box.just = "center",
#       legend.justification = "center",
#       panel.spacing = unit(0.2, "lines")
#     )
#   )
# 
# # ----------------------------
# # Save final combined figure
# # ----------------------------
# ragg::agg_png(
#   OUT_FINAL_PNG,
#   width = WIDTH_PANELS,
#   height = HEIGHT_PANELS,
#   res = RES_MAIN
# )
# print(p_panels)
# dev.off()
# 
# .log("Wrote final combined figure: %s", OUT_FINAL_PNG)
# .log("Done.")

# Study area ----

# scripts/fig_01_wwa_ecoregions_main_with_conus_inset.R
# =============================================================================
# Main map:
#   Southeast bbox region, filled by ecoregion, with state outlines and
#   state abbreviations
#
# Inset:
#   Full CONUS in bottom-left with red Southeast bbox outline
#
# No panel tag
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(maps)
  library(scales)
  library(ragg)
  library(RColorBrewer)
})

sf::sf_use_s2(FALSE)

# ----------------------------
# USER SETTINGS
# ----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
OUT_DIR <- file.path(ROOT, "figures/fig_02_wwa")

CRS_PLOT <- 4326
RES_MAIN <- 300

WIDTH_MAIN  <- 3200
HEIGHT_MAIN <- 2400

BBOX_LONLAT <- c(xmin = -95, xmax = -75, ymin = 25, ymax = 40)

ECO_SHP <- "V:/users/hquintal/phd1_cluster_southeast/data/input/regional_aggregation/ecoregion/Aggr_Ecoregions_2015/Aggr_Ecoregions_2015.shp"

OUT_PNG <- file.path(OUT_DIR, "wwa_ecoregions_main_with_conus_inset.png")

# ----------------------------
# Helpers
# ----------------------------
.log <- function(fmt, ...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))
}

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ----------------------------
# Bounding box
# ----------------------------
bbox_extent <- st_bbox(
  c(
    xmin = as.numeric(BBOX_LONLAT["xmin"]),
    ymin = as.numeric(BBOX_LONLAT["ymin"]),
    xmax = as.numeric(BBOX_LONLAT["xmax"]),
    ymax = as.numeric(BBOX_LONLAT["ymax"])
  ),
  crs = st_crs(4326)
)

bbox_poly <- st_as_sfc(bbox_extent)

# ----------------------------
# Read and prepare ecoregions
# ----------------------------
if (!file.exists(ECO_SHP)) stop("Missing ECO_SHP: ", ECO_SHP)

ecoregions_all <- st_read(ECO_SHP, quiet = TRUE) %>%
  st_transform(CRS_PLOT) %>%
  st_make_valid()

ecoregions_clip <- st_crop(ecoregions_all, bbox_extent)

# ----------------------------
# State outlines
# ----------------------------
states_sf <- maps::map("state", fill = TRUE, plot = FALSE) %>%
  st_as_sf() %>%
  st_transform(CRS_PLOT) %>%
  st_make_valid() %>%
  mutate(abb = state.abb[match(ID, tolower(state.name))])

states_clip <- suppressWarnings(st_intersection(states_sf, bbox_poly))
state_centroids <- suppressWarnings(st_centroid(states_clip))

# ----------------------------
# Ecoregion palette
# ----------------------------
eco_levels <- sort(unique(ecoregions_clip$WSA9_NAME))
eco_pal <- brewer.pal(max(3, min(8, length(eco_levels))), "Accent")
if (length(eco_levels) > length(eco_pal)) {
  eco_pal <- colorRampPalette(brewer.pal(8, "Accent"))(length(eco_levels))
}
names(eco_pal) <- eco_levels

# ----------------------------
# Inset map: full CONUS + red bbox
# ----------------------------
p_inset <- ggplot() +
  geom_sf(
    data = states_sf,
    fill = "white",
    color = "grey55",
    linewidth = 0.30
  ) +
  geom_rect(
    aes(
      xmin = bbox_extent["xmin"],
      xmax = bbox_extent["xmax"],
      ymin = bbox_extent["ymin"],
      ymax = bbox_extent["ymax"]
    ),
    fill = NA,
    color = "red",
    linewidth = 0.9
  ) +
  coord_sf(
    crs = st_crs(CRS_PLOT),
    datum = NA,
    xlim = c(-125, -66),
    ylim = c(24, 50),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
    plot.background  = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

# ----------------------------
# Main map: ecoregions + black states + abbreviations
# ----------------------------
p_main <- ggplot() +
  geom_sf(
    data = ecoregions_clip,
    aes(fill = WSA9_NAME),
    color = "grey35",
    linewidth = 0.25
  ) +
  geom_sf(
    data = states_clip,
    fill = NA,
    color = "black",
    linewidth = 0.45
  ) +
  geom_sf_text(
    data = state_centroids,
    aes(label = abb),
    size = 4.2,
    colour = "grey20"
  ) +
  coord_sf(
    xlim = c(bbox_extent["xmin"], bbox_extent["xmax"]),
    ylim = c(bbox_extent["ymin"], bbox_extent["ymax"]),
    expand = FALSE,
    datum = st_crs(4326),
    label_axes = list(bottom = "E", left = "N")
  ) +
  scale_fill_manual(
    values = eco_pal,
    name = "Ecoregion",
    drop = FALSE
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  ) +
  theme_bw(base_size = 20) +
  theme(
    text = element_text(size = 20),
    panel.grid.major = element_line(color = "grey82", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_blank(),
    axis.text = element_text(size = 16, color = "grey25"),
    axis.ticks = element_line(color = "grey25", linewidth = 0.3),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 11),
    legend.background = element_rect(
      fill = alpha("white", 0.95),
      color = "black",
      linewidth = 0.25
    ),
    plot.margin = margin(12, 12, 12, 12)
  )

# ----------------------------
# Combine main + inset
#   move inset upward/right so it does not cover x/y axes
# ----------------------------
combined <- ggdraw() +
  draw_plot(p_main) +
  draw_plot(
    p_inset,
    x = 0.208,
    y = 0.18,
    width = 0.255,
    height = 0.255
  )

# ----------------------------
# Save
# ----------------------------
ragg::agg_png(
  OUT_PNG,
  width = WIDTH_MAIN,
  height = HEIGHT_MAIN,
  res = RES_MAIN
)
print(combined)
dev.off()

.log("Wrote: %s", OUT_PNG)
.log("Done.")

# Warning counts ----

# scripts/fig_01_wwa_heat_precip_two_panel.R
# =============================================================================
# Build 2-panel Southeast WWA frequency figure:
#   (a) heat union
#   (b) precip union
#
# Standalone Southeast count-map workflow only.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(ragg)
  library(tigris)
  library(fs)
  library(qs)
  library(patchwork)
  library(scales)
  library(grid)
  library(viridis)
})

options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)

# ----------------------------
# USER SETTINGS
# ----------------------------
ROOT <- "/proj/sebastianlab/users/hquintal/cisc"

OUT_DIR <- file.path(ROOT, "figures/fig_02_wwa")

ZONE_SHP <- file.path(ROOT, "data/in/shapefile/nws_public_forecast_zone/z_18mr25.shp")
WFO_SHP  <- "V:/users/hquintal/phd1_cluster_southeast/data/input/regional_aggregation/nws_wfo_cws/w_05mr24/w_05mr24.shp"

BBOX_LONLAT <- c(xmin = -95, xmax = -75, ymin = 25, ymax = 40)
CRS_PLOT <- 4326

RES_MAIN <- 200
WIDTH_PANELS  <- 3800
HEIGHT_PANELS <- 1900

REBUILD_COMBOS <- TRUE

# Heat (left)
HEAT_PHENOM <- c("EH")
HEAT_SIG    <- c("W")

# Precip (right)
PREC_PHENOM <- c("FF")
PREC_SIG    <- c("W")

# Outputs
OUT_HEAT_QS    <- file.path(OUT_DIR, "wwa_freq_combo_heat_AW.qs")
OUT_PREC_QS    <- file.path(OUT_DIR, "wwa_freq_combo_precip_AW.qs")
OUT_PANELS_PNG <- file.path(OUT_DIR, "wwa_freq_panels_heat_vs_precip.png")

# ----------------------------
# Helpers
# ----------------------------
.log <- function(fmt, ...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))
}

.auto_find_out_dir <- function(root, out_dir) {
  if (dir_exists(out_dir)) {
    f0 <- dir_ls(out_dir, recurse = TRUE, type = "file", regexp = "wwa_freq_.*\\.qs$")
    if (length(f0)) return(out_dir)
  }
  
  candidates <- c(
    file.path(root, "figures"),
    file.path(root, "figures/fig_01_wwa"),
    file.path(root, "figures/fig_01_wwa_events"),
    file.path(root, "figures/fig_01_noaa_events")
  )
  candidates <- candidates[dir_exists(candidates)]
  
  for (d in candidates) {
    f <- dir_ls(d, recurse = TRUE, type = "file", regexp = "wwa_freq_.*\\.qs$")
    if (length(f)) {
      .log("[FIX] Found wwa_freq_*.qs under: %s", d)
      return(d)
    }
  }
  
  stop("Could not find any wwa_freq_*.qs under ROOT/figures. Check drive mapping or output location.")
}

idx_qs <- function(out_dir) {
  qs_files <- dir_ls(out_dir, recurse = TRUE, type = "file", regexp = "wwa_freq_.*\\.qs$")
  if (!length(qs_files)) {
    return(data.table(file = character(), fname = character(), PHENOM = character(), SIG = character()))
  }
  
  dt <- data.table(file = qs_files)
  dt[, fname := path_file(file)]
  
  toks <- strsplit(dt$fname, "_", fixed = TRUE)
  n_tok <- vapply(toks, length, integer(1))
  dt <- dt[n_tok >= 4]
  toks <- toks[n_tok >= 4]
  
  dt[, PHENOM := toupper(vapply(toks, `[`, character(1), 3))]
  dt[, SIG    := toupper(vapply(toks, `[`, character(1), 4))]
  
  dt <- dt[nchar(PHENOM) == 2 & nchar(SIG) == 1]
  dt[, .(file, fname, PHENOM, SIG)]
}

.read_union_ts <- function(index_dt, phenoms, sigs) {
  phenoms <- toupper(phenoms)
  sigs <- toupper(sigs)
  wanted <- CJ(PHENOM = phenoms, SIG = sigs, unique = TRUE)
  
  comp <- index_dt[PHENOM %in% phenoms & SIG %in% sigs]
  have <- unique(comp[, .(PHENOM, SIG)])
  missing <- wanted[!have, on = .(PHENOM, SIG)]
  
  if (nrow(missing)) {
    .log(
      "[WARN] Missing %d component TS files: %s",
      nrow(missing),
      paste(sprintf("%s.%s", missing$PHENOM, missing$SIG), collapse = ", ")
    )
  }
  
  if (!nrow(comp)) {
    .log("[WARN] No component TS files found for requested union; returning EMPTY ts.")
    return(data.table(STATE_ZONE = character(), date = as.Date(character()), present = integer()))
  }
  
  ts_list <- lapply(comp$file, function(f) {
    x <- tryCatch(as.data.table(qread(f)), error = function(e) NULL)
    if (is.null(x)) return(NULL)
    if (!all(c("STATE_ZONE", "date") %in% names(x))) return(NULL)
    x <- x[, .(STATE_ZONE, date)]
    x[, STATE_ZONE := toupper(trimws(as.character(STATE_ZONE)))]
    x[, date := as.Date(date)]
    x
  })
  
  ts_list <- Filter(Negate(is.null), ts_list)
  if (!length(ts_list)) {
    return(data.table(STATE_ZONE = character(), date = as.Date(character()), present = integer()))
  }
  
  ts <- rbindlist(ts_list, use.names = TRUE, fill = TRUE)
  if (!nrow(ts)) {
    return(data.table(STATE_ZONE = character(), date = as.Date(character()), present = integer()))
  }
  
  setkey(ts, STATE_ZONE, date)
  ts <- unique(ts)
  ts[, present := 1L]
  ts
}

.freq_from_ts <- function(ts) {
  if (is.null(ts) || !nrow(ts)) {
    return(data.table(STATE_ZONE = character(), freq_days = integer()))
  }
  ts[, .(freq_days = .N), by = STATE_ZONE]
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
      plot.tag = element_text(size = 20, family = "Arial", face = "bold"),
      plot.tag.position = c(0.02, 0.95),
      plot.margin = margin(1, 1, 1, 1)
    )
}

.make_plot <- function(zones_sf, states_sf, wfo_sf, bbox_poly, fill_limits,
                       legend_title = "Count", panel_tag = NULL) {
  bb <- st_bbox(bbox_poly)
  se_land <- st_union(states_sf)
  
  ggplot() +
    geom_sf(data = se_land, fill = "white", color = NA) +
    geom_sf(data = zones_sf, aes(fill = freq_days), color = NA) +
    geom_sf(data = wfo_sf, fill = NA, color = "grey55", linewidth = 0.25) +
    geom_sf(data = states_sf, fill = NA, color = "black", linewidth = 0.35) +
    scale_fill_viridis_b(
      option = "viridis",
      limits = fill_limits,
      breaks = seq(fill_limits[1], fill_limits[2], length.out = 8),
      labels = scales::label_number(accuracy = 1),
      oob = scales::squish,
      na.value = NA,
      name = legend_title,
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
}

# ----------------------------
# Resolve OUT_DIR
# ----------------------------
OUT_DIR <- .auto_find_out_dir(ROOT, OUT_DIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_HEAT_QS    <- file.path(OUT_DIR, "wwa_freq_combo_heat_AW.qs")
OUT_PREC_QS    <- file.path(OUT_DIR, "wwa_freq_combo_precip_AW.qs")
OUT_PANELS_PNG <- file.path(OUT_DIR, "wwa_freq_panels_heat_vs_precip.png")

# ----------------------------
# BBOX mask polygon in lon/lat
# ----------------------------
bbox_poly <- .make_bbox_poly(BBOX_LONLAT, CRS_PLOT)

# ----------------------------
# Geometries
# ----------------------------
if (!file.exists(ZONE_SHP)) stop("Missing ZONE_SHP: ", ZONE_SHP)
zones <- st_read(ZONE_SHP, quiet = TRUE) |> st_make_valid()
zones$STATE_ZONE <- toupper(trimws(as.character(zones$STATE_ZONE)))
zones <- st_transform(zones, CRS_PLOT)
zones_crop <- .clip_to_mask(zones, bbox_poly)

STATES <- tigris::states(cb = TRUE, year = 2022, class = "sf")
drop_statefp <- c("02", "15", "72", "60", "66", "69", "78")
STATES <- STATES[!(STATES$STATEFP %in% drop_statefp), ] |> st_make_valid()
STATES <- st_transform(STATES, CRS_PLOT)
STATES_crop <- .clip_to_mask(STATES, bbox_poly)

if (!file.exists(WFO_SHP)) stop("Missing WFO_SHP: ", WFO_SHP)
WFO <- st_read(WFO_SHP, quiet = TRUE) |> st_make_valid()
WFO <- st_transform(WFO, CRS_PLOT)
WFO_crop <- .clip_to_mask(WFO, bbox_poly)

# ----------------------------
# Index per-code TS .qs files
# ----------------------------
index_dt <- idx_qs(OUT_DIR)
.log("OUT_DIR: %s", OUT_DIR)
.log("Indexed per-code TS files: %d", nrow(index_dt))

# ----------------------------
# Heat union
# ----------------------------
if (REBUILD_COMBOS || !file.exists(OUT_HEAT_QS)) {
  ts_heat <- .read_union_ts(index_dt, HEAT_PHENOM, HEAT_SIG)
  qsave(ts_heat, OUT_HEAT_QS)
  .log("Wrote: %s  (rows=%d)", OUT_HEAT_QS, nrow(ts_heat))
} else {
  ts_heat <- as.data.table(qread(OUT_HEAT_QS))
  .log("Loaded: %s  (rows=%d)", OUT_HEAT_QS, nrow(ts_heat))
}

freq_heat <- .freq_from_ts(ts_heat)
zones_heat <- zones_crop %>%
  left_join(as.data.frame(freq_heat), by = c("STATE_ZONE" = "STATE_ZONE"))

# ----------------------------
# Precip union
# ----------------------------
if (REBUILD_COMBOS || !file.exists(OUT_PREC_QS)) {
  ts_prec <- .read_union_ts(index_dt, PREC_PHENOM, PREC_SIG)
  qsave(ts_prec, OUT_PREC_QS)
  .log("Wrote: %s  (rows=%d)", OUT_PREC_QS, nrow(ts_prec))
} else {
  ts_prec <- as.data.table(qread(OUT_PREC_QS))
  .log("Loaded: %s  (rows=%d)", OUT_PREC_QS, nrow(ts_prec))
}

freq_prec <- .freq_from_ts(ts_prec)
zones_prec <- zones_crop %>%
  left_join(as.data.frame(freq_prec), by = c("STATE_ZONE" = "STATE_ZONE"))

# ----------------------------
# Shared fill limits
# ----------------------------
max_heat <- if ("freq_days" %in% names(zones_heat)) {
  suppressWarnings(max(zones_heat$freq_days, na.rm = TRUE))
} else {
  0
}
max_prec <- if ("freq_days" %in% names(zones_prec)) {
  suppressWarnings(max(zones_prec$freq_days, na.rm = TRUE))
} else {
  0
}
fill_limits <- c(0, max(max_heat, max_prec, 1, na.rm = TRUE))

# ----------------------------
# Build two panels
# ----------------------------
p_heat <- .make_plot(
  zones_heat,
  STATES_crop,
  WFO_crop,
  bbox_poly,
  fill_limits,
  legend_title = "Count",
  panel_tag = "(a)"
)

p_prec <- .make_plot(
  zones_prec,
  STATES_crop,
  WFO_crop,
  bbox_poly,
  fill_limits,
  legend_title = "Count",
  panel_tag = "(b)"
)

p_panels <- p_heat + p_prec +
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) +
  plot_annotation(
    theme = theme(
      text = element_text(family = "Arial", size = 20),
      legend.position = "bottom",
      legend.box.just = "center",
      legend.justification = "center",
      panel.spacing = unit(0.2, "lines")
    )
  )

# ----------------------------
# Save
# ----------------------------
ragg::agg_png(
  OUT_PANELS_PNG,
  width = WIDTH_PANELS,
  height = HEIGHT_PANELS,
  res = RES_MAIN
)
print(p_panels)
dev.off()

.log("Wrote two-panel figure: %s", OUT_PANELS_PNG)
.log("Done.")