# #!/usr/bin/env Rscript
# 
# # scripts/dissertation_01_examples.R
# # =============================================================================
# # Dissertation example-event figures and county confusion matrices
# #
# # Native example events only:
# #   Heat index:
# #     heatindex-heatindex_0.25-2023-06-15-001
# #
# #   Precipitation:
# #     precipitation-precipitation_0.25-2023-08-30-001
# #
# # For each event, write:
# #
# #   1. Event-only figure
# #      - native 0.25-degree return-period raster
# #      - black outline around event footprint
# #      - state boundaries
# #      - event-duration inset at lower left
# #      - NO county X symbols
# #      - NO legend
# #
# #   2. County-confusion figure
# #      - same event map
# #      - one X for every county in the Southeast study domain
# #
# #        True positive  = cluster + warning       = red
# #        True negative  = no cluster + no warning = yellow
# #        False positive = cluster + no warning    = blue
# #        False negative = warning + no cluster    = black
# #
# #   3. Event-specific CSV containing:
# #      - TP
# #      - TN
# #      - FP
# #      - FN
# #      - recall
# #      - precision
# #      - F1
# #
# # Warning definitions:
# #   Heat          = Excessive Heat Warning (EH, W)
# #   Precipitation = Flash Flood Warning    (FF, W)
# #
# # County classification:
# #   cluster = county intersects the nonzero cluster footprint
# #   warning = county intersects >= 1 applicable warning polygon whose warning
# #             interval overlaps the cluster event interval
# #
# # Output:
# #   /proj/sebastianlab/users/hquintal/cisc/figures/dissertation/
# #
# # =============================================================================
# 
# suppressPackageStartupMessages({
#   library(terra)
#   library(sf)
#   library(ggplot2)
#   library(readr)
#   library(dplyr)
#   library(stringr)
#   library(tigris)
#   library(ragg)
#   library(RColorBrewer)
# })
# 
# options(
#   tigris_use_cache = TRUE,
#   tigris_class = "sf"
# )
# 
# sf::sf_use_s2(FALSE)
# 
# # =============================================================================
# # USER SETTINGS
# # =============================================================================
# 
# ROOT <- Sys.getenv(
#   "CISC_ROOT",
#   "/proj/sebastianlab/users/hquintal/cisc"
# )
# 
# GEOTIFF_ROOT <- file.path(
#   ROOT,
#   "data",
#   "out",
#   "08_tifs_pngs",
#   "geotiffs"
# )
# 
# CHAR_ROOT <- file.path(
#   ROOT,
#   "data",
#   "out",
#   "07_characteristics"
# )
# 
# MESONET_ROOT <- file.path(
#   ROOT,
#   "data",
#   "in",
#   "mesonet"
# )
# 
# OUT_DIR <- file.path(
#   ROOT,
#   "figures",
#   "dissertation"
# )
# 
# dir.create(
#   OUT_DIR,
#   recursive = TRUE,
#   showWarnings = FALSE
# )
# 
# # -----------------------------------------------------------------------------
# # Native example events
# # -----------------------------------------------------------------------------
# 
# EVENTS <- tibble::tribble(
#   ~variable,       ~dataset_tag,         ~cluster_id,                                      ~phenom, ~sig,
#   "heatindex",     "heatindex_0.25",     "heatindex-heatindex_0.25-2023-06-15-001",         "EH",    "W",
#   "precipitation", "precipitation_0.25", "precipitation-precipitation_0.25-2018-09-13-001", "FF", "W"
# )
# 
# # These example events come from the RP5 cluster workflow.
# FILTER_RP <- 5L
# 
# # -----------------------------------------------------------------------------
# # Southeast plotting/evaluation domain
# # -----------------------------------------------------------------------------
# 
# BBOX_LONLAT <- c(
#   xmin = -95,
#   xmax = -75,
#   ymin = 25,
#   ymax = 40
# )
# 
# CRS_PLOT <- 4326
# 
# FIG_WIDTH_IN <- 12
# FIG_HEIGHT_IN <- 8.1
# DPI <- 300
# 
# BASE_SIZE <- 18
# 
# # X plotting
# COUNTY_X_SIZE <- 3.0
# COUNTY_X_STROKE <- 1.05
# 
# CONFUSION_COLORS <- c(
#   "TP" = "red",
#   "TN" = "blue",
#   "FP" = "black",
#   "FN" = "grey"
# )
# 
# # =============================================================================
# # RETURN-PERIOD CLASSES
# # =============================================================================
# 
# RP_LEVELS <- c(
#   "2",
#   "5",
#   "10",
#   "25",
#   "50",
#   "100",
#   "200",
#   "500",
#   "1000"
# )
# 
# RP_CUT_BREAKS <- c(
#   2,
#   5,
#   10,
#   25,
#   50,
#   100,
#   200,
#   500,
#   1000,
#   Inf
# )
# 
# HEAT_COLORS <- stats::setNames(
#   c(
#     "#d9d9d9",
#     RColorBrewer::brewer.pal(
#       9,
#       "Reds"
#     )[2:9]
#   ),
#   RP_LEVELS
# )
# 
# PRECIP_COLORS <- stats::setNames(
#   c(
#     "#d9d9d9",
#     RColorBrewer::brewer.pal(
#       9,
#       "Blues"
#     )[2:9]
#   ),
#   RP_LEVELS
# )
# 
# # =============================================================================
# # HELPERS
# # =============================================================================
# 
# .log <- function(fmt, ...) {
#   message(
#     sprintf(
#       "[%s] %s",
#       format(
#         Sys.time(),
#         "%Y-%m-%d %H:%M:%S"
#       ),
#       sprintf(
#         fmt,
#         ...
#       )
#     )
#   )
# }
# 
# .safe_divide <- function(
#     numerator,
#     denominator
# ) {
#   if (
#     is.na(denominator) ||
#     denominator == 0
#   ) {
#     return(NA_real_)
#   }
#   
#   numerator / denominator
# }
# 
# .make_bbox <- function() {
#   sf::st_as_sfc(
#     sf::st_bbox(
#       c(
#         xmin = BBOX_LONLAT[["xmin"]],
#         ymin = BBOX_LONLAT[["ymin"]],
#         xmax = BBOX_LONLAT[["xmax"]],
#         ymax = BBOX_LONLAT[["ymax"]]
#       ),
#       crs = sf::st_crs(4326)
#     )
#   )
# }
# 
# # =============================================================================
# # INPUT PATHS
# # =============================================================================
# 
# .cluster_tif_path <- function(
#     dataset_tag,
#     cluster_id
# ) {
#   file.path(
#     GEOTIFF_ROOT,
#     dataset_tag,
#     paste0(
#       cluster_id,
#       ".tif"
#     )
#   )
# }
# 
# .characteristics_path <- function(
#     dataset_tag
# ) {
#   file.path(
#     CHAR_ROOT,
#     paste0(
#       "filtered_rp",
#       FILTER_RP
#     ),
#     dataset_tag,
#     "cluster_characteristics.csv"
#   )
# }
# 
# # =============================================================================
# # EVENT CHARACTERISTICS
# # =============================================================================
# 
# .read_event_characteristics <- function(
#     dataset_tag,
#     cluster_id
# ) {
#   
#   path <- .characteristics_path(
#     dataset_tag
#   )
#   
#   if (!file.exists(path)) {
#     stop(
#       "Missing characteristics file: ",
#       path,
#       call. = FALSE
#     )
#   }
#   
#   .log(
#     "Reading characteristics: %s",
#     path
#   )
#   
#   x <- readr::read_csv(
#     path,
#     show_col_types = FALSE
#   )
#   
#   if (!("cluster_id" %in% names(x))) {
#     stop(
#       "Characteristics file does not contain cluster_id: ",
#       path,
#       call. = FALSE
#     )
#   }
#   
#   event <- x |>
#     dplyr::filter(
#       cluster_id == !!cluster_id
#     )
#   
#   if (nrow(event) == 0L) {
#     stop(
#       "Cluster not found in characteristics file: ",
#       cluster_id,
#       call. = FALSE
#     )
#   }
#   
#   if (nrow(event) > 1L) {
#     warning(
#       "More than one characteristics row found for ",
#       cluster_id,
#       "; using the first."
#     )
#     
#     event <- event[1, ]
#   }
#   
#   required <- c(
#     "start_time",
#     "end_time"
#   )
#   
#   missing_cols <- setdiff(
#     required,
#     names(event)
#   )
#   
#   if (length(missing_cols)) {
#     stop(
#       "Characteristics file is missing: ",
#       paste(
#         missing_cols,
#         collapse = ", "
#       ),
#       call. = FALSE
#     )
#   }
#   
#   event_start <- as.POSIXct(
#     event$start_time[[1]],
#     tz = "UTC"
#   )
#   
#   event_end <- as.POSIXct(
#     event$end_time[[1]],
#     tz = "UTC"
#   )
#   
#   if (
#     is.na(event_start) ||
#     is.na(event_end)
#   ) {
#     stop(
#       "Could not parse event start/end time for ",
#       cluster_id,
#       call. = FALSE
#     )
#   }
#   
#   list(
#     row = event,
#     start = event_start,
#     end = event_end
#   )
# }
# 
# # =============================================================================
# # RASTER
# # =============================================================================
# 
# .read_cluster_raster <- function(
#     dataset_tag,
#     cluster_id
# ) {
#   
#   path <- .cluster_tif_path(
#     dataset_tag,
#     cluster_id
#   )
#   
#   if (!file.exists(path)) {
#     stop(
#       "Cluster GeoTIFF not found: ",
#       path,
#       call. = FALSE
#     )
#   }
#   
#   .log(
#     "Reading cluster raster: %s",
#     path
#   )
#   
#   r <- terra::rast(
#     path
#   )[[1]]
#   
#   names(r) <- "rp"
#   
#   r
# }
# 
# .raster_to_df <- function(r) {
#   
#   df <- as.data.frame(
#     r,
#     xy = TRUE,
#     na.rm = TRUE
#   )
#   
#   names(df) <- c(
#     "x",
#     "y",
#     "rp"
#   )
#   
#   df |>
#     dplyr::filter(
#       is.finite(rp),
#       rp > 0
#     ) |>
#     dplyr::mutate(
#       rp_class = cut(
#         rp,
#         breaks = RP_CUT_BREAKS,
#         labels = RP_LEVELS,
#         right = FALSE,
#         include.lowest = TRUE
#       ),
#       rp_class = factor(
#         rp_class,
#         levels = RP_LEVELS,
#         ordered = TRUE
#       )
#     )
# }
# 
# .cluster_footprint <- function(r) {
#   
#   mask <- terra::ifel(
#     is.finite(r) & r > 0,
#     1,
#     NA
#   )
#   
#   if (
#     !any(
#       !is.na(
#         terra::values(
#           mask
#         )
#       )
#     )
#   ) {
#     stop(
#       "No nonzero cells were found in cluster raster.",
#       call. = FALSE
#     )
#   }
#   
#   footprint <- terra::as.polygons(
#     mask,
#     dissolve = TRUE,
#     na.rm = TRUE,
#     values = FALSE
#   )
#   
#   footprint <- sf::st_as_sf(
#     footprint
#   )
#   
#   footprint <- sf::st_make_valid(
#     footprint
#   )
#   
#   sf::st_transform(
#     footprint,
#     CRS_PLOT
#   )
# }
# 
# # =============================================================================
# # STATES AND COUNTIES
# # =============================================================================
# 
# .load_states <- function() {
#   
#   states <- tigris::states(
#     cb = TRUE,
#     year = 2023,
#     class = "sf"
#   )
#   
#   states |>
#     dplyr::filter(
#       !STUSPS %in% c(
#         "AK",
#         "HI",
#         "PR",
#         "VI",
#         "GU",
#         "MP",
#         "AS"
#       )
#     ) |>
#     sf::st_make_valid() |>
#     sf::st_transform(
#       CRS_PLOT
#     )
# }
# 
# .load_domain_counties <- function(
#     bbox_poly
# ) {
#   
#   counties <- tigris::counties(
#     cb = TRUE,
#     year = 2023,
#     class = "sf"
#   )
#   
#   counties <- counties |>
#     dplyr::filter(
#       !STATEFP %in% c(
#         "02",
#         "15",
#         "60",
#         "66",
#         "69",
#         "72",
#         "78"
#       )
#     ) |>
#     sf::st_make_valid() |>
#     sf::st_transform(
#       CRS_PLOT
#     )
#   
#   keep <- lengths(
#     sf::st_intersects(
#       counties,
#       bbox_poly
#     )
#   ) > 0L
#   
#   counties <- counties[
#     keep,
#   ]
#   
#   counties$GEOID <- stringr::str_pad(
#     as.character(
#       counties$GEOID
#     ),
#     width = 5,
#     side = "left",
#     pad = "0"
#   )
#   
#   counties
# }
# 
# .make_county_points <- function(
#     counties
# ) {
#   
#   # Generate points in an equal-area CRS so the X remains inside each county.
#   points <- counties |>
#     sf::st_transform(
#       5070
#     ) |>
#     sf::st_point_on_surface() |>
#     sf::st_transform(
#       CRS_PLOT
#     )
#   
#   points
# }
# 
# # =============================================================================
# # MESONET WARNING DATA
# # =============================================================================
# 
# .parse_wwa_time <- function(x) {
#   
#   x <- trimws(
#     as.character(x)
#   )
#   
#   x[x %in% c(
#     "",
#     "NA",
#     "NULL"
#   )] <- NA_character_
#   
#   as.POSIXct(
#     x,
#     format = "%Y%m%d%H%M",
#     tz = "UTC"
#   )
# }
# 
# .find_column <- function(
#     x,
#     choices,
#     label
# ) {
#   
#   hit <- choices[
#     choices %in% names(x)
#   ]
#   
#   if (!length(hit)) {
#     stop(
#       "Could not find ",
#       label,
#       " column. Tried: ",
#       paste(
#         choices,
#         collapse = ", "
#       ),
#       "\nAvailable fields: ",
#       paste(
#         names(x),
#         collapse = ", "
#       ),
#       call. = FALSE
#     )
#   }
#   
#   hit[[1]]
# }
# 
# .read_warning_data <- function(
#     event_start,
#     event_end,
#     phenom,
#     sig
# ) {
#   
#   years <- seq.int(
#     as.integer(
#       format(
#         event_start,
#         "%Y"
#       )
#     ),
#     as.integer(
#       format(
#         event_end,
#         "%Y"
#       )
#     )
#   )
#   
#   pieces <- list()
#   
#   for (yr in years) {
#     
#     year_dir <- file.path(
#       MESONET_ROOT,
#       paste0(
#         yr,
#         "_all"
#       )
#     )
#     
#     shp_files <- list.files(
#       year_dir,
#       pattern = "^wwa_.*\\.shp$",
#       full.names = TRUE
#     )
#     
#     if (!length(shp_files)) {
#       stop(
#         "No Mesonet WWA shapefile found in: ",
#         year_dir,
#         call. = FALSE
#       )
#     }
#     
#     # Normally there is one yearly WWA shapefile.
#     for (shp in shp_files) {
#       
#       .log(
#         "Reading warnings: %s",
#         shp
#       )
#       
#       pieces[[length(pieces) + 1L]] <- sf::st_read(
#         shp,
#         quiet = TRUE
#       )
#     }
#   }
#   
#   warnings <- do.call(
#     rbind,
#     pieces
#   )
#   
#   warnings <- sf::st_make_valid(
#     warnings
#   )
#   
#   phenom_col <- .find_column(
#     warnings,
#     c(
#       "WWA_PHENOM",
#       "PHENOM"
#     ),
#     "phenomenon"
#   )
#   
#   sig_col <- .find_column(
#     warnings,
#     c(
#       "WWA_SIG",
#       "SIG"
#     ),
#     "significance"
#   )
#   
#   issued_col <- .find_column(
#     warnings,
#     c(
#       "ISSUED",
#       "ISSUE"
#     ),
#     "ISSUED"
#   )
#   
#   expired_col <- .find_column(
#     warnings,
#     c(
#       "EXPIRED",
#       "EXPIRE"
#     ),
#     "EXPIRED"
#   )
#   
#   warnings$.phenom <- trimws(
#     as.character(
#       warnings[[phenom_col]]
#     )
#   )
#   
#   warnings$.sig <- trimws(
#     as.character(
#       warnings[[sig_col]]
#     )
#   )
#   
#   warnings$.issued <- .parse_wwa_time(
#     warnings[[issued_col]]
#   )
#   
#   warnings$.expired <- .parse_wwa_time(
#     warnings[[expired_col]]
#   )
#   
#   # Interval-overlap rule:
#   #
#   # warning starts before cluster ends
#   # AND
#   # warning ends after cluster starts
#   
#   warnings <- warnings[
#     warnings$.phenom == phenom &
#       warnings$.sig == sig &
#       !is.na(warnings$.issued) &
#       !is.na(warnings$.expired) &
#       warnings$.issued <= event_end &
#       warnings$.expired >= event_start,
#   ]
#   
#   if (nrow(warnings) == 0L) {
#     
#     .log(
#       "No %s_%s warnings overlap event interval.",
#       phenom,
#       sig
#     )
#     
#     return(
#       warnings
#     )
#   }
#   
#   warnings <- sf::st_transform(
#     warnings,
#     CRS_PLOT
#   )
#   
#   .log(
#     "Selected %d overlapping %s_%s warning polygons.",
#     nrow(warnings),
#     phenom,
#     sig
#   )
#   
#   warnings
# }
# 
# # =============================================================================
# # COUNTY CONFUSION CLASSIFICATION
# # =============================================================================
# 
# .classify_counties <- function(
#     counties,
#     cluster_footprint,
#     warnings
# ) {
#   
#   cluster_present <- lengths(
#     sf::st_intersects(
#       counties,
#       cluster_footprint
#     )
#   ) > 0L
#   
#   if (nrow(warnings) > 0L) {
#     
#     warning_present <- lengths(
#       sf::st_intersects(
#         counties,
#         warnings
#       )
#     ) > 0L
#     
#   } else {
#     
#     warning_present <- rep(
#       FALSE,
#       nrow(counties)
#     )
#   }
#   
#   class <- dplyr::case_when(
#     cluster_present & warning_present  ~ "TP",
#     !cluster_present & !warning_present ~ "TN",
#     cluster_present & !warning_present ~ "FP",
#     !cluster_present & warning_present ~ "FN"
#   )
#   
#   counties$cluster_present <- cluster_present
#   counties$warning_present <- warning_present
#   counties$confusion_class <- factor(
#     class,
#     levels = c(
#       "TP",
#       "TN",
#       "FP",
#       "FN"
#     )
#   )
#   
#   counties
# }
# 
# .calculate_metrics <- function(
#     counties,
#     event
# ) {
#   
#   class <- as.character(
#     counties$confusion_class
#   )
#   
#   tp <- sum(
#     class == "TP",
#     na.rm = TRUE
#   )
#   
#   tn <- sum(
#     class == "TN",
#     na.rm = TRUE
#   )
#   
#   fp <- sum(
#     class == "FP",
#     na.rm = TRUE
#   )
#   
#   fn <- sum(
#     class == "FN",
#     na.rm = TRUE
#   )
#   
#   recall <- .safe_divide(
#     tp,
#     tp + fn
#   )
#   
#   precision <- .safe_divide(
#     tp,
#     tp + fp
#   )
#   
#   f1 <- if (
#     is.na(recall) ||
#     is.na(precision) ||
#     recall + precision == 0
#   ) {
#     NA_real_
#   } else {
#     2 * precision * recall /
#       (
#         precision + recall
#       )
#   }
#   
#   tibble::tibble(
#     variable = event$variable,
#     dataset_tag = event$dataset_tag,
#     cluster_id = event$cluster_id,
#     warning_phenom = event$phenom,
#     warning_sig = event$sig,
#     event_start = format(
#       event$event_start,
#       "%Y-%m-%d %H:%M:%S",
#       tz = "UTC"
#     ),
#     event_end = format(
#       event$event_end,
#       "%Y-%m-%d %H:%M:%S",
#       tz = "UTC"
#     ),
#     n_counties = nrow(counties),
#     true_positive = tp,
#     true_negative = tn,
#     false_positive = fp,
#     false_negative = fn,
#     recall = recall,
#     precision = precision,
#     f1 = f1
#   )
# }
# 
# # =============================================================================
# # DURATION LABEL
# # =============================================================================
# 
# .duration_label <- function(
#     start,
#     end
# ) {
#   
#   duration_days <- as.numeric(
#     difftime(
#       end,
#       start,
#       units = "days"
#     )
#   )
#   
#   start_text <- format(
#     start,
#     "%b %d, %Y %H:%M",
#     tz = "UTC"
#   )
#   
#   end_text <- format(
#     end,
#     "%b %d, %Y %H:%M",
#     tz = "UTC"
#   )
#   
#   sprintf(
#     "%s UTC – %s UTC\nDuration: %.1f days",
#     start_text,
#     end_text,
#     duration_days
#   )
# }
# 
# # =============================================================================
# # PLOTTING
# # =============================================================================
# 
# .event_palette <- function(
#     variable
# ) {
#   
#   if (variable == "heatindex") {
#     return(
#       HEAT_COLORS
#     )
#   }
#   
#   if (variable == "precipitation") {
#     return(
#       PRECIP_COLORS
#     )
#   }
#   
#   stop(
#     "Unsupported variable: ",
#     variable,
#     call. = FALSE
#   )
# }
# 
# .build_event_plot <- function(
#     variable,
#     raster,
#     footprint,
#     states,
#     county_points = NULL,
#     event_start,
#     event_end
# ) {
#   
#   raster_df <- .raster_to_df(
#     raster
#   )
#   
#   dx <- terra::res(
#     raster
#   )[1]
#   
#   dy <- terra::res(
#     raster
#   )[2]
#   
#   palette <- .event_palette(
#     variable
#   )
#   
#   duration_text <- .duration_label(
#     event_start,
#     event_end
#   )
#   
#   p <- ggplot2::ggplot() +
#     
#     # -------------------------------------------------------------------------
#   # States first: white land background
#   # -------------------------------------------------------------------------
#   
#   ggplot2::geom_sf(
#     data = states,
#     fill = "white",
#     color = "black",
#     linewidth = 0.45,
#     inherit.aes = FALSE,
#     show.legend = FALSE
#   ) +
#     
#     # -------------------------------------------------------------------------
#   # Cluster return-period raster plotted over states
#   # -------------------------------------------------------------------------
#   
#   ggplot2::geom_tile(
#     data = raster_df,
#     mapping = ggplot2::aes(
#       x = x,
#       y = y,
#       fill = rp_class
#     ),
#     width = dx,
#     height = dy,
#     show.legend = FALSE
#   ) +
#     
#     ggplot2::scale_fill_manual(
#       values = palette,
#       limits = RP_LEVELS,
#       drop = FALSE,
#       na.value = "white",
#       guide = "none"
#     ) +
#     
#     # -------------------------------------------------------------------------
#   # Cluster footprint
#   # -------------------------------------------------------------------------
#   
#   ggplot2::geom_sf(
#     data = footprint,
#     fill = NA,
#     color = "black",
#     linewidth = 0.65,
#     linejoin = "mitre",
#     inherit.aes = FALSE,
#     show.legend = FALSE
#   )
#   
#   # ---------------------------------------------------------------------------
#   # Optional county confusion X symbols
#   # ---------------------------------------------------------------------------
#   
#   if (!is.null(county_points)) {
#     
#     p <- p +
#       ggplot2::geom_sf(
#         data = county_points,
#         mapping = ggplot2::aes(
#           color = confusion_class
#         ),
#         shape = 4,
#         size = COUNTY_X_SIZE,
#         stroke = COUNTY_X_STROKE,
#         inherit.aes = FALSE,
#         show.legend = FALSE
#       ) +
#       ggplot2::scale_color_manual(
#         values = CONFUSION_COLORS,
#         limits = names(
#           CONFUSION_COLORS
#         ),
#         guide = "none"
#       )
#   }
#   
#   p +
#     
#     # -------------------------------------------------------------------------
#   # Southeast dissertation domain
#   # -------------------------------------------------------------------------
#   
#   ggplot2::coord_sf(
#     crs = sf::st_crs(
#       CRS_PLOT
#     ),
#     datum = sf::st_crs(
#       CRS_PLOT
#     ),
#     xlim = c(
#       BBOX_LONLAT[["xmin"]],
#       BBOX_LONLAT[["xmax"]]
#     ),
#     ylim = c(
#       BBOX_LONLAT[["ymin"]],
#       BBOX_LONLAT[["ymax"]]
#     ),
#     expand = FALSE,
#     label_axes = list(
#       bottom = "E",
#       left = "N"
#     )
#   ) +
#     
#     # -------------------------------------------------------------------------
#   # Event-duration inset
#   # -------------------------------------------------------------------------
#   
#   ggplot2::annotate(
#     "label",
#     x = BBOX_LONLAT[["xmin"]] + 0.55,
#     y = BBOX_LONLAT[["ymin"]] + 0.55,
#     label = duration_text,
#     hjust = 0,
#     vjust = 0,
#     size = 4.5,
#     family = "sans",
#     lineheight = 1.05,
#     label.size = 0.35,
#     label.padding = grid::unit(
#       0.20,
#       "lines"
#     ),
#     fill = scales::alpha(
#       "white",
#       0.88
#     ),
#     color = "black"
#   ) +
#     
#     ggplot2::labs(
#       x = NULL,
#       y = NULL
#     ) +
#     
#     ggplot2::theme_bw(
#       base_size = BASE_SIZE,
#       base_family = "sans"
#     ) +
#     
#     ggplot2::theme(
#       text = ggplot2::element_text(
#         family = "sans",
#         size = BASE_SIZE
#       ),
#       panel.grid.major = ggplot2::element_line(
#         color = "grey82",
#         linewidth = 0.35
#       ),
#       panel.grid.minor = ggplot2::element_blank(),
#       panel.background = ggplot2::element_rect(
#         fill = "white",
#         color = NA
#       ),
#       plot.background = ggplot2::element_rect(
#         fill = "white",
#         color = NA
#       ),
#       axis.title = ggplot2::element_blank(),
#       axis.text = ggplot2::element_text(
#         size = 14,
#         color = "grey25"
#       ),
#       axis.ticks = ggplot2::element_line(
#         color = "grey25",
#         linewidth = 0.3
#       ),
#       panel.border = ggplot2::element_rect(
#         color = "black",
#         fill = NA,
#         linewidth = 0.55
#       ),
#       
#       # Explicitly eliminate every legend.
#       legend.position = "none",
#       
#       plot.margin = ggplot2::margin(
#         8,
#         8,
#         8,
#         8
#       )
#     )
# }
# 
# # =============================================================================
# # SAVE FIGURE
# # =============================================================================
# 
# .save_png <- function(
#     plot,
#     path
# ) {
#   
#   ragg::agg_png(
#     filename = path,
#     width = FIG_WIDTH_IN,
#     height = FIG_HEIGHT_IN,
#     units = "in",
#     res = DPI,
#     background = "white"
#   )
#   
#   print(
#     plot
#   )
#   
#   grDevices::dev.off()
#   
#   .log(
#     "Wrote: %s",
#     path
#   )
# }
# 
# # =============================================================================
# # PROCESS ONE EVENT
# # =============================================================================
# 
# .process_event <- function(
#     variable,
#     dataset_tag,
#     cluster_id,
#     phenom,
#     sig,
#     states,
#     counties,
#     bbox_poly
# ) {
#   
#   .log(
#     "============================================================"
#   )
#   
#   .log(
#     "Processing %s",
#     cluster_id
#   )
#   
#   # ---------------------------------------------------------------------------
#   # Event timing
#   # ---------------------------------------------------------------------------
#   
#   char <- .read_event_characteristics(
#     dataset_tag = dataset_tag,
#     cluster_id = cluster_id
#   )
#   
#   event_start <- char$start
#   event_end <- char$end
#   
#   .log(
#     "Event interval: %s to %s UTC",
#     format(
#       event_start,
#       "%Y-%m-%d %H:%M"
#     ),
#     format(
#       event_end,
#       "%Y-%m-%d %H:%M"
#     )
#   )
#   
#   # ---------------------------------------------------------------------------
#   # Raster / footprint
#   # ---------------------------------------------------------------------------
#   
#   raster <- .read_cluster_raster(
#     dataset_tag = dataset_tag,
#     cluster_id = cluster_id
#   )
#   
#   footprint <- .cluster_footprint(
#     raster
#   )
#   
#   # ---------------------------------------------------------------------------
#   # Warnings overlapping event duration
#   # ---------------------------------------------------------------------------
#   
#   warnings <- .read_warning_data(
#     event_start = event_start,
#     event_end = event_end,
#     phenom = phenom,
#     sig = sig
#   )
#   
#   # ---------------------------------------------------------------------------
#   # County classifications
#   # ---------------------------------------------------------------------------
#   
#   classified_counties <- .classify_counties(
#     counties = counties,
#     cluster_footprint = footprint,
#     warnings = warnings
#   )
#   
#   county_points <- .make_county_points(
#     classified_counties
#   )
#   
#   counts <- table(
#     factor(
#       classified_counties$confusion_class,
#       levels = c(
#         "TP",
#         "TN",
#         "FP",
#         "FN"
#       )
#     )
#   )
#   
#   .log(
#     "County confusion counts: TP=%d | TN=%d | FP=%d | FN=%d",
#     counts[["TP"]],
#     counts[["TN"]],
#     counts[["FP"]],
#     counts[["FN"]]
#   )
#   
#   # ---------------------------------------------------------------------------
#   # Metrics
#   # ---------------------------------------------------------------------------
#   
#   event_info <- list(
#     variable = variable,
#     dataset_tag = dataset_tag,
#     cluster_id = cluster_id,
#     phenom = phenom,
#     sig = sig,
#     event_start = event_start,
#     event_end = event_end
#   )
#   
#   metrics <- .calculate_metrics(
#     counties = classified_counties,
#     event = event_info
#   )
#   
#   .log(
#     "Recall=%.4f | Precision=%.4f | F1=%.4f",
#     metrics$recall,
#     metrics$precision,
#     metrics$f1
#   )
#   
#   # ---------------------------------------------------------------------------
#   # Output-name stem
#   # ---------------------------------------------------------------------------
#   
#   short_name <- if (
#     variable == "heatindex"
#   ) {
#     "heatindex"
#   } else {
#     "precipitation"
#   }
#   
#   # ---------------------------------------------------------------------------
#   # Figure 1: event only
#   # ---------------------------------------------------------------------------
#   
#   p_event <- .build_event_plot(
#     variable = variable,
#     raster = raster,
#     footprint = footprint,
#     states = states,
#     county_points = NULL,
#     event_start = event_start,
#     event_end = event_end
#   )
#   
#   event_png <- file.path(
#     OUT_DIR,
#     paste0(
#       "dissertation_01_examples_",
#       short_name,
#       "_event.png"
#     )
#   )
#   
#   .save_png(
#     p_event,
#     event_png
#   )
#   
#   # ---------------------------------------------------------------------------
#   # Figure 2: all county X symbols classified TP/TN/FP/FN
#   # ---------------------------------------------------------------------------
#   
#   p_confusion <- .build_event_plot(
#     variable = variable,
#     raster = raster,
#     footprint = footprint,
#     states = states,
#     county_points = county_points,
#     event_start = event_start,
#     event_end = event_end
#   )
#   
#   confusion_png <- file.path(
#     OUT_DIR,
#     paste0(
#       "dissertation_01_examples_",
#       short_name,
#       "_confusion.png"
#     )
#   )
#   
#   .save_png(
#     p_confusion,
#     confusion_png
#   )
#   
#   # ---------------------------------------------------------------------------
#   # CSV
#   # ---------------------------------------------------------------------------
#   
#   metrics_csv <- file.path(
#     OUT_DIR,
#     paste0(
#       "dissertation_01_examples_",
#       short_name,
#       "_metrics.csv"
#     )
#   )
#   
#   readr::write_csv(
#     metrics,
#     metrics_csv,
#     na = ""
#   )
#   
#   .log(
#     "Wrote: %s",
#     metrics_csv
#   )
#   
#   invisible(
#     list(
#       metrics = metrics,
#       counties = classified_counties
#     )
#   )
# }
# 
# # =============================================================================
# # LOAD COMMON SPATIAL DATA
# # =============================================================================
# 
# bbox_poly <- .make_bbox()
# 
# .log(
#   "Loading states."
# )
# 
# states <- .load_states()
# 
# .log(
#   "Loading counties."
# )
# 
# counties <- .load_domain_counties(
#   bbox_poly
# )
# 
# .log(
#   "Evaluation domain contains %d counties.",
#   nrow(counties)
# )
# 
# # =============================================================================
# # RUN BOTH EXAMPLE EVENTS
# # =============================================================================
# 
# results <- vector(
#   "list",
#   nrow(EVENTS)
# )
# 
# for (i in seq_len(
#   nrow(EVENTS)
# )) {
#   
#   event <- EVENTS[
#     i,
#   ]
#   
#   results[[i]] <- .process_event(
#     variable = event$variable[[1]],
#     dataset_tag = event$dataset_tag[[1]],
#     cluster_id = event$cluster_id[[1]],
#     phenom = event$phenom[[1]],
#     sig = event$sig[[1]],
#     states = states,
#     counties = counties,
#     bbox_poly = bbox_poly
#   )
# }
# 
# names(results) <- EVENTS$variable
# 
# .log(
#   "============================================================"
# )
# 
# .log(
#   "Complete."
# )
# 
# .log(
#   "Outputs written to: %s",
#   OUT_DIR
# )

# tn only ----
#!/usr/bin/env Rscript

# scripts/dissertation_01_examples.R
# =============================================================================
# Dissertation example-event figures and county confusion matrices
#
# Native example events only:
#
#   Heat index:
#     heatindex-heatindex_0.25-2023-06-15-001
#
#   Precipitation:
#     precipitation-precipitation_0.25-2018-09-13-001
#
# For each event, write:
#
#   1. Event-only figure
#      - native 0.25-degree return-period raster
#      - raster colored by return period
#      - black outline around event footprint
#      - state boundaries
#      - event-duration inset at lower left
#      - NO county X symbols
#      - NO legend
#
#   2. County-confusion figure
#      - NO return-period raster coloring
#      - black outline around event footprint
#      - state boundaries
#      - county X symbols only for:
#
#          True positive  = cluster + warning    = red
#          False positive = cluster + no warning = blue
#          False negative = warning + no cluster = black
#
#      - True negatives are NOT plotted
#      - event-duration inset at lower left
#      - NO legend
#
#   3. Event-specific CSV containing:
#      - TP
#      - TN
#      - FP
#      - FN
#      - recall
#      - precision
#      - F1
#
# Warning definitions:
#   Heat          = Excessive Heat Warning (EH, W)
#   Precipitation = Flash Flood Warning    (FF, W)
#
# County classification:
#   cluster = county intersects the nonzero cluster footprint
#   warning = county intersects >= 1 applicable warning polygon whose warning
#             interval overlaps the cluster event interval
#
# Output:
#   /proj/sebastianlab/users/hquintal/cisc/figures/dissertation/
#
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tigris)
  library(ragg)
  library(RColorBrewer)
})

options(
  tigris_use_cache = TRUE,
  tigris_class = "sf"
)

sf::sf_use_s2(FALSE)

# =============================================================================
# USER SETTINGS
# =============================================================================

ROOT <- Sys.getenv(
  "CISC_ROOT",
  "/proj/sebastianlab/users/hquintal/cisc"
)

GEOTIFF_ROOT <- file.path(
  ROOT,
  "data",
  "out",
  "08_tifs_pngs",
  "geotiffs"
)

CHAR_ROOT <- file.path(
  ROOT,
  "data",
  "out",
  "07_characteristics"
)

MESONET_ROOT <- file.path(
  ROOT,
  "data",
  "in",
  "mesonet"
)

OUT_DIR <- file.path(
  ROOT,
  "figures",
  "dissertation"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# -----------------------------------------------------------------------------
# Native example events
# -----------------------------------------------------------------------------

EVENTS <- tibble::tribble(
  ~variable,       ~dataset_tag,         ~cluster_id,                                      ~phenom, ~sig,
  "heatindex",     "heatindex_0.25",     "heatindex-heatindex_0.25-2023-06-15-001",         "EH",    "W",
  "precipitation", "precipitation_0.25", "precipitation-precipitation_0.25-2018-09-13-001", "FF",    "W"
)

# These example events come from the RP5 cluster workflow.
FILTER_RP <- 5L

# -----------------------------------------------------------------------------
# Southeast plotting/evaluation domain
# -----------------------------------------------------------------------------

BBOX_LONLAT <- c(
  xmin = -95,
  xmax = -75,
  ymin = 25,
  ymax = 40
)

CRS_PLOT <- 4326

FIG_WIDTH_IN <- 12
FIG_HEIGHT_IN <- 8.1
DPI <- 300

BASE_SIZE <- 18

# County X plotting
COUNTY_X_SIZE <- 3.0
COUNTY_X_STROKE <- 1.05

CONFUSION_COLORS <- c(
  "TP" = "red",
  "FP" = "blue",
  "FN" = "black"
)

# =============================================================================
# RETURN-PERIOD CLASSES
# =============================================================================

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

HEAT_COLORS <- stats::setNames(
  c(
    "#d9d9d9",
    RColorBrewer::brewer.pal(
      9,
      "Reds"
    )[2:9]
  ),
  RP_LEVELS
)

PRECIP_COLORS <- stats::setNames(
  c(
    "#d9d9d9",
    RColorBrewer::brewer.pal(
      9,
      "Blues"
    )[2:9]
  ),
  RP_LEVELS
)

# =============================================================================
# HELPERS
# =============================================================================

.log <- function(fmt, ...) {
  message(
    sprintf(
      "[%s] %s",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      ),
      sprintf(
        fmt,
        ...
      )
    )
  )
}

.safe_divide <- function(
    numerator,
    denominator
) {
  if (
    is.na(denominator) ||
    denominator == 0
  ) {
    return(NA_real_)
  }
  
  numerator / denominator
}

.make_bbox <- function() {
  sf::st_as_sfc(
    sf::st_bbox(
      c(
        xmin = BBOX_LONLAT[["xmin"]],
        ymin = BBOX_LONLAT[["ymin"]],
        xmax = BBOX_LONLAT[["xmax"]],
        ymax = BBOX_LONLAT[["ymax"]]
      ),
      crs = sf::st_crs(4326)
    )
  )
}

# =============================================================================
# INPUT PATHS
# =============================================================================

.cluster_tif_path <- function(
    dataset_tag,
    cluster_id
) {
  file.path(
    GEOTIFF_ROOT,
    dataset_tag,
    paste0(
      cluster_id,
      ".tif"
    )
  )
}

.characteristics_path <- function(
    dataset_tag
) {
  file.path(
    CHAR_ROOT,
    paste0(
      "filtered_rp",
      FILTER_RP
    ),
    dataset_tag,
    "cluster_characteristics.csv"
  )
}

# =============================================================================
# EVENT CHARACTERISTICS
# =============================================================================

.read_event_characteristics <- function(
    dataset_tag,
    cluster_id
) {
  
  path <- .characteristics_path(
    dataset_tag
  )
  
  if (!file.exists(path)) {
    stop(
      "Missing characteristics file: ",
      path,
      call. = FALSE
    )
  }
  
  .log(
    "Reading characteristics: %s",
    path
  )
  
  x <- readr::read_csv(
    path,
    show_col_types = FALSE
  )
  
  if (!("cluster_id" %in% names(x))) {
    stop(
      "Characteristics file does not contain cluster_id: ",
      path,
      call. = FALSE
    )
  }
  
  event <- x |>
    dplyr::filter(
      cluster_id == !!cluster_id
    )
  
  if (nrow(event) == 0L) {
    stop(
      "Cluster not found in characteristics file: ",
      cluster_id,
      call. = FALSE
    )
  }
  
  if (nrow(event) > 1L) {
    warning(
      "More than one characteristics row found for ",
      cluster_id,
      "; using the first."
    )
    
    event <- event[1, ]
  }
  
  required <- c(
    "start_time",
    "end_time"
  )
  
  missing_cols <- setdiff(
    required,
    names(event)
  )
  
  if (length(missing_cols)) {
    stop(
      "Characteristics file is missing: ",
      paste(
        missing_cols,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  event_start <- as.POSIXct(
    event$start_time[[1]],
    tz = "UTC"
  )
  
  event_end <- as.POSIXct(
    event$end_time[[1]],
    tz = "UTC"
  )
  
  if (
    is.na(event_start) ||
    is.na(event_end)
  ) {
    stop(
      "Could not parse event start/end time for ",
      cluster_id,
      call. = FALSE
    )
  }
  
  list(
    row = event,
    start = event_start,
    end = event_end
  )
}

# =============================================================================
# RASTER
# =============================================================================

.read_cluster_raster <- function(
    dataset_tag,
    cluster_id
) {
  
  path <- .cluster_tif_path(
    dataset_tag,
    cluster_id
  )
  
  if (!file.exists(path)) {
    stop(
      "Cluster GeoTIFF not found: ",
      path,
      call. = FALSE
    )
  }
  
  .log(
    "Reading cluster raster: %s",
    path
  )
  
  r <- terra::rast(
    path
  )[[1]]
  
  names(r) <- "rp"
  
  r
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
  
  df |>
    dplyr::filter(
      is.finite(rp),
      rp > 0
    ) |>
    dplyr::mutate(
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
}

.cluster_footprint <- function(r) {
  
  mask <- terra::ifel(
    is.finite(r) & r > 0,
    1,
    NA
  )
  
  if (
    !any(
      !is.na(
        terra::values(
          mask
        )
      )
    )
  ) {
    stop(
      "No nonzero cells were found in cluster raster.",
      call. = FALSE
    )
  }
  
  footprint <- terra::as.polygons(
    mask,
    dissolve = TRUE,
    na.rm = TRUE,
    values = FALSE
  )
  
  footprint <- sf::st_as_sf(
    footprint
  )
  
  footprint <- sf::st_make_valid(
    footprint
  )
  
  sf::st_transform(
    footprint,
    CRS_PLOT
  )
}

# =============================================================================
# STATES AND COUNTIES
# =============================================================================

.load_states <- function() {
  
  states <- tigris::states(
    cb = TRUE,
    year = 2023,
    class = "sf"
  )
  
  states |>
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
    ) |>
    sf::st_make_valid() |>
    sf::st_transform(
      CRS_PLOT
    )
}

.load_domain_counties <- function(
    bbox_poly
) {
  
  counties <- tigris::counties(
    cb = TRUE,
    year = 2023,
    class = "sf"
  )
  
  counties <- counties |>
    dplyr::filter(
      !STATEFP %in% c(
        "02",
        "15",
        "60",
        "66",
        "69",
        "72",
        "78"
      )
    ) |>
    sf::st_make_valid() |>
    sf::st_transform(
      CRS_PLOT
    )
  
  keep <- lengths(
    sf::st_intersects(
      counties,
      bbox_poly
    )
  ) > 0L
  
  counties <- counties[
    keep,
  ]
  
  counties$GEOID <- stringr::str_pad(
    as.character(
      counties$GEOID
    ),
    width = 5,
    side = "left",
    pad = "0"
  )
  
  counties
}

.make_county_points <- function(
    counties
) {
  
  counties |>
    sf::st_transform(
      5070
    ) |>
    sf::st_point_on_surface() |>
    sf::st_transform(
      CRS_PLOT
    )
}

# =============================================================================
# MESONET WARNING DATA
# =============================================================================

.parse_wwa_time <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  x[x %in% c(
    "",
    "NA",
    "NULL"
  )] <- NA_character_
  
  as.POSIXct(
    x,
    format = "%Y%m%d%H%M",
    tz = "UTC"
  )
}

.find_column <- function(
    x,
    choices,
    label
) {
  
  hit <- choices[
    choices %in% names(x)
  ]
  
  if (!length(hit)) {
    stop(
      "Could not find ",
      label,
      " column. Tried: ",
      paste(
        choices,
        collapse = ", "
      ),
      "\nAvailable fields: ",
      paste(
        names(x),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  hit[[1]]
}

.read_warning_data <- function(
    event_start,
    event_end,
    phenom,
    sig
) {
  
  years <- seq.int(
    as.integer(
      format(
        event_start,
        "%Y"
      )
    ),
    as.integer(
      format(
        event_end,
        "%Y"
      )
    )
  )
  
  pieces <- list()
  
  for (yr in years) {
    
    year_dir <- file.path(
      MESONET_ROOT,
      paste0(
        yr,
        "_all"
      )
    )
    
    shp_files <- list.files(
      year_dir,
      pattern = "^wwa_.*\\.shp$",
      full.names = TRUE
    )
    
    if (!length(shp_files)) {
      stop(
        "No Mesonet WWA shapefile found in: ",
        year_dir,
        call. = FALSE
      )
    }
    
    for (shp in shp_files) {
      
      .log(
        "Reading warnings: %s",
        shp
      )
      
      pieces[[length(pieces) + 1L]] <- sf::st_read(
        shp,
        quiet = TRUE
      )
    }
  }
  
  warnings <- do.call(
    rbind,
    pieces
  )
  
  warnings <- sf::st_make_valid(
    warnings
  )
  
  phenom_col <- .find_column(
    warnings,
    c(
      "WWA_PHENOM",
      "PHENOM"
    ),
    "phenomenon"
  )
  
  sig_col <- .find_column(
    warnings,
    c(
      "WWA_SIG",
      "SIG"
    ),
    "significance"
  )
  
  issued_col <- .find_column(
    warnings,
    c(
      "ISSUED",
      "ISSUE"
    ),
    "ISSUED"
  )
  
  expired_col <- .find_column(
    warnings,
    c(
      "EXPIRED",
      "EXPIRE"
    ),
    "EXPIRED"
  )
  
  warnings$.phenom <- trimws(
    as.character(
      warnings[[phenom_col]]
    )
  )
  
  warnings$.sig <- trimws(
    as.character(
      warnings[[sig_col]]
    )
  )
  
  warnings$.issued <- .parse_wwa_time(
    warnings[[issued_col]]
  )
  
  warnings$.expired <- .parse_wwa_time(
    warnings[[expired_col]]
  )
  
  # Warning interval overlaps event interval when:
  #
  #   warning begins before event ends
  #   AND
  #   warning ends after event begins
  
  warnings <- warnings[
    warnings$.phenom == phenom &
      warnings$.sig == sig &
      !is.na(warnings$.issued) &
      !is.na(warnings$.expired) &
      warnings$.issued <= event_end &
      warnings$.expired >= event_start,
  ]
  
  if (nrow(warnings) == 0L) {
    
    .log(
      "No %s_%s warnings overlap event interval.",
      phenom,
      sig
    )
    
    return(
      warnings
    )
  }
  
  warnings <- sf::st_transform(
    warnings,
    CRS_PLOT
  )
  
  .log(
    "Selected %d overlapping %s_%s warning polygons.",
    nrow(warnings),
    phenom,
    sig
  )
  
  warnings
}

# =============================================================================
# COUNTY CONFUSION CLASSIFICATION
# =============================================================================

.classify_counties <- function(
    counties,
    cluster_footprint,
    warnings
) {
  
  cluster_present <- lengths(
    sf::st_intersects(
      counties,
      cluster_footprint
    )
  ) > 0L
  
  if (nrow(warnings) > 0L) {
    
    warning_present <- lengths(
      sf::st_intersects(
        counties,
        warnings
      )
    ) > 0L
    
  } else {
    
    warning_present <- rep(
      FALSE,
      nrow(counties)
    )
  }
  
  class <- dplyr::case_when(
    cluster_present & warning_present   ~ "TP",
    !cluster_present & !warning_present ~ "TN",
    cluster_present & !warning_present  ~ "FP",
    !cluster_present & warning_present  ~ "FN"
  )
  
  counties$cluster_present <- cluster_present
  counties$warning_present <- warning_present
  
  counties$confusion_class <- factor(
    class,
    levels = c(
      "TP",
      "TN",
      "FP",
      "FN"
    )
  )
  
  counties
}

.calculate_metrics <- function(
    counties,
    event
) {
  
  class <- as.character(
    counties$confusion_class
  )
  
  tp <- sum(
    class == "TP",
    na.rm = TRUE
  )
  
  tn <- sum(
    class == "TN",
    na.rm = TRUE
  )
  
  fp <- sum(
    class == "FP",
    na.rm = TRUE
  )
  
  fn <- sum(
    class == "FN",
    na.rm = TRUE
  )
  
  recall <- .safe_divide(
    tp,
    tp + fn
  )
  
  precision <- .safe_divide(
    tp,
    tp + fp
  )
  
  f1 <- if (
    is.na(recall) ||
    is.na(precision) ||
    recall + precision == 0
  ) {
    NA_real_
  } else {
    2 * precision * recall /
      (
        precision + recall
      )
  }
  
  tibble::tibble(
    variable = event$variable,
    dataset_tag = event$dataset_tag,
    cluster_id = event$cluster_id,
    warning_phenom = event$phenom,
    warning_sig = event$sig,
    event_start = format(
      event$event_start,
      "%Y-%m-%d %H:%M:%S",
      tz = "UTC"
    ),
    event_end = format(
      event$event_end,
      "%Y-%m-%d %H:%M:%S",
      tz = "UTC"
    ),
    n_counties = nrow(counties),
    true_positive = tp,
    true_negative = tn,
    false_positive = fp,
    false_negative = fn,
    recall = recall,
    precision = precision,
    f1 = f1
  )
}

# =============================================================================
# DURATION LABEL
# =============================================================================

.duration_label <- function(
    start,
    end
) {
  
  duration_days <- as.numeric(
    difftime(
      end,
      start,
      units = "days"
    )
  )
  
  start_text <- format(
    start,
    "%b %d, %Y %H:%M",
    tz = "UTC"
  )
  
  end_text <- format(
    end,
    "%b %d, %Y %H:%M",
    tz = "UTC"
  )
  
  sprintf(
    "%s UTC – %s UTC\nDuration: %.1f days",
    start_text,
    end_text,
    duration_days
  )
}

# =============================================================================
# PLOTTING HELPERS
# =============================================================================

.event_palette <- function(
    variable
) {
  
  if (variable == "heatindex") {
    return(
      HEAT_COLORS
    )
  }
  
  if (variable == "precipitation") {
    return(
      PRECIP_COLORS
    )
  }
  
  stop(
    "Unsupported variable: ",
    variable,
    call. = FALSE
  )
}

.add_common_map_elements <- function(
    p,
    event_start,
    event_end
) {
  
  duration_text <- .duration_label(
    event_start,
    event_end
  )
  
  p +
    ggplot2::coord_sf(
      crs = sf::st_crs(
        CRS_PLOT
      ),
      datum = sf::st_crs(
        CRS_PLOT
      ),
      xlim = c(
        BBOX_LONLAT[["xmin"]],
        BBOX_LONLAT[["xmax"]]
      ),
      ylim = c(
        BBOX_LONLAT[["ymin"]],
        BBOX_LONLAT[["ymax"]]
      ),
      expand = FALSE,
      label_axes = list(
        bottom = "E",
        left = "N"
      )
    ) +
    
    ggplot2::annotate(
      "label",
      x = BBOX_LONLAT[["xmin"]] + 0.55,
      y = BBOX_LONLAT[["ymin"]] + 0.55,
      label = duration_text,
      hjust = 0,
      vjust = 0,
      size = 4.5,
      family = "sans",
      lineheight = 1.05,
      label.size = 0.35,
      label.padding = grid::unit(
        0.20,
        "lines"
      ),
      fill = scales::alpha(
        "white",
        0.88
      ),
      color = "black"
    ) +
    
    ggplot2::labs(
      x = NULL,
      y = NULL
    ) +
    
    ggplot2::theme_bw(
      base_size = BASE_SIZE,
      base_family = "sans"
    ) +
    
    ggplot2::theme(
      text = ggplot2::element_text(
        family = "sans",
        size = BASE_SIZE
      ),
      panel.grid.major = ggplot2::element_line(
        color = "grey82",
        linewidth = 0.35
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(
        size = 14,
        color = "grey25"
      ),
      axis.ticks = ggplot2::element_line(
        color = "grey25",
        linewidth = 0.3
      ),
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.55
      ),
      legend.position = "none",
      plot.margin = ggplot2::margin(
        8,
        8,
        8,
        8
      )
    )
}

# =============================================================================
# EVENT-ONLY FIGURE
# =============================================================================

.build_event_plot <- function(
    variable,
    raster,
    footprint,
    states,
    event_start,
    event_end
) {
  
  raster_df <- .raster_to_df(
    raster
  )
  
  dx <- terra::res(
    raster
  )[1]
  
  dy <- terra::res(
    raster
  )[2]
  
  palette <- .event_palette(
    variable
  )
  
  p <- ggplot2::ggplot() +
    
    # -------------------------------------------------------------------------
  # State background
  # -------------------------------------------------------------------------
  
  ggplot2::geom_sf(
    data = states,
    fill = "white",
    color = "black",
    linewidth = 0.45,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
    
    # -------------------------------------------------------------------------
  # Return-period raster
  # -------------------------------------------------------------------------
  
  ggplot2::geom_tile(
    data = raster_df,
    mapping = ggplot2::aes(
      x = x,
      y = y,
      fill = rp_class
    ),
    width = dx,
    height = dy,
    show.legend = FALSE
  ) +
    
    ggplot2::scale_fill_manual(
      values = palette,
      limits = RP_LEVELS,
      drop = FALSE,
      na.value = "white",
      guide = "none"
    ) +
    
    # -------------------------------------------------------------------------
  # Event footprint
  # -------------------------------------------------------------------------
  
  ggplot2::geom_sf(
    data = footprint,
    fill = NA,
    color = "black",
    linewidth = 0.65,
    linejoin = "mitre",
    inherit.aes = FALSE,
    show.legend = FALSE
  )
  
  .add_common_map_elements(
    p = p,
    event_start = event_start,
    event_end = event_end
  )
}

# =============================================================================
# COUNTY-CONFUSION FIGURE
# =============================================================================

.build_confusion_plot <- function(
    footprint,
    states,
    county_points,
    event_start,
    event_end
) {
  
  # Do not plot TN counties.
  plotted_points <- county_points |>
    dplyr::filter(
      confusion_class %in% c(
        "TP",
        "FP",
        "FN"
      )
    ) |>
    dplyr::mutate(
      confusion_class = factor(
        as.character(
          confusion_class
        ),
        levels = c(
          "TP",
          "FP",
          "FN"
        )
      )
    )
  
  p <- ggplot2::ggplot() +
    
    # -------------------------------------------------------------------------
  # State background
  # -------------------------------------------------------------------------
  
  ggplot2::geom_sf(
    data = states,
    fill = "white",
    color = "black",
    linewidth = 0.45,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
    
    # -------------------------------------------------------------------------
  # Event footprint only.
  #
  # No return-period raster is plotted in the confusion figure.
  # -------------------------------------------------------------------------
  
  ggplot2::geom_sf(
    data = footprint,
    fill = NA,
    color = "black",
    linewidth = 0.65,
    linejoin = "mitre",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
    
    # -------------------------------------------------------------------------
  # TP / FP / FN county X symbols
  # -------------------------------------------------------------------------
  
  ggplot2::geom_sf(
    data = plotted_points,
    mapping = ggplot2::aes(
      color = confusion_class
    ),
    shape = 4,
    size = COUNTY_X_SIZE,
    stroke = COUNTY_X_STROKE,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
    
    ggplot2::scale_color_manual(
      values = CONFUSION_COLORS,
      limits = names(
        CONFUSION_COLORS
      ),
      drop = FALSE,
      guide = "none"
    )
  
  .add_common_map_elements(
    p = p,
    event_start = event_start,
    event_end = event_end
  )
}

# =============================================================================
# SAVE FIGURE
# =============================================================================

.save_png <- function(
    plot,
    path
) {
  
  ragg::agg_png(
    filename = path,
    width = FIG_WIDTH_IN,
    height = FIG_HEIGHT_IN,
    units = "in",
    res = DPI,
    background = "white"
  )
  
  print(
    plot
  )
  
  grDevices::dev.off()
  
  .log(
    "Wrote: %s",
    path
  )
}

# =============================================================================
# PROCESS ONE EVENT
# =============================================================================

.process_event <- function(
    variable,
    dataset_tag,
    cluster_id,
    phenom,
    sig,
    states,
    counties
) {
  
  .log(
    "============================================================"
  )
  
  .log(
    "Processing %s",
    cluster_id
  )
  
  # ---------------------------------------------------------------------------
  # Event timing
  # ---------------------------------------------------------------------------
  
  char <- .read_event_characteristics(
    dataset_tag = dataset_tag,
    cluster_id = cluster_id
  )
  
  event_start <- char$start
  event_end <- char$end
  
  .log(
    "Event interval: %s to %s UTC",
    format(
      event_start,
      "%Y-%m-%d %H:%M"
    ),
    format(
      event_end,
      "%Y-%m-%d %H:%M"
    )
  )
  
  # ---------------------------------------------------------------------------
  # Raster and footprint
  # ---------------------------------------------------------------------------
  
  raster <- .read_cluster_raster(
    dataset_tag = dataset_tag,
    cluster_id = cluster_id
  )
  
  footprint <- .cluster_footprint(
    raster
  )
  
  # ---------------------------------------------------------------------------
  # Warnings overlapping event duration
  # ---------------------------------------------------------------------------
  
  warnings <- .read_warning_data(
    event_start = event_start,
    event_end = event_end,
    phenom = phenom,
    sig = sig
  )
  
  # ---------------------------------------------------------------------------
  # County classifications
  # ---------------------------------------------------------------------------
  
  classified_counties <- .classify_counties(
    counties = counties,
    cluster_footprint = footprint,
    warnings = warnings
  )
  
  county_points <- .make_county_points(
    classified_counties
  )
  
  counts <- table(
    factor(
      classified_counties$confusion_class,
      levels = c(
        "TP",
        "TN",
        "FP",
        "FN"
      )
    )
  )
  
  .log(
    "County confusion counts: TP=%d | TN=%d | FP=%d | FN=%d",
    counts[["TP"]],
    counts[["TN"]],
    counts[["FP"]],
    counts[["FN"]]
  )
  
  # ---------------------------------------------------------------------------
  # Metrics
  # ---------------------------------------------------------------------------
  
  event_info <- list(
    variable = variable,
    dataset_tag = dataset_tag,
    cluster_id = cluster_id,
    phenom = phenom,
    sig = sig,
    event_start = event_start,
    event_end = event_end
  )
  
  metrics <- .calculate_metrics(
    counties = classified_counties,
    event = event_info
  )
  
  .log(
    "Recall=%.4f | Precision=%.4f | F1=%.4f",
    metrics$recall,
    metrics$precision,
    metrics$f1
  )
  
  # ---------------------------------------------------------------------------
  # Output-name stem
  # ---------------------------------------------------------------------------
  
  short_name <- if (
    variable == "heatindex"
  ) {
    "heatindex"
  } else {
    "precipitation"
  }
  
  # ---------------------------------------------------------------------------
  # Figure 1: RP-colored example event
  # ---------------------------------------------------------------------------
  
  p_event <- .build_event_plot(
    variable = variable,
    raster = raster,
    footprint = footprint,
    states = states,
    event_start = event_start,
    event_end = event_end
  )
  
  event_png <- file.path(
    OUT_DIR,
    paste0(
      "dissertation_01_examples_",
      short_name,
      "_event.png"
    )
  )
  
  .save_png(
    p_event,
    event_png
  )
  
  # ---------------------------------------------------------------------------
  # Figure 2: confusion X symbols
  #
  # No RP raster.
  # TP = red
  # FP = blue
  # FN = black
  # TN = not plotted
  # ---------------------------------------------------------------------------
  
  p_confusion <- .build_confusion_plot(
    footprint = footprint,
    states = states,
    county_points = county_points,
    event_start = event_start,
    event_end = event_end
  )
  
  confusion_png <- file.path(
    OUT_DIR,
    paste0(
      "dissertation_01_examples_",
      short_name,
      "_confusion.png"
    )
  )
  
  .save_png(
    p_confusion,
    confusion_png
  )
  
  # ---------------------------------------------------------------------------
  # Metrics CSV
  # ---------------------------------------------------------------------------
  
  metrics_csv <- file.path(
    OUT_DIR,
    paste0(
      "dissertation_01_examples_",
      short_name,
      "_metrics.csv"
    )
  )
  
  readr::write_csv(
    metrics,
    metrics_csv,
    na = ""
  )
  
  .log(
    "Wrote: %s",
    metrics_csv
  )
  
  invisible(
    list(
      metrics = metrics,
      counties = classified_counties
    )
  )
}

# =============================================================================
# LOAD COMMON SPATIAL DATA
# =============================================================================

bbox_poly <- .make_bbox()

.log(
  "Loading states."
)

states <- .load_states()

.log(
  "Loading counties."
)

counties <- .load_domain_counties(
  bbox_poly
)

.log(
  "Evaluation domain contains %d counties.",
  nrow(counties)
)

# =============================================================================
# RUN BOTH EXAMPLE EVENTS
# =============================================================================

results <- vector(
  "list",
  nrow(EVENTS)
)

for (i in seq_len(
  nrow(EVENTS)
)) {
  
  event <- EVENTS[
    i,
  ]
  
  results[[i]] <- .process_event(
    variable = event$variable[[1]],
    dataset_tag = event$dataset_tag[[1]],
    cluster_id = event$cluster_id[[1]],
    phenom = event$phenom[[1]],
    sig = event$sig[[1]],
    states = states,
    counties = counties
  )
}

names(results) <- EVENTS$variable

.log(
  "============================================================"
)

.log(
  "Complete."
)

.log(
  "Outputs written to: %s",
  OUT_DIR
)