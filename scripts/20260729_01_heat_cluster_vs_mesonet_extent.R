#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(terra)
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  library(tigris)
  library(grid)
})

options(
  tigris_use_cache = TRUE,
  tigris_class = "sf"
)

# IEM VTEC shapefiles can contain self-crossing rings. This workflow uses
# planar GEOS operations and performs county topology tests in EPSG:5070.
suppressMessages(
  sf::sf_use_s2(FALSE)
)

# =============================================================================
# Heat-index cluster resampling comparison with county-centered IEM
# warning overlays.
#
# Behavior:
#   1. Read five selected heat-index clusters.
#   2. Determine the duration of each cluster from the filtered-rp5 Parquet.
#   3. Compute the strict temporal intersection shared by all five clusters.
#   4. Read each required IEM annual WWA shapefile only once.
#   5. Retain heat-related warnings that overlap the strict shared cluster interval.
#   6. Loop over every requested PHENOM x SIG combination.
#   7. For each combination, identify counties intersected by at least one
#      matching warning polygon.
#   8. Plot the same affected-county blue x markers on all five panels.
#   9. Write one five-panel, one-row figure per PHENOM x SIG combination.
#
# Layout:
#
#   Native | Bilinear gamma2 | Maximum gamma2 | Bilinear gamma3 | Maximum gamma3
#
# Inputs:
#   data/out/08_tifs_pngs/geotiffs/<dataset_tag>/<cluster_id>.tif
#   data/out/06_clusters/filtered_rp5/<dataset_tag>/cluster_rp5_<year>.parquet
#   data/in/mesonet/<year>_all/wwa_<year>01010000_<year>12312359.shp
#
# Outputs:
#   figures/08_tifs_pngs_panels/
#     fig_heatindex_cluster_resampling_<PHENOM>_<SIG>.png
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

CLUSTER_ROOT <- file.path(
  ROOT,
  "data/out/06_clusters/filtered_rp5"
)

MESONET_ROOT <- file.path(
  ROOT,
  "data/in/mesonet"
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

# Each combination receives a separate output figure.
#
# Heat-related PHENOM codes:
#   EH = Excessive Heat
#   HT = Heat
#   XH = Extreme Heat
WWA_PHENOM_VALUES <- c(
  "EH",
  "HT",
  "XH"
)

# Significance codes:
#   W = Warning
#   A = Watch
#   Y = Advisory
#   S = Statement
WWA_SIG_VALUES <- c(
  "W",
  "A",
  "Y",
  "S"
)

PHENOM_NAMES <- c(
  EH = "Excessive Heat",
  HT = "Heat",
  XH = "Extreme Heat"
)

SIG_NAMES <- c(
  W = "Warning",
  A = "Watch",
  Y = "Advisory",
  S = "Statement"
)

# If TRUE, combinations with zero matching warning records are logged and
# skipped. If FALSE, a figure is still written with no blue x markers.
SKIP_EMPTY_COMBINATIONS <- TRUE
SKIP_EXISTING_FIGURES <- FALSE

# Restrict warning and county processing to the union extent of all five
# selected rasters.
LIMIT_WARNINGS_TO_COMBINED_EXTENT <- TRUE

# Buffer around the combined raster extent. Expected units are degrees because
# the current rasters use longitude/latitude.
COMBINED_EXTENT_BUFFER <- 0.50

WARNING_X_COLOUR <- "#2166ac"
WARNING_X_SIZE <- 1.10
WARNING_X_STROKE <- 0.42

# Optional CSV override.
#
# Required columns:
#   dataset_tag,cluster_id
SELECTION_CSV <- Sys.getenv(
  "CISC_HEAT_PANEL_SELECTION_CSV",
  ""
)

FIG_WIDTH_IN <- 15.0
FIG_HEIGHT_IN <- 3.8
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
    expected_tags
) {
  missing_tags <- setdiff(
    expected_tags,
    names(selection)
  )
  
  if (length(missing_tags)) {
    stop(
      "Selection is missing dataset tags: ",
      paste(missing_tags, collapse = ", "),
      call. = FALSE
    )
  }
  
  selected <- trimws(
    unname(selection[expected_tags])
  )
  
  bad <- is.na(selected) | !nzchar(selected)
  
  if (any(bad)) {
    stop(
      "Selection contains missing or empty cluster IDs for: ",
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
    transmute(
      dataset_tag = trimws(as.character(.data$dataset_tag)),
      cluster_id = trimws(as.character(.data$cluster_id))
    ) |>
    filter(
      nzchar(.data$dataset_tag),
      nzchar(.data$cluster_id)
    )
  
  setNames(
    tab$cluster_id,
    tab$dataset_tag
  )
}

.cluster_tif_path <- function(
    dataset_tag,
    cluster_id
) {
  file.path(
    GEOTIFF_ROOT,
    dataset_tag,
    paste0(cluster_id, ".tif")
  )
}

.cluster_year_from_id <- function(cluster_id) {
  year_text <- stringr::str_extract(
    cluster_id,
    "(?<=-)\\d{4}(?=-\\d{2}-\\d{2}-\\d{3}$)"
  )
  
  if (is.na(year_text)) {
    stop(
      "Could not extract year from cluster ID: ",
      cluster_id,
      call. = FALSE
    )
  }
  
  as.integer(year_text)
}

.cluster_parquet_path <- function(
    dataset_tag,
    year
) {
  file.path(
    CLUSTER_ROOT,
    dataset_tag,
    sprintf(
      "cluster_rp5_%d.parquet",
      year
    )
  )
}

.read_cluster_interval <- function(
    dataset_tag,
    cluster_id
) {
  year <- .cluster_year_from_id(cluster_id)
  path <- .cluster_parquet_path(dataset_tag, year)
  
  if (!file.exists(path)) {
    stop("Cluster Parquet not found: ", path, call. = FALSE)
  }
  
  cluster_tab <- arrow::read_parquet(
    path,
    col_select = dplyr::any_of(
      c(
        "cluster_id",
        "time_utc"
      )
    )
  ) |>
    dplyr::filter(
      .data$cluster_id == .env$cluster_id
    ) |>
    dplyr::mutate(
      time_utc = as.POSIXct(
        .data$time_utc,
        tz = "UTC"
      )
    ) |>
    dplyr::filter(
      !is.na(.data$time_utc)
    )
  
  if (!nrow(cluster_tab)) {
    stop(
      "Cluster ID was not found, or has no valid time_utc values, in ",
      path,
      ": ",
      cluster_id,
      call. = FALSE
    )
  }
  
  retained_ids <- unique(
    as.character(cluster_tab$cluster_id)
  )
  
  if (
    length(retained_ids) != 1L ||
    !identical(retained_ids[[1]], cluster_id)
  ) {
    stop(
      "Cluster filtering failed for dataset ",
      dataset_tag,
      ". Requested cluster ID: ",
      cluster_id,
      "; retained cluster IDs: ",
      paste(retained_ids, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Use the observed cluster timestamps exactly. Do not add time_step_hours:
  # the warning-evaluation window is the temporal intersection shared by all
  # five selected clusters, matching range(cluster_rows$time_utc).
  tibble::tibble(
    dataset_tag = dataset_tag,
    cluster_id = cluster_id,
    start_utc = min(cluster_tab$time_utc),
    end_utc = max(cluster_tab$time_utc),
    n_unique_times = dplyr::n_distinct(cluster_tab$time_utc)
  )
}

.intersect_cluster_intervals <- function(
    selection,
    dataset_tags
) {
  individual <- lapply(
    dataset_tags,
    function(tag) {
      .read_cluster_interval(
        dataset_tag = tag,
        cluster_id = selection[[tag]]
      )
    }
  ) |>
    bind_rows() |>
    arrange(.data$dataset_tag)
  
  common_start <- max(individual$start_utc)
  common_end <- min(individual$end_utc)
  
  if (common_start > common_end) {
    details <- paste0(
      individual$dataset_tag,
      ": ",
      format(individual$start_utc, "%Y-%m-%d %H:%M UTC"),
      " to ",
      format(individual$end_utc, "%Y-%m-%d %H:%M UTC"),
      collapse = "
"
    )
    
    stop(
      "The five selected clusters have no shared temporal intersection.
",
      details,
      call. = FALSE
    )
  }
  
  intersection <- tibble::tibble(
    start_utc = common_start,
    end_utc = common_end
  )
  
  list(
    individual = individual,
    intersection = intersection,
    years = seq(lubridate::year(common_start), lubridate::year(common_end))
  )
}

.mesonet_shapefile_path <- function(year) {
  file.path(
    MESONET_ROOT,
    sprintf("%d_all", year),
    sprintf(
      "wwa_%d01010000_%d12312359.shp",
      year,
      year
    )
  )
}

.parse_wwa_datetime <- function(x) {
  if (inherits(x, "POSIXt")) {
    return(as.POSIXct(x, tz = "UTC"))
  }
  
  text <- trimws(as.character(x))
  text[!nzchar(text) | text %in% c("NA", "<NA>")] <- NA_character_
  
  suppressWarnings(
    as.POSIXct(
      lubridate::parse_date_time(
        text,
        orders = c(
          "ymd HMS",
          "ymd HM",
          "Ymd HMS",
          "Ymd HM",
          "YmdHMS",
          "YmdHM"
        ),
        tz = "UTC",
        quiet = TRUE
      ),
      tz = "UTC"
    )
  )
}

.warning_overlaps_window <- function(
    issued_utc,
    expired_utc,
    evaluation_window
) {
  issued_utc <= evaluation_window$end_utc[[1]] &
    expired_utc >= evaluation_window$start_utc[[1]]
}

.read_mesonet_warning_cache <- function(
    evaluation_intervals,
    combined_bbox,
    target_crs
) {
  evaluation_window <- evaluation_intervals$intersection
  pieces <- vector("list", length(evaluation_intervals$years))
  
  .log(
    "Warning evaluation uses strict shared interval: %s to %s UTC.",
    format(evaluation_window$start_utc[[1]], "%Y-%m-%d %H:%M"),
    format(evaluation_window$end_utc[[1]], "%Y-%m-%d %H:%M")
  )
  
  for (i in seq_along(evaluation_intervals$years)) {
    year <- evaluation_intervals$years[[i]]
    path <- .mesonet_shapefile_path(year)
    
    if (!file.exists(path)) {
      stop("Mesonet warning shapefile not found: ", path, call. = FALSE)
    }
    
    .log("Reading Mesonet warning file once for %d.", year)
    
    wwa <- sf::st_read(
      path,
      quiet = TRUE,
      stringsAsFactors = FALSE
    )
    
    required <- c("ISSUED", "EXPIRED", "PHENOM", "SIG")
    missing_cols <- setdiff(required, names(wwa))
    
    if (length(missing_cols)) {
      stop(
        "Mesonet file is missing required fields: ",
        paste(missing_cols, collapse = ", "),
        "
File: ",
        path,
        call. = FALSE
      )
    }
    
    wwa <- wwa |>
      mutate(
        PHENOM = toupper(trimws(as.character(.data$PHENOM))),
        SIG = toupper(trimws(as.character(.data$SIG))),
        issued_utc = .parse_wwa_datetime(.data$ISSUED),
        expired_utc = .parse_wwa_datetime(.data$EXPIRED)
      ) |>
      filter(
        .data$PHENOM %in% WWA_PHENOM_VALUES,
        .data$SIG %in% WWA_SIG_VALUES,
        !is.na(.data$issued_utc),
        !is.na(.data$expired_utc),
        .data$expired_utc >= .data$issued_utc
      )
    
    if (!nrow(wwa)) {
      pieces[[i]] <- NULL
      next
    }
    
    overlap <- .warning_overlaps_window(
      issued_utc = wwa$issued_utc,
      expired_utc = wwa$expired_utc,
      evaluation_window = evaluation_window
    )
    
    wwa <- wwa[overlap, ]
    
    if (!nrow(wwa)) {
      pieces[[i]] <- NULL
      next
    }
    
    # Hard validation: no record can survive unless it overlaps the strict
    # shared interval.
    still_overlaps <- .warning_overlaps_window(
      issued_utc = wwa$issued_utc,
      expired_utc = wwa$expired_utc,
      evaluation_window = evaluation_window
    )
    
    if (!all(still_overlaps)) {
      stop("Internal error: non-overlapping warnings survived time filtering.", call. = FALSE)
    }
    
    wwa <- suppressWarnings(sf::st_make_valid(wwa))
    wwa <- wwa[!sf::st_is_empty(wwa), ]
    
    if (!nrow(wwa)) {
      pieces[[i]] <- NULL
      next
    }
    
    wwa <- sf::st_transform(wwa, target_crs)
    wwa <- suppressWarnings(sf::st_make_valid(wwa))
    wwa <- wwa[!sf::st_is_empty(wwa), ]
    
    if (LIMIT_WARNINGS_TO_COMBINED_EXTENT && nrow(wwa)) {
      wwa <- suppressWarnings(sf::st_crop(wwa, combined_bbox))
      wwa <- wwa[!sf::st_is_empty(wwa), ]
    }
    
    pieces[[i]] <- if (nrow(wwa)) wwa else NULL
  }
  
  pieces <- Filter(Negate(is.null), pieces)
  
  if (!length(pieces)) {
    return(NULL)
  }
  
  warnings <- do.call(rbind, pieces)
  warnings <- warnings[!sf::st_is_empty(warnings), ]
  
  if (!nrow(warnings)) {
    return(NULL)
  }
  
  .log(
    "Cached warnings span issued times %s to %s UTC and expiration times %s to %s UTC.",
    format(min(warnings$issued_utc), "%Y-%m-%d %H:%M"),
    format(max(warnings$issued_utc), "%Y-%m-%d %H:%M"),
    format(min(warnings$expired_utc), "%Y-%m-%d %H:%M"),
    format(max(warnings$expired_utc), "%Y-%m-%d %H:%M")
  )
  
  warnings
}

.load_states_and_counties <- function(target_crs) {
  excluded_state_abbr <- c(
    "AK",
    "HI",
    "PR",
    "VI",
    "GU",
    "MP",
    "AS"
  )
  
  excluded_state_fips <- c(
    "02",
    "15",
    "60",
    "66",
    "69",
    "72",
    "78"
  )
  
  states_sf <- tigris::states(
    cb = TRUE,
    year = 2023,
    class = "sf"
  ) |>
    filter(
      !.data$STUSPS %in% excluded_state_abbr
    ) |>
    sf::st_transform(
      target_crs
    )
  
  counties_sf <- tigris::counties(
    cb = TRUE,
    year = 2023,
    class = "sf"
  ) |>
    filter(
      !.data$STATEFP %in% excluded_state_fips
    ) |>
    select(
      GEOID,
      NAME,
      STATEFP,
      geometry
    ) |>
    sf::st_transform(
      target_crs
    )
  
  list(
    states = states_sf,
    counties = counties_sf
  )
}

.warning_county_points <- function(
    warnings,
    counties,
    combined_bbox,
    target_crs
) {
  if (is.null(warnings) || !nrow(warnings)) {
    return(NULL)
  }
  
  candidate_counties <- suppressWarnings(
    sf::st_crop(counties, combined_bbox)
  )
  
  if (!nrow(candidate_counties)) {
    return(NULL)
  }
  
  area_crs <- "EPSG:5070"
  
  warnings_area <- warnings |>
    sf::st_transform(area_crs) |>
    sf::st_make_valid()
  
  counties_area <- candidate_counties |>
    sf::st_transform(area_crs) |>
    sf::st_make_valid()
  
  warnings_area <- warnings_area[!sf::st_is_empty(warnings_area), ]
  counties_area <- counties_area[!sf::st_is_empty(counties_area), ]
  
  if (!nrow(warnings_area) || !nrow(counties_area)) {
    return(NULL)
  }
  
  # Binary county classification only. Avoid st_union(), st_intersection(),
  # and st_area(), which are unnecessary and fragile for malformed WWA rings.
  affected_index <- suppressWarnings(
    lengths(
      sf::st_intersects(
        counties_area,
        warnings_area,
        sparse = TRUE
      )
    ) > 0L
  )
  
  affected <- counties_area[affected_index, ]
  
  if (!nrow(affected)) {
    return(NULL)
  }
  
  suppressWarnings(sf::st_point_on_surface(affected)) |>
    sf::st_transform(target_crs)
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
  
  df |>
    mutate(
      rp_class = cut(
        .data$rp,
        breaks = RP_CUT_BREAKS,
        labels = RP_LEVELS,
        right = FALSE,
        include.lowest = TRUE
      ),
      rp_class = factor(
        .data$rp_class,
        levels = RP_LEVELS,
        ordered = TRUE
      )
    )
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
  
  sf::st_as_sf(
    footprint
  ) |>
    sf::st_make_valid()
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

.prepare_panel_data <- function(
    dataset_tag,
    cluster_id
) {
  path <- .cluster_tif_path(
    dataset_tag,
    cluster_id
  )
  
  obj <- .read_cluster_raster(
    path
  )
  
  list(
    dataset_tag = dataset_tag,
    cluster_id = cluster_id,
    raster_path = path,
    raster = obj$raster,
    extent = obj$extent
  )
}

.combined_bbox_from_panels <- function(
    panel_data,
    target_crs
) {
  xmin <- min(
    vapply(
      panel_data,
      function(x) terra::xmin(x$extent),
      numeric(1)
    )
  )
  
  xmax <- max(
    vapply(
      panel_data,
      function(x) terra::xmax(x$extent),
      numeric(1)
    )
  )
  
  ymin <- min(
    vapply(
      panel_data,
      function(x) terra::ymin(x$extent),
      numeric(1)
    )
  )
  
  ymax <- max(
    vapply(
      panel_data,
      function(x) terra::ymax(x$extent),
      numeric(1)
    )
  )
  
  sf::st_bbox(
    c(
      xmin = xmin - COMBINED_EXTENT_BUFFER,
      ymin = ymin - COMBINED_EXTENT_BUFFER,
      xmax = xmax + COMBINED_EXTENT_BUFFER,
      ymax = ymax + COMBINED_EXTENT_BUFFER
    ),
    crs = sf::st_crs(target_crs)
  )
}

.panel_plot <- function(
    panel_data,
    states,
    warning_points,
    panel_title,
    show_x = TRUE,
    show_y = TRUE
) {
  r <- panel_data$raster
  ex <- panel_data$extent
  
  df <- .raster_to_df(
    r
  )
  
  nonzero_outline <- .nonzero_outline_sf(
    r
  )
  
  dx <- terra::res(r)[1]
  dy <- terra::res(r)[2]
  
  p <- ggplot() +
    geom_tile(
      data = df,
      mapping = aes(
        x = .data$x,
        y = .data$y,
        fill = .data$rp_class
      ),
      width = dx,
      height = dy,
      inherit.aes = FALSE,
      show.legend = TRUE
    )
  
  if (!is.null(nonzero_outline)) {
    p <- p +
      geom_sf(
        data = nonzero_outline,
        fill = NA,
        colour = "black",
        linewidth = 0.45,
        linejoin = "mitre",
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
  
  if (
    !is.null(warning_points) &&
    nrow(warning_points)
  ) {
    p <- p +
      geom_sf(
        data = warning_points,
        shape = 4,
        colour = WARNING_X_COLOUR,
        size = WARNING_X_SIZE,
        stroke = WARNING_X_STROKE,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
  
  if (!is.null(states)) {
    p <- p +
      geom_sf(
        data = states,
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

.build_combination_figure <- function(
    panel_data,
    states,
    warning_points,
    evaluation_intervals,
    selected_phenom,
    selected_sig,
    out_file
) {
  native_tag <- "heatindex_0.25"
  gamma2_bilinear_tag <- "heatindex_0.30_bilinear"
  gamma2_max_tag <- "heatindex_0.30_max"
  gamma3_bilinear_tag <- "heatindex_0.36_bilinear"
  gamma3_max_tag <- "heatindex_0.36_max"
  
  p_native <- .panel_plot(
    panel_data = panel_data[[native_tag]],
    states = states,
    warning_points = warning_points,
    panel_title = "Native (0.25°)",
    show_x = TRUE,
    show_y = TRUE
  )
  
  p_g2_bilinear <- .panel_plot(
    panel_data = panel_data[[gamma2_bilinear_tag]],
    states = states,
    warning_points = warning_points,
    panel_title = "Bilinear γ₂ (0.30°)",
    show_x = TRUE,
    show_y = FALSE
  )
  
  p_g2_max <- .panel_plot(
    panel_data = panel_data[[gamma2_max_tag]],
    states = states,
    warning_points = warning_points,
    panel_title = "Maximum γ₂ (0.30°)",
    show_x = TRUE,
    show_y = FALSE
  )
  
  p_g3_bilinear <- .panel_plot(
    panel_data = panel_data[[gamma3_bilinear_tag]],
    states = states,
    warning_points = warning_points,
    panel_title = "Bilinear γ₃ (0.36°)",
    show_x = TRUE,
    show_y = FALSE
  )
  
  p_g3_max <- .panel_plot(
    panel_data = panel_data[[gamma3_max_tag]],
    states = states,
    warning_points = warning_points,
    panel_title = "Maximum γ₃ (0.36°)",
    show_x = TRUE,
    show_y = FALSE
  )
  
  design <- "
ABCDE
FFFFF
"
  
  phenom_name <- unname(
    PHENOM_NAMES[selected_phenom]
  )
  
  sig_name <- unname(
    SIG_NAMES[selected_sig]
  )
  
  if (
    is.na(phenom_name) ||
    !nzchar(phenom_name)
  ) {
    phenom_name <- selected_phenom
  }
  
  if (
    is.na(sig_name) ||
    !nzchar(sig_name)
  ) {
    sig_name <- selected_sig
  }
  
  date_start <- format(
    evaluation_intervals$intersection$start_utc[[1]],
    "%Y-%m-%d"
  )
  
  date_end <- format(
    evaluation_intervals$intersection$end_utc[[1]],
    "%Y-%m-%d"
  )
  
  figure_title <- sprintf(
    "%s (%s) %s (%s) for %s to %s",
    phenom_name,
    selected_phenom,
    sig_name,
    selected_sig,
    date_start,
    date_end
  )
  
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
        widths = rep(1, 5),
        heights = c(
          1,
          0.13
        )
      ) +
      plot_annotation(
        title = figure_title,
        theme = theme(
          plot.title = element_text(
            size = BASE_SIZE + 1,
            face = "bold",
            hjust = 0.5,
            margin = margin(
              b = 2
            )
          )
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

native_tag <- "heatindex_0.25"
gamma2_bilinear_tag <- "heatindex_0.30_bilinear"
gamma2_max_tag <- "heatindex_0.30_max"
gamma3_bilinear_tag <- "heatindex_0.36_bilinear"
gamma3_max_tag <- "heatindex_0.36_max"

expected_tags <- c(
  native_tag,
  gamma2_bilinear_tag,
  gamma2_max_tag,
  gamma3_bilinear_tag,
  gamma3_max_tag
)

heat_selection <- HEAT_CLUSTER_IDS

if (nzchar(SELECTION_CSV)) {
  if (!file.exists(SELECTION_CSV)) {
    stop(
      "Selection CSV does not exist: ",
      SELECTION_CSV,
      call. = FALSE
    )
  }
  
  heat_selection <- .read_selection_csv(
    SELECTION_CSV
  )
}

.validate_selection(
  heat_selection,
  expected_tags
)

panel_data <- lapply(
  expected_tags,
  function(tag) {
    .prepare_panel_data(
      dataset_tag = tag,
      cluster_id = heat_selection[[tag]]
    )
  }
)

names(panel_data) <- expected_tags

native_crs <- terra::crs(
  panel_data[[native_tag]]$raster
)

boundaries <- .load_states_and_counties(
  native_crs
)

combined_bbox <- .combined_bbox_from_panels(
  panel_data = panel_data,
  target_crs = native_crs
)

evaluation_intervals <- .intersect_cluster_intervals(
  selection = heat_selection,
  dataset_tags = expected_tags
)

.log("Individual cluster intervals:")

for (i in seq_len(nrow(evaluation_intervals$individual))) {
  .log(
    "%s | %s to %s UTC | %d unique time(s)",
    evaluation_intervals$individual$cluster_id[[i]],
    format(
      evaluation_intervals$individual$start_utc[[i]],
      "%Y-%m-%d %H:%M"
    ),
    format(
      evaluation_intervals$individual$end_utc[[i]],
      "%Y-%m-%d %H:%M"
    ),
    evaluation_intervals$individual$n_unique_times[[i]]
  )
}

.log(
  "Strict shared interval across all five clusters: %s to %s UTC.",
  format(evaluation_intervals$intersection$start_utc[[1]], "%Y-%m-%d %H:%M"),
  format(evaluation_intervals$intersection$end_utc[[1]], "%Y-%m-%d %H:%M")
)

warning_cache <- .read_mesonet_warning_cache(
  evaluation_intervals = evaluation_intervals,
  combined_bbox = combined_bbox,
  target_crs = native_crs
)

if (is.null(warning_cache)) {
  .log(
    paste0(
      "No Mesonet warning records matched the requested PHENOM/SIG values ",
      "during the strict shared cluster interval."
    )
  )
} else {
  .log(
    "Cached %s warning records across all requested combinations.",
    format(
      nrow(warning_cache),
      big.mark = ","
    )
  )
}

warning_combinations <- tidyr::crossing(
  PHENOM = WWA_PHENOM_VALUES,
  SIG = WWA_SIG_VALUES
)

n_written <- 0L
n_skipped <- 0L

for (i in seq_len(nrow(warning_combinations))) {
  selected_phenom <- warning_combinations$PHENOM[[i]]
  selected_sig <- warning_combinations$SIG[[i]]
  
  if (is.null(warning_cache)) {
    combination_warnings <- NULL
  } else {
    combination_warnings <- warning_cache |>
      filter(
        .data$PHENOM == selected_phenom,
        .data$SIG == selected_sig
      )
  }
  
  n_warning_records <- if (is.null(combination_warnings)) {
    0L
  } else {
    nrow(combination_warnings)
  }
  
  if (
    SKIP_EMPTY_COMBINATIONS &&
    n_warning_records == 0L
  ) {
    .log(
      "%s_%s | no matching warning records; skipping.",
      selected_phenom,
      selected_sig
    )
    
    n_skipped <- n_skipped + 1L
    next
  }
  
  warning_points <- .warning_county_points(
    warnings = combination_warnings,
    counties = boundaries$counties,
    combined_bbox = combined_bbox,
    target_crs = native_crs
  )
  
  n_affected_counties <- if (is.null(warning_points)) {
    0L
  } else {
    nrow(warning_points)
  }
  
  .log(
    "%s_%s | %s warning records | %s affected counties",
    selected_phenom,
    selected_sig,
    format(
      n_warning_records,
      big.mark = ","
    ),
    format(
      n_affected_counties,
      big.mark = ","
    )
  )
  
  out_file <- file.path(
    OUT_DIR,
    sprintf(
      "fig_heatindex_cluster_resampling_%s_%s.png",
      selected_phenom,
      selected_sig
    )
  )
  
  if (SKIP_EXISTING_FIGURES && file.exists(out_file)) {
    .log(
      "%s_%s | output already exists; skipping.",
      selected_phenom,
      selected_sig
    )
    n_skipped <- n_skipped + 1L
    next
  }
  
  .build_combination_figure(
    panel_data = panel_data,
    states = boundaries$states,
    warning_points = warning_points,
    evaluation_intervals = evaluation_intervals,
    selected_phenom = selected_phenom,
    selected_sig = selected_sig,
    out_file = out_file
  )
  
  n_written <- n_written + 1L
}

.log(
  "Complete: wrote %d figure(s) and skipped %d empty combination(s).",
  n_written,
  n_skipped
)