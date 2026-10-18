#!/usr/bin/env Rscript

# =============================================================================
# Per-cluster characteristics for year-chunked clustered exceedance Parquet files
#
# INPUTS:
#   /proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/
#     filtered_rp{5,10,25,50}/
#       <cluster_set>/
#         cluster_rp<filter>_<year>.parquet
#
# OUTPUTS:
#   /proj/sebastianlab/users/hquintal/cisc/data/out/07_characteristics/
#     filtered_rp{5,10,25,50}/
#       <cluster_set>/
#         cluster_characteristics.parquet
#         cluster_characteristics.csv
#
# Also writes one summary per event-filter directory:
#   filtered_rp{5,10,25,50}/cluster_summary_by_resolution.parquet
#   filtered_rp{5,10,25,50}/cluster_summary_by_resolution.csv
#
# Important:
#   - timestep_hr is never inferred from exceedance spacing.
#   - timestep_hr comes from time_step_hours, otherwise defaults to 24.
#   - n_steps and duration_hr are based on integer timestep indices.
# =============================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(fs)
  library(future)
  library(future.apply)
  library(geosphere)
  library(lubridate)
  library(progressr)
})

# =============================================================================
# USER SETTINGS
# =============================================================================

ROOT <- Sys.getenv(
  "CISC_ROOT",
  "/proj/sebastianlab/users/hquintal/cisc"
)

BASE_IN <- fs::path(
  ROOT,
  "data/out/06_clusters"
)

BASE_OUT <- fs::path(
  ROOT,
  "data/out/07_characteristics"
)

FILTER_RP_YEARS <- c(5, 10, 25, 50)

CLUSTER_SETS <- c(
  "heatindex_0.25",
  "heatindex_0.30_bilinear",
  "heatindex_0.30_max",
  "heatindex_0.36_bilinear",
  "heatindex_0.36_max",
  "precipitation_0.25",
  "precipitation_0.11_bilinear",
  "precipitation_0.11_max",
  "precipitation_0.14_bilinear",
  "precipitation_0.14_max"
)

N_WORKERS <- as.integer(
  Sys.getenv("CISC_CHARACTERISTICS_WORKERS", "8")
)

DEFAULT_TIMESTEP_HR <- 24

# =============================================================================
# LOGGING
# =============================================================================

.msg <- function(...) {
  message(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf(...)
  )
}

# =============================================================================
# GEO HELPERS
# =============================================================================

km_per_deg_lat <- function(lat) {
  rep(111.132, length(lat))
}

km_per_deg_lon <- function(lat) {
  111.320 * cospi(lat / 180)
}

cell_area_km2_vec <- function(lat_vec, resx_deg, resy_deg) {
  km_y <- km_per_deg_lat(lat_vec) * resy_deg
  km_x <- km_per_deg_lon(lat_vec) * resx_deg
  km_x * km_y
}

grid_step <- function(v) {
  u <- sort(unique(v[is.finite(v)]))
  
  if (length(u) < 2L) {
    return(NA_real_)
  }
  
  d <- diff(u)
  d <- d[is.finite(d) & d > 0]
  
  if (!length(d)) {
    return(NA_real_)
  }
  
  min(d)
}

as_posix_utc <- function(tt) {
  if (inherits(tt, "POSIXct")) {
    return(as.POSIXct(tt, tz = "UTC"))
  }
  
  if (inherits(tt, "Date")) {
    return(as.POSIXct(tt, tz = "UTC"))
  }
  
  suppressWarnings(as.POSIXct(tt, tz = "UTC"))
}

gc_km <- function(p1, p2) {
  geosphere::distHaversine(p1, p2) / 1000
}

centroid_timestep <- function(x, y, w) {
  keep <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  
  if (!any(keep)) {
    return(c(lon = NA_real_, lat = NA_real_))
  }
  
  c(
    lon = weighted.mean(x[keep], w[keep]),
    lat = weighted.mean(y[keep], w[keep])
  )
}

elongation_major_minor <- function(
    xs,
    ys,
    resx_deg,
    resy_deg
) {
  keep <- is.finite(xs) & is.finite(ys)
  xs <- xs[keep]
  ys <- ys[keep]
  
  n <- length(xs)
  
  if (!n) {
    return(list(
      major_km = NA_real_,
      minor_km = NA_real_,
      elongation_ratio = NA_real_
    ))
  }
  
  lat0 <- mean(ys)
  
  cell_width_km <- km_per_deg_lon(lat0) * resx_deg
  cell_height_km <- km_per_deg_lat(lat0) * resy_deg
  
  base_major <- suppressWarnings(
    max(cell_width_km, cell_height_km, na.rm = TRUE)
  )
  base_minor <- suppressWarnings(
    min(cell_width_km, cell_height_km, na.rm = TRUE)
  )
  
  if (!is.finite(base_major) || base_major < 0) {
    base_major <- 0
  }
  
  if (!is.finite(base_minor) || base_minor <= 0) {
    base_minor <- 1e-6
  }
  
  if (n < 2L) {
    return(list(
      major_km = base_major,
      minor_km = base_minor,
      elongation_ratio = base_major / base_minor
    ))
  }
  
  coords <- cbind(xs, ys)
  
  coords_hull <- if (n >= 3L) {
    hull_index <- chull(coords)
    coords[hull_index, , drop = FALSE]
  } else {
    coords
  }
  
  distance_matrix <- geosphere::distm(coords_hull)
  
  if (
    !length(distance_matrix) ||
    !any(is.finite(distance_matrix))
  ) {
    return(list(
      major_km = base_major,
      minor_km = base_minor,
      elongation_ratio = base_major / base_minor
    ))
  }
  
  distance_matrix[!is.finite(distance_matrix)] <- -Inf
  ij <- which(
    distance_matrix == max(distance_matrix),
    arr.ind = TRUE
  )[1L, ]
  
  p1 <- coords_hull[ij[1L], ]
  p2 <- coords_hull[ij[2L], ]
  
  major_km <- max(distance_matrix) / 1000
  
  if (!is.finite(major_km) || major_km <= 0) {
    major_km <- base_major
  }
  
  theta <- atan2(
    p2[2L] - p1[2L],
    p2[1L] - p1[1L]
  ) + pi / 2
  
  perpendicular_unit <- c(
    cos(theta),
    sin(theta)
  )
  
  x_km <- (xs - p1[1L]) * km_per_deg_lon(lat0)
  y_km <- (ys - p1[2L]) * km_per_deg_lat(lat0)
  
  projection <- (
    x_km * perpendicular_unit[1L] +
      y_km * perpendicular_unit[2L]
  )
  
  minor_km <- diff(range(projection, na.rm = TRUE))
  
  if (!is.finite(minor_km) || minor_km <= 0) {
    minor_km <- base_minor
  }
  
  minor_km <- max(minor_km, 1e-6)
  
  list(
    major_km = as.numeric(major_km),
    minor_km = as.numeric(minor_km),
    elongation_ratio = as.numeric(major_km / minor_km)
  )
}

# =============================================================================
# COLUMN HELPERS
# =============================================================================

normalize_xy <- function(dt) {
  if (all(c("x", "y") %in% names(dt))) {
    return(invisible(NULL))
  }
  
  if (all(c("lon", "lat") %in% names(dt))) {
    data.table::setnames(
      dt,
      c("lon", "lat"),
      c("x", "y")
    )
    return(invisible(NULL))
  }
  
  stop(
    "Missing coordinate columns: expected x/y or lon/lat."
  )
}

guess_time_col <- function(dt) {
  candidates <- c(
    "time_utc",
    "time",
    "datetime",
    "timestamp",
    "day",
    "date"
  )
  
  found <- candidates[candidates %in% names(dt)]
  
  if (!length(found)) {
    stop(
      "Missing time column. Expected one of: ",
      paste(candidates, collapse = ", ")
    )
  }
  
  found[1L]
}

guess_intensity_col <- function(dt, cluster_set) {
  if (
    grepl(
      "^heatindex",
      cluster_set,
      ignore.case = TRUE
    )
  ) {
    candidates <- c(
      "intensity",
      "hi3_intensity",
      "hi_intensity",
      "heatindex",
      "HeatIndex",
      "value"
    )
  } else {
    candidates <- c(
      "intensity",
      "pr24_intensity_mm",
      "intensity_mm",
      "mm",
      "precip_mm",
      "prcp_mm",
      "precip",
      "pr",
      "ppt",
      "value"
    )
  }
  
  found <- candidates[candidates %in% names(dt)]
  
  if (length(found)) {
    return(found[1L])
  }
  
  ignored <- c(
    "x",
    "y",
    "lon",
    "lat",
    "time",
    "time_utc",
    "day",
    "date",
    "cluster",
    "cluster_id",
    "eps",
    "mu",
    "min_samples",
    "rp",
    "rp_max_years",
    "rp_max_layer",
    "cluster_max_rp_years",
    "resolution_deg",
    "time_step_hours",
    "t_ix",
    "row_ix",
    "col_ix",
    "cluster_floor_years",
    "event_filter_rp_years"
  )
  
  numeric_columns <- names(dt)[
    vapply(dt, is.numeric, logical(1L))
  ]
  
  numeric_columns <- setdiff(
    numeric_columns,
    ignored
  )
  
  if (length(numeric_columns) == 1L) {
    return(numeric_columns[1L])
  }
  
  stop(
    "Could not determine intensity column for ",
    cluster_set,
    "."
  )
}

guess_rp_col <- function(dt) {
  candidates <- c(
    "cluster_max_rp_years",
    "rp_max_years",
    "rp",
    "max_rp_exceeded",
    "rp_years",
    "rp_max"
  )
  
  found <- candidates[candidates %in% names(dt)]
  
  if (!length(found)) {
    return(NA_character_)
  }
  
  found[1L]
}

first_finite <- function(x, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  
  if (!length(x)) {
    return(default)
  }
  
  x[1L]
}

first_character <- function(x, default = NA_character_) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  
  if (!length(x)) {
    return(default)
  }
  
  x[1L]
}

# =============================================================================
# NOISE HELPERS
# =============================================================================

is_noise_point <- function(
    cluster_id,
    cluster = NULL
) {
  cluster_id_character <- trimws(
    as.character(cluster_id)
  )
  
  cluster_id_lower <- tolower(
    cluster_id_character
  )
  
  out <- (
    is.na(cluster_id) |
      cluster_id_character == "" |
      cluster_id_lower %in% c(
        "na",
        "nan",
        "noise",
        "-1",
        "0"
      )
  )
  
  if (!is.null(cluster)) {
    cluster_numeric <- suppressWarnings(
      as.numeric(cluster)
    )
    
    out <- (
      out |
        is.na(cluster_numeric) |
        cluster_numeric %in% c(-1, 0)
    )
  }
  
  out
}

# =============================================================================
# SAFE OUTPUT HELPERS
# =============================================================================

write_parquet_atomic <- function(
    x,
    path
) {
  fs::dir_create(
    fs::path_dir(path),
    recurse = TRUE
  )
  
  temporary_path <- paste0(
    path,
    ".tmp"
  )
  
  if (fs::file_exists(temporary_path)) {
    fs::file_delete(temporary_path)
  }
  
  arrow::write_parquet(
    as.data.frame(x),
    sink = temporary_path,
    compression = "zstd"
  )
  
  if (fs::file_exists(path)) {
    fs::file_delete(path)
  }
  
  fs::file_move(
    temporary_path,
    path
  )
}

write_csv_atomic <- function(
    x,
    path
) {
  fs::dir_create(
    fs::path_dir(path),
    recurse = TRUE
  )
  
  temporary_path <- paste0(
    path,
    ".tmp"
  )
  
  if (fs::file_exists(temporary_path)) {
    fs::file_delete(temporary_path)
  }
  
  data.table::fwrite(
    x,
    temporary_path
  )
  
  if (fs::file_exists(path)) {
    fs::file_delete(path)
  }
  
  fs::file_move(
    temporary_path,
    path
  )
}

# =============================================================================
# PER-FILE PROCESSOR
# =============================================================================

process_one_file <- function(
    path,
    cluster_set,
    filter_rp_years
) {
  source_file <- fs::path_file(path)
  
  dt0 <- tryCatch(
    data.table::as.data.table(
      arrow::read_parquet(
        path,
        as_data_frame = TRUE
      )
    ),
    error = function(e) {
      warning(
        "Failed to read ",
        path,
        ": ",
        conditionMessage(e)
      )
      data.table()
    }
  )
  
  if (!nrow(dt0)) {
    return(NULL)
  }
  
  normalize_xy(dt0)
  
  if (!("cluster_id" %in% names(dt0))) {
    warning(
      "Skipping file without cluster_id: ",
      source_file
    )
    return(NULL)
  }
  
  if (!("cluster" %in% names(dt0))) {
    dt0[, cluster := NA_integer_]
  }
  
  n_points_total <- nrow(dt0)
  
  noise_index <- is_noise_point(
    dt0$cluster_id,
    dt0$cluster
  )
  
  n_noise_total <- sum(
    noise_index,
    na.rm = TRUE
  )
  
  n_cluster_points <- (
    n_points_total -
      n_noise_total
  )
  
  dt <- data.table::copy(
    dt0[!noise_index]
  )
  
  if (!nrow(dt)) {
    return(NULL)
  }
  
  time_column <- guess_time_col(dt)
  intensity_column <- guess_intensity_col(
    dt,
    cluster_set
  )
  
  if (time_column != "time_raw") {
    data.table::setnames(
      dt,
      time_column,
      "time_raw"
    )
  }
  
  if (intensity_column != "intensity") {
    data.table::setnames(
      dt,
      intensity_column,
      "intensity"
    )
  }
  
  dt[, time := {
    if (inherits(time_raw, "Date")) {
      as_posix_utc(time_raw)
    } else if (inherits(time_raw, "POSIXct")) {
      as_posix_utc(time_raw)
    } else {
      parsed <- suppressWarnings(
        lubridate::ymd_hms(
          time_raw,
          tz = "UTC",
          quiet = TRUE
        )
      )
      
      if (all(is.na(parsed))) {
        parsed <- suppressWarnings(
          lubridate::ymd(
            time_raw,
            tz = "UTC",
            quiet = TRUE
          )
        )
      }
      
      if (all(is.na(parsed))) {
        parsed <- as_posix_utc(time_raw)
      }
      
      parsed
    }
  }]
  
  dt[, intensity := suppressWarnings(
    as.numeric(intensity)
  )]
  
  dt <- dt[
    is.finite(x) &
      is.finite(y) &
      is.finite(intensity) &
      !is.na(time)
  ]
  
  if (!nrow(dt)) {
    return(NULL)
  }
  
  resolution_from_data <- if (
    "resolution_deg" %in% names(dt)
  ) {
    first_finite(dt$resolution_deg)
  } else {
    NA_real_
  }
  
  resx_deg <- resolution_from_data
  resy_deg <- resolution_from_data
  
  if (!is.finite(resx_deg) || resx_deg <= 0) {
    resx_deg <- grid_step(dt$x)
  }
  
  if (!is.finite(resy_deg) || resy_deg <= 0) {
    resy_deg <- grid_step(dt$y)
  }
  
  dt[, area_km2 := {
    if (
      is.finite(resx_deg) &&
      resx_deg > 0 &&
      is.finite(resy_deg) &&
      resy_deg > 0
    ) {
      cell_area_km2_vec(
        y,
        resx_deg,
        resy_deg
      )
    } else {
      NA_real_
    }
  }]
  
  timestep_hr <- if (
    "time_step_hours" %in% names(dt)
  ) {
    first_finite(
      dt$time_step_hours,
      default = DEFAULT_TIMESTEP_HR
    )
  } else {
    DEFAULT_TIMESTEP_HR
  }
  
  if (
    !is.finite(timestep_hr) ||
    timestep_hr <= 0
  ) {
    timestep_hr <- DEFAULT_TIMESTEP_HR
  }
  
  if (!("t_ix" %in% names(dt))) {
    dt[, cluster_start_time := min(
      time,
      na.rm = TRUE
    ), by = cluster_id]
    
    dt[, t_ix := as.integer(
      round(
        as.numeric(
          difftime(
            time,
            cluster_start_time,
            units = "hours"
          )
        ) / timestep_hr
      )
    ) + 1L]
    
    dt[, cluster_start_time := NULL]
  } else {
    dt[, t_ix := suppressWarnings(
      as.integer(t_ix)
    )]
  }
  
  dt <- dt[
    !is.na(t_ix) &
      t_ix >= 1L
  ]
  
  if (!nrow(dt)) {
    return(NULL)
  }
  
  time_bounds <- dt[, .(
    start_time = min(
      time,
      na.rm = TRUE
    ),
    end_time = max(
      time,
      na.rm = TRUE
    )
  ), by = cluster_id]
  
  duration_table <- dt[, .(
    t_min = min(
      t_ix,
      na.rm = TRUE
    ),
    t_max = max(
      t_ix,
      na.rm = TRUE
    )
  ), by = cluster_id]
  
  duration_table[, `:=`(
    n_steps = as.integer(
      t_max - t_min + 1L
    ),
    duration_hr = as.numeric(
      t_max - t_min + 1L
    ) * timestep_hr
  )]
  
  duration_table[, c(
    "t_min",
    "t_max"
  ) := NULL]
  
  intensity_stats <- dt[, .(
    peak_intensity = max(
      intensity,
      na.rm = TRUE
    ),
    mean_intensity = mean(
      intensity,
      na.rm = TRUE
    ),
    median_intensity = as.numeric(
      stats::median(
        intensity,
        na.rm = TRUE
      )
    )
  ), by = cluster_id]
  
  extent_by_timestep <- dt[, .(
    timestep_extent_km2 = sum(
      area_km2,
      na.rm = TRUE
    )
  ), by = .(
    cluster_id,
    t_ix
  )]
  
  maximum_timestep_extent <- extent_by_timestep[, .(
    max_t_extent_km2 = max(
      timestep_extent_km2,
      na.rm = TRUE
    )
  ), by = cluster_id]
  
  union_extent_cells <- unique(
    dt[, .(
      cluster_id,
      x,
      y,
      area_km2
    )]
  )
  
  total_extent <- union_extent_cells[, .(
    total_extent_km2 = sum(
      area_km2,
      na.rm = TRUE
    )
  ), by = cluster_id]
  
  is_precipitation <- grepl(
    "^precipitation",
    cluster_set,
    ignore.case = TRUE
  )
  
  if (is_precipitation) {
    volume_table <- dt[, .(
      total_volume_km3 = sum(
        (
          intensity / 1000
        ) * (
          area_km2 * 1e6
        ),
        na.rm = TRUE
      ) / 1e9
    ), by = cluster_id]
  } else {
    volume_table <- unique(
      dt[, .(
        cluster_id,
        total_volume_km3 = NA_real_
      )]
    )
  }
  
  centroids <- dt[, {
    centroid <- centroid_timestep(
      x,
      y,
      intensity
    )
    
    .(
      lon = centroid[["lon"]],
      lat = centroid[["lat"]]
    )
  }, by = .(
    cluster_id,
    t_ix
  )]
  
  data.table::setorder(
    centroids,
    cluster_id,
    t_ix
  )
  
  speeds <- centroids[, {
    if (.N < 2L) {
      .(
        max_speed_kmhr = 0,
        mean_speed_kmhr = 0,
        median_speed_kmhr = 0
      )
    } else {
      lon1 <- head(lon, -1L)
      lat1 <- head(lat, -1L)
      lon2 <- tail(lon, -1L)
      lat2 <- tail(lat, -1L)
      
      timestep_differences <- pmax(
        1,
        diff(t_ix)
      )
      
      distances_km <- mapply(
        FUN = function(
    lon_start,
    lat_start,
    lon_end,
    lat_end
        ) {
          gc_km(
            c(lon_start, lat_start),
            c(lon_end, lat_end)
          )
        },
    lon1,
    lat1,
    lon2,
    lat2
      )
      
      elapsed_hours <- (
        timestep_differences *
          timestep_hr
      )
      
      speed_values <- (
        distances_km /
          elapsed_hours
      )
      
      speed_values <- speed_values[
        is.finite(speed_values)
      ]
      
      if (!length(speed_values)) {
        speed_values <- 0
      }
      
      .(
        max_speed_kmhr = max(
          speed_values,
          na.rm = TRUE
        ),
        mean_speed_kmhr = mean(
          speed_values,
          na.rm = TRUE
        ),
        median_speed_kmhr = as.numeric(
          stats::median(
            speed_values,
            na.rm = TRUE
          )
        )
      )
    }
  }, by = cluster_id]
  
  elongation_table <- union_extent_cells[, {
    elongation <- elongation_major_minor(
      x,
      y,
      resx_deg,
      resy_deg
    )
    
    .(
      elongation_ratio = elongation$elongation_ratio,
      major_axis_km = elongation$major_km,
      minor_axis_km = elongation$minor_km
    )
  }, by = cluster_id]
  
  rp_column <- guess_rp_col(dt)
  
  if (is.na(rp_column)) {
    rp_table <- unique(
      dt[, .(
        cluster_id,
        cluster_max_rp_years = NA_real_
      )]
    )
  } else {
    rp_table <- dt[
      is.finite(
        suppressWarnings(
          as.numeric(
            get(rp_column)
          )
        )
      ),
      .(
        cluster_max_rp_years = max(
          suppressWarnings(
            as.numeric(
              get(rp_column)
            )
          ),
          na.rm = TRUE
        )
      ),
      by = cluster_id
    ]
    
    missing_cluster_ids <- setdiff(
      unique(dt$cluster_id),
      rp_table$cluster_id
    )
    
    if (length(missing_cluster_ids)) {
      rp_table <- data.table::rbindlist(
        list(
          rp_table,
          data.table(
            cluster_id = missing_cluster_ids,
            cluster_max_rp_years = NA_real_
          )
        ),
        use.names = TRUE,
        fill = TRUE
      )
    }
  }
  
  cluster_map <- unique(
    dt[, .(
      cluster_id,
      cluster
    )]
  )
  
  source_table <- data.table(
    cluster_id = unique(dt$cluster_id),
    source_file = source_file
  )
  
  file_counts <- data.table(
    cluster_id = unique(dt$cluster_id),
    n_points_file = n_points_total,
    n_noise_file = n_noise_total,
    n_clusterpts_file = n_cluster_points
  )
  
  metadata_columns <- c(
    "variable",
    "dataset_tag",
    "resample_method",
    "cluster_floor_years",
    "event_filter_rp_years",
    "eps",
    "min_samples",
    "time_step_hours"
  )
  
  metadata_columns <- metadata_columns[
    metadata_columns %in% names(dt)
  ]
  
  metadata_values <- lapply(
    metadata_columns,
    function(column_name) {
      values <- dt[[column_name]]
      
      if (is.numeric(values)) {
        first_finite(values)
      } else {
        first_character(values)
      }
    }
  )
  
  names(metadata_values) <- metadata_columns
  
  metadata_values$filter_rp_years <- filter_rp_years
  
  if (
    is.null(
      metadata_values$event_filter_rp_years
    ) ||
    !is.finite(
      suppressWarnings(
        as.numeric(
          metadata_values$event_filter_rp_years
        )
      )
    )
  ) {
    metadata_values$event_filter_rp_years <- filter_rp_years
  }
  
  metadata_table <- data.table(
    cluster_id = unique(dt$cluster_id)
  )
  
  for (metadata_name in names(metadata_values)) {
    metadata_table[, (metadata_name) := metadata_values[[metadata_name]]]
  }
  
  output <- Reduce(
    function(left, right) {
      merge(
        left,
        right,
        by = "cluster_id",
        all = TRUE
      )
    },
    list(
      data.table(
        cluster_id = unique(dt$cluster_id),
        cluster_set = cluster_set
      ),
      metadata_table,
      cluster_map,
      source_table,
      file_counts,
      rp_table,
      time_bounds,
      duration_table,
      intensity_stats,
      maximum_timestep_extent,
      total_extent,
      volume_table,
      speeds,
      elongation_table
    )
  )
  
  output[, resolution_deg := {
    if (
      is.finite(resx_deg) &&
      is.finite(resy_deg) &&
      abs(resx_deg - resy_deg) < 1e-10
    ) {
      resx_deg
    } else if (is.finite(resx_deg)) {
      resx_deg
    } else {
      NA_real_
    }
  }]
  
  output[, timestep_hr := timestep_hr]
  output[, duration_days := duration_hr / 24]
  
  preferred_column_order <- c(
    "cluster_set",
    "variable",
    "dataset_tag",
    "filter_rp_years",
    "event_filter_rp_years",
    "cluster_floor_years",
    "resample_method",
    "resolution_deg",
    "cluster_id",
    "cluster",
    "cluster_max_rp_years",
    "source_file",
    "n_points_file",
    "n_noise_file",
    "n_clusterpts_file",
    "start_time",
    "end_time",
    "n_steps",
    "timestep_hr",
    "duration_hr",
    "duration_days",
    "peak_intensity",
    "mean_intensity",
    "median_intensity",
    "max_t_extent_km2",
    "total_extent_km2",
    "total_volume_km3",
    "max_speed_kmhr",
    "mean_speed_kmhr",
    "median_speed_kmhr",
    "elongation_ratio",
    "major_axis_km",
    "minor_axis_km",
    "eps",
    "min_samples"
  )
  
  preferred_column_order <- preferred_column_order[
    preferred_column_order %in% names(output)
  ]
  
  data.table::setcolorder(
    output,
    c(
      preferred_column_order,
      setdiff(
        names(output),
        preferred_column_order
      )
    )
  )
  
  output[]
}

# =============================================================================
# PER-FILTER PROCESSOR
# =============================================================================

process_filter <- function(
    filter_rp_years
) {
  filter_name <- paste0(
    "filtered_rp",
    filter_rp_years
  )
  
  filter_input_directory <- fs::path(
    BASE_IN,
    filter_name
  )
  
  filter_output_directory <- fs::path(
    BASE_OUT,
    filter_name
  )
  
  fs::dir_create(
    filter_output_directory,
    recurse = TRUE
  )
  
  if (!fs::dir_exists(filter_input_directory)) {
    warning(
      "Missing input filter directory: ",
      filter_input_directory
    )
    return(NULL)
  }
  
  .msg(
    "Starting %s",
    filter_name
  )
  
  filter_summary_list <- list()
  
  for (cluster_set in CLUSTER_SETS) {
    input_directory <- fs::path(
      filter_input_directory,
      cluster_set
    )
    
    output_directory <- fs::path(
      filter_output_directory,
      cluster_set
    )
    
    if (!fs::dir_exists(input_directory)) {
      warning(
        "Missing cluster-set directory; skipping: ",
        input_directory
      )
      next
    }
    
    parquet_files <- fs::dir_ls(
      input_directory,
      recurse = TRUE,
      type = "file",
      glob = "*.parquet"
    )
    
    parquet_files <- parquet_files[
      !grepl(
        "characteristics",
        fs::path_file(parquet_files),
        ignore.case = TRUE
      )
    ]
    
    parquet_files <- sort(
      as.character(parquet_files)
    )
    
    if (!length(parquet_files)) {
      warning(
        "No input Parquet files under: ",
        input_directory
      )
      next
    }
    
    .msg(
      "%s | %s | %s input files",
      filter_name,
      cluster_set,
      format(
        length(parquet_files),
        big.mark = ","
      )
    )
    
    result <- progressr::with_progress({
      progress <- progressr::progressor(
        steps = length(parquet_files)
      )
      
      parts <- future.apply::future_lapply(
        parquet_files,
        function(file_path) {
          progress(
            fs::path_file(file_path)
          )
          
          process_one_file(
            path = file_path,
            cluster_set = cluster_set,
            filter_rp_years = filter_rp_years
          )
        },
        future.seed = TRUE
      )
      
      data.table::rbindlist(
        Filter(
          Negate(is.null),
          parts
        ),
        use.names = TRUE,
        fill = TRUE
      )
    })
    
    if (!nrow(result)) {
      warning(
        "No characteristic rows for ",
        filter_name,
        "/",
        cluster_set
      )
      next
    }
    
    data.table::setDT(result)
    
    result <- result[
      !is.na(cluster_id) &
        nzchar(
          as.character(cluster_id)
        )
    ]
    
    duplicate_keys <- result[
      ,
      .N,
      by = .(
        cluster_set,
        cluster_id
      )
    ][N > 1L]
    
    if (nrow(duplicate_keys)) {
      warning(
        filter_name,
        "/",
        cluster_set,
        " contains ",
        format(
          nrow(duplicate_keys),
          big.mark = ","
        ),
        " duplicated cluster_id values across source files. ",
        "Keeping the first occurrence after source-file sorting."
      )
      
      data.table::setorder(
        result,
        source_file
      )
      
      result <- unique(
        result,
        by = c(
          "cluster_set",
          "cluster_id"
        )
      )
    }
    
    fs::dir_create(
      output_directory,
      recurse = TRUE
    )
    
    output_parquet <- fs::path(
      output_directory,
      "cluster_characteristics.parquet"
    )
    
    output_csv <- fs::path(
      output_directory,
      "cluster_characteristics.csv"
    )
    
    write_parquet_atomic(
      result,
      output_parquet
    )
    
    write_csv_atomic(
      result,
      output_csv
    )
    
    .msg(
      "Wrote %s rows and %s unique clusters to %s",
      format(
        nrow(result),
        big.mark = ","
      ),
      format(
        data.table::uniqueN(
          result$cluster_id
        ),
        big.mark = ","
      ),
      output_directory
    )
    
    file_counts <- unique(
      result[, .(
        cluster_set,
        filter_rp_years,
        resolution_deg,
        source_file,
        n_points_file,
        n_noise_file,
        n_clusterpts_file
      )]
    )
    
    count_summary <- file_counts[, .(
      n_input_files = data.table::uniqueN(
        source_file
      ),
      n_points = sum(
        n_points_file,
        na.rm = TRUE
      ),
      n_noise = sum(
        n_noise_file,
        na.rm = TRUE
      ),
      n_clusterpoints = sum(
        n_clusterpts_file,
        na.rm = TRUE
      )
    ), by = .(
      cluster_set,
      filter_rp_years,
      resolution_deg
    )]
    
    characteristic_summary <- result[, .(
      n_clusters = data.table::uniqueN(
        cluster_id
      ),
      median_duration_days = as.numeric(
        stats::median(
          duration_days,
          na.rm = TRUE
        )
      ),
      median_max_t_extent_km2 = as.numeric(
        stats::median(
          max_t_extent_km2,
          na.rm = TRUE
        )
      ),
      median_total_extent_km2 = as.numeric(
        stats::median(
          total_extent_km2,
          na.rm = TRUE
        )
      )
    ), by = .(
      cluster_set,
      filter_rp_years,
      resolution_deg
    )]
    
    filter_summary_list[[cluster_set]] <- merge(
      count_summary,
      characteristic_summary,
      by = c(
        "cluster_set",
        "filter_rp_years",
        "resolution_deg"
      ),
      all = TRUE
    )
  }
  
  if (!length(filter_summary_list)) {
    warning(
      "No successful cluster sets for ",
      filter_name
    )
    return(NULL)
  }
  
  filter_summary <- data.table::rbindlist(
    filter_summary_list,
    use.names = TRUE,
    fill = TRUE
  )
  
  data.table::setorder(
    filter_summary,
    cluster_set
  )
  
  summary_parquet <- fs::path(
    filter_output_directory,
    "cluster_summary_by_resolution.parquet"
  )
  
  summary_csv <- fs::path(
    filter_output_directory,
    "cluster_summary_by_resolution.csv"
  )
  
  write_parquet_atomic(
    filter_summary,
    summary_parquet
  )
  
  write_csv_atomic(
    filter_summary,
    summary_csv
  )
  
  .msg(
    "Completed %s",
    filter_name
  )
  
  filter_summary[]
}

# =============================================================================
# RUN
# =============================================================================

if (
  !is.finite(N_WORKERS) ||
  N_WORKERS < 1L
) {
  stop(
    "CISC_CHARACTERISTICS_WORKERS must be a positive integer."
  )
}

fs::dir_create(
  BASE_OUT,
  recurse = TRUE
)

future::plan(
  future::multisession,
  workers = N_WORKERS
)

progressr::handlers(
  global = TRUE
)

all_summaries <- lapply(
  FILTER_RP_YEARS,
  process_filter
)

all_summaries <- Filter(
  Negate(is.null),
  all_summaries
)

if (length(all_summaries)) {
  combined_summary <- data.table::rbindlist(
    all_summaries,
    use.names = TRUE,
    fill = TRUE
  )
  
  data.table::setorder(
    combined_summary,
    filter_rp_years,
    cluster_set
  )
  
  combined_summary_parquet <- fs::path(
    BASE_OUT,
    "cluster_summary_by_resolution.parquet"
  )
  
  combined_summary_csv <- fs::path(
    BASE_OUT,
    "cluster_summary_by_resolution.csv"
  )
  
  write_parquet_atomic(
    combined_summary,
    combined_summary_parquet
  )
  
  write_csv_atomic(
    combined_summary,
    combined_summary_csv
  )
  
  .msg(
    "Wrote combined summaries to %s",
    BASE_OUT
  )
  
  print(combined_summary)
} else {
  warning(
    "No characteristic outputs were created."
  )
}

future::plan(
  future::sequential
)