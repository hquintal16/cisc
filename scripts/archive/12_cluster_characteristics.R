# # =============================================================================
# # Per-cluster characteristics for YEAR-CHUNKED clustered exceedance .qs files
# #
# # INPUTS (new format):
# #   V:/users/hquintal/phd1_cluster_southeast_v2/data/out/06_clusters/hi{H}_pr{P}/{set}/*.qs
# #
# # OUTPUTS:
# #   data/out/07_characteristics/hi{H}_pr{P}/cluster_characteristics_{set}.qs (+ .csv)
# #
# # IMPORTANT FIX:
# #   - timestep_hr is NEVER inferred from exceedance spacing
# #   - timestep_hr is taken from dt$time_step_hours (written during clustering) else defaults to 24
# #   - n_steps & duration_hr are computed from integer timestep indices t_ix
# # =============================================================================
# 
# suppressPackageStartupMessages({
#   library(data.table)
#   library(qs)
#   library(fs)
#   library(stringr)
#   library(lubridate)
#   library(future)
#   library(future.apply)
#   library(progressr)
#   library(geosphere)
# })
# 
# # ----------------------------
# # USER SETTINGS
# # ----------------------------
# BASE <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out"
# 
# HI_FLOOR_YR <- 25
# PR_FLOOR_YR <- 25
# 
# BASE_IN <- fs::path(BASE, sprintf("06_clusters/hi%d_pr%d", HI_FLOOR_YR, PR_FLOOR_YR))
# 
# OUT_DIR <- fs::path(BASE, sprintf("07_characteristics/hi%d_pr%d", HI_FLOOR_YR, PR_FLOOR_YR))
# dir_create(OUT_DIR, recurse = TRUE)
# 
# N_WORKERS <- 30
# 
# CLUSTER_SETS <- c(
#   "heatindex_0.25",
#   "heatindex_0.30_bilinear",
#   "heatindex_0.30_max",
#   "heatindex_0.36_bilinear",
#   "heatindex_0.36_max",
#   "precipitation_0.25",
#   "precipitation_0.11_bilinear",
#   "precipitation_0.11_max",
#   "precipitation_0.14_bilinear",
#   "precipitation_0.14_max"
# )
# 
# # ----------------------------
# # GEO HELPERS
# # ----------------------------
# km_per_deg_lat <- function(lat) 111.132
# km_per_deg_lon <- function(lat) 111.320 * cospi(lat/180)
# 
# cell_area_km2_vec <- function(lat_vec, resx_deg, resy_deg) {
#   km_y <- km_per_deg_lat(lat_vec) * resy_deg
#   km_x <- km_per_deg_lon(lat_vec) * resx_deg
#   km_x * km_y
# }
# 
# grid_step <- function(v) {
#   u <- sort(unique(v))
#   if (length(u) < 2) return(NA_real_)
#   suppressWarnings(min(diff(u)))
# }
# 
# as_posix_utc <- function(tt) {
#   if (inherits(tt, "POSIXct")) return(as.POSIXct(tt, tz = "UTC"))
#   if (inherits(tt, "Date"))   return(as.POSIXct(tt, tz = "UTC"))
#   suppressWarnings(as.POSIXct(tt, tz = "UTC"))
# }
# 
# gc_km <- function(p1, p2) geosphere::distHaversine(p1, p2) / 1000
# 
# centroid_timestep <- function(x, y, w) {
#   pos <- is.finite(w) & w > 0
#   if (!any(pos)) return(c(lon = NA_real_, lat = NA_real_))
#   c(lon = mean(x[pos]), lat = mean(y[pos]))
# }
# 
# elongation_major_minor <- function(xs, ys, resx_deg, resy_deg) {
#   n <- length(xs)
#   lat0 <- mean(ys)
#   
#   base_major <- max(km_per_deg_lon(lat0) * resx_deg, km_per_deg_lat(lat0) * resy_deg, na.rm = TRUE)
#   base_minor <- max(min(km_per_deg_lon(lat0) * resx_deg, km_per_deg_lat(lat0) * resy_deg, na.rm = TRUE), 1e-6)
#   
#   if (!is.finite(base_major)) base_major <- 0
#   if (!is.finite(base_minor) || base_minor <= 0) base_minor <- 1e-6
#   
#   if (n < 2) {
#     return(list(major_km = base_major, minor_km = base_minor, elongation_ratio = base_major / base_minor))
#   }
#   
#   coords <- cbind(xs, ys)
#   coords_h <- if (n >= 6) {
#     idx <- chull(coords)
#     coords[c(idx, idx[1]), , drop = FALSE]
#   } else coords
#   
#   D <- geosphere::distm(coords_h)
#   if (!length(D) || !is.finite(max(D))) {
#     return(list(major_km = base_major, minor_km = base_minor, elongation_ratio = base_major / base_minor))
#   }
#   
#   ij <- which(D == max(D), arr.ind = TRUE)[1, ]
#   p1 <- coords_h[ij[1], ]
#   p2 <- coords_h[ij[2], ]
#   major_km <- max(D) / 1000
#   
#   theta <- atan2(p2[2] - p1[2], p2[1] - p1[1]) + pi/2
#   u_perp <- c(cos(theta), sin(theta))
#   
#   xkm <- (xs - p1[1]) * km_per_deg_lon(lat0)
#   ykm <- (ys - p1[2]) * km_per_deg_lat(lat0)
#   proj <- xkm * u_perp[1] + ykm * u_perp[2]
#   
#   minor_km <- diff(range(proj))
#   if (!is.finite(minor_km) || minor_km <= 0) minor_km <- base_minor
#   minor_km <- max(minor_km, 1e-6)
#   
#   list(
#     major_km = as.numeric(major_km),
#     minor_km = as.numeric(minor_km),
#     elongation_ratio = as.numeric(major_km / minor_km)
#   )
# }
# 
# # ----------------------------
# # COLUMN GUESSERS
# # ----------------------------
# normalize_xy <- function(dt) {
#   if (all(c("x","y") %in% names(dt))) return(invisible(NULL))
#   if (all(c("lon","lat") %in% names(dt))) {
#     setnames(dt, c("lon","lat"), c("x","y"))
#     return(invisible(NULL))
#   }
#   stop("Missing coordinate columns: need x,y (or lon,lat).")
# }
# 
# guess_time_col <- function(dt) {
#   for (cc in c("time_utc","time","datetime","timestamp","day","date")) {
#     if (cc %in% names(dt)) return(cc)
#   }
#   stop("Missing time column (expected one of time_utc, time, day/date).")
# }
# 
# guess_intensity_col <- function(dt, cluster_set) {
#   if (grepl("^heatindex", cluster_set, ignore.case = TRUE)) {
#     for (cc in c("hi3_intensity","hi_intensity","heatindex","HeatIndex","value","intensity")) {
#       if (cc %in% names(dt)) return(cc)
#     }
#   }
#   
#   # precipitation (explicit priority)
#   for (cc in c("pr24_intensity_mm","intensity_mm","mm","precip_mm","prcp_mm","precip","pr","ppt","value","intensity")) {
#     if (cc %in% names(dt)) return(cc)
#   }
#   
#   ignore <- c("x","y","lon","lat","time","time_utc","day","date",
#               "cluster","cluster_id","eps","mu",
#               "rp","rp_max_years","rp_max_layer","resolution_deg",
#               "time_step_hours","t_ix","row_ix","col_ix","floor_yr","dataset_tag")
#   num_cols <- names(dt)[vapply(dt, is.numeric, TRUE)]
#   num_cols <- setdiff(num_cols, ignore)
#   if (length(num_cols) == 1) return(num_cols[1])
#   
#   stop("Could not determine intensity column for set=", cluster_set,
#        ". Rename intensity column or extend guess_intensity_col().")
# }
# 
# guess_rp_col <- function(dt) {
#   for (cc in c("rp_max_years","rp","max_rp_exceeded","rp_years","rp_max")) {
#     if (cc %in% names(dt)) return(cc)
#   }
#   return(NA_character_)
# }
# 
# # ----------------------------
# # PER-FILE PROCESSOR
# # ----------------------------
# process_one_file <- function(path, cluster_set) {
#   base <- path_file(path)
#   dt <- tryCatch(as.data.table(qread(path)), error = function(e) data.table())
#   if (!nrow(dt)) return(NULL)
#   
#   normalize_xy(dt)
#   
#   if (!("cluster_id" %in% names(dt))) {
#     warning("Skipping (no cluster_id): ", base)
#     return(NULL)
#   }
#   if (!("cluster" %in% names(dt))) dt[, cluster := NA_integer_]
#   
#   # keep only assigned clusters
#   dt <- dt[!is.na(cluster_id) & nzchar(as.character(cluster_id))]
#   if (!nrow(dt)) return(NULL)
#   
#   time_col  <- guess_time_col(dt)
#   inten_col <- guess_intensity_col(dt, cluster_set)
#   
#   setnames(dt, time_col,  "time_raw")
#   setnames(dt, inten_col, "intensity")
#   
#   # parse time to POSIXct UTC
#   dt[, time := {
#     if (inherits(time_raw, "Date")) as_posix_utc(time_raw)
#     else if (inherits(time_raw, "POSIXct")) as_posix_utc(time_raw)
#     else {
#       tt <- suppressWarnings(ymd_hms(time_raw, tz = "UTC"))
#       if (all(is.na(tt))) tt <- suppressWarnings(ymd(time_raw, tz = "UTC"))
#       if (all(is.na(tt))) tt <- as_posix_utc(time_raw)
#       tt
#     }
#   }]
#   
#   dt <- dt[is.finite(x) & is.finite(y) & is.finite(intensity) & !is.na(time)]
#   if (!nrow(dt)) return(NULL)
#   
#   # resolution
#   resx_deg <- if ("resolution_deg" %in% names(dt)) unique(dt$resolution_deg)[1] else grid_step(dt$x)
#   resy_deg <- if ("resolution_deg" %in% names(dt)) unique(dt$resolution_deg)[1] else grid_step(dt$y)
#   if (!is.finite(resx_deg) || resx_deg <= 0) resx_deg <- grid_step(dt$x)
#   if (!is.finite(resy_deg) || resy_deg <= 0) resy_deg <- grid_step(dt$y)
#   
#   dt[, area_km2 := if (is.finite(resx_deg) && is.finite(resy_deg)) cell_area_km2_vec(y, resx_deg, resy_deg) else NA_real_]
#   
#   # ---- FIXED timestep: NEVER infer from exceedance spacing ----
#   timestep_hr <- if ("time_step_hours" %in% names(dt)) as.numeric(unique(dt$time_step_hours)[1]) else 24
#   if (!is.finite(timestep_hr) || timestep_hr <= 0) timestep_hr <- 24
#   
#   # t_ix: use existing if present, else compute from cluster start (keeps gaps)
#   if (!("t_ix" %in% names(dt))) {
#     dt[, start0 := min(time, na.rm = TRUE), by = cluster_id]
#     dt[, t_ix := as.integer(round(as.numeric(difftime(time, start0, units = "hours")) / timestep_hr)) + 1L]
#     dt[, start0 := NULL]
#   }
#   
#   # bounds
#   time_bounds <- dt[, .(
#     start_time = min(time, na.rm = TRUE),
#     end_time   = max(time, na.rm = TRUE)
#   ), by = cluster_id]
#   
#   # duration based on t_ix (inclusive)
#   duration_tbl <- dt[, .(
#     t_min = min(t_ix, na.rm = TRUE),
#     t_max = max(t_ix, na.rm = TRUE)
#   ), by = cluster_id][, .(
#     n_steps = as.integer(t_max - t_min + 1L),
#     duration_hr = as.numeric(t_max - t_min + 1L) * timestep_hr
#   ), by = cluster_id]
#   
#   # intensity stats
#   intensity_stats <- dt[, .(
#     peak_intensity   = max(intensity, na.rm = TRUE),
#     mean_intensity   = mean(intensity, na.rm = TRUE),
#     median_intensity = as.numeric(stats::median(intensity, na.rm = TRUE))
#   ), by = cluster_id]
#   
#   # extents by timestep (use t_ix, not unique times)
#   extent_by_t <- dt[, .(t_extent_km2 = sum(area_km2, na.rm = TRUE)),
#                     by = .(cluster_id, t_ix)]
#   max_t_extent <- extent_by_t[, .(max_t_extent_km2 = max(t_extent_km2, na.rm = TRUE)), by = cluster_id]
#   
#   union_extent <- unique(dt[, .(cluster_id, x, y, area_km2)])
#   total_extent <- union_extent[, .(total_extent_km2 = sum(area_km2, na.rm = TRUE)), by = cluster_id]
#   
#   # volume only for precipitation sets
#   is_precip <- grepl("^precipitation", cluster_set, ignore.case = TRUE)
#   if (is_precip) {
#     vol_tbl <- dt[, .(
#       total_volume_km3 = sum((intensity / 1000) * (area_km2 * 1e6), na.rm = TRUE) / 1e9
#     ), by = cluster_id]
#   } else {
#     vol_tbl <- unique(dt[, .(cluster_id, total_volume_km3 = as.numeric(NA))])
#   }
#   
#   # centroid speeds (per t_ix), using fixed timestep_hr
#   centroids <- dt[, {
#     cxy <- centroid_timestep(x, y, intensity)
#     .(lon = cxy[["lon"]], lat = cxy[["lat"]])
#   }, by = .(cluster_id, t_ix)]
#   setorder(centroids, cluster_id, t_ix)
#   
#   speeds <- centroids[, {
#     if (.N < 2) {
#       data.table(max_speed_kmhr = 0, mean_speed_kmhr = 0, median_speed_kmhr = 0)
#     } else {
#       lon1 <- head(lon, -1); lat1 <- head(lat, -1)
#       lon2 <- tail(lon, -1); lat2 <- tail(lat, -1)
#       dstep <- pmax(1, diff(t_ix))
#       dkm <- mapply(function(a,b,c,d) gc_km(c(a,b), c(c,d)), lon1, lat1, lon2, lat2)
#       dh <- dstep * timestep_hr
#       spd <- dkm / dh
#       spd <- spd[is.finite(spd)]
#       if (!length(spd)) spd <- 0
#       data.table(
#         max_speed_kmhr    = max(spd, na.rm = TRUE),
#         mean_speed_kmhr   = mean(spd, na.rm = TRUE),
#         median_speed_kmhr = as.numeric(stats::median(spd, na.rm = TRUE))
#       )
#     }
#   }, by = cluster_id]
#   
#   # elongation (spatial footprint union)
#   elong_tbl <- union_extent[, {
#     el <- elongation_major_minor(x, y, resx_deg, resy_deg)
#     data.table(
#       elongation_ratio = el$elongation_ratio,
#       major_axis_km    = el$major_km,
#       minor_axis_km    = el$minor_km
#     )
#   }, by = cluster_id]
#   
#   # RP per cluster (max)
#   rp_col <- guess_rp_col(dt)
#   if (is.na(rp_col)) {
#     rp_tbl <- unique(dt[, .(cluster_id, rp = as.numeric(NA))])
#   } else {
#     rp_tbl <- dt[is.finite(get(rp_col)), .(rp = max(get(rp_col), na.rm = TRUE)), by = cluster_id]
#   }
#   
#   cluster_map <- unique(dt[, .(cluster_id, cluster)])
#   file_col    <- data.table(cluster_id = unique(dt$cluster_id), source_file = base)
#   
#   # floors (if present)
#   hi_floor <- if ("hi_floor_yr" %in% names(dt)) unique(dt$hi_floor_yr)[1] else NA_real_
#   pr_floor <- if ("pr_floor_yr" %in% names(dt)) unique(dt$pr_floor_yr)[1] else NA_real_
#   # new clustering script writes floor_yr; map it to hi/pr using set name
#   if ("floor_yr" %in% names(dt)) {
#     if (grepl("^heatindex", cluster_set, ignore.case=TRUE)) hi_floor <- unique(dt$floor_yr)[1]
#     if (grepl("^precipitation", cluster_set, ignore.case=TRUE)) pr_floor <- unique(dt$floor_yr)[1]
#   }
#   floor_tbl <- unique(data.table(cluster_id = unique(dt$cluster_id),
#                                  hi_floor_yr = hi_floor,
#                                  pr_floor_yr = pr_floor))
#   
#   out <- Reduce(function(a,b) merge(a, b, by = "cluster_id", all = TRUE), list(
#     data.table(cluster_id = unique(dt$cluster_id), cluster_set = cluster_set),
#     floor_tbl, cluster_map, file_col, rp_tbl, time_bounds, duration_tbl,
#     intensity_stats, max_t_extent, total_extent, vol_tbl,
#     speeds, elong_tbl
#   ))
#   
#   out[, resolution_deg := if (is.finite(resx_deg)) resx_deg else NA_real_]
#   out[, timestep_hr := timestep_hr]
#   
#   setcolorder(out, c(
#     "cluster_set","hi_floor_yr","pr_floor_yr",
#     "cluster_id","cluster","rp","source_file",
#     "start_time","end_time","n_steps","timestep_hr","duration_hr",
#     "peak_intensity","mean_intensity","median_intensity",
#     "max_t_extent_km2","total_extent_km2","total_volume_km3",
#     "max_speed_kmhr","mean_speed_kmhr","median_speed_kmhr",
#     "elongation_ratio","major_axis_km","minor_axis_km",
#     "resolution_deg"
#   ))
#   
#   out[]
# }
# 
# # ----------------------------
# # RUN
# # ----------------------------
# plan(multisession, workers = N_WORKERS)
# handlers(global = TRUE)
# 
# for (cluster_set in CLUSTER_SETS) {
#   in_dir <- fs::path(BASE_IN, cluster_set)
#   if (!dir_exists(in_dir)) {
#     warning("Missing cluster set folder (skipping): ", in_dir)
#     next
#   }
#   
#   qs_files <- dir_ls(in_dir, recurse = TRUE, type = "file", glob = "*.qs")
#   qs_files <- qs_files[!grepl("characteristics", qs_files, ignore.case = TRUE)]
#   qs_files <- sort(as.character(qs_files))
#   
#   if (!length(qs_files)) {
#     warning("No .qs found under: ", in_dir)
#     next
#   }
#   
#   message("============================================================")
#   message("Cluster set: ", cluster_set)
#   message("Inputs: ", length(qs_files))
#   message("============================================================")
#   
#   res <- with_progress({
#     p <- progressor(steps = length(qs_files))
#     parts <- future_lapply(qs_files, function(f) {
#       p(path_file(f))
#       process_one_file(f, cluster_set = cluster_set)
#     })
#     rbindlist(Filter(Negate(is.null), parts), use.names = TRUE, fill = TRUE)
#   })
#   
#   if (!nrow(res)) {
#     warning("No output rows for cluster_set=", cluster_set)
#     next
#   }
#   
#   setDT(res)
#   res <- res[!is.na(cluster_id) & nzchar(as.character(cluster_id))]
#   res <- unique(res, by = c("cluster_set","cluster_id","hi_floor_yr","pr_floor_yr"))
#   
#   out_qs  <- fs::path(OUT_DIR, paste0("cluster_characteristics_", cluster_set, ".qs"))
#   out_csv <- fs::path(OUT_DIR, paste0("cluster_characteristics_", cluster_set, ".csv"))
#   
#   tmp <- paste0(out_qs, ".tmp")
#   if (file_exists(tmp)) file_delete(tmp)
#   qsave(res, tmp); file_move(tmp, out_qs)
#   fwrite(res, out_csv)
#   
#   message("✅ Wrote: ", out_qs)
#   message("✅ Wrote: ", out_csv)
#   message("Rows: ", nrow(res), " | Unique clusters: ", uniqueN(res$cluster_id))
#   message("Unique timestep_hr: ", paste(sort(unique(res$timestep_hr)), collapse = ", "))
# }

# =============================================================================
# Per-cluster characteristics for YEAR-CHUNKED clustered exceedance .qs files
#
# INPUTS (new format):
#   V:/users/hquintal/phd1_cluster_southeast_v2/data/out/06_clusters/hi{H}_pr{P}/{set}/*.qs
#
# OUTPUTS:
#   data/out/07_characteristics/hi{H}_pr{P}/cluster_characteristics_{set}.qs (+ .csv)
#   data/out/07_characteristics/hi{H}_pr{P}/cluster_summary_by_resolution.csv
#
# IMPORTANT FIX:
#   - timestep_hr is NEVER inferred from exceedance spacing
#   - timestep_hr is taken from dt$time_step_hours (written during clustering) else defaults to 24
#   - n_steps & duration_hr are computed from integer timestep indices t_ix
#
# NEW:
#   - also calculates per-resolution summary metrics:
#       # points
#       # noise
#       # clusters
#       median duration (days)
#       median extent (km2)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(fs)
  library(stringr)
  library(lubridate)
  library(future)
  library(future.apply)
  library(progressr)
  library(geosphere)
})

# ----------------------------
# USER SETTINGS
# ----------------------------
BASE <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out"

HI_FLOOR_YR <- 25
PR_FLOOR_YR <- 25

BASE_IN <- fs::path(BASE, sprintf("06_clusters/hi%d_pr%d", HI_FLOOR_YR, PR_FLOOR_YR))

OUT_DIR <- fs::path(BASE, sprintf("07_characteristics/hi%d_pr%d", HI_FLOOR_YR, PR_FLOOR_YR))
dir_create(OUT_DIR, recurse = TRUE)

N_WORKERS <- 30

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

# ----------------------------
# GEO HELPERS
# ----------------------------
km_per_deg_lat <- function(lat) 111.132
km_per_deg_lon <- function(lat) 111.320 * cospi(lat/180)

cell_area_km2_vec <- function(lat_vec, resx_deg, resy_deg) {
  km_y <- km_per_deg_lat(lat_vec) * resy_deg
  km_x <- km_per_deg_lon(lat_vec) * resx_deg
  km_x * km_y
}

grid_step <- function(v) {
  u <- sort(unique(v))
  if (length(u) < 2) return(NA_real_)
  suppressWarnings(min(diff(u)))
}

as_posix_utc <- function(tt) {
  if (inherits(tt, "POSIXct")) return(as.POSIXct(tt, tz = "UTC"))
  if (inherits(tt, "Date"))   return(as.POSIXct(tt, tz = "UTC"))
  suppressWarnings(as.POSIXct(tt, tz = "UTC"))
}

gc_km <- function(p1, p2) geosphere::distHaversine(p1, p2) / 1000

centroid_timestep <- function(x, y, w) {
  pos <- is.finite(w) & w > 0
  if (!any(pos)) return(c(lon = NA_real_, lat = NA_real_))
  c(lon = mean(x[pos]), lat = mean(y[pos]))
}

elongation_major_minor <- function(xs, ys, resx_deg, resy_deg) {
  n <- length(xs)
  lat0 <- mean(ys)
  
  base_major <- max(km_per_deg_lon(lat0) * resx_deg, km_per_deg_lat(lat0) * resy_deg, na.rm = TRUE)
  base_minor <- max(min(km_per_deg_lon(lat0) * resx_deg, km_per_deg_lat(lat0) * resy_deg, na.rm = TRUE), 1e-6)
  
  if (!is.finite(base_major)) base_major <- 0
  if (!is.finite(base_minor) || base_minor <= 0) base_minor <- 1e-6
  
  if (n < 2) {
    return(list(major_km = base_major, minor_km = base_minor, elongation_ratio = base_major / base_minor))
  }
  
  coords <- cbind(xs, ys)
  coords_h <- if (n >= 6) {
    idx <- chull(coords)
    coords[c(idx, idx[1]), , drop = FALSE]
  } else coords
  
  D <- geosphere::distm(coords_h)
  if (!length(D) || !is.finite(max(D))) {
    return(list(major_km = base_major, minor_km = base_minor, elongation_ratio = base_major / base_minor))
  }
  
  ij <- which(D == max(D), arr.ind = TRUE)[1, ]
  p1 <- coords_h[ij[1], ]
  p2 <- coords_h[ij[2], ]
  major_km <- max(D) / 1000
  
  theta <- atan2(p2[2] - p1[2], p2[1] - p1[1]) + pi/2
  u_perp <- c(cos(theta), sin(theta))
  
  xkm <- (xs - p1[1]) * km_per_deg_lon(lat0)
  ykm <- (ys - p1[2]) * km_per_deg_lat(lat0)
  proj <- xkm * u_perp[1] + ykm * u_perp[2]
  
  minor_km <- diff(range(proj))
  if (!is.finite(minor_km) || minor_km <= 0) minor_km <- base_minor
  minor_km <- max(minor_km, 1e-6)
  
  list(
    major_km = as.numeric(major_km),
    minor_km = as.numeric(minor_km),
    elongation_ratio = as.numeric(major_km / minor_km)
  )
}

# ----------------------------
# COLUMN GUESSERS
# ----------------------------
normalize_xy <- function(dt) {
  if (all(c("x","y") %in% names(dt))) return(invisible(NULL))
  if (all(c("lon","lat") %in% names(dt))) {
    setnames(dt, c("lon","lat"), c("x","y"))
    return(invisible(NULL))
  }
  stop("Missing coordinate columns: need x,y (or lon,lat).")
}

guess_time_col <- function(dt) {
  for (cc in c("time_utc","time","datetime","timestamp","day","date")) {
    if (cc %in% names(dt)) return(cc)
  }
  stop("Missing time column (expected one of time_utc, time, day/date).")
}

guess_intensity_col <- function(dt, cluster_set) {
  if (grepl("^heatindex", cluster_set, ignore.case = TRUE)) {
    for (cc in c("hi3_intensity","hi_intensity","heatindex","HeatIndex","value","intensity")) {
      if (cc %in% names(dt)) return(cc)
    }
  }
  
  for (cc in c("pr24_intensity_mm","intensity_mm","mm","precip_mm","prcp_mm","precip","pr","ppt","value","intensity")) {
    if (cc %in% names(dt)) return(cc)
  }
  
  ignore <- c("x","y","lon","lat","time","time_utc","day","date",
              "cluster","cluster_id","eps","mu",
              "rp","rp_max_years","rp_max_layer","resolution_deg",
              "time_step_hours","t_ix","row_ix","col_ix","floor_yr","dataset_tag")
  num_cols <- names(dt)[vapply(dt, is.numeric, TRUE)]
  num_cols <- setdiff(num_cols, ignore)
  if (length(num_cols) == 1) return(num_cols[1])
  
  stop("Could not determine intensity column for set=", cluster_set,
       ". Rename intensity column or extend guess_intensity_col().")
}

guess_rp_col <- function(dt) {
  for (cc in c("rp_max_years","rp","max_rp_exceeded","rp_years","rp_max")) {
    if (cc %in% names(dt)) return(cc)
  }
  return(NA_character_)
}

# ----------------------------
# NOISE / POINT HELPERS
# ----------------------------
is_noise_point <- function(cluster_id, cluster = NULL) {
  cid_chr <- trimws(as.character(cluster_id))
  cid_low <- tolower(cid_chr)
  
  out <- is.na(cluster_id) |
    cid_chr == "" |
    cid_low %in% c("na", "nan", "noise", "-1", "0")
  
  if (!is.null(cluster)) {
    out <- out | is.na(cluster) | cluster %in% c(-1, 0)
  }
  
  out
}

# ----------------------------
# PER-FILE PROCESSOR
# ----------------------------
process_one_file <- function(path, cluster_set) {
  base <- path_file(path)
  dt0 <- tryCatch(as.data.table(qread(path)), error = function(e) data.table())
  if (!nrow(dt0)) return(NULL)
  
  normalize_xy(dt0)
  
  if (!("cluster_id" %in% names(dt0))) {
    warning("Skipping (no cluster_id): ", base)
    return(NULL)
  }
  if (!("cluster" %in% names(dt0))) dt0[, cluster := NA_integer_]
  
  # count points / noise BEFORE filtering to assigned clusters
  n_points_total <- nrow(dt0)
  noise_idx <- is_noise_point(dt0$cluster_id, dt0$cluster)
  n_noise_total <- sum(noise_idx, na.rm = TRUE)
  n_cluster_points <- n_points_total - n_noise_total
  
  # keep only assigned clusters
  dt <- copy(dt0[!noise_idx])
  if (!nrow(dt)) return(NULL)
  
  time_col  <- guess_time_col(dt)
  inten_col <- guess_intensity_col(dt, cluster_set)
  
  setnames(dt, time_col,  "time_raw")
  setnames(dt, inten_col, "intensity")
  
  dt[, time := {
    if (inherits(time_raw, "Date")) as_posix_utc(time_raw)
    else if (inherits(time_raw, "POSIXct")) as_posix_utc(time_raw)
    else {
      tt <- suppressWarnings(ymd_hms(time_raw, tz = "UTC"))
      if (all(is.na(tt))) tt <- suppressWarnings(ymd(time_raw, tz = "UTC"))
      if (all(is.na(tt))) tt <- as_posix_utc(time_raw)
      tt
    }
  }]
  
  dt <- dt[is.finite(x) & is.finite(y) & is.finite(intensity) & !is.na(time)]
  if (!nrow(dt)) return(NULL)
  
  # resolution
  resx_deg <- if ("resolution_deg" %in% names(dt)) unique(dt$resolution_deg)[1] else grid_step(dt$x)
  resy_deg <- if ("resolution_deg" %in% names(dt)) unique(dt$resolution_deg)[1] else grid_step(dt$y)
  if (!is.finite(resx_deg) || resx_deg <= 0) resx_deg <- grid_step(dt$x)
  if (!is.finite(resy_deg) || resy_deg <= 0) resy_deg <- grid_step(dt$y)
  
  dt[, area_km2 := if (is.finite(resx_deg) && is.finite(resy_deg)) cell_area_km2_vec(y, resx_deg, resy_deg) else NA_real_]
  
  # ---- FIXED timestep: NEVER infer from exceedance spacing ----
  timestep_hr <- if ("time_step_hours" %in% names(dt)) as.numeric(unique(dt$time_step_hours)[1]) else 24
  if (!is.finite(timestep_hr) || timestep_hr <= 0) timestep_hr <- 24
  
  # t_ix: use existing if present, else compute from cluster start
  if (!("t_ix" %in% names(dt))) {
    dt[, start0 := min(time, na.rm = TRUE), by = cluster_id]
    dt[, t_ix := as.integer(round(as.numeric(difftime(time, start0, units = "hours")) / timestep_hr)) + 1L]
    dt[, start0 := NULL]
  }
  
  # bounds
  time_bounds <- dt[, .(
    start_time = min(time, na.rm = TRUE),
    end_time   = max(time, na.rm = TRUE)
  ), by = cluster_id]
  
  # duration based on t_ix (inclusive)
  duration_tbl <- dt[, .(
    t_min = min(t_ix, na.rm = TRUE),
    t_max = max(t_ix, na.rm = TRUE)
  ), by = cluster_id][, .(
    n_steps = as.integer(t_max - t_min + 1L),
    duration_hr = as.numeric(t_max - t_min + 1L) * timestep_hr
  ), by = cluster_id]
  
  # intensity stats
  intensity_stats <- dt[, .(
    peak_intensity   = max(intensity, na.rm = TRUE),
    mean_intensity   = mean(intensity, na.rm = TRUE),
    median_intensity = as.numeric(stats::median(intensity, na.rm = TRUE))
  ), by = cluster_id]
  
  # extents by timestep
  extent_by_t <- dt[, .(t_extent_km2 = sum(area_km2, na.rm = TRUE)),
                    by = .(cluster_id, t_ix)]
  max_t_extent <- extent_by_t[, .(max_t_extent_km2 = max(t_extent_km2, na.rm = TRUE)), by = cluster_id]
  
  union_extent <- unique(dt[, .(cluster_id, x, y, area_km2)])
  total_extent <- union_extent[, .(total_extent_km2 = sum(area_km2, na.rm = TRUE)), by = cluster_id]
  
  # volume only for precipitation sets
  is_precip <- grepl("^precipitation", cluster_set, ignore.case = TRUE)
  if (is_precip) {
    vol_tbl <- dt[, .(
      total_volume_km3 = sum((intensity / 1000) * (area_km2 * 1e6), na.rm = TRUE) / 1e9
    ), by = cluster_id]
  } else {
    vol_tbl <- unique(dt[, .(cluster_id, total_volume_km3 = as.numeric(NA))])
  }
  
  # centroid speeds
  centroids <- dt[, {
    cxy <- centroid_timestep(x, y, intensity)
    .(lon = cxy[["lon"]], lat = cxy[["lat"]])
  }, by = .(cluster_id, t_ix)]
  setorder(centroids, cluster_id, t_ix)
  
  speeds <- centroids[, {
    if (.N < 2) {
      data.table(max_speed_kmhr = 0, mean_speed_kmhr = 0, median_speed_kmhr = 0)
    } else {
      lon1 <- head(lon, -1); lat1 <- head(lat, -1)
      lon2 <- tail(lon, -1); lat2 <- tail(lat, -1)
      dstep <- pmax(1, diff(t_ix))
      dkm <- mapply(function(a,b,c,d) gc_km(c(a,b), c(c,d)), lon1, lat1, lon2, lat2)
      dh <- dstep * timestep_hr
      spd <- dkm / dh
      spd <- spd[is.finite(spd)]
      if (!length(spd)) spd <- 0
      data.table(
        max_speed_kmhr    = max(spd, na.rm = TRUE),
        mean_speed_kmhr   = mean(spd, na.rm = TRUE),
        median_speed_kmhr = as.numeric(stats::median(spd, na.rm = TRUE))
      )
    }
  }, by = cluster_id]
  
  # elongation
  elong_tbl <- union_extent[, {
    el <- elongation_major_minor(x, y, resx_deg, resy_deg)
    data.table(
      elongation_ratio = el$elongation_ratio,
      major_axis_km    = el$major_km,
      minor_axis_km    = el$minor_km
    )
  }, by = cluster_id]
  
  # RP per cluster
  rp_col <- guess_rp_col(dt)
  if (is.na(rp_col)) {
    rp_tbl <- unique(dt[, .(cluster_id, rp = as.numeric(NA))])
  } else {
    rp_tbl <- dt[is.finite(get(rp_col)), .(rp = max(get(rp_col), na.rm = TRUE)), by = cluster_id]
  }
  
  cluster_map <- unique(dt[, .(cluster_id, cluster)])
  file_col    <- data.table(cluster_id = unique(dt$cluster_id), source_file = base)
  
  # floors
  hi_floor <- if ("hi_floor_yr" %in% names(dt)) unique(dt$hi_floor_yr)[1] else NA_real_
  pr_floor <- if ("pr_floor_yr" %in% names(dt)) unique(dt$pr_floor_yr)[1] else NA_real_
  if ("floor_yr" %in% names(dt)) {
    if (grepl("^heatindex", cluster_set, ignore.case = TRUE)) hi_floor <- unique(dt$floor_yr)[1]
    if (grepl("^precipitation", cluster_set, ignore.case = TRUE)) pr_floor <- unique(dt$floor_yr)[1]
  }
  floor_tbl <- unique(data.table(cluster_id = unique(dt$cluster_id),
                                 hi_floor_yr = hi_floor,
                                 pr_floor_yr = pr_floor))
  
  # file-level summary repeated per cluster row so it can be recovered later
  counts_tbl <- unique(data.table(
    cluster_id        = unique(dt$cluster_id),
    n_points_file     = n_points_total,
    n_noise_file      = n_noise_total,
    n_clusterpts_file = n_cluster_points
  ))
  
  out <- Reduce(function(a,b) merge(a, b, by = "cluster_id", all = TRUE), list(
    data.table(cluster_id = unique(dt$cluster_id), cluster_set = cluster_set),
    floor_tbl, cluster_map, file_col, counts_tbl, rp_tbl, time_bounds, duration_tbl,
    intensity_stats, max_t_extent, total_extent, vol_tbl,
    speeds, elong_tbl
  ))
  
  out[, resolution_deg := if (is.finite(resx_deg)) resx_deg else NA_real_]
  out[, timestep_hr := timestep_hr]
  
  setcolorder(out, c(
    "cluster_set","hi_floor_yr","pr_floor_yr",
    "cluster_id","cluster","rp","source_file",
    "n_points_file","n_noise_file","n_clusterpts_file",
    "start_time","end_time","n_steps","timestep_hr","duration_hr",
    "peak_intensity","mean_intensity","median_intensity",
    "max_t_extent_km2","total_extent_km2","total_volume_km3",
    "max_speed_kmhr","mean_speed_kmhr","median_speed_kmhr",
    "elongation_ratio","major_axis_km","minor_axis_km",
    "resolution_deg"
  ))
  
  out[]
}

# ----------------------------
# RUN
# ----------------------------
plan(multisession, workers = N_WORKERS)
handlers(global = TRUE)

summary_list <- list()

for (cluster_set in CLUSTER_SETS) {
  in_dir <- fs::path(BASE_IN, cluster_set)
  if (!dir_exists(in_dir)) {
    warning("Missing cluster set folder (skipping): ", in_dir)
    next
  }
  
  qs_files <- dir_ls(in_dir, recurse = TRUE, type = "file", glob = "*.qs")
  qs_files <- qs_files[!grepl("characteristics", qs_files, ignore.case = TRUE)]
  qs_files <- sort(as.character(qs_files))
  
  if (!length(qs_files)) {
    warning("No .qs found under: ", in_dir)
    next
  }
  
  message("============================================================")
  message("Cluster set: ", cluster_set)
  message("Inputs: ", length(qs_files))
  message("============================================================")
  
  res <- with_progress({
    p <- progressor(steps = length(qs_files))
    parts <- future_lapply(qs_files, function(f) {
      p(path_file(f))
      process_one_file(f, cluster_set = cluster_set)
    })
    rbindlist(Filter(Negate(is.null), parts), use.names = TRUE, fill = TRUE)
  })
  
  if (!nrow(res)) {
    warning("No output rows for cluster_set=", cluster_set)
    next
  }
  
  setDT(res)
  res <- res[!is.na(cluster_id) & nzchar(as.character(cluster_id))]
  res <- unique(res, by = c("cluster_set","cluster_id","hi_floor_yr","pr_floor_yr"))
  
  # convert to days for convenience
  res[, duration_days := duration_hr / 24]
  
  out_qs  <- fs::path(OUT_DIR, paste0("cluster_characteristics_", cluster_set, ".qs"))
  out_csv <- fs::path(OUT_DIR, paste0("cluster_characteristics_", cluster_set, ".csv"))
  
  tmp <- paste0(out_qs, ".tmp")
  if (file_exists(tmp)) file_delete(tmp)
  qsave(res, tmp); file_move(tmp, out_qs)
  fwrite(res, out_csv)
  
  message("✅ Wrote: ", out_qs)
  message("✅ Wrote: ", out_csv)
  message("Rows: ", nrow(res), " | Unique clusters: ", uniqueN(res$cluster_id))
  message("Unique timestep_hr: ", paste(sort(unique(res$timestep_hr)), collapse = ", "))
  
  # ----------------------------
  # per-resolution summary
  # ----------------------------
  # sum file-level counts once per source file
  file_counts <- unique(
    res[, .(cluster_set, resolution_deg, source_file, n_points_file, n_noise_file, n_clusterpts_file)]
  )
  
  sum_counts <- file_counts[, .(
    n_points        = sum(n_points_file, na.rm = TRUE),
    n_noise         = sum(n_noise_file, na.rm = TRUE),
    n_clusterpoints = sum(n_clusterpts_file, na.rm = TRUE)
  ), by = .(cluster_set, resolution_deg)]
  
  cluster_stats <- res[, .(
    n_clusters            = uniqueN(cluster_id),
    median_duration_days  = as.numeric(stats::median(duration_days, na.rm = TRUE)),
    median_extent_km2     = as.numeric(stats::median(total_extent_km2, na.rm = TRUE))
  ), by = .(cluster_set, resolution_deg)]
  
  summary_dt <- merge(sum_counts, cluster_stats,
                      by = c("cluster_set", "resolution_deg"), all = TRUE)
  
  summary_dt[, `:=`(
    hi_floor_yr = unique(na.omit(res$hi_floor_yr))[1],
    pr_floor_yr = unique(na.omit(res$pr_floor_yr))[1]
  )]
  
  setcolorder(summary_dt, c(
    "cluster_set", "resolution_deg", "hi_floor_yr", "pr_floor_yr",
    "n_points", "n_noise", "n_clusterpoints", "n_clusters",
    "median_duration_days", "median_extent_km2"
  ))
  
  summary_list[[cluster_set]] <- summary_dt
}

# ----------------------------
# WRITE COMBINED SUMMARY
# ----------------------------
if (length(summary_list)) {
  summary_all <- rbindlist(summary_list, use.names = TRUE, fill = TRUE)
  
  out_summary_csv <- fs::path(OUT_DIR, "cluster_summary_by_resolution.csv")
  out_summary_qs  <- fs::path(OUT_DIR, "cluster_summary_by_resolution.qs")
  
  tmp_sum <- paste0(out_summary_qs, ".tmp")
  if (file_exists(tmp_sum)) file_delete(tmp_sum)
  qsave(summary_all, tmp_sum); file_move(tmp_sum, out_summary_qs)
  fwrite(summary_all, out_summary_csv)
  
  message("============================================================")
  message("✅ Wrote combined summary: ", out_summary_qs)
  message("✅ Wrote combined summary: ", out_summary_csv)
  print(summary_all)
}
