# scripts/07_experimental_covariance_daily_NCAregions_parallel_relic.R
# =============================================================================
# DAILY experimental covariance (RELIC estimator) by NCA Region — PARALLEL
#
# Computes (for each day × region):
#   - space-only covariance vs spatial lag (tLag = 0)
#   - time-only  covariance vs temporal lag (rLag = 0)
# Trims AFTER first negative, but KEEPS first negative (inclusive)
#
# Supports BOTH:
#   - Heat Index: uses hourly heat_index within each day
#   - Precipitation: uses hourly intensity within each day, thresholding < 1 mm/hr to NA (configurable)
#
# Outputs:
#   V:/users/hquintal/phd1_cluster_southeast_v2/data/out/03_cov_model/nca/
#     heatindex/daily/   heatindex_<REGION>_day_<YYYY-MM-DD>_{cov_space,cov_time}.txt
#     precipitation/daily/ precip_<REGION>_day_<YYYY-MM-DD>_{cov_space,cov_time}.txt
#   plus run logs and the (cached) cell->region lookup table(s).
# =============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(data.table)
  library(lubridate)
  library(sf)
  library(future)
  library(future.apply)
  library(progressr)
  library(fs)
})

# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------
IN_DIR_HI <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/heatindex"
IN_DIR_PR <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation"

NCA_SHP <- "V:/users/hquintal/phd2_cluster_conus/data/in/shapefiles/NCA_Regions/NCA_Regions.shp"

OUT_BASE <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/03_cov_model/nca"
OUT_HI   <- file.path(OUT_BASE, "heatindex", "daily")
OUT_PR   <- file.path(OUT_BASE, "precipitation", "daily")

dir.create(OUT_HI, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_PR, recursive = TRUE, showWarnings = FALSE)

WORKERS <- 16
plan(multisession, workers = WORKERS)

# Performance / safety (controls kron vs kronloop)
max_pairs <- 2e7

# Precip hourly threshold (set 0 to keep all non-NA hourly values)
PRECIP_HOURLY_MIN_MM <- 1

# -----------------------------------------------------------------------------
# RELIC CORE (from your scripts, unchanged)
# -----------------------------------------------------------------------------
coord2dist <- function(c1, c2) {
  if (is.data.frame(c1)) c1 <- as.matrix(c1)
  if (is.data.frame(c2)) c2 <- as.matrix(c2)
  if (is.list(c1)) { c1 <- c1[[1]]; c2 <- c2[[1]] }
  x1_sq <- rowSums(c1^2); x2_sq <- rowSums(c2^2)
  D2 <- outer(x1_sq, rep(1, nrow(c2))) +
    outer(rep(1, nrow(c1)), x2_sq) -
    2 * tcrossprod(c1, c2)
  D2[D2 < 0] <- 0
  sqrt(D2)
}

pairsindex <- function(c1, c2, rLag, rLagTol, dist=list("coord2dist"), method=NULL) {
  if (is.data.frame(c1)) c1 <- as.matrix(c1)
  if (is.data.frame(c2)) c2 <- as.matrix(c2)
  n1 <- nrow(c1); n2 <- nrow(c2); nr <- length(rLag)
  if (!is.list(dist)) stop("dist must be a list")
  
  if (is.null(method)) {
    method <- if (n1 * n2 <= max_pairs) "kron" else "kronloop"
  }
  
  idxpairs <- vector("list", nr)
  
  call_dist <- function(A, B) {
    if (length(dist) == 1) {
      do.call(dist[[1]], list(A, B))
    } else {
      do.call(dist[[1]], list(A, B, dist[[2]]))
    }
  }
  
  switch(method,
         kron = {
           D <- call_dist(c1, c2)
           for (ir in seq_len(nr)) {
             mask <- (D >= rLag[ir] - rLagTol[ir]) & (D <= rLag[ir] + rLagTol[ir])
             idxpairs[[ir]] <- which(mask, arr.ind = TRUE)
           }
         },
         kronloop = {
           for (ir in seq_len(nr)) idxpairs[[ir]] <- matrix(numeric(0), ncol = 2)
           for (i in seq_len(n1)) {
             d_i <- call_dist(matrix(c1[i,], nrow = 1), c2)
             for (ir in seq_len(nr)) {
               js <- which(d_i >= rLag[ir] - rLagTol[ir] & d_i <= rLag[ir] + rLagTol[ir])
               if (length(js) > 0) idxpairs[[ir]] <- rbind(idxpairs[[ir]], cbind(i, js))
             }
           }
         },
         stop("method must be 'kron' or 'kronloop'")
  )
  
  for (ir in seq_len(nr)) {
    if (nrow(idxpairs[[ir]]) > 0) {
      idxpairs[[ir]] <- idxpairs[[ir]][order(idxpairs[[ir]][,1], idxpairs[[ir]][,2]), , drop = FALSE]
    }
  }
  idxpairs
}

stcov <- function(
    Zi, cMSi, tMEi,
    Zj, cMSj, tMEj,
    rLag, rLagTol,
    tLag, tLagTol,
    dist = list("coord2dist"),
    method = NULL
) {
  nMSi <- nrow(Zi);  nMEi <- ncol(Zi)
  if (!all(dim(cMSi) == c(nMSi, 2))) stop("cMSi must be nMSi × 2")
  if (length(tMEi) != nMEi) stop("tMEi must be length nMEi")
  
  nMSj <- nrow(Zj);  nMEj <- ncol(Zj)
  if (!all(dim(cMSj) == c(nMSj, 2))) stop("cMSj must be nMSj × 2")
  if (length(tMEj) != nMEj) stop("tMEj must be length nMEj")
  
  nr <- length(rLag); if (length(rLagTol) != nr) stop("rLag and rLagTol length mismatch")
  nt <- length(tLag); if (length(tLagTol) != nt) stop("tLag and tLagTol length mismatch")
  
  if (is.null(method)) {
    method <- if (nMSi * nMSj <= max_pairs) "kron" else "kronloop"
  }
  
  idxpairsMS <- pairsindex(cMSi, cMSj, rLag, rLagTol, dist, method)
  tMat <- outer(tMEi, tMEj, function(ti, tj) tj - ti)
  
  C  <- matrix(NA_real_, nrow = nr, ncol = nt)
  np <- matrix(0L,       nrow = nr, ncol = nt)
  
  for (ir in seq_len(nr)) {
    pairs_ij <- idxpairsMS[[ir]]
    if (nrow(pairs_ij) == 0L) next
    
    for (it in seq_len(nt)) {
      mask <- (tLag[it] - tLagTol[it] <= tMat) & (tMat <= tLag[it] + tLagTol[it])
      if (!any(mask)) next
      ev_idx <- which(mask, arr.ind = TRUE)
      
      nSP <- nrow(pairs_ij); nEP <- nrow(ev_idx)
      sp_rep <- pairs_ij[rep(seq_len(nSP), times = nEP), , drop = FALSE]
      ep_rep <- ev_idx[rep(seq_len(nEP), each = nSP), , drop = FALSE]
      
      head_vals <- Zi[cbind(sp_rep[,1], ep_rep[,1])]
      tail_vals <- Zj[cbind(sp_rep[,2], ep_rep[,2])]
      
      valid <- !is.na(head_vals) & !is.na(tail_vals)
      np[ir,it] <- sum(valid)
      if (np[ir,it] == 0L) next
      
      h <- head_vals[valid]
      t <- tail_vals[valid]
      C[ir,it] <- mean(h * t) - mean(h) * mean(t)
    }
  }
  
  list(C = C, np = np, method = method)
}

trim_after_negative_inclusive <- function(x, lags, np = NULL) {
  keep <- !is.na(x) & !is.na(lags)
  x <- x[keep]; lags <- lags[keep]
  if (!is.null(np)) np <- np[keep]
  
  if (length(x) == 0L) return(list(cov = numeric(0), lags = numeric(0), np = integer(0)))
  
  neg_idx <- which(x < 0)[1]
  if (!is.na(neg_idx)) {
    x <- x[1:neg_idx]
    lags <- lags[1:neg_idx]
    if (!is.null(np)) np <- np[1:neg_idx]
  }
  
  if (is.null(np)) np <- rep(NA_integer_, length(x))
  list(cov = x, lags = lags, np = np)
}

estimate_base_dist <- function(coords) {
  ux <- sort(unique(coords[,1]))
  uy <- sort(unique(coords[,2]))
  dx <- diff(ux); dy <- diff(uy)
  dx <- dx[dx > 0]; dy <- dy[dy > 0]
  base <- suppressWarnings(min(median(dx, na.rm = TRUE), median(dy, na.rm = TRUE), na.rm = TRUE))
  if (!is.finite(base) || base <= 0) base <- 0.5
  base
}

make_lags <- function(max_dist, max_time_hr, base_dist) {
  r_raw <- base_dist * c(0:4, 6, 8, 12)
  rLag <- sort(unique(r_raw[r_raw <= max_dist]))
  if (length(rLag) < 2) rLag <- c(0, min(base_dist, max_dist))
  rLagTol <- c(0, diff(rLag)/2)
  
  t_raw <- c(0:4, 6, 8, 12)
  tLag <- sort(unique(t_raw[t_raw <= max_time_hr]))
  if (length(tLag) < 2) tLag <- c(0, min(1, max_time_hr))
  tLagTol <- c(0, diff(tLag)/2)
  
  list(rLag = rLag, rLagTol = rLagTol, tLag = tLag, tLagTol = tLagTol)
}

mb <- function(x) as.numeric(object.size(x)) / (1024^2)

# -----------------------------------------------------------------------------
# GENERIC HELPERS
# -----------------------------------------------------------------------------
ensure_datetime_utc <- function(dt, qs_file) {
  if ("datetime" %in% names(dt) && inherits(dt$datetime, "POSIXct")) {
    attr(dt$datetime, "tzone") <- "UTC"
    return(dt)
  }
  if ("datetime_utc" %in% names(dt)) {
    dt[, datetime := as.POSIXct(datetime_utc, tz = "UTC")]
    return(dt)
  }
  stop("No datetime or datetime_utc column found in: ", qs_file)
}

replace_below_threshold_with_na <- function(x, thresh) {
  if (!is.finite(thresh) || thresh <= 0) return(x)
  x[!is.na(x) & x < thresh] <- NA_real_
  x
}

# Build a day matrix from already-filtered dt with columns: lon, lat, datetime, value
day_dt_to_matrix <- function(dt, value_col) {
  dt <- dt[!is.na(get(value_col)), .(lon, lat, datetime, val = as.numeric(get(value_col)))]
  if (nrow(dt) == 0L) return(NULL)
  
  dt[, cell_id := paste0(sprintf("%.6f", lon), "_", sprintf("%.6f", lat))]
  
  times <- sort(unique(dt$datetime))
  if (length(times) < 2L) return(NULL)
  
  dt[, t_idx := match(datetime, times)]
  
  setkey(dt, cell_id, t_idx)
  wide <- dcast(
    dt,
    cell_id + lon + lat ~ t_idx,
    value.var = "val",
    fun.aggregate = mean
  )
  if (nrow(wide) < 2L) return(NULL)
  
  coords <- as.matrix(wide[, .(lon, lat)])
  time_cols <- setdiff(names(wide), c("cell_id", "lon", "lat"))
  time_cols_ord <- time_cols[order(as.integer(time_cols))]
  
  values <- as.matrix(wide[, ..time_cols_ord])
  if (ncol(values) < 2L) return(NULL)
  
  t_hr <- as.numeric(difftime(times, times[1], units = "hours"))
  
  list(coords = coords, values = values, t_hr = t_hr)
}

# -----------------------------------------------------------------------------
# NCA REGION LOOKUP (cache cell_id -> RegionName)
# -----------------------------------------------------------------------------
# We build separate caches for HI and PR grids (in case grids differ).
build_cell_to_region_lookup <- function(example_qs, out_qs) {
  if (file.exists(out_qs)) {
    message("[lookup] using cached: ", out_qs)
    return(qread(out_qs))
  }
  
  message("[lookup] building from example: ", example_qs)
  
  # read example day and extract unique lon/lat
  dt0 <- as.data.table(qread(example_qs))
  stopifnot(all(c("lon","lat") %in% names(dt0)))
  dt0 <- unique(dt0[, .(lon, lat)])
  
  # NCA polygons (your shapefile is in EPSG:3857). Transform to WGS84 lon/lat.
  nca <- suppressWarnings(st_read(NCA_SHP, quiet = TRUE))
  nca <- st_make_valid(nca)
  nca_ll <- st_transform(nca, 4326)
  
  # points in lon/lat WGS84
  pts <- st_as_sf(dt0, coords = c("lon","lat"), crs = 4326, remove = FALSE)
  
  # spatial join: keep only RegionName
  joined <- st_join(pts, nca_ll[, c("RegionName")], left = TRUE, join = st_within)
  
  out <- as.data.table(joined)
  out[, cell_id := paste0(sprintf("%.6f", lon), "_", sprintf("%.6f", lat))]
  out <- out[, .(cell_id, RegionName)]
  out[, RegionName := as.character(RegionName)]
  
  # Some points may fall outside (NA). Keep them as NA and we’ll drop during processing.
  qsave(out, out_qs)
  message("[lookup] wrote: ", out_qs)
  out
}

# -----------------------------------------------------------------------------
# PROCESS ONE (day × region) FOR A GIVEN VARIABLE
# -----------------------------------------------------------------------------
process_one_day_region <- function(qs_file, region_name, lookup_dt,
                                   var_tag = c("heatindex","precip"),
                                   value_col,
                                   out_dir,
                                   precip_hourly_min_mm = 0) {
  
  var_tag <- match.arg(var_tag)
  day_str <- gsub("\\.qs$", "", basename(qs_file))
  t0 <- Sys.time()
  
  # read + datetime
  dt <- as.data.table(qread(qs_file))
  dt <- ensure_datetime_utc(dt, qs_file)
  
  # value col
  if (!value_col %in% names(dt)) {
    return(list(day = day_str, region = region_name, ok = FALSE,
                reason = paste0("missing_col:", value_col),
                dt_sec = as.numeric(difftime(Sys.time(), t0, units="secs")),
                n_space = 0L, n_time = 0L, method = NA_character_,
                mem_coords_mb = NA_real_, mem_values_mb = NA_real_,
                n_points_in = nrow(dt), n_points_kept = 0L))
  }
  
  # stable cell id for join
  dt[, cell_id := paste0(sprintf("%.6f", lon), "_", sprintf("%.6f", lat))]
  
  # attach region
  dt <- merge(dt, lookup_dt, by = "cell_id", all.x = TRUE, sort = FALSE)
  
  # filter to this region (drop NA RegionName)
  dt <- dt[!is.na(RegionName) & RegionName == region_name]
  
  if (nrow(dt) == 0L) {
    return(list(day = day_str, region = region_name, ok = FALSE,
                reason = "no_points_in_region",
                dt_sec = as.numeric(difftime(Sys.time(), t0, units="secs")),
                n_space = 0L, n_time = 0L, method = NA_character_,
                mem_coords_mb = NA_real_, mem_values_mb = NA_real_,
                n_points_in = NA_integer_, n_points_kept = 0L))
  }
  
  # precip thresholding (hourly)
  if (var_tag == "precip") {
    dt[, (value_col) := replace_below_threshold_with_na(get(value_col), precip_hourly_min_mm)]
  }
  
  # matrix
  mat <- day_dt_to_matrix(dt, value_col = value_col)
  if (is.null(mat)) {
    return(list(day = day_str, region = region_name, ok = FALSE,
                reason = "NULL_matrix",
                dt_sec = as.numeric(difftime(Sys.time(), t0, units="secs")),
                n_space = 0L, n_time = 0L, method = NA_character_,
                mem_coords_mb = NA_real_, mem_values_mb = NA_real_,
                n_points_in = nrow(dt), n_points_kept = sum(!is.na(dt[[value_col]]))))
  }
  
  coords <- mat$coords
  values <- mat$values
  t_hr   <- mat$t_hr
  
  base_dist <- estimate_base_dist(coords)
  max_dist  <- max(coord2dist(coords, coords), na.rm = TRUE)
  max_time  <- max(t_hr, na.rm = TRUE)
  lags <- make_lags(max_dist, max_time, base_dist)
  
  mem_coords_mb <- mb(coords)
  mem_values_mb <- mb(values)
  
  # space-only
  out_space <- stcov(values, coords, t_hr, values, coords, t_hr,
                     rLag = lags$rLag, rLagTol = lags$rLagTol,
                     tLag = 0,         tLagTol = 0,
                     method = NULL)
  Cr  <- out_space$C[, 1]
  npr <- out_space$np[, 1]
  
  # time-only
  out_time <- stcov(values, coords, t_hr, values, coords, t_hr,
                    rLag = 0,         rLagTol = 0,
                    tLag = lags$tLag, tLagTol = lags$tLagTol,
                    method = out_space$method)
  
  Ct  <- out_time$C[1, ]
  npt <- out_time$np[1, ]
  
  res_space <- trim_after_negative_inclusive(Cr, lags$rLag, np = npr)
  res_time  <- trim_after_negative_inclusive(Ct, lags$tLag, np = npt)
  
  dt_space <- unique(data.table(rLag = res_space$lags, Cov = res_space$cov, np = res_space$np), by = "rLag")
  dt_time  <- unique(data.table(tLag = res_time$lags,  Cov = res_time$cov,  np = res_time$np),  by = "tLag")
  
  # safe region tag for filenames
  region_tag <- gsub("[^A-Za-z0-9]+", "_", region_name)
  
  tag <- sprintf("%s_%s_day_%s", var_tag, region_tag, day_str)
  fwrite(dt_space, file.path(out_dir, paste0(tag, "_cov_space.txt")), sep = "\t")
  fwrite(dt_time,  file.path(out_dir, paste0(tag, "_cov_time.txt")),  sep = "\t")
  
  t1 <- Sys.time()
  list(day = day_str, region = region_name, ok = TRUE,
       reason = NA_character_,
       dt_sec = as.numeric(difftime(t1, t0, units="secs")),
       n_space = nrow(dt_space), n_time = nrow(dt_time),
       method = out_space$method,
       mem_coords_mb = mem_coords_mb, mem_values_mb = mem_values_mb,
       n_points_in = nrow(dt), n_points_kept = sum(!is.na(dt[[value_col]])))
}

# -----------------------------------------------------------------------------
# RUNNER (daily, per region)
# -----------------------------------------------------------------------------
run_daily_by_region <- function(in_dir, out_dir, var_tag,
                                value_col,
                                lookup_cache_qs,
                                precip_hourly_min_mm = 0) {
  
  files <- sort(list.files(in_dir, pattern = "^\\d{4}-\\d{2}-\\d{2}\\.qs$", full.names = TRUE))
  stopifnot(length(files) > 0)
  
  # Build / load lookup from an example file
  lookup_dt <- build_cell_to_region_lookup(example_qs = files[1], out_qs = lookup_cache_qs)
  
  regions <- sort(unique(na.omit(lookup_dt$RegionName)))
  if (length(regions) == 0L) stop("No regions found in lookup (all NA RegionName). Check CRS/join.")
  
  # Task table = day × region
  day_strs <- gsub("\\.qs$", "", basename(files))
  map_dt <- CJ(day = day_strs, region = regions, unique = TRUE)
  map_dt[, file := files[match(day, day_strs)]]
  
  # Output paths for skip check
  region_tag <- gsub("[^A-Za-z0-9]+", "_", map_dt$region)
  tags <- sprintf("%s_%s_day_%s", var_tag, region_tag, map_dt$day)
  
  out_space_paths <- file.path(out_dir, paste0(tags, "_cov_space.txt"))
  out_time_paths  <- file.path(out_dir, paste0(tags, "_cov_time.txt"))
  
  have_space <- file.exists(out_space_paths)
  have_time  <- file.exists(out_time_paths)
  is_done <- have_space & have_time
  
  message(sprintf(
    "[%s NCA DAILY SKIP CHECK] total=%d | done=%d | todo=%d | space_missing=%d | time_missing=%d",
    var_tag, nrow(map_dt), sum(is_done), sum(!is_done), sum(!have_space), sum(!have_time)
  ))
  
  todo_dt <- map_dt[!is_done]
  handlers(global = TRUE)
  
  t_all0 <- Sys.time()
  
  if (nrow(todo_dt) == 0L) {
    message("✅ Nothing to do — all outputs already exist.")
    
    log_dt <- copy(map_dt)
    log_dt[, `:=`(
      ok = TRUE, reason = NA_character_, dt_sec = NA_real_,
      n_space = NA_integer_, n_time = NA_integer_,
      method = NA_character_,
      mem_coords_mb = NA_real_, mem_values_mb = NA_real_,
      n_points_in = NA_integer_, n_points_kept = NA_integer_,
      had_cov_space = have_space, had_cov_time = have_time,
      skipped_already_done = TRUE
    )]
    
    fwrite(log_dt, file.path(out_dir, sprintf("run_log_parallel_relic_daily_%s_NCA.csv", var_tag)))
    message(sprintf("✅ DONE | workers=%d | total_elapsed=%s",
                    WORKERS, format(difftime(Sys.time(), t_all0), digits = 3)))
    return(invisible(NULL))
  }
  
  # Parallel execution over (day, region)
  with_progress({
    p <- progressor(along = seq_len(nrow(todo_dt)))
    results <- future_lapply(seq_len(nrow(todo_dt)), function(i) {
      p(sprintf("%s | %s", todo_dt$day[i], todo_dt$region[i]))
      process_one_day_region(
        qs_file = todo_dt$file[i],
        region_name = todo_dt$region[i],
        lookup_dt = lookup_dt,
        var_tag = var_tag,
        value_col = value_col,
        out_dir = out_dir,
        precip_hourly_min_mm = precip_hourly_min_mm
      )
    })
  })
  
  log_done <- rbindlist(lapply(results, as.data.table), fill = TRUE)
  log_done[, `:=`(had_cov_space = TRUE, had_cov_time = TRUE, skipped_already_done = FALSE)]
  
  # skipped rows
  if (any(is_done)) {
    skip_dt <- map_dt[is_done]
    skip_dt[, `:=`(
      ok = TRUE, reason = NA_character_, dt_sec = NA_real_,
      n_space = NA_integer_, n_time = NA_integer_,
      method = NA_character_,
      mem_coords_mb = NA_real_, mem_values_mb = NA_real_,
      n_points_in = NA_integer_, n_points_kept = NA_integer_,
      had_cov_space = have_space[is_done], had_cov_time = have_time[is_done],
      skipped_already_done = TRUE
    )]
    log_dt <- rbindlist(list(log_done, skip_dt), use.names = TRUE, fill = TRUE)
  } else {
    log_dt <- log_done
  }
  
  # Final filesystem truth (robust if job failed mid-write)
  region_tag_all <- gsub("[^A-Za-z0-9]+", "_", log_dt$region)
  tags_all <- sprintf("%s_%s_day_%s", var_tag, region_tag_all, log_dt$day)
  log_dt[, `:=`(
    had_cov_space = file.exists(file.path(out_dir, paste0(tags_all, "_cov_space.txt"))),
    had_cov_time  = file.exists(file.path(out_dir, paste0(tags_all, "_cov_time.txt")))
  )]
  
  setorder(log_dt, day, region)
  fwrite(log_dt, file.path(out_dir, sprintf("run_log_parallel_relic_daily_%s_NCA.csv", var_tag)))
  
  message(sprintf("✅ DONE | workers=%d | processed=%d | skipped=%d | total_elapsed=%s",
                  WORKERS, nrow(todo_dt), sum(is_done),
                  format(difftime(Sys.time(), t_all0), digits = 3)))
}

# -----------------------------------------------------------------------------
# RUN: HEAT INDEX + PRECIP (daily, NCA regions)
# # -----------------------------------------------------------------------------
# # Heat index value column: "heat_index" (change here if your files differ)
# run_daily_by_region(
#   in_dir = IN_DIR_HI,
#   out_dir = OUT_HI,
#   var_tag = "heatindex",
#   value_col = "heat_index",
#   lookup_cache_qs = file.path(OUT_BASE, "lookup_cell_to_NCAregion_heatindex.qs"),
#   precip_hourly_min_mm = 0
# )

# Precip value column: your daily precip script uses "intensity" (change if needed)
run_daily_by_region(
  in_dir = IN_DIR_PR,
  out_dir = OUT_PR,
  var_tag = "precip",
  value_col = "intensity",
  lookup_cache_qs = file.path(OUT_BASE, "lookup_cell_to_NCAregion_precip.qs"),
  precip_hourly_min_mm = PRECIP_HOURLY_MIN_MM
)