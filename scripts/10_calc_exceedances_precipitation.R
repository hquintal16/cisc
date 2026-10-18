# # scripts/05_pr24_exceedances_native_vs_resampled_MAXRP_reanalysis_all_hourly_v1.R
# # ------------------------------------------------------------------------------
# # PR24 exceedances (24-hr daily summed precipitation) vs PR24 return periods
# #
# # Inputs:
# #   - Hourly precip daily files (qs): data/out/01_vars/precipitation/YYYY-MM-DD.qs
# #   - Return period raster (nc):      data/out/04_rp/precip_idf_pr24_return_periods.nc
# #   - Reference grids (tif):          data/out/09_performance/ref_grids/<tag>/ref_grid_<tag>.tif
# #
# # Outputs (year-chunked .qs):
# #   - Native (0.25°):                 data/out/05_exceedances/precipitation/pr24/0.25
# #   - Resampled (target; bilinear):   data/out/05_exceedances/precipitation/pr24/<target>_bilinear
# #   - Resampled (target; neigh-max):  data/out/05_exceedances/precipitation/pr24/<target>_max
# #
# # Targets included:
# #   - 0.17° (bilinear + max)
# #   - 0.12° (bilinear + max)
# # ------------------------------------------------------------------------------
# 
# suppressPackageStartupMessages({
#   library(terra)
#   library(data.table)
#   library(qs)
#   library(matrixStats)
#   library(future)
#   library(future.apply)
# })
# 
# # -----------------------------
# # USER SETTINGS
# # -----------------------------
# ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
# 
# IN_HOURLY_DIR <- file.path(ROOT, "data/out/01_vars/precipitation")
# RP_NC         <- file.path(ROOT, "data/out/04_rp/precip_idf_pr24_return_periods.nc")
# 
# REF_025 <- file.path(ROOT, "data/out/09_performance/ref_grids/precipitation_0.25/ref_grid_precipitation_0.25.tif")
# REF_017 <- file.path(ROOT, "data/out/09_performance/ref_grids/precipitation_0.17_bilinear/ref_grid_precipitation_0.17_bilinear.tif")
# REF_012 <- file.path(ROOT, "data/out/09_performance/ref_grids/precipitation_0.12_bilinear/ref_grid_precipitation_0.12_bilinear.tif")
# 
# # keep outputs compatible with legacy structure (NOT the earlier *_reanalysis folder)
# OUT_025     <- file.path(ROOT, "data/out/05_exceedances/precipitation/pr24/0.25")
# 
# OUT_017_BIL <- file.path(ROOT, "data/out/05_exceedances/precipitation/pr24/0.17_bilinear")
# OUT_017_MAX <- file.path(ROOT, "data/out/05_exceedances/precipitation/pr24/0.17_max")
# 
# OUT_012_BIL <- file.path(ROOT, "data/out/05_exceedances/precipitation/pr24/0.12_bilinear")
# OUT_012_MAX <- file.path(ROOT, "data/out/05_exceedances/precipitation/pr24/0.12_max")
# 
# dir.create(OUT_025, recursive = TRUE, showWarnings = FALSE)
# dir.create(OUT_017_BIL, recursive = TRUE, showWarnings = FALSE)
# dir.create(OUT_017_MAX, recursive = TRUE, showWarnings = FALSE)
# dir.create(OUT_012_BIL, recursive = TRUE, showWarnings = FALSE)
# dir.create(OUT_012_MAX, recursive = TRUE, showWarnings = FALSE)
# 
# RP_YEARS <- c(1,2,5,10,25,50,100,200,500,1000)
# 
# WORKERS <- 30L
# OVERWRITE_EMPTY <- TRUE
# MONTHS_KEEP <- 1:12
# 
# # TARGET GRIDS (0.25 is native; these are targets derived from 0.25)
# TARGETS <- list(
#   `0.17` = list(ref = REF_017, out_bil = OUT_017_BIL, out_max = OUT_017_MAX, resolution_deg = 0.17),
#   `0.12` = list(ref = REF_012, out_bil = OUT_012_BIL, out_max = OUT_012_MAX, resolution_deg = 0.12)
# )
# 
# # -----------------------------
# # HELPERS
# # -----------------------------
# .is_empty_qs <- function(path) {
#   if (!file.exists(path)) return(FALSE)
#   fi <- file.info(path)
#   if (!is.finite(fi$size) || fi$size <= 2048) return(TRUE)
#   x <- tryCatch(qread(path), error = function(e) NULL)
#   if (is.null(x)) return(TRUE)
#   if (is.data.frame(x) || data.table::is.data.table(x)) return(nrow(x) == 0 || ncol(x) == 0)
#   FALSE
# }
# 
# .should_skip <- function(path) {
#   if (!file.exists(path)) return(FALSE)
#   if (OVERWRITE_EMPTY && .is_empty_qs(path)) return(FALSE)
#   TRUE
# }
# 
# .read_hourly_day <- function(path) {
#   dt <- as.data.table(qread(path))
#   
#   lon_col <- if ("x" %in% names(dt)) "x" else if ("lon" %in% names(dt)) "lon" else NA_character_
#   lat_col <- if ("y" %in% names(dt)) "y" else if ("lat" %in% names(dt)) "lat" else NA_character_
#   if (is.na(lon_col) || is.na(lat_col)) stop("Missing x/y or lon/lat in: ", path, call. = FALSE)
#   
#   dt_col <- if ("datetime_utc" %in% names(dt)) "datetime_utc" else if ("datetime" %in% names(dt)) "datetime" else NA_character_
#   if (is.na(dt_col)) stop("Missing datetime/datetime_utc in: ", path, call. = FALSE)
#   
#   # precipitation intensity column: try common names
#   val_col <- NA_character_
#   for (cand in c("precipitation","precip","pr","rain","tp","intensity","value","val")) {
#     if (cand %in% names(dt)) { val_col <- cand; break }
#   }
#   if (is.na(val_col)) stop("Missing precipitation value column in: ", path, call. = FALSE)
#   
#   data.table(
#     lon = as.numeric(dt[[lon_col]]),
#     lat = as.numeric(dt[[lat_col]]),
#     dt  = dt[[dt_col]],
#     val = as.numeric(dt[[val_col]])
#   )
# }
# 
# # daily 24h sum per 0.25 cell (sum hourly, per cell)
# .daily_sum_vec_025 <- function(dt_hourly, ref_025, ncell025) {
#   cell <- terra::cellFromXY(ref_025, cbind(dt_hourly$lon, dt_hourly$lat))
#   dt_hourly[, cell := cell]
#   dt_hourly <- dt_hourly[is.finite(val) & !is.na(cell)]
#   out <- rep(NA_real_, ncell025)
#   if (nrow(dt_hourly) == 0L) return(out)
#   
#   s <- dt_hourly[, .(psum = sum(val, na.rm = TRUE)), by = cell]
#   out[s$cell] <- s$psum
#   out
# }
# 
# .max_rp_layer <- function(x, thr_mat) {
#   rp_layer <- integer(length(x))
#   for (k in seq_len(ncol(thr_mat))) {
#     thr <- thr_mat[, k]
#     idx <- which(!is.na(x) & !is.na(thr) & (x >= thr))
#     if (length(idx)) rp_layer[idx] <- k
#   }
#   rp_layer
# }
# 
# .exceed_dt_one_day <- function(p_vec, thr_mat, xy_cells, day_chr) {
#   rp_layer <- .max_rp_layer(p_vec, thr_mat)
#   keep <- which(rp_layer >= 1L)
#   if (!length(keep)) {
#     return(data.table(
#       x = numeric(), y = numeric(), day = character(),
#       pr24_intensity_mm = numeric(), rp_max_years = integer(), rp_max_layer = integer()
#     ))
#   }
#   data.table(
#     x = xy_cells[keep, 1],
#     y = xy_cells[keep, 2],
#     day = day_chr,
#     pr24_intensity_mm = p_vec[keep],
#     rp_max_years = RP_YEARS[rp_layer[keep]],
#     rp_max_layer = rp_layer[keep]
#   )
# }
# 
# # -----------------------------
# # RESAMPLING MAPS (0.25 -> TARGET)
# #   Bilinear (renormalized) and "max of 4 bilinear neighbors"
# #   NOTE: This is the SAME logic pattern as the HI script.
# # -----------------------------
# .build_bilinear_map_target_to_025 <- function(ref_025, ref_tgt) {
#   xy_tgt <- terra::xyFromCell(ref_tgt, seq_len(terra::ncell(ref_tgt)))
#   dx <- terra::res(ref_025)[1]
#   dy <- terra::res(ref_025)[2]
#   e <- terra::ext(ref_025)
#   xmin <- e[1]; xmax <- e[2]; ymin <- e[3]; ymax <- e[4]
#   
#   cf <- (xy_tgt[, 1] - xmin) / dx + 0.5
#   rf <- (ymax - xy_tgt[, 2]) / dy + 0.5
#   
#   c0 <- floor(cf); r0 <- floor(rf)
#   tx <- cf - c0; ty <- rf - r0
#   
#   r_t <- r0
#   r_b <- r0 + 1L
#   c_l <- c0
#   c_r <- c0 + 1L
#   
#   ncol025 <- terra::ncol(ref_025)
#   nrow025 <- terra::nrow(ref_025)
#   inb <- function(r, c) (r >= 1L & r <= nrow025 & c >= 1L & c <= ncol025)
#   
#   tl_ok <- inb(r_t, c_l)
#   tr_ok <- inb(r_t, c_r)
#   bl_ok <- inb(r_b, c_l)
#   br_ok <- inb(r_b, c_r)
#   
#   idx <- matrix(NA_integer_, nrow = nrow(xy_tgt), ncol = 4)
#   idx[tl_ok, 1] <- terra::cellFromRowCol(ref_025, r_t[tl_ok], c_l[tl_ok])
#   idx[tr_ok, 2] <- terra::cellFromRowCol(ref_025, r_t[tr_ok], c_r[tr_ok])
#   idx[bl_ok, 3] <- terra::cellFromRowCol(ref_025, r_b[bl_ok], c_l[bl_ok])
#   idx[br_ok, 4] <- terra::cellFromRowCol(ref_025, r_b[br_ok], c_r[br_ok])
#   
#   w <- cbind(
#     (1 - tx) * (1 - ty),
#     (tx)     * (1 - ty),
#     (1 - tx) * (ty),
#     (tx)     * (ty)
#   )
#   
#   list(xy_tgt = xy_tgt, idx = idx, w = w)
# }
# 
# .apply_bilinear_renorm <- function(src_vec_025, idx_4, w_4) {
#   n <- nrow(idx_4)
#   v <- matrix(NA_real_, nrow = n, ncol = 4)
#   
#   for (j in 1:4) {
#     ii <- idx_4[, j]
#     ok <- !is.na(ii)
#     if (any(ok)) v[ok, j] <- src_vec_025[ii[ok]]
#   }
#   
#   ok <- is.finite(v) & is.finite(w_4)
#   w_use <- w_4
#   w_use[!ok] <- 0
#   
#   num <- rowSums(v * w_use, na.rm = TRUE)
#   den <- rowSums(w_use, na.rm = TRUE)
#   
#   out <- num / den
#   out[den == 0] <- NA_real_
#   out
# }
# 
# .apply_neighbor_max <- function(src_vec_025, idx_4) {
#   n <- nrow(idx_4)
#   v <- matrix(NA_real_, nrow = n, ncol = 4)
#   for (j in 1:4) {
#     ii <- idx_4[, j]
#     ok <- !is.na(ii)
#     if (any(ok)) v[ok, j] <- src_vec_025[ii[ok]]
#   }
#   out <- matrixStats::rowMaxs(v, na.rm = TRUE)
#   out[!is.finite(out)] <- NA_real_
#   out
# }
# 
# # -----------------------------
# # DISCOVER ALL DATES + YEARS
# # -----------------------------
# all_files <- list.files(IN_HOURLY_DIR, pattern = "\\.qs$", full.names = TRUE)
# stopifnot(length(all_files) > 0L)
# 
# bn <- basename(all_files)
# date_chr_all <- sub("\\.qs$", "", bn)
# ok <- grepl("^\\d{4}-\\d{2}-\\d{2}$", date_chr_all)
# all_files <- all_files[ok]
# date_chr_all <- date_chr_all[ok]
# 
# mo <- as.integer(substr(date_chr_all, 6, 7))
# keep_mo <- mo %in% MONTHS_KEEP
# all_files <- all_files[keep_mo]
# date_chr_all <- date_chr_all[keep_mo]
# 
# years_all <- as.integer(substr(date_chr_all, 1, 4))
# years <- sort(unique(years_all))
# 
# date_to_file <- setNames(all_files, date_chr_all)
# 
# # -----------------------------
# # WORKER-LOCAL CACHE (terra pointers built inside worker)
# # -----------------------------
# .WORKER_CACHE <- new.env(parent = emptyenv())
# 
# .get_worker_cache <- function() {
#   if (isTRUE(.WORKER_CACHE$ready)) return(.WORKER_CACHE)
#   
#   ref_025 <- terra::rast(REF_025)
#   rp <- terra::rast(RP_NC)
#   stopifnot(terra::nlyr(rp) == 10L)
#   
#   n025 <- terra::ncell(ref_025)
#   xy_025 <- terra::xyFromCell(ref_025, seq_len(n025))
#   
#   thr_025_df <- terra::extract(rp, xy_025)
#   if ("ID" %in% names(thr_025_df)) thr_025_df$ID <- NULL
#   thr_025 <- as.matrix(thr_025_df)
#   rm(thr_025_df)
#   
#   targets_cache <- list()
#   
#   for (tag in names(TARGETS)) {
#     ref_t <- terra::rast(TARGETS[[tag]]$ref)
#     nt   <- terra::ncell(ref_t)
#     xy_t <- terra::xyFromCell(ref_t, seq_len(nt))
#     
#     bil_map <- .build_bilinear_map_target_to_025(ref_025, ref_t)
#     
#     thr_t_bil <- matrix(NA_real_, nrow = nt, ncol = 10)
#     thr_t_max <- matrix(NA_real_, nrow = nt, ncol = 10)
#     for (k in 1:10) {
#       thr_t_bil[, k] <- .apply_bilinear_renorm(thr_025[, k], bil_map$idx, bil_map$w)
#       thr_t_max[, k] <- .apply_neighbor_max(thr_025[, k], bil_map$idx)
#     }
#     
#     targets_cache[[tag]] <- list(
#       ref = ref_t,
#       n   = nt,
#       xy  = xy_t,
#       bil_map = bil_map,
#       thr_bil = thr_t_bil,
#       thr_max = thr_t_max
#     )
#   }
#   
#   .WORKER_CACHE$ready <- TRUE
#   .WORKER_CACHE$ref_025 <- ref_025
#   .WORKER_CACHE$n025 <- n025
#   .WORKER_CACHE$xy_025 <- xy_025
#   .WORKER_CACHE$thr_025 <- thr_025
#   .WORKER_CACHE$targets <- targets_cache
#   
#   .WORKER_CACHE
# }
# 
# # -----------------------------
# # YEAR PROCESSOR
# # -----------------------------
# .process_year <- function(yy) {
#   cache <- .get_worker_cache()
#   
#   ref_025 <- cache$ref_025
#   n025    <- cache$n025
#   xy_025  <- cache$xy_025
#   thr_025 <- cache$thr_025
#   
#   out_025 <- file.path(OUT_025, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
#   
#   # early skip if all outputs are already done
#   all_skip <- .should_skip(out_025)
#   for (tag in names(TARGETS)) {
#     tinfo <- TARGETS[[tag]]
#     out_bil <- file.path(tinfo$out_bil, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
#     out_max <- file.path(tinfo$out_max, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
#     all_skip <- all_skip && .should_skip(out_bil) && .should_skip(out_max)
#   }
#   if (all_skip) return(invisible(TRUE))
#   
#   d_this <- sort(date_chr_all[years_all == yy])
#   
#   # write empties consistently if no dates for the year
#   if (!length(d_this)) {
#     if (!file.exists(out_025)) qsave(data.table(), out_025, preset = "high")
#     for (tag in names(TARGETS)) {
#       tinfo <- TARGETS[[tag]]
#       out_bil <- file.path(tinfo$out_bil, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
#       out_max <- file.path(tinfo$out_max, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
#       if (!file.exists(out_bil)) qsave(data.table(), out_bil, preset = "high")
#       if (!file.exists(out_max)) qsave(data.table(), out_max, preset = "high")
#     }
#     return(invisible(TRUE))
#   }
#   
#   nd <- length(d_this)
#   
#   # build daily 24h sums on native 0.25 grid
#   P_025 <- matrix(NA_real_, nrow = n025, ncol = nd)
#   for (i in seq_len(nd)) {
#     f <- date_to_file[[d_this[i]]]
#     dt_hourly <- .read_hourly_day(f)
#     P_025[, i] <- .daily_sum_vec_025(dt_hourly, ref_025, n025)
#   }
#   
#   # -----------------------------
#   # NATIVE 0.25 OUTPUT
#   # -----------------------------
#   if (!.should_skip(out_025)) {
#     rows <- vector("list", nd)
#     for (i in seq_len(nd)) {
#       d <- .exceed_dt_one_day(P_025[, i], thr_025, xy_025, d_this[i])
#       if (nrow(d)) d[, `:=`(resolution_deg = 0.25, resample_method = "native")]
#       rows[[i]] <- d
#     }
#     out_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
#     setorder(out_dt, day, y, x)
#     qsave(out_dt, out_025, preset = "high")
#   }
#   
#   # -----------------------------
#   # TARGET GRIDS OUTPUTS (bilinear + max)
#   # -----------------------------
#   for (tag in names(TARGETS)) {
#     tinfo <- TARGETS[[tag]]
#     tcache <- cache$targets[[tag]]
#     
#     out_bil <- file.path(tinfo$out_bil, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
#     out_max <- file.path(tinfo$out_max, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
#     if (.should_skip(out_bil) && .should_skip(out_max)) next
#     
#     nt <- tcache$n
#     xy_t <- tcache$xy
#     bil_map <- tcache$bil_map
#     thr_bil <- tcache$thr_bil
#     thr_max <- tcache$thr_max
#     
#     P_t_bil <- matrix(NA_real_, nrow = nt, ncol = nd)
#     P_t_max <- matrix(NA_real_, nrow = nt, ncol = nd)
#     
#     for (i in seq_len(nd)) {
#       v <- P_025[, i]
#       P_t_bil[, i] <- .apply_bilinear_renorm(v, bil_map$idx, bil_map$w)
#       P_t_max[, i] <- .apply_neighbor_max(v, bil_map$idx)
#     }
#     
#     if (!.should_skip(out_bil)) {
#       rows <- vector("list", nd)
#       for (i in seq_len(nd)) {
#         d <- .exceed_dt_one_day(P_t_bil[, i], thr_bil, xy_t, d_this[i])
#         if (nrow(d)) d[, `:=`(resolution_deg = tinfo$resolution_deg, resample_method = "bilinear")]
#         rows[[i]] <- d
#       }
#       out_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
#       setorder(out_dt, day, y, x)
#       qsave(out_dt, out_bil, preset = "high")
#     }
#     
#     if (!.should_skip(out_max)) {
#       rows <- vector("list", nd)
#       for (i in seq_len(nd)) {
#         d <- .exceed_dt_one_day(P_t_max[, i], thr_max, xy_t, d_this[i])
#         if (nrow(d)) d[, `:=`(resolution_deg = tinfo$resolution_deg, resample_method = "max")]
#         rows[[i]] <- d
#       }
#       out_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
#       setorder(out_dt, day, y, x)
#       qsave(out_dt, out_max, preset = "high")
#     }
#   }
#   
#   invisible(TRUE)
# }
# 
# # -----------------------------
# # RUN
# # -----------------------------
# plan(multisession, workers = WORKERS)
# invisible(future_lapply(years, .process_year))
# plan(sequential)

# Southeast grid ----
# scripts/00_build_ref_grids_from_stm_precipitation_southeast.R
# =============================================================================
# Build precipitation reference grids from STM1/STM2 ALL-MONTH MEDIAN values
# (NO hard-coded 0.17/0.12; resolutions come from STM diagnostics CSV)
#
# Inputs (ONLY):
#   1) A CSV with columns: variable, domain, day, stm1, stm2  (at minimum)
#   2) A reference GeoTIFF mask/grid (land mask) defining CRS + extent (+ NA mask)
#
# Outputs:
#   data/out/03_cov_model/precipitation/
#     precipitation_0.25/
#       ref_grid_precipitation_0.25.tif
#       ref_grid_cells_precipitation_0.25.qs
#     precipitation_<R1>_bilinear/
#       ref_grid_precipitation_<R1>_bilinear.tif
#       ref_grid_cells_precipitation_<R1>_bilinear.qs
#     precipitation_<R1>_max/
#       ref_grid_precipitation_<R1>_max.tif
#       ref_grid_cells_precipitation_<R1>_max.qs
#     precipitation_<R2>_bilinear/
#     precipitation_<R2>_max/
#
# Also writes:
#   ref_map_by_resolution_precipitation.qs
#   dataset_to_refgrid_precipitation.qs
# =============================================================================

suppressPackageStartupMessages({
  library(fs)
  library(data.table)
  library(lubridate)
  library(qs)
  library(terra)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

# STM diagnostics CSV (SOUTHEAST)
PR_COV_CSV <- file.path(
  ROOT,
  "data/out/03_cov_model/precipitation/precipitation_daily_optimal_cov_model_southeast.csv"
)

# Land mask GeoTIFF from experimental covariance scripts (lon/lat)
MASK_TIF <- file.path(ROOT, "data/out/01_vars/precipitation/_mask_target.tif")

# Where to write ref grids
REF_DIR <- file.path(ROOT, "data/out/03_cov_model/precipitation")
dir_create(REF_DIR, recurse = TRUE)

DOMAIN_NAME <- "southeast"
MONTHS_KEEP <- 1:12         # ALL months for precipitation
BASE_RES_DEG <- 0.25
RES_TAG_DIGITS <- 2
OVERWRITE_REF <- FALSE

.log <- function(fmt, ...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))

stopifnot(file_exists(PR_COV_CSV), file_exists(MASK_TIF))

# -----------------------------
# Helpers
# -----------------------------
.fmt_res_tag <- function(res_deg, digits = 2) sprintf(paste0("%.", digits, "f"), res_deg)

.parse_day <- function(x) {
  d <- as.IDate(x, format = "%Y-%m-%d")
  if (all(is.na(d))) d <- as.IDate(x, format = "%m/%d/%Y")
  if (all(is.na(d))) d <- as.IDate(x, format = "%m/%d/%y")
  d
}

.get_stm_resolutions <- function(cov_csv, domain, months_keep, base_res_deg) {
  dt <- fread(cov_csv)
  
  req <- c("variable", "domain", "day", "stm1", "stm2")
  miss <- setdiff(req, names(dt))
  if (length(miss)) stop("Missing columns in cov CSV: ", paste(miss, collapse = ", "))
  
  dt <- dt[variable == "precipitation" & domain == domain]
  if (!nrow(dt)) stop("No rows for variable=precipitation and domain=", domain)
  
  dt[, day_date := .parse_day(day)]
  if (all(is.na(dt$day_date))) stop("Could not parse 'day' in: ", cov_csv)
  
  dt <- dt[month(day_date) %in% months_keep]
  if (!nrow(dt)) stop("No rows after filtering to months=", paste(months_keep, collapse=","))
  
  stm1_med <- median(as.numeric(dt$stm1), na.rm = TRUE)
  stm2_med <- median(as.numeric(dt$stm2), na.rm = TRUE)
  
  if (!is.finite(stm1_med) || !is.finite(stm2_med)) {
    stop("Non-finite STM medians; check stm1/stm2 values.")
  }
  
  res1 <- base_res_deg * stm1_med
  res2 <- base_res_deg * stm2_med
  
  if (!is.finite(res1) || !is.finite(res2) || res1 <= 0 || res2 <= 0) {
    stop("Derived resolutions invalid: res1=", res1, " res2=", res2)
  }
  
  list(stm1_median = stm1_med, stm2_median = stm2_med, res_deg = c(res1, res2))
}

.build_ref_from_mask <- function(mask_tif, res_deg, out_tif, overwrite = FALSE) {
  if (file_exists(out_tif) && !overwrite) return(out_tif)
  
  r_mask <- terra::rast(mask_tif)[[1]]
  
  r_tmpl <- terra::rast(ext = terra::ext(r_mask), resolution = res_deg, crs = terra::crs(r_mask))
  try(terra::origin(r_tmpl) <- terra::origin(r_mask), silent = TRUE)
  
  # nearest-neighbor preserves NA mask pattern
  r_ref <- terra::resample(r_mask, r_tmpl, method = "near")
  
  dir_create(path_dir(out_tif), recurse = TRUE)
  terra::writeRaster(r_ref, out_tif, overwrite = TRUE)
  
  out_tif
}

.write_ref_cells_qs <- function(ref_tif, out_qs, overwrite = FALSE) {
  if (file_exists(out_qs) && !overwrite) return(out_qs)
  
  r <- terra::rast(ref_tif)[[1]]
  df <- as.data.frame(r, xy = TRUE, cells = TRUE, na.rm = FALSE)
  dt <- as.data.table(df)
  
  dt_out <- dt[, .(
    cell_id = as.integer(cell),
    x = as.numeric(x),
    y = as.numeric(y)
  )]
  
  dir_create(path_dir(out_qs), recurse = TRUE)
  qs::qsave(dt_out, out_qs, preset = "high")
  out_qs
}

# -----------------------------
# Compute STM-based resolutions
# -----------------------------
stm <- .get_stm_resolutions(PR_COV_CSV, DOMAIN_NAME, MONTHS_KEEP, BASE_RES_DEG)

res_raw <- unique(stm$res_deg)
res_tag <- vapply(res_raw, .fmt_res_tag, character(1), digits = RES_TAG_DIGITS)

.log("STM medians (all months) | stm1=%.6f | stm2=%.6f", stm$stm1_median, stm$stm2_median)
.log("Derived resolutions (deg) | raw: %s | tags: %s",
     paste(sprintf("%.6f", res_raw), collapse = ", "),
     paste(res_tag, collapse = ", "))

needed <- data.table(
  res_deg = c(BASE_RES_DEG, res_raw),
  res_tag = c(.fmt_res_tag(BASE_RES_DEG, RES_TAG_DIGITS), res_tag)
)
needed <- unique(needed, by = "res_tag")
setorder(needed, res_deg)

.log("Unique resolutions to build (n=%d): %s", nrow(needed), paste(needed$res_tag, collapse = ", "))

# -----------------------------
# Build ref grids
# -----------------------------
ref_map <- rbindlist(lapply(seq_len(nrow(needed)), function(i) {
  res_deg <- needed$res_deg[i]
  tag     <- needed$res_tag[i]
  
  if (abs(res_deg - BASE_RES_DEG) < 1e-12) {
    ref_grid_tag <- paste0("precipitation_", tag)
    ref_tif <- file.path(REF_DIR, ref_grid_tag, paste0("ref_grid_", ref_grid_tag, ".tif"))
    ref_qs  <- file.path(REF_DIR, ref_grid_tag, paste0("ref_grid_cells_", ref_grid_tag, ".qs"))
    
    .build_ref_from_mask(MASK_TIF, res_deg, ref_tif, overwrite = OVERWRITE_REF)
    .write_ref_cells_qs(ref_tif, ref_qs, overwrite = OVERWRITE_REF)
    
    return(data.table(
      dataset_tag   = ref_grid_tag,
      ref_grid_tag  = ref_grid_tag,
      res_deg       = res_deg,
      ref_tif       = ref_tif,
      ref_cells_qs  = ref_qs
    ))
  } else {
    out <- rbindlist(lapply(c("bilinear", "max"), function(method) {
      ref_grid_tag <- paste0("precipitation_", tag, "_", method)
      ref_tif <- file.path(REF_DIR, ref_grid_tag, paste0("ref_grid_", ref_grid_tag, ".tif"))
      ref_qs  <- file.path(REF_DIR, ref_grid_tag, paste0("ref_grid_cells_", ref_grid_tag, ".qs"))
      
      .build_ref_from_mask(MASK_TIF, res_deg, ref_tif, overwrite = OVERWRITE_REF)
      .write_ref_cells_qs(ref_tif, ref_qs, overwrite = OVERWRITE_REF)
      
      data.table(
        dataset_tag   = ref_grid_tag,
        ref_grid_tag  = ref_grid_tag,
        res_deg       = res_deg,
        ref_tif       = ref_tif,
        ref_cells_qs  = ref_qs
      )
    }))
    return(out)
  }
}), use.names = TRUE, fill = TRUE)

setkey(ref_map, dataset_tag)

dataset_map <- ref_map[, .(
  dataset_tag,
  ref_grid_tag,
  res_deg,
  ref_tif,
  ref_cells_qs
)]
setkey(dataset_map, dataset_tag)

# -----------------------------
# Write mappings
# -----------------------------
OUT_REF_MAP     <- file.path(REF_DIR, "ref_map_by_resolution_precipitation.qs")
OUT_DATASET_MAP <- file.path(REF_DIR, "dataset_to_refgrid_precipitation.qs")

qsave(ref_map,     OUT_REF_MAP, preset = "high")
qsave(dataset_map, OUT_DATASET_MAP, preset = "high")

.log("Wrote: %s (rows=%d)", OUT_REF_MAP, nrow(ref_map))
.log("Wrote: %s (rows=%d)", OUT_DATASET_MAP, nrow(dataset_map))
.log("REF GRID OUTPUTS:")
print(dataset_map[])

# Southeast exceedances ----
# scripts/05_pr24_exceedances_native_vs_resampled_MAXRP_reanalysis_all_hourly_v1.R
# ------------------------------------------------------------------------------
# PR24 exceedances (24-hr daily summed precipitation) vs PR24 return periods
#
# Updates:
#   1) NO hard-coded target resolutions (no fixed 0.17 / 0.12 assumptions)
#   2) Targets are read from:
#        data/out/03_cov_model/precipitation/ref_grids/dataset_to_refgrid_precipitation.qs
#      created by the STM-driven precipitation ref-grid builder
#   3) Native remains 0.25
#
# Inputs:
#   - Hourly precip daily files (qs): data/out/01_vars/precipitation/YYYY-MM-DD.qs
#   - Return period raster (nc):      data/out/04_rp/precip_idf_pr24_return_periods.nc
#   - Ref-grid map (qs):              data/out/03_cov_model/precipitation/ref_grids/dataset_to_refgrid_precipitation.qs
#
# Outputs (year-chunked .qs):
#   - Native (0.25°):                 data/out/05_exceedances/precipitation/pr24/0.25
#   - Resampled (targets):            data/out/05_exceedances/precipitation/pr24/<res>_bilinear
#                                     data/out/05_exceedances/precipitation/pr24/<res>_max
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(data.table)
  library(qs)
  library(matrixStats)
  library(future)
  library(future.apply)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

IN_HOURLY_DIR <- file.path(ROOT, "data/out/01_vars/precipitation")
RP_NC         <- file.path(ROOT, "data/out/04_rp/precip_idf_pr24_return_periods_bbox.nc")

# STM-driven ref-grid mapping
REF_MAP_QS <- file.path(ROOT, "data/out/03_cov_model/precipitation/dataset_to_refgrid_precipitation.qs")

# keep outputs compatible with legacy structure
OUT_BASE <- file.path(ROOT, "data/out/05_exceedances/precipitation/pr24")
OUT_025  <- file.path(OUT_BASE, "0.25")
dir.create(OUT_025, recursive = TRUE, showWarnings = FALSE)

RP_YEARS <- c(2,5,10,25,50,100,200,500,1000)

WORKERS <- 30L
OVERWRITE_EMPTY <- TRUE
MONTHS_KEEP <- 1:12

# -----------------------------
# HELPERS
# -----------------------------
.is_empty_qs <- function(path) {
  if (!file.exists(path)) return(FALSE)
  fi <- file.info(path)
  if (!is.finite(fi$size) || fi$size <= 2048) return(TRUE)
  x <- tryCatch(qread(path), error = function(e) NULL)
  if (is.null(x)) return(TRUE)
  if (is.data.frame(x) || data.table::is.data.table(x)) return(nrow(x) == 0 || ncol(x) == 0)
  FALSE
}

.should_skip <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (OVERWRITE_EMPTY && .is_empty_qs(path)) return(FALSE)
  TRUE
}

.read_hourly_day <- function(path) {
  dt <- as.data.table(qread(path))
  
  lon_col <- if ("x" %in% names(dt)) "x" else if ("lon" %in% names(dt)) "lon" else NA_character_
  lat_col <- if ("y" %in% names(dt)) "y" else if ("lat" %in% names(dt)) "lat" else NA_character_
  if (is.na(lon_col) || is.na(lat_col)) stop("Missing x/y or lon/lat in: ", path, call. = FALSE)
  
  dt_col <- if ("datetime_utc" %in% names(dt)) "datetime_utc" else if ("datetime" %in% names(dt)) "datetime" else NA_character_
  if (is.na(dt_col)) stop("Missing datetime/datetime_utc in: ", path, call. = FALSE)
  
  val_col <- NA_character_
  for (cand in c("precipitation","precip","pr","rain","tp","intensity","value","val")) {
    if (cand %in% names(dt)) { val_col <- cand; break }
  }
  if (is.na(val_col)) stop("Missing precipitation value column in: ", path, call. = FALSE)
  
  data.table(
    lon = as.numeric(dt[[lon_col]]),
    lat = as.numeric(dt[[lat_col]]),
    dt  = dt[[dt_col]],
    val = as.numeric(dt[[val_col]])
  )
}

.daily_sum_vec_025 <- function(dt_hourly, ref_025, ncell025) {
  cell <- terra::cellFromXY(ref_025, cbind(dt_hourly$lon, dt_hourly$lat))
  dt_hourly[, cell := cell]
  dt_hourly <- dt_hourly[is.finite(val) & !is.na(cell)]
  out <- rep(NA_real_, ncell025)
  if (nrow(dt_hourly) == 0L) return(out)
  
  s <- dt_hourly[, .(psum = sum(val, na.rm = TRUE)), by = cell]
  out[s$cell] <- s$psum
  out
}

.max_rp_layer <- function(x, thr_mat) {
  rp_layer <- integer(length(x))
  for (k in seq_len(ncol(thr_mat))) {
    thr <- thr_mat[, k]
    idx <- which(!is.na(x) & !is.na(thr) & (x >= thr))
    if (length(idx)) rp_layer[idx] <- k
  }
  rp_layer
}

.exceed_dt_one_day <- function(p_vec, thr_mat, xy_cells, day_chr) {
  rp_layer <- .max_rp_layer(p_vec, thr_mat)
  keep <- which(rp_layer >= 1L)
  if (!length(keep)) {
    return(data.table(
      x = numeric(), y = numeric(), day = character(),
      pr24_intensity_mm = numeric(), rp_max_years = integer(), rp_max_layer = integer()
    ))
  }
  data.table(
    x = xy_cells[keep, 1],
    y = xy_cells[keep, 2],
    day = day_chr,
    pr24_intensity_mm = p_vec[keep],
    rp_max_years = RP_YEARS[rp_layer[keep]],
    rp_max_layer = rp_layer[keep]
  )
}

# -----------------------------
# RESAMPLING MAPS (0.25 -> TARGET)
#   Bilinear (renormalized) and "max of 4 bilinear neighbors"
# -----------------------------
.build_bilinear_map_target_to_025 <- function(ref_025, ref_tgt) {
  xy_tgt <- terra::xyFromCell(ref_tgt, seq_len(terra::ncell(ref_tgt)))
  dx <- terra::res(ref_025)[1]
  dy <- terra::res(ref_025)[2]
  e <- terra::ext(ref_025)
  xmin <- e[1]; xmax <- e[2]; ymin <- e[3]; ymax <- e[4]
  
  cf <- (xy_tgt[, 1] - xmin) / dx + 0.5
  rf <- (ymax - xy_tgt[, 2]) / dy + 0.5
  
  c0 <- floor(cf); r0 <- floor(rf)
  tx <- cf - c0; ty <- rf - r0
  
  r_t <- r0
  r_b <- r0 + 1L
  c_l <- c0
  c_r <- c0 + 1L
  
  ncol025 <- terra::ncol(ref_025)
  nrow025 <- terra::nrow(ref_025)
  inb <- function(r, c) (r >= 1L & r <= nrow025 & c >= 1L & c <= ncol025)
  
  tl_ok <- inb(r_t, c_l)
  tr_ok <- inb(r_t, c_r)
  bl_ok <- inb(r_b, c_l)
  br_ok <- inb(r_b, c_r)
  
  idx <- matrix(NA_integer_, nrow = nrow(xy_tgt), ncol = 4)
  idx[tl_ok, 1] <- terra::cellFromRowCol(ref_025, r_t[tl_ok], c_l[tl_ok])
  idx[tr_ok, 2] <- terra::cellFromRowCol(ref_025, r_t[tr_ok], c_r[tr_ok])
  idx[bl_ok, 3] <- terra::cellFromRowCol(ref_025, r_b[bl_ok], c_l[bl_ok])
  idx[br_ok, 4] <- terra::cellFromRowCol(ref_025, r_b[br_ok], c_r[br_ok])
  
  w <- cbind(
    (1 - tx) * (1 - ty),
    (tx)     * (1 - ty),
    (1 - tx) * (ty),
    (tx)     * (ty)
  )
  
  list(xy_tgt = xy_tgt, idx = idx, w = w)
}

.apply_bilinear_renorm <- function(src_vec_025, idx_4, w_4) {
  n <- nrow(idx_4)
  v <- matrix(NA_real_, nrow = n, ncol = 4)
  
  for (j in 1:4) {
    ii <- idx_4[, j]
    ok <- !is.na(ii)
    if (any(ok)) v[ok, j] <- src_vec_025[ii[ok]]
  }
  
  ok <- is.finite(v) & is.finite(w_4)
  w_use <- w_4
  w_use[!ok] <- 0
  
  num <- rowSums(v * w_use, na.rm = TRUE)
  den <- rowSums(w_use, na.rm = TRUE)
  
  out <- num / den
  out[den == 0] <- NA_real_
  out
}

.apply_neighbor_max <- function(src_vec_025, idx_4) {
  n <- nrow(idx_4)
  v <- matrix(NA_real_, nrow = n, ncol = 4)
  for (j in 1:4) {
    ii <- idx_4[, j]
    ok <- !is.na(ii)
    if (any(ok)) v[ok, j] <- src_vec_025[ii[ok]]
  }
  out <- matrixStats::rowMaxs(v, na.rm = TRUE)
  out[!is.finite(out)] <- NA_real_
  out
}

# -----------------------------
# LOAD TARGETS from REF_MAP_QS
# -----------------------------
stopifnot(file.exists(RP_NC), file.exists(REF_MAP_QS))

ref_map <- as.data.table(qread(REF_MAP_QS))

if (!all(c("dataset_tag","ref_tif","res_deg") %in% names(ref_map))) {
  stop("REF_MAP_QS missing required columns (dataset_tag, ref_tif, res_deg): ", REF_MAP_QS)
}

# Native 0.25 ref
row_025 <- ref_map[dataset_tag == "precipitation_0.25"]
if (nrow(row_025) != 1L) stop("Could not find unique 'precipitation_0.25' in: ", REF_MAP_QS)
REF_025 <- row_025$ref_tif[1]

# Derived targets
targets_dt <- ref_map[dataset_tag != "precipitation_0.25"]
targets_dt <- targets_dt[grepl("^precipitation_[0-9]+\\.[0-9]{2}_(bilinear|max)$", dataset_tag)]

if (!nrow(targets_dt)) stop("No derived precipitation targets found in REF_MAP_QS: ", REF_MAP_QS)

targets_dt[, res_tag := sub("^precipitation_([0-9]+\\.[0-9]{2})_.*$", "\\1", dataset_tag)]
targets_dt[, method  := sub("^precipitation_[0-9]+\\.[0-9]{2}_(bilinear|max)$", "\\1", dataset_tag)]

ok_res <- targets_dt[, .(has_bil = any(method == "bilinear"), has_max = any(method == "max")), by = res_tag][has_bil & has_max, res_tag]
targets_dt <- targets_dt[res_tag %in% ok_res]

TARGETS <- lapply(ok_res, function(rt) {
  bil <- targets_dt[res_tag == rt & method == "bilinear"][1]
  mx  <- targets_dt[res_tag == rt & method == "max"][1]
  
  list(
    res_tag = rt,
    resolution_deg = as.numeric(bil$res_deg),
    ref_bil = bil$ref_tif,
    ref_max = mx$ref_tif,
    out_bil = file.path(OUT_BASE, paste0(rt, "_bilinear")),
    out_max = file.path(OUT_BASE, paste0(rt, "_max"))
  )
})
names(TARGETS) <- ok_res

for (rt in names(TARGETS)) {
  dir.create(TARGETS[[rt]]$out_bil, recursive = TRUE, showWarnings = FALSE)
  dir.create(TARGETS[[rt]]$out_max, recursive = TRUE, showWarnings = FALSE)
}

message(sprintf("[TARGETS] Loaded %d derived resolution(s): %s", length(TARGETS), paste(names(TARGETS), collapse = ", ")))
message("[REF_025] ", REF_025)
message("[RP_NC]   ", RP_NC)

# -----------------------------
# DISCOVER ALL DATES + YEARS
# -----------------------------
all_files <- list.files(IN_HOURLY_DIR, pattern = "\\.qs$", full.names = TRUE)
stopifnot(length(all_files) > 0L)

bn <- basename(all_files)
date_chr_all <- sub("\\.qs$", "", bn)
ok <- grepl("^\\d{4}-\\d{2}-\\d{2}$", date_chr_all)
all_files <- all_files[ok]
date_chr_all <- date_chr_all[ok]

mo <- as.integer(substr(date_chr_all, 6, 7))
keep_mo <- mo %in% MONTHS_KEEP
all_files <- all_files[keep_mo]
date_chr_all <- date_chr_all[keep_mo]

years_all <- as.integer(substr(date_chr_all, 1, 4))
years <- sort(unique(years_all))

date_to_file <- setNames(all_files, date_chr_all)

# -----------------------------
# WORKER-LOCAL CACHE
# -----------------------------
.WORKER_CACHE <- new.env(parent = emptyenv())

.get_worker_cache <- function() {
  if (isTRUE(.WORKER_CACHE$ready)) return(.WORKER_CACHE)
  
  ref_025 <- terra::rast(REF_025)
  rp <- terra::rast(RP_NC)
  stopifnot(terra::nlyr(rp) == 10L)
  
  n025 <- terra::ncell(ref_025)
  xy_025 <- terra::xyFromCell(ref_025, seq_len(n025))
  
  thr_025_df <- terra::extract(rp, xy_025)
  if ("ID" %in% names(thr_025_df)) thr_025_df$ID <- NULL
  thr_025 <- as.matrix(thr_025_df)
  rm(thr_025_df)
  
  targets_cache <- list()
  
  for (rt in names(TARGETS)) {
    tinfo <- TARGETS[[rt]]
    
    ref_bil <- terra::rast(tinfo$ref_bil)
    ref_max <- terra::rast(tinfo$ref_max)
    
    nt_bil <- terra::ncell(ref_bil)
    nt_max <- terra::ncell(ref_max)
    if (nt_bil != nt_max) stop("Target bilinear/max grids differ in cell count for res_tag=", rt)
    
    xy_t <- terra::xyFromCell(ref_bil, seq_len(nt_bil))
    bil_map <- .build_bilinear_map_target_to_025(ref_025, ref_bil)
    
    thr_t_bil <- matrix(NA_real_, nrow = nt_bil, ncol = 10)
    thr_t_max <- matrix(NA_real_, nrow = nt_bil, ncol = 10)
    for (k in 1:10) {
      thr_t_bil[, k] <- .apply_bilinear_renorm(thr_025[, k], bil_map$idx, bil_map$w)
      thr_t_max[, k] <- .apply_neighbor_max(thr_025[, k], bil_map$idx)
    }
    
    targets_cache[[rt]] <- list(
      n = nt_bil,
      xy = xy_t,
      bil_map = bil_map,
      thr_bil = thr_t_bil,
      thr_max = thr_t_max,
      resolution_deg = tinfo$resolution_deg
    )
  }
  
  .WORKER_CACHE$ready <- TRUE
  .WORKER_CACHE$ref_025 <- ref_025
  .WORKER_CACHE$n025 <- n025
  .WORKER_CACHE$xy_025 <- xy_025
  .WORKER_CACHE$thr_025 <- thr_025
  .WORKER_CACHE$targets <- targets_cache
  
  .WORKER_CACHE
}

# -----------------------------
# YEAR PROCESSOR
# -----------------------------
.process_year <- function(yy) {
  cache <- .get_worker_cache()
  
  ref_025 <- cache$ref_025
  n025    <- cache$n025
  xy_025  <- cache$xy_025
  thr_025 <- cache$thr_025
  
  out_025 <- file.path(OUT_025, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
  
  all_skip <- .should_skip(out_025)
  for (rt in names(TARGETS)) {
    tinfo <- TARGETS[[rt]]
    out_bil <- file.path(tinfo$out_bil, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
    out_max <- file.path(tinfo$out_max, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
    all_skip <- all_skip && .should_skip(out_bil) && .should_skip(out_max)
  }
  if (all_skip) return(invisible(TRUE))
  
  d_this <- sort(date_chr_all[years_all == yy])
  
  if (!length(d_this)) {
    if (!file.exists(out_025)) qsave(data.table(), out_025, preset = "high")
    for (rt in names(TARGETS)) {
      tinfo <- TARGETS[[rt]]
      out_bil <- file.path(tinfo$out_bil, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
      out_max <- file.path(tinfo$out_max, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
      if (!file.exists(out_bil)) qsave(data.table(), out_bil, preset = "high")
      if (!file.exists(out_max)) qsave(data.table(), out_max, preset = "high")
    }
    return(invisible(TRUE))
  }
  
  nd <- length(d_this)
  
  # Build daily 24h sums on native 0.25 grid
  P_025 <- matrix(NA_real_, nrow = n025, ncol = nd)
  for (i in seq_len(nd)) {
    f <- date_to_file[[d_this[i]]]
    dt_hourly <- .read_hourly_day(f)
    P_025[, i] <- .daily_sum_vec_025(dt_hourly, ref_025, n025)
  }
  
  # -----------------------------
  # NATIVE 0.25 OUTPUT
  # -----------------------------
  if (!.should_skip(out_025)) {
    rows <- vector("list", nd)
    for (i in seq_len(nd)) {
      d <- .exceed_dt_one_day(P_025[, i], thr_025, xy_025, d_this[i])
      if (nrow(d)) d[, `:=`(resolution_deg = 0.25, resample_method = "native")]
      rows[[i]] <- d
    }
    out_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
    setorder(out_dt, day, y, x)
    qsave(out_dt, out_025, preset = "high")
  }
  
  # -----------------------------
  # TARGET GRIDS OUTPUTS
  # -----------------------------
  for (rt in names(TARGETS)) {
    tinfo <- TARGETS[[rt]]
    tcache <- cache$targets[[rt]]
    
    out_bil <- file.path(tinfo$out_bil, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
    out_max <- file.path(tinfo$out_max, sprintf("pr24_exceedances_maxrp_%d.qs", yy))
    if (.should_skip(out_bil) && .should_skip(out_max)) next
    
    nt <- tcache$n
    xy_t <- tcache$xy
    bil_map <- tcache$bil_map
    thr_bil <- tcache$thr_bil
    thr_max <- tcache$thr_max
    res_deg <- tcache$resolution_deg
    
    P_t_bil <- matrix(NA_real_, nrow = nt, ncol = nd)
    P_t_max <- matrix(NA_real_, nrow = nt, ncol = nd)
    
    for (i in seq_len(nd)) {
      v <- P_025[, i]
      P_t_bil[, i] <- .apply_bilinear_renorm(v, bil_map$idx, bil_map$w)
      P_t_max[, i] <- .apply_neighbor_max(v, bil_map$idx)
    }
    
    if (!.should_skip(out_bil)) {
      rows <- vector("list", nd)
      for (i in seq_len(nd)) {
        d <- .exceed_dt_one_day(P_t_bil[, i], thr_bil, xy_t, d_this[i])
        if (nrow(d)) d[, `:=`(resolution_deg = res_deg, resample_method = "bilinear")]
        rows[[i]] <- d
      }
      out_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
      setorder(out_dt, day, y, x)
      qsave(out_dt, out_bil, preset = "high")
    }
    
    if (!.should_skip(out_max)) {
      rows <- vector("list", nd)
      for (i in seq_len(nd)) {
        d <- .exceed_dt_one_day(P_t_max[, i], thr_max, xy_t, d_this[i])
        if (nrow(d)) d[, `:=`(resolution_deg = res_deg, resample_method = "max")]
        rows[[i]] <- d
      }
      out_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
      setorder(out_dt, day, y, x)
      qsave(out_dt, out_max, preset = "high")
    }
  }
  
  invisible(TRUE)
}

# -----------------------------
# RUN
# -----------------------------
plan(multisession, workers = WORKERS)
invisible(future_lapply(years, .process_year))
plan(sequential)
