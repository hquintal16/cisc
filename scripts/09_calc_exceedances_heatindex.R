# # scripts/05_hi3_exceedances_native_vs_resampled_MAXRP_reanalysis_all_hourly_v2.R
# # ------------------------------------------------------------------------------
# # HI3 exceedances (3-day rolling max of 3-day mean of DAILY MAX HI) vs HI3 return periods
# #
# # Inputs:
# #   - Hourly HI daily files (qs):  data/out/01_vars/heatindex/YYYY-MM-DD.qs
# #   - Return period raster (nc):   data/out/04_rp/heatindex_idf_hi3_return_periods.nc
# #   - Reference grids (tif):       data/out/09_performance/ref_grids/<tag>/ref_grid_<tag>.tif
# #
# # Outputs (year-chunked .qs):
# #   - Native (0.25°):                 data/out/05_exceedances/heatindex/0.25
# #   - Resampled (target; bilinear):   data/out/05_exceedances/heatindex/<target>_bilinear
# #   - Resampled (target; neigh-max):  data/out/05_exceedances/heatindex/<target>_max
# #
# # Targets included:
# #   - 0.44° (bilinear + max)
# #   - 0.34° (bilinear + max)
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
# IN_HOURLY_DIR <- file.path(ROOT, "data/out/01_vars/heatindex")
# RP_NC         <- file.path(ROOT, "data/out/04_rp/heatindex_idf_hi3_return_periods.nc")
# 
# REF_025 <- file.path(ROOT, "data/out/09_performance/ref_grids/heatindex_0.25/ref_grid_heatindex_0.25.tif")
# REF_044 <- file.path(ROOT, "data/out/09_performance/ref_grids/heatindex_0.44_bilinear/ref_grid_heatindex_0.44_bilinear.tif")
# REF_034 <- file.path(ROOT, "data/out/09_performance/ref_grids/heatindex_0.34_bilinear/ref_grid_heatindex_0.34_bilinear.tif")
# 
# # keep outputs compatible with legacy layout
# OUT_025     <- file.path(ROOT, "data/out/05_exceedances/heatindex/0.25")
# OUT_044_BIL <- file.path(ROOT, "data/out/05_exceedances/heatindex/0.44_bilinear")
# OUT_044_MAX <- file.path(ROOT, "data/out/05_exceedances/heatindex/0.44_max")
# OUT_034_BIL <- file.path(ROOT, "data/out/05_exceedances/heatindex/0.34_bilinear")
# OUT_034_MAX <- file.path(ROOT, "data/out/05_exceedances/heatindex/0.34_max")
# 
# dir.create(OUT_025,     recursive = TRUE, showWarnings = FALSE)
# dir.create(OUT_044_BIL, recursive = TRUE, showWarnings = FALSE)
# dir.create(OUT_044_MAX, recursive = TRUE, showWarnings = FALSE)
# dir.create(OUT_034_BIL, recursive = TRUE, showWarnings = FALSE)
# dir.create(OUT_034_MAX, recursive = TRUE, showWarnings = FALSE)
# 
# RP_YEARS <- c(1,2,5,10,25,50,100,200,500,1000)
# 
# WORKERS <- 30L
# OVERWRITE_EMPTY <- TRUE
# MONTHS_KEEP <- 1:12
# 
# # TARGET GRIDS: add more here if needed
# TARGETS <- list(
#   `0.44` = list(ref = REF_044, out_bil = OUT_044_BIL, out_max = OUT_044_MAX, resolution_deg = 0.44),
#   `0.34` = list(ref = REF_034, out_bil = OUT_034_BIL, out_max = OUT_034_MAX, resolution_deg = 0.34)
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
#   val_col <- if ("heat_index" %in% names(dt)) "heat_index" else if ("intensity" %in% names(dt)) "intensity" else NA_character_
#   if (is.na(val_col)) stop("Missing heat_index/intensity in: ", path, call. = FALSE)
#   
#   dt[, .(
#     lon = as.numeric(get(lon_col)),
#     lat = as.numeric(get(lat_col)),
#     dt  = get(dt_col),
#     val = as.numeric(get(val_col))
#   )]
# }
# 
# .daily_max_vec_025 <- function(dt_hourly, ref_025, ncell025) {
#   cell <- terra::cellFromXY(ref_025, cbind(dt_hourly$lon, dt_hourly$lat))
#   dt_hourly[, cell := cell]
#   dt_hourly <- dt_hourly[is.finite(val) & !is.na(cell)]
#   out <- rep(NA_real_, ncell025)
#   if (nrow(dt_hourly) == 0L) return(out)
#   m <- dt_hourly[, .(dmax = max(val, na.rm = TRUE)), by = cell]
#   out[m$cell] <- m$dmax
#   out
# }
# 
# .hi3_from_dailymax <- function(H) {
#   H <- as.matrix(H)
#   storage.mode(H) <- "double"
#   ncell <- nrow(H)
#   nd <- ncol(H)
#   hi3 <- matrix(NA_real_, nrow = ncell, ncol = nd)
#   if (nd == 0) return(hi3)
#   
#   m_end <- matrix(NA_real_, nrow = ncell, ncol = nd)
#   if (nd >= 3) {
#     for (t in 3:nd) {
#       m_end[, t] <- matrixStats::rowMeans2(H, cols = (t-2):t, na.rm = TRUE)
#     }
#   }
#   
#   for (d in 1:nd) {
#     idx <- d:min(nd, d + 2L)
#     hi3[, d] <- matrixStats::rowMaxs(m_end[, idx, drop = FALSE], na.rm = TRUE)
#   }
#   hi3
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
# .exceed_dt_one_day <- function(hi3_vec, thr_mat, xy_cells, day_chr) {
#   rp_layer <- .max_rp_layer(hi3_vec, thr_mat)
#   keep <- which(rp_layer >= 1L)
#   if (!length(keep)) {
#     return(data.table(
#       x = numeric(), y = numeric(), day = character(),
#       hi3_intensity = numeric(), rp_max_years = integer(), rp_max_layer = integer()
#     ))
#   }
#   data.table(
#     x = xy_cells[keep, 1],
#     y = xy_cells[keep, 2],
#     day = day_chr,
#     hi3_intensity = hi3_vec[keep],
#     rp_max_years = RP_YEARS[rp_layer[keep]],
#     rp_max_layer = rp_layer[keep]
#   )
# }
# 
# # -----------------------------
# # RESAMPLING MAPS (0.25 -> TARGET)
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
# .build_max_agg_map_025_to_target <- function(ref_025, ref_tgt) {
#   xy_025 <- terra::xyFromCell(ref_025, seq_len(terra::ncell(ref_025)))
#   g <- terra::cellFromXY(ref_tgt, xy_025)
#   
#   ord <- order(g, na.last = TRUE)
#   g_ord <- g[ord]
#   ok <- !is.na(g_ord)
#   ord <- ord[ok]
#   g_ord <- g_ord[ok]
#   
#   r <- rle(g_ord)
#   starts <- c(1L, cumsum(r$lengths)[-length(r$lengths)] + 1L)
#   
#   list(
#     xy_tgt = terra::xyFromCell(ref_tgt, seq_len(terra::ncell(ref_tgt))),
#     ord = ord,
#     groups = r$values,
#     starts = starts,
#     lens = r$lengths,
#     ntgt = terra::ncell(ref_tgt)
#   )
# }
# 
# .apply_group_max <- function(src_vec_025, map) {
#   out <- rep(NA_real_, map$ntgt)
#   v <- src_vec_025[map$ord]
#   for (i in seq_along(map$groups)) {
#     s <- map$starts[i]
#     e <- s + map$lens[i] - 1L
#     m <- suppressWarnings(max(v[s:e], na.rm = TRUE))
#     out[map$groups[i]] <- if (is.finite(m)) m else NA_real_
#   }
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
# # WORKER-LOCAL CACHE (terra pointers must be built inside each worker)
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
#     max_map <- .build_max_agg_map_025_to_target(ref_025, ref_t)
#     
#     thr_t_bil <- matrix(NA_real_, nrow = nt, ncol = 10)
#     thr_t_max <- matrix(NA_real_, nrow = nt, ncol = 10)
#     
#     for (k in 1:10) {
#       thr_t_bil[, k] <- .apply_bilinear_renorm(thr_025[, k], bil_map$idx, bil_map$w)
#       thr_t_max[, k] <- .apply_group_max(thr_025[, k], max_map)
#     }
#     
#     targets_cache[[tag]] <- list(
#       ref = ref_t,
#       n   = nt,
#       xy  = xy_t,
#       bil_map = bil_map,
#       max_map = max_map,
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
#   out_025 <- file.path(OUT_025, sprintf("hi3_exceedances_maxrp_%d.qs", yy))
#   
#   # early skip if all outputs are already done
#   all_skip <- .should_skip(out_025)
#   for (tag in names(TARGETS)) {
#     tinfo <- TARGETS[[tag]]
#     out_bil <- file.path(tinfo$out_bil, sprintf("hi3_exceedances_maxrp_%d.qs", yy))
#     out_max <- file.path(tinfo$out_max, sprintf("hi3_exceedances_maxrp_%d.qs", yy))
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
#       out_bil <- file.path(tinfo$out_bil, sprintf("hi3_exceedances_maxrp_%d.qs", yy))
#       out_max <- file.path(tinfo$out_max, sprintf("hi3_exceedances_maxrp_%d.qs", yy))
#       if (!file.exists(out_bil)) qsave(data.table(), out_bil, preset = "high")
#       if (!file.exists(out_max)) qsave(data.table(), out_max, preset = "high")
#     }
#     return(invisible(TRUE))
#   }
#   
#   dmin <- as.Date(min(d_this))
#   dmax <- as.Date(max(d_this))
#   buf_dates <- seq(dmin - 2, dmax + 2, by = "day")
#   buf_chr <- format(buf_dates, "%Y-%m-%d")
#   buf_chr <- buf_chr[buf_chr %in% names(date_to_file)]
#   buf_dates <- as.Date(buf_chr)
#   
#   o <- order(buf_dates)
#   buf_dates <- buf_dates[o]
#   buf_chr <- buf_chr[o]
#   
#   nd_buf <- length(buf_chr)
#   if (nd_buf == 0L) {
#     if (!file.exists(out_025)) qsave(data.table(), out_025, preset = "high")
#     for (tag in names(TARGETS)) {
#       tinfo <- TARGETS[[tag]]
#       out_bil <- file.path(tinfo$out_bil, sprintf("hi3_exceedances_maxrp_%d.qs", yy))
#       out_max <- file.path(tinfo$out_max, sprintf("hi3_exceedances_maxrp_%d.qs", yy))
#       if (!file.exists(out_bil)) qsave(data.table(), out_bil, preset = "high")
#       if (!file.exists(out_max)) qsave(data.table(), out_max, preset = "high")
#     }
#     return(invisible(TRUE))
#   }
#   
#   # build daily maxima on native 0.25 grid
#   H_buf <- matrix(NA_real_, nrow = n025, ncol = nd_buf)
#   for (i in seq_len(nd_buf)) {
#     f <- date_to_file[[buf_chr[i]]]
#     dt_hourly <- .read_hourly_day(f)
#     H_buf[, i] <- .daily_max_vec_025(dt_hourly, ref_025, n025)
#   }
#   
#   hi3_buf_025 <- .hi3_from_dailymax(H_buf)
#   
#   idx_year <- which(buf_chr %in% d_this)
#   days_year_chr <- buf_chr[idx_year]
#   nd_year <- length(idx_year)
#   
#   # -----------------------------
#   # NATIVE 0.25 OUTPUT
#   # -----------------------------
#   if (!.should_skip(out_025)) {
#     rows <- vector("list", nd_year)
#     for (j in seq_along(idx_year)) {
#       i <- idx_year[j]
#       d <- .exceed_dt_one_day(hi3_buf_025[, i], thr_025, xy_025, days_year_chr[j])
#       if (nrow(d)) d[, `:=`(resolution_deg = 0.25, resample_method = "native")]
#       rows[[j]] <- d
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
#     out_bil <- file.path(tinfo$out_bil, sprintf("hi3_exceedances_maxrp_%d.qs", yy))
#     out_max <- file.path(tinfo$out_max, sprintf("hi3_exceedances_maxrp_%d.qs", yy))
#     
#     if (.should_skip(out_bil) && .should_skip(out_max)) next
#     
#     nt <- tcache$n
#     xy_t <- tcache$xy
#     bil_map <- tcache$bil_map
#     max_map <- tcache$max_map
#     thr_bil <- tcache$thr_bil
#     thr_max <- tcache$thr_max
#     
#     hi3_t_bil <- matrix(NA_real_, nrow = nt, ncol = nd_year)
#     hi3_t_max <- matrix(NA_real_, nrow = nt, ncol = nd_year)
#     
#     for (j in seq_along(idx_year)) {
#       i <- idx_year[j]
#       v <- hi3_buf_025[, i]
#       hi3_t_bil[, j] <- .apply_bilinear_renorm(v, bil_map$idx, bil_map$w)
#       hi3_t_max[, j] <- .apply_group_max(v, max_map)
#     }
#     
#     if (!.should_skip(out_bil)) {
#       rows <- vector("list", nd_year)
#       for (j in seq_len(nd_year)) {
#         d <- .exceed_dt_one_day(hi3_t_bil[, j], thr_bil, xy_t, days_year_chr[j])
#         if (nrow(d)) d[, `:=`(resolution_deg = tinfo$resolution_deg, resample_method = "bilinear")]
#         rows[[j]] <- d
#       }
#       out_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
#       setorder(out_dt, day, y, x)
#       qsave(out_dt, out_bil, preset = "high")
#     }
#     
#     if (!.should_skip(out_max)) {
#       rows <- vector("list", nd_year)
#       for (j in seq_len(nd_year)) {
#         d <- .exceed_dt_one_day(hi3_t_max[, j], thr_max, xy_t, days_year_chr[j])
#         if (nrow(d)) d[, `:=`(resolution_deg = tinfo$resolution_deg, resample_method = "max")]
#         rows[[j]] <- d
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
# scripts/00_build_ref_grids_from_stm_heatindex_southeast.R
# =============================================================================
# Build heatindex reference grids from STM1/STM2 summertime MEDIAN values
# (NO hard-coded 0.34/0.44; resolutions come from STM diagnostics CSV)
#
# Inputs (ONLY):
#   1) A CSV with columns: variable, domain, day, stm1, stm2  (at minimum)
#   2) A reference GeoTIFF mask/grid (land mask) defining CRS + extent (+ NA mask)
#
# Outputs:
#   data/out/09_performance/ref_grids/
#     heatindex_0.25/
#       ref_grid_heatindex_0.25.tif
#       ref_grid_cells_heatindex_0.25.qs
#     heatindex_<R1>_bilinear/
#       ref_grid_heatindex_<R1>_bilinear.tif
#       ref_grid_cells_heatindex_<R1>_bilinear.qs
#     heatindex_<R1>_max/
#       ref_grid_heatindex_<R1>_max.tif
#       ref_grid_cells_heatindex_<R1>_max.qs
#     heatindex_<R2>_bilinear/
#     heatindex_<R2>_max/
#
# Also writes:
#   ref_grids/ref_map_by_resolution_heatindex.qs
#   ref_grids/dataset_to_refgrid_heatindex.qs
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
HI_COV_CSV <- file.path(
  ROOT,
  "data/out/03_cov_model/heatindex/heatindex_daily_optimal_cov_model_southeast.csv"
)

# Land mask GeoTIFF from experimental covariance scripts (lon/lat)
MASK_TIF <- file.path(ROOT, "data/out/01_vars/precipitation/_mask_target.tif")

# Where to write ref grids
REF_DIR <- file.path(ROOT, "data/out/03_cov_model/heatindex")
dir_create(REF_DIR, recurse = TRUE)

DOMAIN_NAME <- "southeast"
MONTHS_KEEP <- 6:8          # summertime (JJA)
BASE_RES_DEG <- 0.25        # native base grid resolution
RES_TAG_DIGITS <- 2         # folder tag precision
OVERWRITE_REF <- FALSE

.log <- function(fmt, ...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))

stopifnot(file_exists(HI_COV_CSV), file_exists(MASK_TIF))

# -----------------------------
# Helpers
# -----------------------------
.fmt_res_tag <- function(res_deg, digits = 2) sprintf(paste0("%.", digits, "f"), res_deg)

# Parse dates robustly (your CSV uses m/d/YYYY in "day")
.parse_day <- function(x) {
  d <- as.IDate(x, format = "%Y-%m-%d")
  if (all(is.na(d))) d <- as.IDate(x, format = "%m/%d/%Y")
  if (all(is.na(d))) d <- as.IDate(x, format = "%m/%d/%y")
  d
}

# Extract STM-based target resolutions (median, months filtered)
.get_stm_resolutions <- function(cov_csv, domain, months_keep, base_res_deg) {
  dt <- fread(cov_csv)
  
  req <- c("variable", "domain", "day", "stm1", "stm2")
  miss <- setdiff(req, names(dt))
  if (length(miss)) stop("Missing columns in cov CSV: ", paste(miss, collapse = ", "))
  
  dt <- dt[variable == "heatindex" & domain == domain]
  if (!nrow(dt)) stop("No rows for variable=heatindex and domain=", domain)
  
  dt[, day_date := .parse_day(day)]
  if (all(is.na(dt$day_date))) stop("Could not parse 'day' in: ", cov_csv)
  
  dt <- dt[month(day_date) %in% months_keep]
  if (!nrow(dt)) stop("No rows after filtering to months=", paste(months_keep, collapse=","))
  
  stm1_med <- median(as.numeric(dt$stm1), na.rm = TRUE)
  stm2_med <- median(as.numeric(dt$stm2), na.rm = TRUE)
  
  if (!is.finite(stm1_med) || !is.finite(stm2_med)) stop("Non-finite STM medians; check stm1/stm2 values.")
  
  res1 <- base_res_deg * stm1_med
  res2 <- base_res_deg * stm2_med
  
  if (!is.finite(res1) || !is.finite(res2) || res1 <= 0 || res2 <= 0) {
    stop("Derived resolutions invalid: res1=", res1, " res2=", res2)
  }
  
  list(stm1_median = stm1_med, stm2_median = stm2_med, res_deg = c(res1, res2))
}

# Build a ref grid by resampling the mask to the new resolution
# This preserves NA mask and guarantees CRS/extent compatibility.
.build_ref_from_mask <- function(mask_tif, res_deg, out_tif, overwrite = FALSE) {
  if (file_exists(out_tif) && !overwrite) return(out_tif)
  
  r_mask <- terra::rast(mask_tif)[[1]]
  
  # Create target template grid with same extent/crs; align origin to mask for stable grid lines
  r_tmpl <- terra::rast(ext = terra::ext(r_mask), resolution = res_deg, crs = terra::crs(r_mask))
  try(terra::origin(r_tmpl) <- terra::origin(r_mask), silent = TRUE)
  
  # Resample mask to template (nearest keeps NA/1 patterns)
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
  
  # Keep cell index + coordinates; keep all cells (including NA mask) so mappings are stable
  # (Downstream scripts can ignore NA cells if needed.)
  dt_out <- dt[, .(cell_id = as.integer(cell), x = as.numeric(x), y = as.numeric(y))]
  
  dir_create(path_dir(out_qs), recurse = TRUE)
  qs::qsave(dt_out, out_qs, preset = "high")
  out_qs
}

# -----------------------------
# Compute STM-based resolutions
# -----------------------------
stm <- .get_stm_resolutions(HI_COV_CSV, DOMAIN_NAME, MONTHS_KEEP, BASE_RES_DEG)

res_raw <- stm$res_deg
res_raw <- unique(res_raw)

# Round to 2 decimals for scenario naming, but keep numeric resolution as-is for raster template
# NOTE: if rounding collapses the two, we still build only one unique res_tag (and both bilinear+max for it).
res_tag <- vapply(res_raw, .fmt_res_tag, character(1), digits = RES_TAG_DIGITS)

.log("STM medians (JJA) | stm1=%.6f | stm2=%.6f", stm$stm1_median, stm$stm2_median)
.log("Derived resolutions (deg) | raw: %s | tags: %s",
     paste(sprintf("%.6f", res_raw), collapse = ", "),
     paste(res_tag, collapse = ", "))

# Always include native 0.25 ref grid
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
  
  # Native gets a simple tag; derived resolutions become two scenario tags (bilinear+max)
  if (abs(res_deg - BASE_RES_DEG) < 1e-12) {
    ref_grid_tag <- paste0("heatindex_", tag)
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
    # Create BOTH bilinear + max folders at this derived resolution
    out <- rbindlist(lapply(c("bilinear", "max"), function(method) {
      ref_grid_tag <- paste0("heatindex_", tag, "_", method)
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

# Convenience mapping: dataset_tag -> ref grid
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
OUT_REF_MAP     <- file.path(REF_DIR, "ref_map_by_resolution_heatindex.qs")
OUT_DATASET_MAP <- file.path(REF_DIR, "dataset_to_refgrid_heatindex.qs")

qsave(ref_map,     OUT_REF_MAP, preset = "high")
qsave(dataset_map, OUT_DATASET_MAP, preset = "high")

.log("Wrote: %s (rows=%d)", OUT_REF_MAP, nrow(ref_map))
.log("Wrote: %s (rows=%d)", OUT_DATASET_MAP, nrow(dataset_map))

# Print what your downstream exceedance script should point to
.log("REF GRID OUTPUTS:")
print(dataset_map[])

# Southeast exceedances ----
# scripts/05_hi1_exceedances_native_vs_resampled_MAXRP_reanalysis_all_hourly_v2.R
# ------------------------------------------------------------------------------
# HI1 exceedances (SINGLE-DAY DAILY MAX HI) vs HI1 return periods
#
# Updates:
#   1) Uses 1-day return periods (HI1) instead of HI3
#   2) Exceedance intensity is the DAILY MAX for that single day (no 3-day averaging)
#   3) NO hard-coded target resolutions (0.34/0.44). Targets are read from:
#        data/out/03_cov_model/heatindex/ref_grids/dataset_to_refgrid_heatindex.qs
#      created by your STM-driven ref-grid builder script.
#
# Inputs:
#   - Hourly HI daily files (qs):  data/out/01_vars/heatindex/YYYY-MM-DD.qs
#   - Return period raster (nc):   data/out/04_rp/heatindex_idf_hi1_return_periods.nc
#   - Ref-grid map (qs):           data/out/03_cov_model/heatindex/ref_grids/dataset_to_refgrid_heatindex.qs
#
# Outputs (year-chunked .qs):
#   - Native (0.25°):                 data/out/05_exceedances/heatindex/0.25
#   - Resampled (targets):            data/out/05_exceedances/heatindex/<res>_bilinear
#                                    data/out/05_exceedances/heatindex/<res>_max
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(data.table)
  library(qs)
  library(future)
  library(future.apply)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

IN_HOURLY_DIR <- file.path(ROOT, "data/out/01_vars/heatindex")

# HI1 return periods (1-day duration)
RP_NC <- file.path(ROOT, "data/out/04_rp/heatindex_idf_hi1_return_periods_bbox.nc")

# Ref-grid mapping built from STM medians
REF_MAP_QS <- file.path(ROOT, "data/out/03_cov_model/heatindex/dataset_to_refgrid_heatindex.qs")

# Outputs (keep legacy layout style)
OUT_BASE <- file.path(ROOT, "data/out/05_exceedances/heatindex")
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
  
  val_col <- if ("heat_index" %in% names(dt)) "heat_index" else if ("intensity" %in% names(dt)) "intensity" else NA_character_
  if (is.na(val_col)) stop("Missing heat_index/intensity in: ", path, call. = FALSE)
  
  dt[, .(
    lon = as.numeric(get(lon_col)),
    lat = as.numeric(get(lat_col)),
    dt  = get(dt_col),
    val = as.numeric(get(val_col))
  )]
}

.daily_max_vec_025 <- function(dt_hourly, ref_025, ncell025) {
  cell <- terra::cellFromXY(ref_025, cbind(dt_hourly$lon, dt_hourly$lat))
  dt_hourly[, cell := cell]
  dt_hourly <- dt_hourly[is.finite(val) & !is.na(cell)]
  out <- rep(NA_real_, ncell025)
  if (nrow(dt_hourly) == 0L) return(out)
  m <- dt_hourly[, .(dmax = max(val, na.rm = TRUE)), by = cell]
  out[m$cell] <- m$dmax
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

.exceed_dt_one_day <- function(hi1_vec, thr_mat, xy_cells, day_chr) {
  rp_layer <- .max_rp_layer(hi1_vec, thr_mat)
  keep <- which(rp_layer >= 1L)
  
  if (!length(keep)) {
    return(data.table(
      x = numeric(), y = numeric(), day = character(),
      hi1_intensity = numeric(), rp_max_years = integer(), rp_max_layer = integer()
    ))
  }
  
  data.table(
    x = xy_cells[keep, 1],
    y = xy_cells[keep, 2],
    day = day_chr,
    hi1_intensity = hi1_vec[keep],
    rp_max_years = RP_YEARS[rp_layer[keep]],
    rp_max_layer = rp_layer[keep]
  )
}

# -----------------------------
# RESAMPLING MAPS (0.25 -> TARGET)
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

.build_max_agg_map_025_to_target <- function(ref_025, ref_tgt) {
  xy_025 <- terra::xyFromCell(ref_025, seq_len(terra::ncell(ref_025)))
  g <- terra::cellFromXY(ref_tgt, xy_025)
  
  ord <- order(g, na.last = TRUE)
  g_ord <- g[ord]
  ok <- !is.na(g_ord)
  ord <- ord[ok]
  g_ord <- g_ord[ok]
  
  r <- rle(g_ord)
  starts <- c(1L, cumsum(r$lengths)[-length(r$lengths)] + 1L)
  
  list(
    xy_tgt = terra::xyFromCell(ref_tgt, seq_len(terra::ncell(ref_tgt))),
    ord = ord,
    groups = r$values,
    starts = starts,
    lens = r$lengths,
    ntgt = terra::ncell(ref_tgt)
  )
}

.apply_group_max <- function(src_vec_025, map) {
  out <- rep(NA_real_, map$ntgt)
  v <- src_vec_025[map$ord]
  for (i in seq_along(map$groups)) {
    s <- map$starts[i]
    e <- s + map$lens[i] - 1L
    m <- suppressWarnings(max(v[s:e], na.rm = TRUE))
    out[map$groups[i]] <- if (is.finite(m)) m else NA_real_
  }
  out
}

# -----------------------------
# LOAD TARGETS from REF_MAP_QS
# -----------------------------
stopifnot(file.exists(RP_NC), file.exists(REF_MAP_QS))

ref_map <- as.data.table(qread(REF_MAP_QS))

# Expect dataset_tag like: heatindex_0.25 OR heatindex_0.16_bilinear / heatindex_0.16_max
if (!all(c("dataset_tag","ref_tif","res_deg") %in% names(ref_map))) {
  stop("REF_MAP_QS missing required columns (dataset_tag, ref_tif, res_deg): ", REF_MAP_QS)
}

# Native 0.25 ref
row_025 <- ref_map[dataset_tag == "heatindex_0.25"]
if (nrow(row_025) != 1L) stop("Could not find unique 'heatindex_0.25' in: ", REF_MAP_QS)
REF_025 <- row_025$ref_tif[1]

# Build targets: for each unique derived resolution tag, require both bilinear and max rows
targets_dt <- ref_map[dataset_tag != "heatindex_0.25"]
targets_dt <- targets_dt[grepl("^heatindex_[0-9]+\\.[0-9]{2}_(bilinear|max)$", dataset_tag)]

if (!nrow(targets_dt)) stop("No derived heatindex targets found in REF_MAP_QS: ", REF_MAP_QS)

# Derive a "res_tag" and method from dataset_tag
targets_dt[, res_tag := sub("^heatindex_([0-9]+\\.[0-9]{2})_.*$", "\\1", dataset_tag)]
targets_dt[, method  := sub("^heatindex_[0-9]+\\.[0-9]{2}_(bilinear|max)$", "\\1", dataset_tag)]

# Keep only res_tags that have BOTH methods present
ok_res <- targets_dt[, .(has_bil = any(method=="bilinear"), has_max = any(method=="max")), by = res_tag][has_bil & has_max, res_tag]
targets_dt <- targets_dt[res_tag %in% ok_res]

# Build TARGETS list keyed by res_tag (folder naming uses res_tag)
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

# Create output dirs
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
# WORKER-LOCAL CACHE (terra pointers must be built inside each worker)
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
    
    # they should share same geometry; but we keep separate refs because you store them separately
    nt   <- terra::ncell(ref_bil)
    xy_t <- terra::xyFromCell(ref_bil, seq_len(nt))
    
    bil_map <- .build_bilinear_map_target_to_025(ref_025, ref_bil)
    max_map <- .build_max_agg_map_025_to_target(ref_025, ref_bil)  # geometry same as bil
    
    thr_t_bil <- matrix(NA_real_, nrow = nt, ncol = 10)
    thr_t_max <- matrix(NA_real_, nrow = nt, ncol = 10)
    
    for (k in 1:10) {
      thr_t_bil[, k] <- .apply_bilinear_renorm(thr_025[, k], bil_map$idx, bil_map$w)
      thr_t_max[, k] <- .apply_group_max(thr_025[, k], max_map)
    }
    
    targets_cache[[rt]] <- list(
      n   = nt,
      xy  = xy_t,
      bil_map = bil_map,
      max_map = max_map,
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
  
  out_025 <- file.path(OUT_025, sprintf("hi1_exceedances_maxrp_%d.qs", yy))
  
  # early skip if all outputs already done
  all_skip <- .should_skip(out_025)
  for (rt in names(TARGETS)) {
    tinfo <- TARGETS[[rt]]
    out_bil <- file.path(tinfo$out_bil, sprintf("hi1_exceedances_maxrp_%d.qs", yy))
    out_max <- file.path(tinfo$out_max, sprintf("hi1_exceedances_maxrp_%d.qs", yy))
    all_skip <- all_skip && .should_skip(out_bil) && .should_skip(out_max)
  }
  if (all_skip) return(invisible(TRUE))
  
  d_this <- sort(date_chr_all[years_all == yy])
  
  # write empties if no dates
  if (!length(d_this)) {
    if (!file.exists(out_025)) qsave(data.table(), out_025, preset = "high")
    for (rt in names(TARGETS)) {
      tinfo <- TARGETS[[rt]]
      out_bil <- file.path(tinfo$out_bil, sprintf("hi1_exceedances_maxrp_%d.qs", yy))
      out_max <- file.path(tinfo$out_max, sprintf("hi1_exceedances_maxrp_%d.qs", yy))
      if (!file.exists(out_bil)) qsave(data.table(), out_bil, preset = "high")
      if (!file.exists(out_max)) qsave(data.table(), out_max, preset = "high")
    }
    return(invisible(TRUE))
  }
  
  nd_year <- length(d_this)
  
  # -----------------------------
  # NATIVE 0.25 OUTPUT
  # -----------------------------
  if (!.should_skip(out_025)) {
    rows <- vector("list", nd_year)
    for (j in seq_len(nd_year)) {
      day_chr <- d_this[j]
      f <- date_to_file[[day_chr]]
      dt_hourly <- .read_hourly_day(f)
      hi1_vec_025 <- .daily_max_vec_025(dt_hourly, ref_025, n025)
      
      d <- .exceed_dt_one_day(hi1_vec_025, thr_025, xy_025, day_chr)
      if (nrow(d)) d[, `:=`(resolution_deg = 0.25, resample_method = "native")]
      rows[[j]] <- d
    }
    out_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
    setorder(out_dt, day, y, x)
    qsave(out_dt, out_025, preset = "high")
  }
  
  # -----------------------------
  # TARGET GRIDS OUTPUTS (bilinear + max)
  # -----------------------------
  for (rt in names(TARGETS)) {
    tinfo <- TARGETS[[rt]]
    tcache <- cache$targets[[rt]]
    
    out_bil <- file.path(tinfo$out_bil, sprintf("hi1_exceedances_maxrp_%d.qs", yy))
    out_max <- file.path(tinfo$out_max, sprintf("hi1_exceedances_maxrp_%d.qs", yy))
    
    if (.should_skip(out_bil) && .should_skip(out_max)) next
    
    nt <- tcache$n
    xy_t <- tcache$xy
    bil_map <- tcache$bil_map
    max_map <- tcache$max_map
    thr_bil <- tcache$thr_bil
    thr_max <- tcache$thr_max
    res_deg <- tcache$resolution_deg
    
    if (!.should_skip(out_bil)) {
      rows <- vector("list", nd_year)
      for (j in seq_len(nd_year)) {
        day_chr <- d_this[j]
        f <- date_to_file[[day_chr]]
        dt_hourly <- .read_hourly_day(f)
        hi1_vec_025 <- .daily_max_vec_025(dt_hourly, ref_025, n025)
        
        hi1_t_bil <- .apply_bilinear_renorm(hi1_vec_025, bil_map$idx, bil_map$w)
        d <- .exceed_dt_one_day(hi1_t_bil, thr_bil, xy_t, day_chr)
        if (nrow(d)) d[, `:=`(resolution_deg = res_deg, resample_method = "bilinear")]
        rows[[j]] <- d
      }
      out_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
      setorder(out_dt, day, y, x)
      qsave(out_dt, out_bil, preset = "high")
    }
    
    if (!.should_skip(out_max)) {
      rows <- vector("list", nd_year)
      for (j in seq_len(nd_year)) {
        day_chr <- d_this[j]
        f <- date_to_file[[day_chr]]
        dt_hourly <- .read_hourly_day(f)
        hi1_vec_025 <- .daily_max_vec_025(dt_hourly, ref_025, n025)
        
        hi1_t_max <- .apply_group_max(hi1_vec_025, max_map)
        d <- .exceed_dt_one_day(hi1_t_max, thr_max, xy_t, day_chr)
        if (nrow(d)) d[, `:=`(resolution_deg = res_deg, resample_method = "max")]
        rows[[j]] <- d
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
