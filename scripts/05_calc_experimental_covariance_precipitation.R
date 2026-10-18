# # scripts/07_experimental_covariance_precip_daily_weekly_parallel_relic.R
# # PRECIPITATION experimental covariance (RELIC estimator)
# #   (A) DAILY: uses HOURLY precipitation within each day  (t axis in hours)
# #   (B) WEEKLY: uses DAILY TOTAL precipitation within each ISO week (UTC) (t axis in days)
# #
# # Key changes vs heat index:
# # - variable is precipitation (auto-detected from common column names)
# # - replace 0 with NA (fast skip for zeros; matches your spec)
# # - weekly uses sum across all hours of each day (daily totals), then covariance over days
# #
# # I/O:
# #   IN  = V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation
# #   OUT = V:/users/hquintal/phd1_cluster_southeast_v2/data/out/02_exp_cov/precipitation/{daily,weekly}
# 
# library(sf)
# test <- sf::st_read('V:/users/hquintal/phd2_cluster_conus/data/in/shapefiles/NCA_Regions/NCA_Regions.shp')
# names
# 
# suppressPackageStartupMessages({
#   library(qs)
#   library(data.table)
#   library(lubridate)
#   library(future)
#   library(future.apply)
#   library(progressr)
# })
# 
# # USER SETTINGS----
# IN_DIR_DAILY  <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation"
# OUT_DIR_DAILY <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/02_exp_cov/precipitation/daily"
# OUT_DIR_WEEK  <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/02_exp_cov/precipitation/weekly"
# 
# dir.create(OUT_DIR_DAILY, recursive = TRUE, showWarnings = FALSE)
# dir.create(OUT_DIR_WEEK,  recursive = TRUE, showWarnings = FALSE)
# 
# WORKERS <- 28
# plan(multisession, workers = WORKERS)
# 
# # Performance / safety (controls kron vs kronloop)
# max_pairs <- 2e7  # n_cells^2 cutoff
# 
# # RELIC CORE (unchanged)----
# coord2dist <- function(c1, c2) {
#   if (is.data.frame(c1)) c1 <- as.matrix(c1)
#   if (is.data.frame(c2)) c2 <- as.matrix(c2)
#   if (is.list(c1)) { c1 <- c1[[1]]; c2 <- c2[[1]] }
#   x1_sq <- rowSums(c1^2); x2_sq <- rowSums(c2^2)
#   D2 <- outer(x1_sq, rep(1, nrow(c2))) +
#     outer(rep(1, nrow(c1)), x2_sq) -
#     2 * tcrossprod(c1, c2)
#   D2[D2 < 0] <- 0
#   sqrt(D2)
# }
# 
# pairsindex <- function(c1, c2, rLag, rLagTol, dist=list("coord2dist"), method=NULL) {
#   if (is.data.frame(c1)) c1 <- as.matrix(c1)
#   if (is.data.frame(c2)) c2 <- as.matrix(c2)
#   n1 <- nrow(c1); n2 <- nrow(c2); nr <- length(rLag)
#   if (!is.list(dist)) stop("dist must be a list")
#   
#   if (is.null(method)) {
#     method <- if (n1 * n2 <= max_pairs) "kron" else "kronloop"
#   }
#   
#   idxpairs <- vector("list", nr)
#   
#   call_dist <- function(A, B) {
#     if (length(dist) == 1) {
#       do.call(dist[[1]], list(A, B))
#     } else {
#       do.call(dist[[1]], list(A, B, dist[[2]]))
#     }
#   }
#   
#   switch(method,
#          kron = {
#            D <- call_dist(c1, c2)
#            for (ir in seq_len(nr)) {
#              mask <- (D >= rLag[ir] - rLagTol[ir]) & (D <= rLag[ir] + rLagTol[ir])
#              idxpairs[[ir]] <- which(mask, arr.ind = TRUE)
#            }
#          },
#          kronloop = {
#            for (ir in seq_len(nr)) idxpairs[[ir]] <- matrix(numeric(0), ncol = 2)
#            for (i in seq_len(n1)) {
#              d_i <- call_dist(matrix(c1[i,], nrow = 1), c2)
#              for (ir in seq_len(nr)) {
#                js <- which(d_i >= rLag[ir] - rLagTol[ir] & d_i <= rLag[ir] + rLagTol[ir])
#                if (length(js) > 0) idxpairs[[ir]] <- rbind(idxpairs[[ir]], cbind(i, js))
#              }
#            }
#          },
#          stop("method must be 'kron' or 'kronloop'")
#   )
#   
#   for (ir in seq_len(nr)) {
#     if (nrow(idxpairs[[ir]]) > 0) {
#       idxpairs[[ir]] <- idxpairs[[ir]][order(idxpairs[[ir]][,1], idxpairs[[ir]][,2]), , drop = FALSE]
#     }
#   }
#   idxpairs
# }
# 
# stcov <- function(
#     Zi, cMSi, tMEi,
#     Zj, cMSj, tMEj,
#     rLag, rLagTol,
#     tLag, tLagTol,
#     dist = list("coord2dist"),
#     method = NULL
# ) {
#   nMSi <- nrow(Zi);  nMEi <- ncol(Zi)
#   if (!all(dim(cMSi) == c(nMSi, 2))) stop("cMSi must be nMSi × 2")
#   if (length(tMEi) != nMEi) stop("tMEi must be length nMEi")
#   
#   nMSj <- nrow(Zj);  nMEj <- ncol(Zj)
#   if (!all(dim(cMSj) == c(nMSj, 2))) stop("cMSj must be nMSj × 2")
#   if (length(tMEj) != nMEj) stop("tMEj must be length nMEj")
#   
#   nr <- length(rLag); if (length(rLagTol) != nr) stop("rLag and rLagTol length mismatch")
#   nt <- length(tLag); if (length(tLagTol) != nt) stop("tLag and tLagTol length mismatch")
#   
#   if (is.null(method)) {
#     method <- if (nMSi * nMSj <= max_pairs) "kron" else "kronloop"
#   }
#   
#   idxpairsMS <- pairsindex(cMSi, cMSj, rLag, rLagTol, dist, method)
#   tMat <- outer(tMEi, tMEj, function(ti, tj) tj - ti)
#   
#   C  <- matrix(NA_real_, nrow = nr, ncol = nt)
#   np <- matrix(0L,       nrow = nr, ncol = nt)
#   
#   for (ir in seq_len(nr)) {
#     pairs_ij <- idxpairsMS[[ir]]
#     if (nrow(pairs_ij) == 0L) next
#     
#     for (it in seq_len(nt)) {
#       mask <- (tLag[it] - tLagTol[it] <= tMat) & (tMat <= tLag[it] + tLagTol[it])
#       if (!any(mask)) next
#       ev_idx <- which(mask, arr.ind = TRUE)
#       
#       nSP <- nrow(pairs_ij); nEP <- nrow(ev_idx)
#       sp_rep <- pairs_ij[rep(seq_len(nSP), times = nEP), , drop = FALSE]
#       ep_rep <- ev_idx[rep(seq_len(nEP), each = nSP), , drop = FALSE]
#       
#       head_vals <- Zi[cbind(sp_rep[,1], ep_rep[,1])]
#       tail_vals <- Zj[cbind(sp_rep[,2], ep_rep[,2])]
#       
#       valid <- !is.na(head_vals) & !is.na(tail_vals)
#       np[ir,it] <- sum(valid)
#       if (np[ir,it] == 0L) next
#       
#       h <- head_vals[valid]
#       t <- tail_vals[valid]
#       C[ir,it] <- mean(h * t) - mean(h) * mean(t)
#     }
#   }
#   
#   list(C = C, np = np, method = method)
# }
# 
# trim_after_negative_inclusive <- function(x, lags, np = NULL) {
#   keep <- !is.na(x) & !is.na(lags)
#   x <- x[keep]; lags <- lags[keep]
#   if (!is.null(np)) np <- np[keep]
#   
#   if (length(x) == 0L) return(list(cov = numeric(0), lags = numeric(0), np = integer(0)))
#   
#   neg_idx <- which(x < 0)[1]
#   if (!is.na(neg_idx)) {
#     x <- x[1:neg_idx]
#     lags <- lags[1:neg_idx]
#     if (!is.null(np)) np <- np[1:neg_idx]
#   }
#   
#   if (is.null(np)) np <- rep(NA_integer_, length(x))
#   list(cov = x, lags = lags, np = np)
# }
# 
# estimate_base_dist <- function(coords) {
#   ux <- sort(unique(coords[,1]))
#   uy <- sort(unique(coords[,2]))
#   dx <- diff(ux); dy <- diff(uy)
#   dx <- dx[dx > 0]; dy <- dy[dy > 0]
#   base <- suppressWarnings(min(median(dx, na.rm = TRUE), median(dy, na.rm = TRUE), na.rm = TRUE))
#   if (!is.finite(base) || base <= 0) base <- 0.5
#   base
# }
# 
# make_lags <- function(max_dist, max_time_hr, base_dist) {
#   r_raw <- base_dist * c(0:4, 6, 8, 12)
#   rLag <- sort(unique(r_raw[r_raw <= max_dist]))
#   if (length(rLag) < 2) rLag <- c(0, min(base_dist, max_dist))
#   rLagTol <- c(0, diff(rLag)/2)
#   
#   t_raw <- c(0:4, 6, 8, 12)
#   tLag <- sort(unique(t_raw[t_raw <= max_time_hr]))
#   if (length(tLag) < 2) tLag <- c(0, min(1, max_time_hr))
#   tLagTol <- c(0, diff(tLag)/2)
#   
#   list(rLag = rLag, rLagTol = rLagTol, tLag = tLag, tLagTol = tLagTol)
# }
# 
# make_lags_weekly_days <- function(max_dist, max_time_day, base_dist) {
#   r_raw <- base_dist * c(0:4, 6, 8, 12)
#   rLag <- sort(unique(r_raw[r_raw <= max_dist]))
#   if (length(rLag) < 2) rLag <- c(0, min(base_dist, max_dist))
#   rLagTol <- c(0, diff(rLag)/2)
#   
#   t_raw <- c(0:6)  # lags in days within a week
#   tLag <- sort(unique(t_raw[t_raw <= max_time_day]))
#   if (length(tLag) < 2) tLag <- c(0, min(1, max_time_day))
#   tLagTol <- c(0, diff(tLag)/2)
#   
#   list(rLag = rLag, rLagTol = rLagTol, tLag = tLag, tLagTol = tLagTol)
# }
# 
# mb <- function(x) as.numeric(object.size(x)) / (1024^2)
# 
# # PRECIP HELPERS----
# detect_precip_col <- function(dt) {
#   # try common names; include your actual column name
#   cand <- c(
#     "intensity",            # <- your files
#     "precip", "precipitation", "prcp", "rain", "tp", "apcp",
#     "precip_mm", "precip_m"
#   )
#   hit <- intersect(cand, names(dt))
#   if (length(hit) == 0L) {
#     stop("Could not detect precipitation column. Available names: ",
#          paste(names(dt), collapse = ", "))
#   }
#   hit[1]
# }
# 
# ensure_datetime_utc <- function(dt, qs_file) {
#   if ("datetime" %in% names(dt) && inherits(dt$datetime, "POSIXct")) {
#     attr(dt$datetime, "tzone") <- "UTC"
#     return(dt)
#   }
#   if ("datetime_utc" %in% names(dt)) {
#     dt[, datetime := as.POSIXct(datetime_utc, tz = "UTC")]
#     return(dt)
#   }
#   stop("No datetime or datetime_utc column found in: ", qs_file)
# }
# 
# replace_below_threshold_with_na <- function(x, thresh = 1) {
#   # keep NA as NA
#   x[!is.na(x) & x < thresh] <- NA_real_
#   x
# }
# 
# # (A) DAILY: HOURLY -> stcov per day (time axis in hours)----
# read_day_to_matrix_precip <- function(qs_file) {
#   dt <- as.data.table(qread(qs_file))
#   dt <- ensure_datetime_utc(dt, qs_file)
#   
#   if (!all(c("lon","lat") %in% names(dt))) {
#     stop("Missing lon/lat columns in: ", qs_file)
#   }
#   
#   pcol <- detect_precip_col(dt)
#   
#   # Keep only needed cols; replace 0 with NA for faster cov compute
#   dt <- dt[, .(lon, lat, datetime, precip = as.numeric(get(pcol)))]
#   dt[, precip := replace_below_threshold_with_na(precip, thresh = 1)]
#   
#   # Drop NA instances (precip now NA for zeros and missing)
#   dt <- dt[!is.na(precip)]
#   if (nrow(dt) == 0L) return(NULL)
#   
#   dt[, cell_id := paste0(sprintf("%.6f", lon), "_", sprintf("%.6f", lat))]
#   
#   # Stable ordered time axis from datetimes present
#   times <- sort(unique(dt$datetime))
#   if (length(times) < 2L) return(NULL
#                                  
#   )
#   dt[, t_idx := match(datetime, times)]
#   
#   setkey(dt, cell_id, t_idx)
#   wide <- dcast(
#     dt,
#     cell_id + lon + lat ~ t_idx,
#     value.var = "precip",
#     fun.aggregate = mean
#   )
#   if (nrow(wide) < 2L) return(NULL)
#   
#   coords <- as.matrix(wide[, .(lon, lat)])
#   time_cols <- setdiff(names(wide), c("cell_id", "lon", "lat"))
#   time_cols_ord <- time_cols[order(as.integer(time_cols))]
#   values <- as.matrix(wide[, ..time_cols_ord])
#   if (ncol(values) < 2L) return(NULL)
#   
#   t_hr <- as.numeric(difftime(times, times[1], units = "hours"))
#   
#   if (!is.finite(max(t_hr)) || max(t_hr) <= 0) {
#     warning(sprintf(
#       "[%s] suspicious time axis: n_times=%d max_time_hr=%g",
#       basename(qs_file), length(times), max(t_hr)
#     ))
#   }
#   
#   list(coords = coords, values = values, times = times, t_hr = t_hr)
# }
# 
# process_one_day_precip <- function(qs_file) {
#   day_str <- gsub("\\.qs$", "", basename(qs_file))
#   t0 <- Sys.time()
#   
#   mat <- read_day_to_matrix_precip(qs_file)
#   if (is.null(mat)) {
#     return(list(day = day_str, ok = FALSE, reason = "NULL matrix",
#                 dt_sec = as.numeric(difftime(Sys.time(), t0, units="secs")),
#                 n_space = 0L, n_time = 0L, method = NA_character_,
#                 mem_coords_mb = NA_real_, mem_values_mb = NA_real_))
#   }
#   
#   coords <- mat$coords
#   values <- mat$values
#   t_hr   <- mat$t_hr
#   
#   base_dist <- estimate_base_dist(coords)
#   max_dist  <- max(coord2dist(coords, coords), na.rm = TRUE)
#   max_time  <- max(t_hr, na.rm = TRUE)
#   lags <- make_lags(max_dist, max_time, base_dist)
#   
#   mem_coords_mb <- mb(coords)
#   mem_values_mb <- mb(values)
#   
#   out_space <- stcov(values, coords, t_hr, values, coords, t_hr,
#                      rLag = lags$rLag, rLagTol = lags$rLagTol,
#                      tLag = 0,         tLagTol = 0,
#                      method = NULL)
#   
#   Cr  <- out_space$C[, 1]
#   npr <- out_space$np[, 1]
#   
#   out_time <- stcov(values, coords, t_hr, values, coords, t_hr,
#                     rLag = 0,         rLagTol = 0,
#                     tLag = lags$tLag, tLagTol = lags$tLagTol,
#                     method = out_space$method)
#   
#   Ct  <- out_time$C[1, ]
#   npt <- out_time$np[1, ]
#   
#   res_space <- trim_after_negative_inclusive(Cr, lags$rLag, np = npr)
#   res_time  <- trim_after_negative_inclusive(Ct, lags$tLag, np = npt)
#   
#   dt_space <- unique(data.table(rLag = res_space$lags, Cov = res_space$cov, np = res_space$np), by = "rLag")
#   dt_time  <- unique(data.table(tLag = res_time$lags,  Cov = res_time$cov,  np = res_time$np),  by = "tLag")
#   
#   tag <- sprintf("precip_CONUS_day_%s", day_str)
#   fwrite(dt_space, file.path(OUT_DIR_DAILY, paste0(tag, "_cov_space.txt")), sep = "\t")
#   fwrite(dt_time,  file.path(OUT_DIR_DAILY, paste0(tag, "_cov_time.txt")),  sep = "\t")
#   
#   t1 <- Sys.time()
#   list(day = day_str, ok = TRUE, reason = NA_character_,
#        dt_sec = as.numeric(difftime(t1, t0, units="secs")),
#        n_space = nrow(dt_space), n_time = nrow(dt_time),
#        method = out_space$method,
#        mem_coords_mb = mem_coords_mb, mem_values_mb = mem_values_mb)
# }
# 
# run_daily_precip <- function() {
#   files <- sort(list.files(IN_DIR_DAILY, pattern = "^\\d{4}-\\d{2}-\\d{2}\\.qs$", full.names = TRUE))
#   stopifnot(length(files) > 0)
#   
#   day_strs <- gsub("\\.qs$", "", basename(files))
#   tags     <- sprintf("precip_CONUS_day_%s", day_strs)
#   
#   out_space_paths <- file.path(OUT_DIR_DAILY, paste0(tags, "_cov_space.txt"))
#   out_time_paths  <- file.path(OUT_DIR_DAILY, paste0(tags, "_cov_time.txt"))
#   have_space <- file.exists(out_space_paths)
#   have_time  <- file.exists(out_time_paths)
#   is_done <- have_space & have_time
#   
#   message(sprintf(
#     "[DAILY SKIP CHECK] total_days=%d | done=%d | todo=%d | space_missing=%d | time_missing=%d",
#     length(files), sum(is_done), sum(!is_done), sum(!have_space), sum(!have_time)
#   ))
#   
#   files_todo <- files[!is_done]
#   
#   handlers(global = TRUE)
#   t_all0 <- Sys.time()
#   
#   if (length(files_todo) == 0L) {
#     message("✅ Nothing to do — all DAILY outputs already exist.")
#     log_dt <- data.table(
#       day = day_strs,
#       ok = TRUE,
#       reason = NA_character_,
#       dt_sec = NA_real_,
#       n_space = NA_integer_,
#       n_time = NA_integer_,
#       method = NA_character_,
#       mem_coords_mb = NA_real_,
#       mem_values_mb = NA_real_,
#       had_cov_space = have_space,
#       had_cov_time  = have_time,
#       skipped_already_done = TRUE
#     )
#     fwrite(log_dt, file.path(OUT_DIR_DAILY, "run_log_parallel_relic_daily_precip.csv"))
#     message(sprintf("✅ DAILY DONE | workers=%d | total_elapsed=%s",
#                     WORKERS, format(difftime(Sys.time(), t_all0), digits = 3)))
#     return(invisible(NULL))
#   }
#   
#   with_progress({
#     p <- progressor(along = files_todo)
#     results <- future_lapply(files_todo, function(f) {
#       p(basename(f))
#       process_one_day_precip(f)
#     })
#   })
#   
#   log_done <- rbindlist(lapply(results, as.data.table), fill = TRUE)
#   log_done[, `:=`(had_cov_space = TRUE, had_cov_time = TRUE, skipped_already_done = FALSE)]
#   
#   if (any(is_done)) {
#     log_skip <- data.table(
#       day = day_strs[is_done],
#       ok = TRUE,
#       reason = NA_character_,
#       dt_sec = NA_real_,
#       n_space = NA_integer_,
#       n_time = NA_integer_,
#       method = NA_character_,
#       mem_coords_mb = NA_real_,
#       mem_values_mb = NA_real_,
#       had_cov_space = have_space[is_done],
#       had_cov_time  = have_time[is_done],
#       skipped_already_done = TRUE
#     )
#     log_dt <- rbindlist(list(log_done, log_skip), use.names = TRUE, fill = TRUE)
#   } else {
#     log_dt <- log_done
#   }
#   
#   # final filesystem truth
#   all_days_logged <- log_dt$day
#   all_tags_logged <- sprintf("precip_CONUS_day_%s", all_days_logged)
#   log_dt[, `:=`(
#     had_cov_space = file.exists(file.path(OUT_DIR_DAILY, paste0(all_tags_logged, "_cov_space.txt"))),
#     had_cov_time  = file.exists(file.path(OUT_DIR_DAILY, paste0(all_tags_logged, "_cov_time.txt")))
#   )]
#   
#   setorder(log_dt, day)
#   fwrite(log_dt, file.path(OUT_DIR_DAILY, "run_log_parallel_relic_daily_precip.csv"))
#   
#   message(sprintf("✅ DAILY DONE | workers=%d | processed=%d | skipped=%d | total_elapsed=%s",
#                   WORKERS, length(files_todo), sum(is_done),
#                   format(difftime(Sys.time(), t_all0), digits = 3)))
# }
# 
# # (B) WEEKLY: DAILY TOTALS -> stcov per ISO week (time axis in days)----
# read_day_to_dailytotal_precip <- function(qs_file, daily_thresh = 1) {
#   
#   dt <- as.data.table(qread(qs_file))
#   dt <- ensure_datetime_utc(dt, qs_file)
#   
#   if (!all(c("lon","lat") %in% names(dt))) {
#     stop("Missing lon/lat columns in: ", qs_file)
#   }
#   
#   pcol <- detect_precip_col(dt)
#   
#   # Use raw hourly values (DO NOT threshold hourly here)
#   dt <- dt[, .(lon, lat, precip = as.numeric(get(pcol)))]
#   
#   # Compute full 24-hour total per cell
#   out <- dt[, .(
#     prcp_day = if (all(is.na(precip))) NA_real_
#     else sum(precip, na.rm = TRUE)
#   ), by = .(lon, lat)]
#   
#   # Apply threshold at DAILY scale
#   out[!is.na(prcp_day) & prcp_day < daily_thresh,
#       prcp_day := NA_real_]
#   
#   out
# }
# 
# week_files_to_matrix_dailytotal <- function(files_week) {
#   day_dates <- as.Date(gsub("\\.qs$", "", basename(files_week)))
#   ord <- order(day_dates)
#   files_week <- files_week[ord]
#   day_dates  <- day_dates[ord]
#   
#   daily_list <- lapply(seq_along(files_week), function(i) {
#     dd <- read_day_to_dailytotal_precip(files_week[i])
#     dd[, day_idx := i]
#     dd
#   })
#   
#   dt <- rbindlist(daily_list, use.names = TRUE, fill = TRUE)
#   if (nrow(dt) == 0L) return(NULL)
#   
#   dt[, cell_id := paste0(sprintf("%.6f", lon), "_", sprintf("%.6f", lat))]
#   
#   setkey(dt, cell_id, day_idx)
#   wide <- dcast(
#     dt,
#     cell_id + lon + lat ~ day_idx,
#     value.var = "prcp_day",
#     fun.aggregate = mean
#   )
#   if (nrow(wide) < 2L) return(NULL)
#   
#   coords <- as.matrix(wide[, .(lon, lat)])
#   time_cols <- setdiff(names(wide), c("cell_id", "lon", "lat"))
#   time_cols_ord <- time_cols[order(as.integer(time_cols))]
#   values <- as.matrix(wide[, ..time_cols_ord])
#   if (ncol(values) < 2L) return(NULL)
#   
#   # time in DAYS: 0..(n_days-1)
#   t_day <- (seq_len(ncol(values)) - 1)
#   
#   list(coords = coords, values = values, t_day = t_day, day_dates = day_dates)
# }
# 
# process_one_week_precip <- function(week_id, files_week) {
#   t0 <- Sys.time()
#   
#   mat <- week_files_to_matrix_dailytotal(files_week)
#   if (is.null(mat)) {
#     return(list(week = week_id, ok = FALSE, reason = "NULL matrix",
#                 dt_sec = as.numeric(difftime(Sys.time(), t0, units="secs")),
#                 n_space = 0L, n_time = 0L, method = NA_character_,
#                 mem_coords_mb = NA_real_, mem_values_mb = NA_real_,
#                 n_days = length(files_week)))
#   }
#   
#   coords <- mat$coords
#   values <- mat$values
#   t_day  <- mat$t_day
#   
#   base_dist <- estimate_base_dist(coords)
#   max_dist  <- max(coord2dist(coords, coords), na.rm = TRUE)
#   max_time  <- max(t_day, na.rm = TRUE)
#   
#   lags <- make_lags_weekly_days(max_dist, max_time, base_dist)
#   
#   mem_coords_mb <- mb(coords)
#   mem_values_mb <- mb(values)
#   
#   out_space <- stcov(values, coords, t_day, values, coords, t_day,
#                      rLag = lags$rLag, rLagTol = lags$rLagTol,
#                      tLag = 0,         tLagTol = 0,
#                      method = NULL)
#   
#   Cr  <- out_space$C[, 1]
#   npr <- out_space$np[, 1]
#   
#   out_time <- stcov(values, coords, t_day, values, coords, t_day,
#                     rLag = 0,         rLagTol = 0,
#                     tLag = lags$tLag, tLagTol = lags$tLagTol,
#                     method = out_space$method)
#   
#   Ct  <- out_time$C[1, ]
#   npt <- out_time$np[1, ]
#   
#   res_space <- trim_after_negative_inclusive(Cr, lags$rLag, np = npr)
#   res_time  <- trim_after_negative_inclusive(Ct, lags$tLag, np = npt)
#   
#   dt_space <- unique(data.table(rLag = res_space$lags, Cov = res_space$cov, np = res_space$np), by = "rLag")
#   dt_time  <- unique(data.table(tLag = res_time$lags,  Cov = res_time$cov,  np = res_time$np),  by = "tLag")
#   
#   tag <- sprintf("precip_CONUS_week_%s", week_id)
#   fwrite(dt_space, file.path(OUT_DIR_WEEK, paste0(tag, "_cov_space.txt")), sep = "\t")
#   fwrite(dt_time,  file.path(OUT_DIR_WEEK, paste0(tag, "_cov_time.txt")),  sep = "\t")
#   
#   t1 <- Sys.time()
#   list(week = week_id, ok = TRUE, reason = NA_character_,
#        dt_sec = as.numeric(difftime(t1, t0, units="secs")),
#        n_space = nrow(dt_space), n_time = nrow(dt_time),
#        method = out_space$method,
#        mem_coords_mb = mem_coords_mb, mem_values_mb = mem_values_mb,
#        n_days = ncol(values))
# }
# 
# run_weekly_precip <- function() {
#   files <- sort(list.files(IN_DIR_DAILY, pattern = "^\\d{4}-\\d{2}-\\d{2}\\.qs$", full.names = TRUE))
#   stopifnot(length(files) > 0)
#   
#   day_dates <- as.Date(gsub("\\.qs$", "", basename(files)))
#   iso_year <- isoyear(day_dates)
#   iso_week <- isoweek(day_dates)
#   week_id  <- sprintf("%04d-W%02d", iso_year, iso_week)
#   
#   map_dt <- data.table(file = files, day = day_dates, week = week_id)
#   setorder(map_dt, week, day)
#   
#   week_list <- split(map_dt$file, map_dt$week)
#   weeks <- names(week_list)
#   
#   tags <- sprintf("precip_CONUS_week_%s", weeks)
#   out_space_paths <- file.path(OUT_DIR_WEEK, paste0(tags, "_cov_space.txt"))
#   out_time_paths  <- file.path(OUT_DIR_WEEK, paste0(tags, "_cov_time.txt"))
#   have_space <- file.exists(out_space_paths)
#   have_time  <- file.exists(out_time_paths)
#   is_done <- have_space & have_time
#   
#   message(sprintf(
#     "[WEEKLY SKIP CHECK] total_weeks=%d | done=%d | todo=%d | space_missing=%d | time_missing=%d",
#     length(weeks), sum(is_done), sum(!is_done), sum(!have_space), sum(!have_time)
#   ))
#   
#   weeks_todo <- weeks[!is_done]
#   
#   handlers(global = TRUE)
#   t_all0 <- Sys.time()
#   
#   if (length(weeks_todo) == 0L) {
#     message("✅ Nothing to do — all WEEKLY outputs already exist.")
#     log_dt <- data.table(
#       week = weeks,
#       ok = TRUE,
#       reason = NA_character_,
#       dt_sec = NA_real_,
#       n_space = NA_integer_,
#       n_time = NA_integer_,
#       method = NA_character_,
#       mem_coords_mb = NA_real_,
#       mem_values_mb = NA_real_,
#       n_days = NA_integer_,
#       had_cov_space = have_space,
#       had_cov_time  = have_time,
#       skipped_already_done = TRUE
#     )
#     fwrite(log_dt, file.path(OUT_DIR_WEEK, "run_log_parallel_relic_weekly_precip.csv"))
#     message(sprintf("✅ WEEKLY DONE | workers=%d | total_elapsed=%s",
#                     WORKERS, format(difftime(Sys.time(), t_all0), digits = 3)))
#     return(invisible(NULL))
#   }
#   
#   with_progress({
#     p <- progressor(along = weeks_todo)
#     results <- future_lapply(weeks_todo, function(wk) {
#       p(wk)
#       process_one_week_precip(wk, week_list[[wk]])
#     })
#   })
#   
#   log_done <- rbindlist(lapply(results, as.data.table), fill = TRUE)
#   log_done[, `:=`(had_cov_space = TRUE, had_cov_time = TRUE, skipped_already_done = FALSE)]
#   
#   if (any(is_done)) {
#     log_skip <- data.table(
#       week = weeks[is_done],
#       ok = TRUE,
#       reason = NA_character_,
#       dt_sec = NA_real_,
#       n_space = NA_integer_,
#       n_time = NA_integer_,
#       method = NA_character_,
#       mem_coords_mb = NA_real_,
#       mem_values_mb = NA_real_,
#       n_days = NA_integer_,
#       had_cov_space = have_space[is_done],
#       had_cov_time  = have_time[is_done],
#       skipped_already_done = TRUE
#     )
#     log_dt <- rbindlist(list(log_done, log_skip), use.names = TRUE, fill = TRUE)
#   } else {
#     log_dt <- log_done
#   }
#   
#   all_weeks_logged <- log_dt$week
#   all_tags_logged  <- sprintf("precip_CONUS_week_%s", all_weeks_logged)
#   log_dt[, `:=`(
#     had_cov_space = file.exists(file.path(OUT_DIR_WEEK, paste0(all_tags_logged, "_cov_space.txt"))),
#     had_cov_time  = file.exists(file.path(OUT_DIR_WEEK, paste0(all_tags_logged, "_cov_time.txt")))
#   )]
#   
#   setorder(log_dt, week)
#   fwrite(log_dt, file.path(OUT_DIR_WEEK, "run_log_parallel_relic_weekly_precip.csv"))
#   
#   message(sprintf("✅ WEEKLY DONE | workers=%d | processed=%d | skipped=%d | total_elapsed=%s",
#                   WORKERS, length(weeks_todo), sum(is_done),
#                   format(difftime(Sys.time(), t_all0), digits = 3)))
# }
# 
# # RUN (choose one or both)----
# 
# # 1) Hourly->Daily covariance per day
# run_daily_precip()
# 
# # Test ----
# # files <- sort(list.files(IN_DIR_DAILY,
# #                          pattern = "^\\d{4}-\\d{2}-\\d{2}\\.qs$",
# #                          full.names = TRUE))
# # 
# # day_dates <- as.Date(gsub("\\.qs$", "", basename(files)))
# # iso_year  <- isoyear(day_dates)
# # iso_week  <- isoweek(day_dates)
# # week_id   <- sprintf("%04d-W%02d", iso_year, iso_week)
# # 
# # map_dt <- data.table(file = files, day = day_dates, week = week_id)
# # setorder(map_dt, week, day)
# # 
# # week_list <- split(map_dt$file, map_dt$week)
# # 
# # plan(sequential)
# # 
# # process_one_week_precip(
# #   "1940-W01",
# #   week_list[["1940-W01"]]
# # )
# 
# # 2) Daily totals->Weekly covariance per ISO week (UTC)
# # run_weekly_precip()

# # Southeast ----
# # scripts/07_experimental_covariance_precip_daily_parallel_relic_SOUTHEAST.R
# # PRECIPITATION experimental covariance (RELIC estimator) — DAILY ONLY
# # - DAILY: uses HOURLY precipitation within each day (t axis in hours)
# # - Uses stcov() + pairsindex() (trusted)
# # - Trims AFTER first negative, but KEEPS the first negative (inclusive)
# # - Parallelize across days with future_lapply + progressr
# # - NEW: Spatial subset to Southeast bbox in lon/lat BEFORE gridding:
# #        BBOX_LONLAT <- c(xmin=-95, xmax=-75, ymin=22, ymax=40)
# # - NEW: Outputs written to .../02_exp_cov/precipitation/southeast
# 
# suppressPackageStartupMessages({
#   library(qs)
#   library(data.table)
#   library(lubridate)
#   library(future)
#   library(future.apply)
#   library(progressr)
#   library(fs)
# })
# 
# # -----------------------------
# # USER SETTINGS
# # -----------------------------
# ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
# 
# IN_DIR_DAILY  <- file.path(ROOT, "data/out/01_vars/precipitation")
# 
# # NEW OUT DIR
# OUT_DIR_DAILY <- file.path(ROOT, "data/out/02_exp_cov/precipitation/southeast")
# dir.create(OUT_DIR_DAILY, recursive = TRUE, showWarnings = FALSE)
# 
# # NEW BBOX (lon/lat)
# BBOX_LONLAT <- c(xmin = -95, xmax = -75, ymin = 22, ymax = 40)
# 
# WORKERS <- 28
# plan(multisession, workers = WORKERS)
# 
# # Performance / safety (controls kron vs kronloop)
# max_pairs <- 2e7  # n_cells^2 cutoff
# 
# # -----------------------------
# # RELIC CORE (UNCHANGED)
# # -----------------------------
# coord2dist <- function(c1, c2) {
#   if (is.data.frame(c1)) c1 <- as.matrix(c1)
#   if (is.data.frame(c2)) c2 <- as.matrix(c2)
#   if (is.list(c1)) { c1 <- c1[[1]]; c2 <- c2[[1]] }
#   x1_sq <- rowSums(c1^2); x2_sq <- rowSums(c2^2)
#   D2 <- outer(x1_sq, rep(1, nrow(c2))) +
#     outer(rep(1, nrow(c1)), x2_sq) -
#     2 * tcrossprod(c1, c2)
#   D2[D2 < 0] <- 0
#   sqrt(D2)
# }
# 
# pairsindex <- function(c1, c2, rLag, rLagTol, dist=list("coord2dist"), method=NULL) {
#   if (is.data.frame(c1)) c1 <- as.matrix(c1)
#   if (is.data.frame(c2)) c2 <- as.matrix(c2)
#   n1 <- nrow(c1); n2 <- nrow(c2); nr <- length(rLag)
#   if (!is.list(dist)) stop("dist must be a list")
#   
#   if (is.null(method)) {
#     method <- if (n1 * n2 <= max_pairs) "kron" else "kronloop"
#   }
#   
#   idxpairs <- vector("list", nr)
#   
#   call_dist <- function(A, B) {
#     if (length(dist) == 1) {
#       do.call(dist[[1]], list(A, B))
#     } else {
#       do.call(dist[[1]], list(A, B, dist[[2]]))
#     }
#   }
#   
#   switch(method,
#          kron = {
#            D <- call_dist(c1, c2)
#            for (ir in seq_len(nr)) {
#              mask <- (D >= rLag[ir] - rLagTol[ir]) & (D <= rLag[ir] + rLagTol[ir])
#              idxpairs[[ir]] <- which(mask, arr.ind = TRUE)
#            }
#          },
#          kronloop = {
#            for (ir in seq_len(nr)) idxpairs[[ir]] <- matrix(numeric(0), ncol = 2)
#            for (i in seq_len(n1)) {
#              d_i <- call_dist(matrix(c1[i,], nrow = 1), c2)
#              for (ir in seq_len(nr)) {
#                js <- which(d_i >= rLag[ir] - rLagTol[ir] & d_i <= rLag[ir] + rLagTol[ir])
#                if (length(js) > 0) idxpairs[[ir]] <- rbind(idxpairs[[ir]], cbind(i, js))
#              }
#            }
#          },
#          stop("method must be 'kron' or 'kronloop'")
#   )
#   
#   for (ir in seq_len(nr)) {
#     if (nrow(idxpairs[[ir]]) > 0) {
#       idxpairs[[ir]] <- idxpairs[[ir]][order(idxpairs[[ir]][,1], idxpairs[[ir]][,2]), , drop = FALSE]
#     }
#   }
#   idxpairs
# }
# 
# stcov <- function(
#     Zi, cMSi, tMEi,
#     Zj, cMSj, tMEj,
#     rLag, rLagTol,
#     tLag, tLagTol,
#     dist = list("coord2dist"),
#     method = NULL
# ) {
#   nMSi <- nrow(Zi);  nMEi <- ncol(Zi)
#   if (!all(dim(cMSi) == c(nMSi, 2))) stop("cMSi must be nMSi × 2")
#   if (length(tMEi) != nMEi) stop("tMEi must be length nMEi")
#   
#   nMSj <- nrow(Zj);  nMEj <- ncol(Zj)
#   if (!all(dim(cMSj) == c(nMSj, 2))) stop("cMSj must be nMSj × 2")
#   if (length(tMEj) != nMEj) stop("tMEj must be length nMEj")
#   
#   nr <- length(rLag); if (length(rLagTol) != nr) stop("rLag and rLagTol length mismatch")
#   nt <- length(tLag); if (length(tLagTol) != nt) stop("tLag and tLagTol length mismatch")
#   
#   if (is.null(method)) {
#     method <- if (nMSi * nMSj <= max_pairs) "kron" else "kronloop"
#   }
#   
#   idxpairsMS <- pairsindex(cMSi, cMSj, rLag, rLagTol, dist, method)
#   tMat <- outer(tMEi, tMEj, function(ti, tj) tj - ti)
#   
#   C  <- matrix(NA_real_, nrow = nr, ncol = nt)
#   np <- matrix(0L,       nrow = nr, ncol = nt)
#   
#   for (ir in seq_len(nr)) {
#     pairs_ij <- idxpairsMS[[ir]]
#     if (nrow(pairs_ij) == 0L) next
#     
#     for (it in seq_len(nt)) {
#       mask <- (tLag[it] - tLagTol[it] <= tMat) & (tMat <= tLag[it] + tLagTol[it])
#       if (!any(mask)) next
#       ev_idx <- which(mask, arr.ind = TRUE)
#       
#       nSP <- nrow(pairs_ij); nEP <- nrow(ev_idx)
#       sp_rep <- pairs_ij[rep(seq_len(nSP), times = nEP), , drop = FALSE]
#       ep_rep <- ev_idx[rep(seq_len(nEP), each = nSP), , drop = FALSE]
#       
#       head_vals <- Zi[cbind(sp_rep[,1], ep_rep[,1])]
#       tail_vals <- Zj[cbind(sp_rep[,2], ep_rep[,2])]
#       
#       valid <- !is.na(head_vals) & !is.na(tail_vals)
#       np[ir,it] <- sum(valid)
#       if (np[ir,it] == 0L) next
#       
#       h <- head_vals[valid]
#       t <- tail_vals[valid]
#       C[ir,it] <- mean(h * t) - mean(h) * mean(t)
#     }
#   }
#   
#   list(C = C, np = np, method = method)
# }
# 
# trim_after_negative_inclusive <- function(x, lags, np = NULL) {
#   keep <- !is.na(x) & !is.na(lags)
#   x <- x[keep]; lags <- lags[keep]
#   if (!is.null(np)) np <- np[keep]
#   
#   if (length(x) == 0L) return(list(cov = numeric(0), lags = numeric(0), np = integer(0)))
#   
#   neg_idx <- which(x < 0)[1]
#   if (!is.na(neg_idx)) {
#     x <- x[1:neg_idx]
#     lags <- lags[1:neg_idx]
#     if (!is.null(np)) np <- np[1:neg_idx]
#   }
#   
#   if (is.null(np)) np <- rep(NA_integer_, length(x))
#   list(cov = x, lags = lags, np = np)
# }
# 
# estimate_base_dist <- function(coords) {
#   ux <- sort(unique(coords[,1]))
#   uy <- sort(unique(coords[,2]))
#   dx <- diff(ux); dy <- diff(uy)
#   dx <- dx[dx > 0]; dy <- dy[dy > 0]
#   base <- suppressWarnings(min(median(dx, na.rm = TRUE), median(dy, na.rm = TRUE), na.rm = TRUE))
#   if (!is.finite(base) || base <= 0) base <- 0.5
#   base
# }
# 
# make_lags <- function(max_dist, max_time_hr, base_dist) {
#   r_raw <- base_dist * c(0:4, 6, 8, 12)
#   rLag <- sort(unique(r_raw[r_raw <= max_dist]))
#   if (length(rLag) < 2) rLag <- c(0, min(base_dist, max_dist))
#   rLagTol <- c(0, diff(rLag)/2)
#   
#   t_raw <- c(0:4, 6, 8, 12)
#   tLag <- sort(unique(t_raw[t_raw <= max_time_hr]))
#   if (length(tLag) < 2) tLag <- c(0, min(1, max_time_hr))
#   tLagTol <- c(0, diff(tLag)/2)
#   
#   list(rLag = rLag, rLagTol = rLagTol, tLag = tLag, tLagTol = tLagTol)
# }
# 
# mb <- function(x) as.numeric(object.size(x)) / (1024^2)
# 
# # -----------------------------
# # PRECIP HELPERS
# # -----------------------------
# detect_precip_col <- function(dt) {
#   cand <- c("intensity", "precip", "precipitation", "prcp", "rain", "tp", "apcp",
#             "precip_mm", "precip_m")
#   hit <- intersect(cand, names(dt))
#   if (length(hit) == 0L) {
#     stop("Could not detect precipitation column. Available names: ",
#          paste(names(dt), collapse = ", "))
#   }
#   hit[1]
# }
# 
# ensure_datetime_utc <- function(dt, qs_file) {
#   if ("datetime" %in% names(dt) && inherits(dt$datetime, "POSIXct")) {
#     attr(dt$datetime, "tzone") <- "UTC"
#     return(dt)
#   }
#   if ("datetime_utc" %in% names(dt)) {
#     dt[, datetime := as.POSIXct(datetime_utc, tz = "UTC")]
#     return(dt)
#   }
#   stop("No datetime or datetime_utc column found in: ", qs_file)
# }
# 
# replace_below_threshold_with_na <- function(x, thresh = 1) {
#   x[!is.na(x) & x < thresh] <- NA_real_
#   x
# }
# 
# # NEW: bbox helper (lon/lat)
# .in_bbox <- function(lon, lat, bb) {
#   lon >= bb["xmin"] & lon <= bb["xmax"] & lat >= bb["ymin"] & lat <= bb["ymax"]
# }
# 
# # -----------------------------
# # DAILY: HOURLY precip within each day + bbox subset
# # -----------------------------
# read_day_to_matrix_precip <- function(qs_file, bbox = BBOX_LONLAT) {
#   dt <- as.data.table(qread(qs_file))
#   dt <- ensure_datetime_utc(dt, qs_file)
#   
#   if (!all(c("lon","lat") %in% names(dt))) stop("Missing lon/lat columns in: ", qs_file)
#   
#   pcol <- detect_precip_col(dt)
#   
#   # Keep only needed cols
#   dt <- dt[, .(lon, lat, datetime, precip = as.numeric(get(pcol)))]
#   n_before <- nrow(dt)
#   if (n_before == 0L) return(NULL)
#   
#   # NEW: bbox subset
#   keep <- .in_bbox(dt$lon, dt$lat, bbox)
#   dt <- dt[keep]
#   n_after <- nrow(dt)
#   if (n_after == 0L) return(NULL)
#   
#   # Replace small values with NA (0/trace skip); then drop NA
#   dt[, precip := replace_below_threshold_with_na(precip, thresh = 0.5)]
#   dt <- dt[!is.na(precip)]
#   if (nrow(dt) == 0L) return(NULL)
#   
#   dt[, cell_id := paste0(sprintf("%.6f", lon), "_", sprintf("%.6f", lat))]
#   
#   times <- sort(unique(dt$datetime))
#   if (length(times) < 2L) return(NULL)
#   
#   dt[, t_idx := match(datetime, times)]
#   
#   setkey(dt, cell_id, t_idx)
#   wide <- dcast(
#     dt,
#     cell_id + lon + lat ~ t_idx,
#     value.var = "precip",
#     fun.aggregate = mean
#   )
#   if (nrow(wide) < 2L) return(NULL)
#   
#   coords <- as.matrix(wide[, .(lon, lat)])
#   
#   time_cols <- setdiff(names(wide), c("cell_id", "lon", "lat"))
#   time_cols_ord <- time_cols[order(as.integer(time_cols))]
#   
#   values <- as.matrix(wide[, ..time_cols_ord])
#   if (ncol(values) < 2L) return(NULL)
#   
#   t_hr <- as.numeric(difftime(times, times[1], units = "hours"))
#   
#   list(
#     coords = coords, values = values, times = times, t_hr = t_hr,
#     n_before = n_before, n_after = n_after
#   )
# }
# 
# process_one_day_precip <- function(qs_file) {
#   day_str <- gsub("\\.qs$", "", basename(qs_file))
#   t0 <- Sys.time()
#   
#   mat <- read_day_to_matrix_precip(qs_file, bbox = BBOX_LONLAT)
#   if (is.null(mat)) {
#     return(list(
#       day = day_str, ok = FALSE, reason = "NULL matrix (after bbox/threshold or insufficient coverage)",
#       dt_sec = as.numeric(difftime(Sys.time(), t0, units="secs")),
#       n_space = 0L, n_time = 0L, method = NA_character_,
#       mem_coords_mb = NA_real_, mem_values_mb = NA_real_,
#       n_pts_before_bbox = NA_integer_, n_pts_after_bbox = NA_integer_, pct_pts_kept = NA_real_
#     ))
#   }
#   
#   coords <- mat$coords
#   values <- mat$values
#   t_hr   <- mat$t_hr
#   
#   n_before <- mat$n_before
#   n_after  <- mat$n_after
#   pct_kept <- if (is.finite(n_before) && n_before > 0) 100 * (n_after / n_before) else NA_real_
#   
#   base_dist <- estimate_base_dist(coords)
#   max_dist  <- max(coord2dist(coords, coords), na.rm = TRUE)
#   max_time  <- max(t_hr, na.rm = TRUE)
#   lags <- make_lags(max_dist, max_time, base_dist)
#   
#   mem_coords_mb <- mb(coords)
#   mem_values_mb <- mb(values)
#   
#   # Space-only
#   out_space <- stcov(values, coords, t_hr, values, coords, t_hr,
#                      rLag = lags$rLag, rLagTol = lags$rLagTol,
#                      tLag = 0,         tLagTol = 0,
#                      method = NULL)
#   Cr  <- out_space$C[, 1]
#   npr <- out_space$np[, 1]
#   
#   # Time-only
#   out_time <- stcov(values, coords, t_hr, values, coords, t_hr,
#                     rLag = 0,         rLagTol = 0,
#                     tLag = lags$tLag, tLagTol = lags$tLagTol,
#                     method = out_space$method)
#   
#   Ct  <- out_time$C[1, ]
#   npt <- out_time$np[1, ]
#   
#   res_space <- trim_after_negative_inclusive(Cr, lags$rLag, np = npr)
#   res_time  <- trim_after_negative_inclusive(Ct, lags$tLag, np = npt)
#   
#   dt_space <- unique(data.table(rLag = res_space$lags, Cov = res_space$cov, np = res_space$np), by = "rLag")
#   dt_time  <- unique(data.table(tLag = res_time$lags,  Cov = res_time$cov,  np = res_time$np),  by = "tLag")
#   
#   # NEW TAG: southeast
#   tag <- sprintf("precip_southeast_day_%s", day_str)
#   
#   fwrite(dt_space, file.path(OUT_DIR_DAILY, paste0(tag, "_cov_space.txt")), sep = "\t")
#   fwrite(dt_time,  file.path(OUT_DIR_DAILY, paste0(tag, "_cov_time.txt")),  sep = "\t")
#   
#   t1 <- Sys.time()
#   list(
#     day = day_str, ok = TRUE, reason = NA_character_,
#     dt_sec = as.numeric(difftime(t1, t0, units="secs")),
#     n_space = nrow(dt_space), n_time = nrow(dt_time),
#     method = out_space$method,
#     mem_coords_mb = mem_coords_mb, mem_values_mb = mem_values_mb,
#     n_pts_before_bbox = n_before, n_pts_after_bbox = n_after, pct_pts_kept = pct_kept
#   )
# }
# 
# # -----------------------------
# # RUN (restart-safe)
# # -----------------------------
# files <- sort(list.files(IN_DIR_DAILY, pattern = "^\\d{4}-\\d{2}-\\d{2}\\.qs$", full.names = TRUE))
# stopifnot(length(files) > 0)
# 
# day_strs <- gsub("\\.qs$", "", basename(files))
# 
# # NEW TAGS: southeast
# tags <- sprintf("precip_southeast_day_%s", day_strs)
# 
# out_space_paths <- file.path(OUT_DIR_DAILY, paste0(tags, "_cov_space.txt"))
# out_time_paths  <- file.path(OUT_DIR_DAILY, paste0(tags, "_cov_time.txt"))
# 
# have_space <- file.exists(out_space_paths)
# have_time  <- file.exists(out_time_paths)
# is_done <- have_space & have_time
# 
# message(sprintf(
#   "[DAILY SKIP CHECK] total_days=%d | done=%d | todo=%d | space_missing=%d | time_missing=%d",
#   length(files), sum(is_done), sum(!is_done), sum(!have_space), sum(!have_time)
# ))
# 
# files_todo <- files[!is_done]
# 
# handlers(global = TRUE)
# t_all0 <- Sys.time()
# 
# if (length(files_todo) == 0L) {
#   message("✅ Nothing to do — all DAILY outputs already exist.")
#   
#   log_dt <- data.table(
#     day = day_strs,
#     ok = TRUE,
#     reason = NA_character_,
#     dt_sec = NA_real_,
#     n_space = NA_integer_,
#     n_time = NA_integer_,
#     method = NA_character_,
#     mem_coords_mb = NA_real_,
#     mem_values_mb = NA_real_,
#     n_pts_before_bbox = NA_integer_,
#     n_pts_after_bbox  = NA_integer_,
#     pct_pts_kept      = NA_real_,
#     had_cov_space = have_space,
#     had_cov_time  = have_time,
#     skipped_already_done = TRUE
#   )
#   fwrite(log_dt, file.path(OUT_DIR_DAILY, "run_log_parallel_relic_daily_precip_southeast.csv"))
#   
#   message(sprintf("✅ DAILY DONE | workers=%d | total_elapsed=%s",
#                   WORKERS, format(difftime(Sys.time(), t_all0), digits = 3)))
# } else {
#   
#   with_progress({
#     p <- progressor(along = files_todo)
#     results <- future_lapply(files_todo, function(f) {
#       p(basename(f))
#       process_one_day_precip(f)
#     })
#   })
#   
#   log_done <- rbindlist(lapply(results, as.data.table), fill = TRUE)
#   log_done[, `:=`(had_cov_space = TRUE, had_cov_time = TRUE, skipped_already_done = FALSE)]
#   
#   if (any(is_done)) {
#     log_skip <- data.table(
#       day = day_strs[is_done],
#       ok = TRUE,
#       reason = NA_character_,
#       dt_sec = NA_real_,
#       n_space = NA_integer_,
#       n_time = NA_integer_,
#       method = NA_character_,
#       mem_coords_mb = NA_real_,
#       mem_values_mb = NA_real_,
#       n_pts_before_bbox = NA_integer_,
#       n_pts_after_bbox  = NA_integer_,
#       pct_pts_kept      = NA_real_,
#       had_cov_space = have_space[is_done],
#       had_cov_time  = have_time[is_done],
#       skipped_already_done = TRUE
#     )
#     log_dt <- rbindlist(list(log_done, log_skip), use.names = TRUE, fill = TRUE)
#   } else {
#     log_dt <- log_done
#   }
#   
#   # Final filesystem truth check
#   all_days_logged <- log_dt$day
#   all_tags_logged <- sprintf("precip_southeast_day_%s", all_days_logged)
#   log_dt[, `:=`(
#     had_cov_space = file.exists(file.path(OUT_DIR_DAILY, paste0(all_tags_logged, "_cov_space.txt"))),
#     had_cov_time  = file.exists(file.path(OUT_DIR_DAILY, paste0(all_tags_logged, "_cov_time.txt")))
#   )]
#   
#   setorder(log_dt, day)
#   fwrite(log_dt, file.path(OUT_DIR_DAILY, "run_log_parallel_relic_daily_precip_southeast.csv"))
#   
#   message(sprintf("✅ DAILY DONE | workers=%d | processed=%d | skipped=%d | total_elapsed=%s",
#                   WORKERS, length(files_todo), sum(is_done),
#                   format(difftime(Sys.time(), t_all0), digits = 3)))
# }

# Southeast 2 ----
# scripts/07_experimental_covariance_precip_daily_parallel_relic_SOUTHEAST.R
# PRECIPITATION experimental covariance (RELIC estimator) — DAILY ONLY
#
# Updates included (mirrors heatindex script):
# 1) Build a LOOKUP TABLE of valid cells (lon/lat) that satisfy:
#    - land mask (REF_MASK_TIF) AND
#    - bbox (BBOX_LONLAT)
#    Lookup is built once in the main process and saved to OUT_DIR for reuse.
# 2) Workers ONLY use the lookup's allowed cell_ids (NO terra objects in workers).
# 3) Robust terra::extract handling (does NOT assume column 2 exists).
# 4) Overwrite controls (OVERWRITE_MODE): "none" | "partial" | "all"
# 5) Atomic writes.
# 6) DAILY ONLY: uses HOURLY precipitation within each day (t axis in hours)
# 7) Applies precip thresholding (values < PRECIP_MIN_MM become NA and are dropped)

suppressPackageStartupMessages({
  library(qs)
  library(data.table)
  library(lubridate)
  library(future)
  library(future.apply)
  library(progressr)
  library(fs)
  library(terra)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

IN_DIR_DAILY  <- file.path(ROOT, "data/out/01_vars/precipitation")

OUT_DIR_DAILY <- file.path(ROOT, "data/out/02_exp_cov/precipitation/southeast")
dir.create(OUT_DIR_DAILY, recursive = TRUE, showWarnings = FALSE)

# Southeast bbox (lon/lat)
BBOX_LONLAT <- c(xmin = -95, xmax = -75, ymin = 22, ymax = 40)

# Land mask tif (lon/lat)
REF_MASK_TIF <- file.path(ROOT, "data/out/01_vars/precipitation/_mask_target.tif")
if (!file.exists(REF_MASK_TIF)) stop("Missing REF_MASK_TIF: ", REF_MASK_TIF)

# Lookup saved to disk so subsequent runs are fast
LOOKUP_QS <- file.path(OUT_DIR_DAILY, "lookup_valid_cells_bbox_land.qs")
REBUILD_LOOKUP <- FALSE
LOOKUP_SEED_FILE <- NULL  # if NULL, uses first file in IN_DIR_DAILY

WORKERS <- 28
options(future.globals.maxSize = 8 * 1024^3)
plan(multisession, workers = WORKERS)

# Performance / safety
max_pairs <- 2e7  # controls kron vs kronloop (n_cells^2 cutoff)

# Overwrite behavior
OVERWRITE_MODE <- "all"  # "none" | "partial" | "all"

# Precip filtering (mm per hour, assuming your hourly precip units are mm)
PRECIP_MIN_MM <- 0.5

# -----------------------------
# Utility helpers
# -----------------------------
.in_bbox <- function(lon, lat, bb) {
  lon >= bb["xmin"] & lon <= bb["xmax"] & lat >= bb["ymin"] & lat <= bb["ymax"]
}

.atomic_fwrite <- function(dt, out_path, ...) {
  tmp <- paste0(out_path, ".tmp_", Sys.getpid(), "_", format(Sys.time(), "%Y%m%d%H%M%OS3"))
  data.table::fwrite(dt, tmp, ...)
  file.rename(tmp, out_path)
}

.extract_vals_robust <- function(r, xy) {
  ex <- terra::extract(r, xy)
  if (is.null(ex)) return(rep(NA_real_, nrow(xy)))
  
  if (is.atomic(ex) && is.vector(ex)) {
    return(as.numeric(ex))
  }
  
  ex <- as.data.frame(ex)
  
  if (ncol(ex) == 0L) return(rep(NA_real_, nrow(xy)))
  if (ncol(ex) == 1L) return(as.numeric(ex[[1]]))
  
  as.numeric(ex[[2]])
}

# -----------------------------
# RELIC CORE (UNCHANGED)
# -----------------------------
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

# -----------------------------
# PRECIP HELPERS
# -----------------------------
detect_precip_col <- function(dt) {
  cand <- c("intensity", "precip", "precipitation", "prcp", "rain", "tp", "apcp",
            "precip_mm", "precip_m")
  hit <- intersect(cand, names(dt))
  if (length(hit) == 0L) {
    stop("Could not detect precipitation column. Available names: ",
         paste(names(dt), collapse = ", "))
  }
  hit[1]
}

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
  x[!is.na(x) & x < thresh] <- NA_real_
  x
}

# -----------------------------
# LOOKUP TABLE (same logic as heatindex)
# -----------------------------
.build_lookup <- function(seed_qs, mask_tif, bbox, out_qs) {
  message("[LOOKUP] Building allowed-cell lookup from: ", seed_qs)
  
  dt <- as.data.table(qread(seed_qs))
  if (!all(c("lon","lat") %in% names(dt))) stop("[LOOKUP] Missing lon/lat columns in: ", seed_qs)
  
  grid <- unique(dt[, .(lon, lat)])
  setDT(grid)
  message(sprintf("[LOOKUP] Seed unique coords: %d", nrow(grid)))
  
  grid <- grid[.in_bbox(lon, lat, bbox)]
  setDT(grid)
  message(sprintf("[LOOKUP] After bbox: %d", nrow(grid)))
  if (nrow(grid) == 0L) stop("[LOOKUP] No cells remain after bbox.")
  
  mask_r <- terra::rast(mask_tif)
  xy <- as.matrix(grid[, c("lon","lat"), drop = FALSE])
  v  <- .extract_vals_robust(mask_r, xy)
  
  keep <- is.finite(v) & (v > 0)
  grid <- grid[keep]
  message(sprintf("[LOOKUP] After land mask: %d", nrow(grid)))
  if (nrow(grid) == 0L) stop("[LOOKUP] No cells remain after land mask.")
  
  grid[, cell_id := paste0(sprintf("%.6f", lon), "_", sprintf("%.6f", lat))]
  
  qs::qsave(grid, out_qs)
  message("[LOOKUP] Saved: ", out_qs)
  
  invisible(grid)
}

.get_allowed_cells <- function(files, lookup_seed_file, lookup_qs, rebuild = FALSE) {
  if (!rebuild && file.exists(lookup_qs)) {
    message("[LOOKUP] Using existing lookup: ", lookup_qs)
    grid <- as.data.table(qs::qread(lookup_qs))
    if (!all(c("lon","lat","cell_id") %in% names(grid))) {
      stop("[LOOKUP] Lookup missing required columns lon/lat/cell_id: ", lookup_qs)
    }
    return(grid)
  }
  
  seed <- lookup_seed_file
  if (is.null(seed)) {
    if (length(files) == 0L) stop("[LOOKUP] No files available to build lookup.")
    seed <- files[1]
  }
  .build_lookup(seed, REF_MASK_TIF, BBOX_LONLAT, lookup_qs)
  as.data.table(qs::qread(lookup_qs))
}

# -----------------------------
# DAILY: HOURLY precip within day -> matrix, filtered by allowed cells
# -----------------------------
read_day_to_matrix_precip <- function(qs_file, allowed_cell_ids) {
  dt <- as.data.table(qread(qs_file))
  dt <- ensure_datetime_utc(dt, qs_file)
  
  if (!all(c("lon","lat") %in% names(dt))) stop("Missing lon/lat columns in: ", qs_file)
  
  pcol <- detect_precip_col(dt)
  dt <- dt[, .(lon, lat, datetime, precip = as.numeric(get(pcol)))]
  n_before <- nrow(dt)
  if (n_before == 0L) return(NULL)
  
  # stable cell_id then filter by allowed lookup
  dt[, cell_id := paste0(sprintf("%.6f", lon), "_", sprintf("%.6f", lat))]
  dt <- dt[cell_id %in% allowed_cell_ids]
  n_after_cells <- nrow(dt)
  if (n_after_cells == 0L) return(NULL)
  
  # threshold + drop
  dt[, precip := replace_below_threshold_with_na(precip, thresh = PRECIP_MIN_MM)]
  dt <- dt[!is.na(precip)]
  n_after_thresh <- nrow(dt)
  if (n_after_thresh == 0L) return(NULL)
  
  times <- sort(unique(dt$datetime))
  if (length(times) < 2L) return(NULL)
  
  dt[, t_idx := match(datetime, times)]
  
  setkey(dt, cell_id, t_idx)
  wide <- dcast(
    dt,
    cell_id + lon + lat ~ t_idx,
    value.var = "precip",
    fun.aggregate = mean
  )
  if (nrow(wide) < 2L) return(NULL)
  
  coords <- as.matrix(wide[, .(lon, lat)])
  
  time_cols <- setdiff(names(wide), c("cell_id", "lon", "lat"))
  time_cols_ord <- time_cols[order(as.integer(time_cols))]
  
  values <- as.matrix(wide[, ..time_cols_ord])
  if (ncol(values) < 2L) return(NULL)
  
  t_hr <- as.numeric(difftime(times, times[1], units = "hours"))
  
  list(
    coords = coords, values = values, times = times, t_hr = t_hr,
    n_before = n_before,
    n_after_cells = n_after_cells,
    n_after_thresh = n_after_thresh
  )
}

# -----------------------------
# process one day
# -----------------------------
process_one_day_precip <- function(qs_file, allowed_cell_ids, overwrite_mode = OVERWRITE_MODE) {
  stopifnot(overwrite_mode %in% c("none", "partial", "all"))
  
  day_str <- gsub("\\.qs$", "", basename(qs_file))
  t0 <- Sys.time()
  
  tag <- sprintf("precip_southeast_day_%s", day_str)
  out_space_file <- file.path(OUT_DIR_DAILY, paste0(tag, "_cov_space.txt"))
  out_time_file  <- file.path(OUT_DIR_DAILY, paste0(tag, "_cov_time.txt"))
  
  if (overwrite_mode == "none") {
    if (file.exists(out_space_file) || file.exists(out_time_file)) {
      return(list(day=day_str, ok=TRUE, reason="SKIP (OVERWRITE_MODE=none; output exists)",
                  dt_sec=as.numeric(difftime(Sys.time(), t0, units="secs")),
                  wrote=FALSE))
    }
  } else if (overwrite_mode == "partial") {
    if (file.exists(out_space_file) && file.exists(out_time_file)) {
      return(list(day=day_str, ok=TRUE, reason="SKIP (already done)",
                  dt_sec=as.numeric(difftime(Sys.time(), t0, units="secs")),
                  wrote=FALSE))
    }
  }
  
  mat <- read_day_to_matrix_precip(qs_file, allowed_cell_ids = allowed_cell_ids)
  if (is.null(mat)) {
    return(list(
      day = day_str, ok = FALSE, reason = "NULL matrix (after allowed-cells+threshold filter or insufficient coverage)",
      dt_sec = as.numeric(difftime(Sys.time(), t0, units="secs")),
      n_space = 0L, n_time = 0L, method = NA_character_,
      mem_coords_mb = NA_real_, mem_values_mb = NA_real_,
      n_pts_before = NA_integer_, n_pts_after_cells = NA_integer_, n_pts_after_thresh = NA_integer_,
      wrote = FALSE
    ))
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
  
  out_space <- stcov(values, coords, t_hr, values, coords, t_hr,
                     rLag = lags$rLag, rLagTol = lags$rLagTol,
                     tLag = 0,         tLagTol = 0,
                     method = NULL)
  Cr  <- out_space$C[, 1]
  npr <- out_space$np[, 1]
  
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
  
  .atomic_fwrite(dt_space, out_space_file, sep = "\t")
  .atomic_fwrite(dt_time,  out_time_file,  sep = "\t")
  
  t1 <- Sys.time()
  list(
    day = day_str, ok = TRUE, reason = NA_character_,
    dt_sec = as.numeric(difftime(t1, t0, units="secs")),
    n_space = nrow(dt_space), n_time = nrow(dt_time),
    method = out_space$method,
    mem_coords_mb = mem_coords_mb, mem_values_mb = mem_values_mb,
    n_pts_before = mat$n_before,
    n_pts_after_cells = mat$n_after_cells,
    n_pts_after_thresh = mat$n_after_thresh,
    wrote = TRUE
  )
}

# -----------------------------
# RUN
# -----------------------------
files <- sort(list.files(IN_DIR_DAILY, pattern = "^\\d{4}-\\d{2}-\\d{2}\\.qs$", full.names = TRUE))
stopifnot(length(files) > 0)

lookup_dt <- .get_allowed_cells(
  files = files,
  lookup_seed_file = LOOKUP_SEED_FILE,
  lookup_qs = LOOKUP_QS,
  rebuild = REBUILD_LOOKUP
)
allowed_cell_ids <- lookup_dt$cell_id
message(sprintf("[LOOKUP] allowed_cells=%d", length(allowed_cell_ids)))

day_strs <- gsub("\\.qs$", "", basename(files))
tags <- sprintf("precip_southeast_day_%s", day_strs)

out_space_paths <- file.path(OUT_DIR_DAILY, paste0(tags, "_cov_space.txt"))
out_time_paths  <- file.path(OUT_DIR_DAILY, paste0(tags, "_cov_time.txt"))
have_space <- file.exists(out_space_paths)
have_time  <- file.exists(out_time_paths)

if (OVERWRITE_MODE == "all") {
  files_todo <- files
} else if (OVERWRITE_MODE == "partial") {
  is_done <- have_space & have_time
  files_todo <- files[!is_done]
} else if (OVERWRITE_MODE == "none") {
  is_clean <- !have_space & !have_time
  files_todo <- files[is_clean]
} else {
  stop("Invalid OVERWRITE_MODE: ", OVERWRITE_MODE)
}

message(sprintf(
  "[MODE=%s] total_days=%d | todo=%d | have_both=%d | have_neither=%d",
  OVERWRITE_MODE, length(files), length(files_todo),
  sum(have_space & have_time), sum(!have_space & !have_time)
))

handlers(global = TRUE)
t_all0 <- Sys.time()

if (length(files_todo) == 0L) {
  message("✅ Nothing to do — per OVERWRITE_MODE, all required outputs already exist.")
  log_dt <- data.table(day = day_strs, had_cov_space = have_space, had_cov_time = have_time)
  fwrite(log_dt, file.path(OUT_DIR_DAILY, "run_log_parallel_relic_daily_precip_southeast.csv"))
  message(sprintf("✅ DAILY DONE | workers=%d | total_elapsed=%s",
                  WORKERS, format(difftime(Sys.time(), t_all0), digits = 3)))
} else {
  
  with_progress({
    p <- progressor(along = files_todo)
    results <- future_lapply(files_todo, function(f) {
      p(basename(f))
      process_one_day_precip(f, allowed_cell_ids = allowed_cell_ids, overwrite_mode = OVERWRITE_MODE)
    })
  })
  
  t_all1 <- Sys.time()
  log_dt <- rbindlist(lapply(results, as.data.table), fill = TRUE)
  
  all_days_logged <- log_dt$day
  all_tags_logged <- sprintf("precip_southeast_day_%s", all_days_logged)
  log_dt[, `:=`(
    had_cov_space = file.exists(file.path(OUT_DIR_DAILY, paste0(all_tags_logged, "_cov_space.txt"))),
    had_cov_time  = file.exists(file.path(OUT_DIR_DAILY, paste0(all_tags_logged, "_cov_time.txt")))
  )]
  
  setorder(log_dt, day)
  fwrite(log_dt, file.path(OUT_DIR_DAILY, "run_log_parallel_relic_daily_precip_southeast.csv"))
  
  message(sprintf(
    "✅ DAILY DONE | workers=%d | processed=%d | total_elapsed=%s | mode=%s",
    WORKERS, length(files_todo),
    format(difftime(t_all1, t_all0), digits = 3),
    OVERWRITE_MODE
  ))
}