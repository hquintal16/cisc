# # =============================================================================
# # 24-hour precipitation return periods (daily SUM then annual maxima)
# #
# # Stage A: yearly matrices of DAILY SUM precipitation on CONUS ref grid
# # Stage B: per-cell distribution fit on ANNUAL MAX of daily sums (24-hr totals)
# # Stage C: compile NetCDF on CONUS ref grid
# # =============================================================================
# 
# suppressPackageStartupMessages({
#   library(qs)
#   library(data.table)
#   library(fs)
#   library(terra)
#   library(float)
#   library(stringr)
#   library(future)
#   library(future.apply)
#   library(progressr)
# })
# 
# # ----------------------------
# # USER SETTINGS
# # ----------------------------
# ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
# 
# # Daily files (one per day) containing hourly/subdaily precip samples
# IN_PR_DIR <- file.path(ROOT, "data/out/01_vars/precipitation")
# 
# # CONUS ref grid (same idea as your HI fix)
# # REF_NC_CONUS <- "V:/users/hquintal/phd1_cluster_southeast/data/input/temperature/1940_01_era5_temperature_conus.nc"
# REF_NC_CONUS <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation/_mask_target.tif"
# 
# OUT_DIR <- file.path(ROOT, "data/out/04_rp")
# OUT_TS_DIR <- file.path(OUT_DIR, "precip_dailySUM_timeseries_by_year")
# OUT_RP_CHUNKS_DIR <- file.path(OUT_DIR, "precip_idf_pr24_return_periods_chunks")
# OUT_NC_CONUS <- file.path(OUT_DIR, "precip_idf_pr24_return_periods.nc")
# dir_create(OUT_TS_DIR); dir_create(OUT_RP_CHUNKS_DIR)
# 
# LOOKUP_QS <- file.path(OUT_DIR, "precip_refgrid_lonlat_to_cell.qs")
# 
# # Controls
# RUN_STAGE_A <- TRUE
# RUN_STAGE_B <- TRUE
# RUN_STAGE_C <- TRUE
# 
# OVERWRITE_LOOKUP <- TRUE
# OVERWRITE_STAGE_A <- TRUE
# OVERWRITE_STAGE_B <- TRUE
# OVERWRITE_STAGE_C <- TRUE
# 
# WORKERS_A <- 18L
# WORKERS_B <- 18L
# CHUNK_SIZE <- 300L
# 
# QS_PRESET <- "high"
# KEY_DIGITS <- 6
# kfmt <- function(x, digits = KEY_DIGITS) sprintf(paste0("%.", digits, "f"), x)
# 
# # Column name in your daily precip .qs
# PR_COL <- "intensity"  # <-- change if your column is named differently (e.g., "precip_mm", "pr", "precipitation")
# test <- qs::qread('V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation/1940-01-01.qs')
# # test
# # Optional: treat daily totals < 1 mm as NA (your earlier preference)
# DROP_DAILY_TOTAL_LT_MM <- NA_real_  # set to NA_real_ / NULL to disable
# 
# # Return periods
# RP_LABELS  <- c(1, 2, 5, 10, 25, 50, 100, 200, 500, 1000)
# RP_COMPUTE <- c(1.01, 2, 5, 10, 25, 50, 100, 200, 500, 1000)
# 
# # Fit settings (same as your HI workflow)
# MIN_YEARS_REQUIRED <- 20L
# PVAL_THRESHOLD     <- 0.4
# CANDIDATE_DISTS    <- c("cosine", "hypsecant", "logistic", "norm", "uniform")
# 
# # ----------------------------
# # Helpers: file parsing + atomic writes
# # ----------------------------
# is_daily_qs <- function(p) grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}\\.qs$", fs::path_file(p))
# day_from_file_chr <- function(p) sub("\\.qs$", "", fs::path_file(p))
# 
# cleanup_tmp_files <- function(final_path) {
#   patt <- paste0(fs::path_file(final_path), ".tmp_*")
#   junk <- fs::dir_ls(dirname(final_path), glob = patt, fail = FALSE, type="file")
#   if (length(junk) > 0) fs::file_delete(junk)
# }
# atomic_qsave <- function(obj, final_path, preset = "high") {
#   tmp <- paste0(final_path, ".tmp_", Sys.getpid())
#   qs::qsave(obj, tmp, preset = preset)
#   fs::file_move(tmp, final_path)
#   invisible(TRUE)
# }
# 
# # ----------------------------
# # Lookup: CONUS ref grid (cell index + lon/lat + join keys)
# # ----------------------------
# build_lookup_from_conus_ref <- function() {
#   if (!OVERWRITE_LOOKUP && fs::file_exists(LOOKUP_QS)) {
#     grid <- as.data.table(qread(LOOKUP_QS))
#     setorder(grid, cell); setkey(grid, lon_key, lat_key)
#     return(grid)
#   }
#   r <- terra::rast(REF_NC_CONUS)[[1]]
#   v <- terra::values(r, mat = FALSE)
#   cells <- which(!is.na(v))
#   if (length(cells) == 0) stop("REF_NC_CONUS has 0 non-NA cells.")
#   xy <- terra::xyFromCell(r, cells)
#   grid <- data.table(
#     cell = as.integer(cells),
#     lon  = as.numeric(xy[,1]),
#     lat  = as.numeric(xy[,2])
#   )
#   grid[, `:=`(lon_key = kfmt(lon), lat_key = kfmt(lat))]
#   setorder(grid, cell); setkey(grid, lon_key, lat_key)
#   qsave(grid, LOOKUP_QS, preset = QS_PRESET)
#   grid
# }
# 
# # ----------------------------
# # Stage A: DAILY SUM precipitation vectors
# # ----------------------------
# daily_sum_vec_one_day <- function(day_file, grid) {
#   dt <- as.data.table(qread(day_file))
#   
#   req <- c("lon","lat",PR_COL)
#   miss <- setdiff(req, names(dt))
#   if (length(miss) > 0) stop("Missing columns in ", fs::path_file(day_file), ": ", paste(miss, collapse=", "))
#   
#   dt[, `:=`(
#     lon = as.numeric(lon),
#     lat = as.numeric(lat)
#   )]
#   dt[, (PR_COL) := as.numeric(get(PR_COL))]
#   
#   # keep only finite precip
#   dt <- dt[is.finite(get(PR_COL))]
#   if (nrow(dt) == 0) return(rep(NA_real_, nrow(grid)))
#   
#   # daily sum per location (sum across hours/subdaily rows)
#   ds <- dt[, .(pr_sum = sum(get(PR_COL), na.rm = TRUE)), by = .(lon, lat)]
#   
#   # optional: drop tiny daily totals
#   if (is.finite(DROP_DAILY_TOTAL_LT_MM)) {
#     ds[pr_sum < DROP_DAILY_TOTAL_LT_MM, pr_sum := NA_real_]
#   }
#   
#   ds[, `:=`(lon_key = kfmt(lon), lat_key = kfmt(lat))]
#   setkey(ds, lon_key, lat_key)
#   
#   # align to grid order
#   m <- ds[grid, on = .(lon_key, lat_key)]
#   as.numeric(m$pr_sum)
# }
# 
# stageA_process_one_year <- function(y, files_year, days_year, grid) {
#   out_qs <- file.path(OUT_TS_DIR, sprintf("pr_dailySUM_timeseries_%d.qs", y))
#   if (!OVERWRITE_STAGE_A && fs::file_exists(out_qs)) return(list(year=y, status="skip_exists", file=out_qs))
#   cleanup_tmp_files(out_qs)
#   
#   n_cells <- nrow(grid); n_days <- length(files_year)
#   mat <- matrix(NA_real_, nrow = n_cells, ncol = n_days)
#   
#   for (j in seq_len(n_days)) {
#     v <- daily_sum_vec_one_day(files_year[j], grid)
#     if (length(v) != n_cells) stop("Day vec length mismatch for ", files_year[j])
#     mat[, j] <- v
#     if (j %% 25 == 0) gc()
#   }
#   
#   obj <- list(
#     meta = list(
#       variable="precipitation",
#       statistic="daily_sum",
#       units="mm",
#       year=y,
#       created_utc = format(Sys.time(), tz="UTC", usetz=TRUE),
#       ref_nc=REF_NC_CONUS,
#       pr_col=PR_COL,
#       drop_daily_total_lt_mm=DROP_DAILY_TOTAL_LT_MM,
#       key_digits=KEY_DIGITS
#     ),
#     grid = grid[, .(cell, lon, lat)],
#     days = as.Date(days_year),
#     pr_sum = float::fl(mat)
#   )
#   atomic_qsave(obj, out_qs, preset=QS_PRESET)
#   rm(obj, mat); gc()
#   list(year=y, status="ok", file=out_qs)
# }
# 
# run_stageA <- function(grid) {
#   files <- sort(as.character(fs::dir_ls(IN_PR_DIR, type="file", regexp="\\.qs$", fail=FALSE)))
#   files <- files[is_daily_qs(files)]
#   if (length(files) == 0) stop("No daily precip .qs files in: ", IN_PR_DIR)
#   
#   day_chr <- vapply(files, day_from_file_chr, character(1))
#   days <- as.Date(day_chr)
#   o <- order(days); files <- files[o]; days <- days[o]
#   
#   years <- as.integer(format(days, "%Y"))
#   by_year <- split(seq_along(files), years)
#   years_u <- sort(unique(years))
#   
#   # resume scan
#   existing <- fs::dir_ls(OUT_TS_DIR, type="file", regexp="^pr_dailySUM_timeseries_[0-9]{4}\\.qs$", fail=FALSE)
#   done_years <- integer(0)
#   if (length(existing) > 0) {
#     done_years <- as.integer(str_match(fs::path_file(existing), "^pr_dailySUM_timeseries_([0-9]{4})\\.qs$")[,2])
#     done_years <- done_years[is.finite(done_years)]
#   }
#   todo <- years_u
#   if (!OVERWRITE_STAGE_A) todo <- setdiff(todo, done_years)
#   
#   message(sprintf("[Stage A] years total=%d | done=%d | todo=%d | n_cells(CONUS)=%d",
#                   length(years_u), length(done_years), length(todo), nrow(grid)))
#   if (length(todo) == 0) return(data.table(year=integer(), status=character(), file=character(), error=character()))
#   
#   # ---- DEBUG FIRST: sequential is much easier to diagnose on V: ----
#   if (WORKERS_A <= 1L) {
#     future::plan(sequential)
#     res <- lapply(todo, function(y) {
#       idx <- by_year[[as.character(y)]]
#       tryCatch(
#         stageA_process_one_year(y, files[idx], days[idx], grid),
#         error = function(e) list(year=y, status="error", file=file.path(OUT_TS_DIR, sprintf("pr_dailySUM_timeseries_%d.qs", y)),
#                                  error=conditionMessage(e))
#       )
#     })
#   } else {
#     future::plan(multisession, workers=WORKERS_A)
#     
#     # NOTE: progressr + multisession on network drives can be noisy; keep it simple
#     res <- future.apply::future_lapply(todo, function(y) {
#       idx <- by_year[[as.character(y)]]
#       tryCatch(
#         stageA_process_one_year(y, files[idx], days[idx], grid),
#         error = function(e) list(year=y, status="error", file=file.path(OUT_TS_DIR, sprintf("pr_dailySUM_timeseries_%d.qs", y)),
#                                  error=conditionMessage(e))
#       )
#     }, future.seed = TRUE)
#   }
#   
#   res_dt <- data.table::rbindlist(res, fill=TRUE)
#   message("\n[Stage A] Summary:")
#   print(res_dt[, .N, by=status][order(-N)])
#   
#   if (res_dt[status == "error", .N] > 0) {
#     message("\n[Stage A] First error examples:")
#     print(res_dt[status=="error", .(year, error)][order(year)][1:min(10, .N)])
#   }
#   
#   invisible(res_dt)
# }
# # ----------------------------
# # Stage B: fit distribution on ANNUAL MAX of daily sums (pr24)
# # ----------------------------
# dcosine_sc <- function(x, mu, s, log = FALSE) {
#   y <- (x - mu) / s
#   dens <- ifelse(abs(y) <= 1, 1 / (2 * s) * (1 + cos(pi * y)), 0)
#   if (log) log(dens) else dens
# }
# pcosine_sc <- function(q, mu, s) {
#   y <- (q - mu) / s
#   ifelse(y <= -1, 0, ifelse(y >= 1, 1, 0.5*(1 + y + sin(pi*y)/pi)))
# }
# qcosine_sc <- function(p, mu, s) {
#   f <- function(y, target) 0.5*(1 + y + sin(pi*y)/pi) - target
#   y <- vapply(p, function(pp) uniroot(f, interval=c(-1,1), target=pp)$root, numeric(1))
#   mu + s*y
# }
# dhypsec_sc <- function(x, mu, sigma, log = FALSE) {
#   z <- (pi/2)*((x-mu)/sigma)
#   dens <- (1/(2*sigma))*(1/cosh(z))
#   if (log) log(dens) else dens
# }
# phypsec_sc <- function(q, mu, sigma) {
#   z <- (pi/2)*((q-mu)/sigma)
#   (2/pi)*atan(exp(z))
# }
# qhypsec_sc <- function(p, mu, sigma) {
#   mu + sigma*(2/pi)*log(tan(pi*p/2))
# }
# 
# fit_dist_mle <- function(x, dist_name) {
#   x <- x[is.finite(x)]
#   if (length(x) < MIN_YEARS_REQUIRED) stop("Too few years.")
#   mu0 <- mean(x); sd0 <- sd(x); if (!is.finite(sd0) || sd0<=0) sd0 <- 1
#   eps <- 1e-12
#   nll <- switch(
#     dist_name,
#     "norm" = function(par) -sum(dnorm(x, mean=par[1], sd=exp(par[2]), log=TRUE)),
#     "logistic" = function(par) -sum(dlogis(x, location=par[1], scale=exp(par[2]), log=TRUE)),
#     "uniform" = function(par) { a<-par[1]; b<-par[2]; if (b<=a) return(Inf); if (any(x<a|x>b)) return(Inf); -sum(dunif(x,a,b,log=TRUE)) },
#     "cosine" = function(par) { mu<-par[1]; s<-exp(par[2]); if (s<=0) return(Inf); if (max(abs((x-mu)/s))>1) return(Inf); d<-pmax(dcosine_sc(x,mu,s,FALSE),eps); -sum(log(d)) },
#     "hypsecant" = function(par) { mu<-par[1]; sg<-exp(par[2]); if (sg<=0) return(Inf); d<-pmax(dhypsec_sc(x,mu,sg,FALSE),eps); -sum(log(d)) },
#     stop("Unknown dist")
#   )
#   start <- switch(
#     dist_name,
#     "norm"=c(mu0,log(sd0)),
#     "logistic"=c(mu0,log(sd0)),
#     "hypsecant"=c(mu0,log(sd0)),
#     "cosine"={ s0<-max(abs(x-mu0)); if(!is.finite(s0)||s0<=0) s0<-sd0; c(mu0,log(s0*1.05)) },
#     "uniform"=c(min(x),max(x))
#   )
#   opt <- optim(start, nll, method="Nelder-Mead", control=list(maxit=5000))
#   params <- switch(
#     dist_name,
#     "norm"=list(mu=opt$par[1], sd=exp(opt$par[2])),
#     "logistic"=list(location=opt$par[1], scale=exp(opt$par[2])),
#     "cosine"=list(mu=opt$par[1], s=exp(opt$par[2])),
#     "hypsecant"=list(mu=opt$par[1], sigma=exp(opt$par[2])),
#     "uniform"=list(min=opt$par[1], max=opt$par[2])
#   )
#   list(dist=dist_name, params=params)
# }
# ks_for_fit <- function(x, fit) {
#   x <- x[is.finite(x)]
#   cdf_fun <- switch(
#     fit$dist,
#     "norm"=function(q) pnorm(q, fit$params$mu, fit$params$sd),
#     "logistic"=function(q) plogis(q, fit$params$location, fit$params$scale),
#     "uniform"=function(q) punif(q, fit$params$min, fit$params$max),
#     "cosine"=function(q) pcosine_sc(q, fit$params$mu, fit$params$s),
#     "hypsecant"=function(q) phypsec_sc(q, fit$params$mu, fit$params$sigma)
#   )
#   ks <- suppressWarnings(stats::ks.test(x, cdf_fun))
#   list(stat=unname(ks$statistic), p=as.numeric(ks$p.value))
# }
# pick_best <- function(x) {
#   res <- rbindlist(lapply(CANDIDATE_DISTS, function(dn) {
#     fit <- fit_dist_mle(x, dn)
#     ks <- ks_for_fit(x, fit)
#     data.table(dist=dn, ks_stat=ks$stat, p_value=ks$p, params=list(fit$params))
#   }))
#   res <- res[order(-p_value, ks_stat, dist)]
#   i <- 1L; while (i <= nrow(res) && res$p_value[i] < PVAL_THRESHOLD) i <- i+1L
#   if (i > nrow(res)) i <- 1L
#   res[i]
# }
# return_levels <- function(fit_row) {
#   dist <- fit_row$dist
#   params <- fit_row$params[[1]]
#   p <- 1 - 1/RP_COMPUTE
#   p <- pmin(pmax(p, 1e-6), 1-1e-6)
#   q_fun <- switch(
#     dist,
#     "norm"=function(pp) qnorm(pp, params$mu, params$sd),
#     "logistic"=function(pp) qlogis(pp, params$location, params$scale),
#     "uniform"=function(pp) qunif(pp, params$min, params$max),
#     "cosine"=function(pp) qcosine_sc(pp, params$mu, params$s),
#     "hypsecant"=function(pp) qhypsec_sc(pp, params$mu, params$sigma)
#   )
#   as.numeric(q_fun(p))
# }
# 
# row_max_na <- function(m) apply(m, 1, function(v) if (all(!is.finite(v))) NA_real_ else max(v, na.rm=TRUE))
# 
# run_stageB <- function() {
#   # discover Stage A files
#   all_qs <- fs::dir_ls(OUT_TS_DIR, type="file", regexp="\\.qs$", fail=FALSE)
#   all_qs <- sort(as.character(all_qs))
#   year_files <- all_qs[grepl("^pr_dailySUM_timeseries_[0-9]{4}\\.qs$", fs::path_file(all_qs))]
#   if (length(year_files) == 0) stop("No Stage A files in: ", OUT_TS_DIR)
#   
#   years <- as.integer(str_match(fs::path_file(year_files), "^pr_dailySUM_timeseries_([0-9]{4})\\.qs$")[,2])
#   ord <- order(years); year_files <- year_files[ord]; years <- years[ord]
#   
#   obj0 <- qread(year_files[1])
#   grid_master <- as.data.table(obj0$grid); rm(obj0); gc()
#   grid_master[, cell := as.integer(cell)]
#   setorder(grid_master, cell)
#   n_cells <- nrow(grid_master)
#   
#   chunk_ids <- seq_len(ceiling(n_cells / CHUNK_SIZE))
#   chunk_rows <- lapply(chunk_ids, function(k) {
#     i1 <- (k-1L) * CHUNK_SIZE + 1L
#     i2 <- min(k * CHUNK_SIZE, n_cells)
#     i1:i2
#   })
#   
#   # resume scan
#   all_chunks <- fs::dir_ls(OUT_RP_CHUNKS_DIR, type="file", regexp="\\.qs$", fail=FALSE)
#   all_chunks <- sort(as.character(all_chunks))
#   done_files <- all_chunks[grepl("^rp_pr24_chunk_[0-9]{4}\\.qs$", fs::path_file(all_chunks))]
#   done_ids <- integer(0)
#   if (length(done_files) > 0) {
#     done_ids <- as.integer(str_match(fs::path_file(done_files), "^rp_pr24_chunk_([0-9]{4})\\.qs$")[,2])
#     done_ids <- done_ids[is.finite(done_ids)]
#   }
#   done_ids <- sort(unique(done_ids))
#   
#   todo <- chunk_ids
#   if (!OVERWRITE_STAGE_B) todo <- setdiff(todo, done_ids)
#   
#   message(sprintf("[Stage B] years=%d (%d–%d) | chunks total=%d | done=%d | todo=%d",
#                   length(years), min(years), max(years),
#                   length(chunk_ids), length(done_ids), length(todo)))
#   if (length(todo) == 0) return(invisible(TRUE))
#   
#   process_chunk <- function(cid) {
#     out_qs <- file.path(OUT_RP_CHUNKS_DIR, sprintf("rp_pr24_chunk_%04d.qs", cid))
#     if (!OVERWRITE_STAGE_B && fs::file_exists(out_qs)) return(list(chunk=cid,status="skip_exists"))
#     cleanup_tmp_files(out_qs)
#     
#     idx <- chunk_rows[[cid]]
#     grid_chunk <- grid_master[idx]
#     n_chunk <- nrow(grid_chunk)
#     
#     ann_max <- matrix(NA_real_, nrow=n_chunk, ncol=length(years))
#     
#     for (j in seq_along(year_files)) {
#       obj <- qread(year_files[j])
#       
#       g <- as.data.table(obj$grid); g[, cell := as.integer(cell)]
#       pos <- match(grid_chunk$cell, g$cell)
#       if (anyNA(pos)) stop("Grid mismatch in year ", years[j])
#       
#       xmat <- obj$pr_sum[pos, , drop=FALSE]   # subset while still float32
#       xmat <- float::dbl(xmat)               # convert to base double matrix
#       ann_max[, j] <- row_max_na(xmat)
#       
#       rm(obj, g, pos, xmat); gc()
#     }
#     
#     out <- data.table(cell=grid_chunk$cell, lon=grid_chunk$lon, lat=grid_chunk$lat)
#     for (rp in RP_LABELS) out[[sprintf("rp_%dyr", rp)]] <- NA_real_
#     out[, `:=`(dist=NA_character_, ks_p=NA_real_, n_years_used=0L)]
#     
#     for (i in seq_len(n_chunk)) {
#       x <- ann_max[i, ]
#       x <- x[is.finite(x)]
#       out$n_years_used[i] <- length(x)
#       if (length(x) < MIN_YEARS_REQUIRED) next
#       
#       pick <- tryCatch(pick_best(x), error=function(e) NULL)
#       if (is.null(pick)) next
#       
#       vals <- tryCatch(return_levels(pick), error=function(e) rep(NA_real_, length(RP_LABELS)))
#       
#       out$dist[i] <- as.character(pick$dist)
#       out$ks_p[i] <- as.numeric(pick$p_value)
#       
#       for (k in seq_along(RP_LABELS)) {
#         out[[sprintf("rp_%dyr", RP_LABELS[k])]][i] <- vals[k]
#       }
#     }
#     
#     atomic_qsave(out, out_qs, preset=QS_PRESET)
#     rm(out, ann_max); gc()
#     list(chunk=cid,status="ok")
#   }
#   
#   future::plan(multisession, workers = WORKERS_B)
#   res <- future.apply::future_lapply(todo, process_chunk, future.seed = TRUE)
#   
#   res_dt <- rbindlist(res, fill=TRUE)
#   print(res_dt[, .N, by=status][order(-N)])
#   invisible(TRUE)
# }
# 
# # ----------------------------
# # Stage C: Compile NetCDF on CONUS ref grid
# # ----------------------------
# run_stageC <- function() {
#   all_files <- fs::dir_ls(OUT_RP_CHUNKS_DIR, type="file", fail=FALSE)
#   all_files <- sort(as.character(all_files))
#   rp_files <- all_files[grepl("^rp_pr24_chunk_[0-9]{4}\\.qs$", fs::path_file(all_files))]
#   
#   message(sprintf("[Stage C] found %d rp chunk files in: %s", length(rp_files), OUT_RP_CHUNKS_DIR))
#   if (length(rp_files) == 0) stop("No rp chunk files in: ", OUT_RP_CHUNKS_DIR)
#   
#   r_ref <- terra::rast(REF_NC_CONUS)[[1]]
#   ncell_ref <- terra::ncell(r_ref)
#   
#   vals_mat <- matrix(NA_real_, nrow = ncell_ref, ncol = length(RP_LABELS))
#   colnames(vals_mat) <- sprintf("rp_%dyr", RP_LABELS)
#   
#   required_cols <- c("cell", sprintf("rp_%dyr", RP_LABELS))
#   
#   for (f in rp_files) {
#     dt <- as.data.table(qread(f))
#     miss <- setdiff(required_cols, names(dt))
#     if (length(miss) > 0) stop("Missing columns in ", fs::path_file(f), ": ", paste(miss, collapse=", "))
#     
#     cell_idx <- as.integer(dt$cell)
#     keep <- is.finite(cell_idx) & !is.na(cell_idx) & cell_idx >= 1 & cell_idx <= ncell_ref
#     if (!any(keep)) next
#     
#     cell_idx <- cell_idx[keep]
#     dt <- dt[keep]
#     
#     for (k in seq_along(RP_LABELS)) {
#       col <- sprintf("rp_%dyr", RP_LABELS[k])
#       v <- as.numeric(dt[[col]])
#       ok <- is.finite(v)
#       if (any(ok)) vals_mat[cell_idx[ok], k] <- v[ok]
#     }
#     
#     rm(dt); gc()
#   }
#   
#   stack <- terra::rast(replicate(length(RP_LABELS), r_ref, simplify = FALSE))
#   names(stack) <- colnames(vals_mat)
#   stack <- terra::setValues(stack, vals_mat)
#   
#   terra::writeCDF(
#     stack,
#     filename  = OUT_NC_CONUS,
#     overwrite = OVERWRITE_STAGE_C,
#     varname   = "precip_idf_pr24",
#     unit      = "mm"
#   )
#   
#   message("[WRITE] ", OUT_NC_CONUS)
#   invisible(OUT_NC_CONUS)
# }
# 
# # ----------------------------
# # RUN ALL
# # ----------------------------
# run_all <- function() {
#   stopifnot(fs::file_exists(REF_NC_CONUS))
#   grid <- build_lookup_from_conus_ref()
#   message("[REF CONUS] ", paste(as.vector(terra::ext(terra::rast(REF_NC_CONUS)[[1]])), collapse=", "))
#   
#   if (RUN_STAGE_A) {
#     resA <- run_stageA(grid)
#     
#     # hard stop if stage A produced no outputs
#     n_ok <- if (nrow(resA) == 0) 0L else resA[status == "ok", .N]
#     if (n_ok == 0L) {
#       stop("Stage A produced 0 yearly files. Fix Stage A error(s) above before Stage B.")
#     }
#   }
#   
#   if (RUN_STAGE_B) run_stageB()
#   if (RUN_STAGE_C) run_stageC()
#   invisible(TRUE)
# }
# 
# # ---- RUN ----
# # For max stability first:
# future::plan(sequential); WORKERS_A <- 6L; WORKERS_B <- 6L
# run_all()
# 
# # future::plan(multisession, workers = 6L)
# # WORKERS_A <- 6L
# # run_stageA()
# # 
# # WORKERS_B <- 18L
# # run_stageB()
# # 
# # run_stageC()
# 
# test <- terra::rast('V:/users/hquintal/phd1_cluster_southeast_v2/data/out/04_rp/precip_idf_pr24_return_periods.nc')
# plot(test[[10]])
# 
# # test_diff <- test$precip_idf_pr24_3-test$precip_idf_pr24_2
# # plot(test_diff)

# Southeast ----
# =============================================================================
# 24-hour precipitation return periods (daily SUM then annual maxima)
#
# Stage A: yearly matrices of DAILY SUM precipitation on BBOX-cropped ref grid
# Stage B: per-cell distribution fit on ANNUAL MAX of daily sums (24-hr totals)
# Stage C: compile NetCDF on BBOX-cropped ref grid
#
# Key requirements (matches your HI + experimental covariance masking logic):
#  - Use REF_NC_CONUS as the NA mask + target grid
#  - Crop that raster to BBOX_LONLAT; DO NOT trim
#  - Build lookup ONLY from non-NA cells in cropped ref raster (land/valid-cell mask)
#  - Stage A filters dt by BBOX first, then joins to lookup grid (implicitly enforces mask)
#  - Stage C writes NetCDF on the cropped ref raster so extent/mask align exactly
#
# Robustness fixes:
#  - Stage B converts float32 matrices with float::dbl() BEFORE subsetting
#    (avoids aperm/as.matrix edge-case errors)
# =============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(data.table)
  library(fs)
  library(terra)
  library(float)
  library(stringr)
  library(future)
  library(future.apply)
  library(progressr)
})

# ----------------------------
# USER SETTINGS
# ----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

# Daily files (one per day) with subdaily/hourly precip samples
IN_PR_DIR <- file.path(ROOT, "data/out/01_vars/precipitation")

# Mask/target raster (lon/lat) defining allowed cells via NA pattern
REF_NC_CONUS <- file.path(ROOT, "data/out/01_vars/precipitation/_mask_target.tif")

# BBOX in lon/lat
BBOX_LONLAT <- c(xmin = -95, xmax = -75, ymin = 22, ymax = 40)

OUT_DIR <- file.path(ROOT, "data/out/04_rp")

# Outputs (bbox)
OUT_TS_DIR        <- file.path(OUT_DIR, "precip_dailySUM_timeseries_by_year_bbox_pr24")
OUT_RP_CHUNKS_DIR <- file.path(OUT_DIR, "precip_idf_pr24_return_periods_chunks_bbox")
OUT_NC_BBOX       <- file.path(OUT_DIR, "precip_idf_pr24_return_periods_bbox.nc")

dir_create(OUT_TS_DIR)
dir_create(OUT_RP_CHUNKS_DIR)

LOOKUP_QS <- file.path(OUT_DIR, "precip_refgrid_lonlat_to_cell_bbox_pr24.qs")

# Controls
RUN_STAGE_A <- TRUE
RUN_STAGE_B <- TRUE
RUN_STAGE_C <- TRUE

OVERWRITE_LOOKUP   <- TRUE
OVERWRITE_STAGE_A  <- TRUE
OVERWRITE_STAGE_B  <- TRUE
OVERWRITE_STAGE_C  <- TRUE

# Parallelism
WORKERS_A  <- 10L   # IO-heavy
WORKERS_B  <- 10L   # compute-heavy
CHUNK_SIZE <- 300L

QS_PRESET  <- "high"
KEY_DIGITS <- 6
options(future.globals.maxSize = 8 * 1024^3)

kfmt <- function(x, digits = KEY_DIGITS) sprintf(paste0("%.", digits, "f"), x)

# Column in daily precip .qs
PR_COL <- "intensity"   # change if needed

# Optional: drop tiny daily totals (disable with NA_real_ or NULL)
DROP_DAILY_TOTAL_LT_MM <- NA_real_

# Return periods
RP_LABELS  <- c(1, 2, 5, 10, 25, 50, 100, 200, 500, 1000)
RP_COMPUTE <- c(1.01, 2, 5, 10, 25, 50, 100, 200, 500, 1000)

# Fit settings
MIN_YEARS_REQUIRED <- 20L
PVAL_THRESHOLD     <- 0.4
CANDIDATE_DISTS    <- c("cosine", "hypsecant", "logistic", "norm", "uniform")

# ----------------------------
# Helpers: file parsing + atomic writes
# ----------------------------
is_daily_qs <- function(p) grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}\\.qs$", fs::path_file(p))
day_from_file_chr <- function(p) sub("\\.qs$", "", fs::path_file(p))

cleanup_tmp_files <- function(final_path) {
  patt <- paste0(fs::path_file(final_path), ".tmp_*")
  junk <- fs::dir_ls(dirname(final_path), glob = patt, fail = FALSE, type = "file")
  if (length(junk) > 0) fs::file_delete(junk)
}
atomic_qsave <- function(obj, final_path, preset = "high") {
  tmp <- paste0(final_path, ".tmp_", Sys.getpid())
  qs::qsave(obj, tmp, preset = preset)
  fs::file_move(tmp, final_path)
  invisible(TRUE)
}

# ----------------------------
# Stage 0: BBOX ref raster + lookup (mask enforced by non-NA cells)
# ----------------------------
get_ref_bbox_raster <- function() {
  stopifnot(fs::file_exists(REF_NC_CONUS))
  r0 <- terra::rast(REF_NC_CONUS)[[1]]
  
  e <- terra::ext(
    BBOX_LONLAT["xmin"], BBOX_LONLAT["xmax"],
    BBOX_LONLAT["ymin"], BBOX_LONLAT["ymax"]
  )
  
  # Crop to bbox; preserve NA mask
  terra::crop(r0, e, snap = "out")
}

build_lookup_from_ref_bbox <- function() {
  if (!OVERWRITE_LOOKUP && fs::file_exists(LOOKUP_QS)) {
    grid <- as.data.table(qread(LOOKUP_QS))
    if (!all(c("cell","lon","lat","lon_key","lat_key") %in% names(grid))) {
      stop("[LOOKUP] Existing lookup missing required columns.")
    }
    setorder(grid, cell)
    setkey(grid, lon_key, lat_key)
    return(grid)
  }
  
  message("[LOOKUP] Building lookup from REF bbox grid: ", REF_NC_CONUS)
  r <- get_ref_bbox_raster()
  
  v <- terra::values(r, mat = FALSE)
  cells <- which(!is.na(v))
  if (length(cells) == 0) stop("REF bbox raster has 0 non-NA cells (mask empty in bbox).")
  
  xy <- terra::xyFromCell(r, cells)
  grid <- data.table(
    cell = as.integer(cells),   # cell index into the *bbox raster*
    lon  = as.numeric(xy[, 1]),
    lat  = as.numeric(xy[, 2])
  )
  grid[, `:=`(lon_key = kfmt(lon), lat_key = kfmt(lat))]
  setorder(grid, cell)
  setkey(grid, lon_key, lat_key)
  
  qsave(grid, LOOKUP_QS, preset = QS_PRESET)
  message("[LOOKUP] Wrote: ", LOOKUP_QS, " | n_cells=", nrow(grid))
  grid
}

# ----------------------------
# Stage A: daily SUM vector for one day (bbox-filter, join to grid)
# ----------------------------
daily_sum_vec_one_day_bbox <- function(day_file, grid) {
  dt <- as.data.table(qread(day_file))
  
  req <- c("lon","lat",PR_COL)
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) stop("Missing columns in ", fs::path_file(day_file), ": ", paste(miss, collapse = ", "))
  
  dt[, `:=`(lon = as.numeric(lon), lat = as.numeric(lat))]
  dt[, (PR_COL) := as.numeric(get(PR_COL))]
  
  # BBOX filter first (fast)
  dt <- dt[
    lon >= BBOX_LONLAT["xmin"] & lon <= BBOX_LONLAT["xmax"] &
      lat >= BBOX_LONLAT["ymin"] & lat <= BBOX_LONLAT["ymax"]
  ]
  
  # finite only
  dt <- dt[is.finite(get(PR_COL))]
  if (nrow(dt) == 0) return(rep(NA_real_, nrow(grid)))
  
  # daily sum per location (sum across hours/subdaily rows)
  ds <- dt[, .(pr_sum = sum(get(PR_COL), na.rm = TRUE)), by = .(lon, lat)]
  
  # optional tiny-total drop
  if (is.finite(DROP_DAILY_TOTAL_LT_MM)) {
    ds[pr_sum < DROP_DAILY_TOTAL_LT_MM, pr_sum := NA_real_]
  }
  
  ds[, `:=`(lon_key = kfmt(lon), lat_key = kfmt(lat))]
  setkey(ds, lon_key, lat_key)
  
  m <- ds[grid, on = .(lon_key, lat_key)]
  as.numeric(m$pr_sum)
}

stageA_process_one_year <- function(y, files_year, days_year, grid) {
  out_qs <- file.path(OUT_TS_DIR, sprintf("pr_dailySUM_timeseries_%d.qs", y))
  if (!OVERWRITE_STAGE_A && fs::file_exists(out_qs)) return(list(year = y, status = "skip_exists", file = out_qs))
  cleanup_tmp_files(out_qs)
  
  n_cells <- nrow(grid)
  n_days  <- length(files_year)
  if (n_days == 0) return(list(year = y, status = "skip_empty", file = out_qs))
  
  mat <- matrix(NA_real_, nrow = n_cells, ncol = n_days)
  for (j in seq_len(n_days)) {
    v <- daily_sum_vec_one_day_bbox(files_year[j], grid)
    if (length(v) != n_cells) stop("Day vec length mismatch for ", files_year[j])
    mat[, j] <- v
    if (j %% 25 == 0) gc()
  }
  
  obj <- list(
    meta = list(
      variable    = "precipitation",
      statistic   = "daily_sum",
      duration    = "24-hr",
      units       = "mm",
      year        = y,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      ref_mask    = REF_NC_CONUS,
      bbox_lonlat = BBOX_LONLAT,
      pr_col      = PR_COL,
      drop_daily_total_lt_mm = DROP_DAILY_TOTAL_LT_MM,
      key_digits  = KEY_DIGITS
    ),
    grid   = grid[, .(cell, lon, lat)],
    days   = as.Date(days_year),
    pr_sum = float::fl(mat)
  )
  
  atomic_qsave(obj, out_qs, preset = QS_PRESET)
  rm(obj, mat); gc()
  list(year = y, status = "ok", file = out_qs)
}

run_stageA <- function(grid) {
  message("\n==================== STAGE A: Daily SUM precip by year (BBOX) ====================")
  
  files <- sort(as.character(fs::dir_ls(IN_PR_DIR, type = "file", glob = "*.qs", fail = FALSE)))
  files <- files[is_daily_qs(files)]
  if (length(files) == 0) stop("No daily precip .qs files in: ", IN_PR_DIR)
  
  day_chr <- vapply(files, day_from_file_chr, character(1))
  days <- as.Date(day_chr)
  o <- order(days); files <- files[o]; days <- days[o]
  
  years <- as.integer(format(days, "%Y"))
  by_year <- split(seq_along(files), years)
  years_u <- sort(unique(years))
  
  # resume scan
  existing <- fs::dir_ls(OUT_TS_DIR, type = "file", regexp = "\\.qs$", fail = FALSE)
  existing <- as.character(existing)
  done_years <- integer(0)
  if (length(existing) > 0) {
    done_years <- as.integer(str_match(fs::path_file(existing), "^pr_dailySUM_timeseries_([0-9]{4})\\.qs$")[, 2])
    done_years <- done_years[is.finite(done_years)]
  }
  todo <- years_u
  if (!OVERWRITE_STAGE_A) todo <- setdiff(todo, done_years)
  
  message(sprintf("[Stage A] years=%d (%d–%d) | done=%d | todo=%d | n_cells(bbox)=%d",
                  length(years_u), min(years_u), max(years_u),
                  length(done_years), length(todo), nrow(grid)))
  if (length(todo) == 0) return(invisible(TRUE))
  
  future::plan(multisession, workers = WORKERS_A)
  progressr::handlers(global = TRUE)
  progressr::handlers("txtprogressbar")
  
  res <- progressr::with_progress({
    p <- progressr::progressor(along = todo)
    future.apply::future_lapply(todo, function(y) {
      on.exit(p(sprintf("year=%d", y)), add = TRUE)
      idx <- by_year[[as.character(y)]]
      tryCatch(
        stageA_process_one_year(y, files[idx], days[idx], grid),
        error = function(e) list(year = y, status = "error",
                                 file = file.path(OUT_TS_DIR, sprintf("pr_dailySUM_timeseries_%d.qs", y)),
                                 error = conditionMessage(e))
      )
    }, future.seed = TRUE)
  })
  
  res_dt <- rbindlist(res, fill = TRUE)
  message("\n[STAGE A] Summary:")
  print(res_dt[, .N, by = status][order(-N)])
  if (res_dt[status == "error", .N] > 0) {
    message("\n[ERROR YEARS]")
    print(res_dt[status == "error", .(year, error)][order(year)])
  }
  
  message("[OUT] ", OUT_TS_DIR)
  invisible(TRUE)
}

# ----------------------------
# Stage B: distribution fit on ANNUAL MAX of DAILY SUM (pr24)
# ----------------------------
dcosine_sc <- function(x, mu, s, log = FALSE) {
  y <- (x - mu) / s
  dens <- ifelse(abs(y) <= 1, 1 / (2 * s) * (1 + cos(pi * y)), 0)
  if (log) log(dens) else dens
}
pcosine_sc <- function(q, mu, s) {
  y <- (q - mu) / s
  ifelse(y <= -1, 0, ifelse(y >= 1, 1, 0.5 * (1 + y + sin(pi * y) / pi)))
}
qcosine_sc <- function(p, mu, s) {
  f <- function(y, target) 0.5 * (1 + y + sin(pi * y) / pi) - target
  y <- vapply(p, function(pp) uniroot(f, interval = c(-1, 1), target = pp)$root, numeric(1))
  mu + s * y
}
dhypsec_sc <- function(x, mu, sigma, log = FALSE) {
  z <- (pi / 2) * ((x - mu) / sigma)
  dens <- (1 / (2 * sigma)) * (1 / cosh(z))
  if (log) log(dens) else dens
}
phypsec_sc <- function(q, mu, sigma) {
  z <- (pi / 2) * ((q - mu) / sigma)
  (2 / pi) * atan(exp(z))
}
qhypsec_sc <- function(p, mu, sigma) {
  mu + sigma * (2 / pi) * log(tan(pi * p / 2))
}

fit_dist_mle <- function(x, dist_name) {
  x <- x[is.finite(x)]
  if (length(x) < MIN_YEARS_REQUIRED) stop("Too few years.")
  mu0 <- mean(x); sd0 <- sd(x); if (!is.finite(sd0) || sd0 <= 0) sd0 <- 1
  eps <- 1e-12
  
  nll <- switch(
    dist_name,
    "norm" = function(par) -sum(dnorm(x, mean = par[1], sd = exp(par[2]), log = TRUE)),
    "logistic" = function(par) -sum(dlogis(x, location = par[1], scale = exp(par[2]), log = TRUE)),
    "uniform" = function(par) {
      a <- par[1]; b <- par[2]
      if (b <= a) return(Inf)
      if (any(x < a | x > b)) return(Inf)
      -sum(dunif(x, a, b, log = TRUE))
    },
    "cosine" = function(par) {
      mu <- par[1]; s <- exp(par[2])
      if (s <= 0) return(Inf)
      if (max(abs((x - mu) / s)) > 1) return(Inf)
      d <- pmax(dcosine_sc(x, mu, s, FALSE), eps)
      -sum(log(d))
    },
    "hypsecant" = function(par) {
      mu <- par[1]; sg <- exp(par[2])
      if (sg <= 0) return(Inf)
      d <- pmax(dhypsec_sc(x, mu, sg, FALSE), eps)
      -sum(log(d))
    },
    stop("Unknown dist")
  )
  
  start <- switch(
    dist_name,
    "norm" = c(mu0, log(sd0)),
    "logistic" = c(mu0, log(sd0)),
    "hypsecant" = c(mu0, log(sd0)),
    "cosine" = {
      s0 <- max(abs(x - mu0)); if (!is.finite(s0) || s0 <= 0) s0 <- sd0
      c(mu0, log(s0 * 1.05))
    },
    "uniform" = c(min(x), max(x))
  )
  
  opt <- optim(start, nll, method = "Nelder-Mead", control = list(maxit = 5000))
  
  params <- switch(
    dist_name,
    "norm" = list(mu = opt$par[1], sd = exp(opt$par[2])),
    "logistic" = list(location = opt$par[1], scale = exp(opt$par[2])),
    "cosine" = list(mu = opt$par[1], s = exp(opt$par[2])),
    "hypsecant" = list(mu = opt$par[1], sigma = exp(opt$par[2])),
    "uniform" = list(min = opt$par[1], max = opt$par[2])
  )
  
  list(dist = dist_name, params = params)
}

ks_for_fit <- function(x, fit) {
  x <- x[is.finite(x)]
  cdf_fun <- switch(
    fit$dist,
    "norm" = function(q) pnorm(q, fit$params$mu, fit$params$sd),
    "logistic" = function(q) plogis(q, fit$params$location, fit$params$scale),
    "uniform" = function(q) punif(q, fit$params$min, fit$params$max),
    "cosine" = function(q) pcosine_sc(q, fit$params$mu, fit$params$s),
    "hypsecant" = function(q) phypsec_sc(q, fit$params$mu, fit$params$sigma)
  )
  ks <- suppressWarnings(stats::ks.test(x, cdf_fun))
  list(stat = unname(ks$statistic), p = as.numeric(ks$p.value))
}

pick_best <- function(x) {
  res <- rbindlist(lapply(CANDIDATE_DISTS, function(dn) {
    fit <- fit_dist_mle(x, dn)
    ks <- ks_for_fit(x, fit)
    data.table(dist = dn, ks_stat = ks$stat, p_value = ks$p, params = list(fit$params))
  }))
  res <- res[order(-p_value, ks_stat, dist)]
  i <- 1L
  while (i <= nrow(res) && res$p_value[i] < PVAL_THRESHOLD) i <- i + 1L
  if (i > nrow(res)) i <- 1L
  res[i]
}

return_levels <- function(fit_row) {
  dist <- fit_row$dist
  params <- fit_row$params[[1]]
  p <- 1 - 1 / RP_COMPUTE
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  q_fun <- switch(
    dist,
    "norm" = function(pp) qnorm(pp, params$mu, params$sd),
    "logistic" = function(pp) qlogis(pp, params$location, params$scale),
    "uniform" = function(pp) qunif(pp, params$min, params$max),
    "cosine" = function(pp) qcosine_sc(pp, params$mu, params$s),
    "hypsecant" = function(pp) qhypsec_sc(pp, params$mu, params$sigma)
  )
  as.numeric(q_fun(p))
}

row_max_na <- function(m) {
  apply(m, 1, function(v) if (all(!is.finite(v))) NA_real_ else max(v, na.rm = TRUE))
}

run_stageB <- function() {
  message("\n==================== STAGE B: Return periods (pr24 annual maxima) ====================")
  
  all_qs <- fs::dir_ls(OUT_TS_DIR, type = "file", regexp = "\\.qs$", fail = FALSE)
  all_qs <- sort(as.character(all_qs))
  year_files <- all_qs[grepl("^pr_dailySUM_timeseries_[0-9]{4}\\.qs$", fs::path_file(all_qs))]
  if (length(year_files) == 0) stop("No Stage A files in: ", OUT_TS_DIR)
  
  years <- as.integer(str_match(fs::path_file(year_files), "^pr_dailySUM_timeseries_([0-9]{4})\\.qs$")[, 2])
  ord <- order(years); years <- years[ord]; year_files <- year_files[ord]
  
  obj0 <- qread(year_files[1])
  grid_master <- as.data.table(obj0$grid); rm(obj0); gc()
  grid_master[, cell := as.integer(cell)]
  setorder(grid_master, cell)
  n_cells <- nrow(grid_master)
  
  chunk_ids <- seq_len(ceiling(n_cells / CHUNK_SIZE))
  chunk_rows <- lapply(chunk_ids, function(k) {
    i1 <- (k - 1L) * CHUNK_SIZE + 1L
    i2 <- min(k * CHUNK_SIZE, n_cells)
    i1:i2
  })
  
  # resume scan
  all_chunks <- fs::dir_ls(OUT_RP_CHUNKS_DIR, type = "file", regexp = "\\.qs$", fail = FALSE)
  all_chunks <- sort(as.character(all_chunks))
  done_files <- all_chunks[grepl("^rp_pr24_chunk_[0-9]{4}\\.qs$", fs::path_file(all_chunks))]
  done_ids <- integer(0)
  if (length(done_files) > 0) {
    done_ids <- as.integer(str_match(fs::path_file(done_files), "^rp_pr24_chunk_([0-9]{4})\\.qs$")[, 2])
    done_ids <- done_ids[is.finite(done_ids)]
  }
  done_ids <- sort(unique(done_ids))
  
  todo <- chunk_ids
  if (!OVERWRITE_STAGE_B) todo <- setdiff(todo, done_ids)
  
  message(sprintf("[Stage B] years=%d (%d–%d) | chunks total=%d | done=%d | todo=%d",
                  length(years), min(years), max(years),
                  length(chunk_ids), length(done_ids), length(todo)))
  if (length(todo) == 0) return(invisible(TRUE))
  
  process_chunk <- function(cid) {
    out_qs <- file.path(OUT_RP_CHUNKS_DIR, sprintf("rp_pr24_chunk_%04d.qs", cid))
    if (!OVERWRITE_STAGE_B && fs::file_exists(out_qs)) return(list(chunk = cid, status = "skip_exists"))
    cleanup_tmp_files(out_qs)
    
    idx <- chunk_rows[[cid]]
    grid_chunk <- grid_master[idx]
    n_chunk <- nrow(grid_chunk)
    
    ann_max <- matrix(NA_real_, nrow = n_chunk, ncol = length(years))
    
    for (j in seq_along(year_files)) {
      obj <- qread(year_files[j])
      
      g <- as.data.table(obj$grid); g[, cell := as.integer(cell)]
      pos <- match(grid_chunk$cell, g$cell)
      if (anyNA(pos)) stop("Grid mismatch in year ", years[j])
      
      # ROBUST: convert float32 -> base matrix FIRST, then subset
      xmat_all <- float::dbl(obj$pr_sum)
      xmat     <- xmat_all[pos, , drop = FALSE]
      ann_max[, j] <- row_max_na(xmat)
      
      rm(obj, g, pos, xmat_all, xmat); gc()
    }
    
    out <- data.table(cell = grid_chunk$cell, lon = grid_chunk$lon, lat = grid_chunk$lat)
    for (rp in RP_LABELS) out[[sprintf("rp_%dyr", rp)]] <- NA_real_
    out[, `:=`(dist = NA_character_, ks_p = NA_real_, n_years_used = 0L)]
    
    for (i in seq_len(n_chunk)) {
      x <- ann_max[i, ]
      x <- x[is.finite(x)]
      out$n_years_used[i] <- length(x)
      if (length(x) < MIN_YEARS_REQUIRED) next
      
      pick <- tryCatch(pick_best(x), error = function(e) NULL)
      if (is.null(pick)) next
      
      vals <- tryCatch(return_levels(pick), error = function(e) rep(NA_real_, length(RP_LABELS)))
      
      out$dist[i] <- as.character(pick$dist)
      out$ks_p[i] <- as.numeric(pick$p_value)
      for (k in seq_along(RP_LABELS)) {
        out[[sprintf("rp_%dyr", RP_LABELS[k])]][i] <- vals[k]
      }
    }
    
    atomic_qsave(out, out_qs, preset = QS_PRESET)
    rm(out, ann_max); gc()
    list(chunk = cid, status = "ok")
  }
  
  future::plan(multisession, workers = WORKERS_B)
  res <- future.apply::future_lapply(todo, process_chunk, future.seed = TRUE)
  
  res_dt <- rbindlist(res, fill = TRUE)
  message("\n[STAGE B] Summary:")
  print(res_dt[, .N, by = status][order(-N)])
  invisible(TRUE)
}

# ----------------------------
# Stage C: compile NetCDF on BBOX ref grid (mask preserved)
# ----------------------------
run_stageC <- function() {
  message("\n==================== STAGE C: NetCDF compile (bbox) ====================")
  
  all_files <- fs::dir_ls(OUT_RP_CHUNKS_DIR, type = "file", fail = FALSE)
  all_files <- sort(as.character(all_files))
  rp_files <- all_files[grepl("^rp_pr24_chunk_[0-9]{4}\\.qs$", fs::path_file(all_files))]
  
  message(sprintf("[Stage C] found %d rp chunk files in: %s", length(rp_files), OUT_RP_CHUNKS_DIR))
  if (length(rp_files) == 0) stop("No rp chunk files in: ", OUT_RP_CHUNKS_DIR)
  
  r_ref <- get_ref_bbox_raster()
  ncell_ref <- terra::ncell(r_ref)
  
  vals_mat <- matrix(NA_real_, nrow = ncell_ref, ncol = length(RP_LABELS))
  colnames(vals_mat) <- sprintf("rp_%dyr", RP_LABELS)
  
  required_cols <- c("cell", sprintf("rp_%dyr", RP_LABELS))
  
  for (f in rp_files) {
    dt <- as.data.table(qread(f))
    miss <- setdiff(required_cols, names(dt))
    if (length(miss) > 0) stop("Missing columns in ", fs::path_file(f), ": ", paste(miss, collapse = ", "))
    
    cell_idx <- as.integer(dt$cell)
    keep <- is.finite(cell_idx) & !is.na(cell_idx) & cell_idx >= 1 & cell_idx <= ncell_ref
    if (!any(keep)) next
    
    cell_idx <- cell_idx[keep]
    dt <- dt[keep]
    
    for (k in seq_along(RP_LABELS)) {
      col <- sprintf("rp_%dyr", RP_LABELS[k])
      v <- as.numeric(dt[[col]])
      ok <- is.finite(v)
      if (any(ok)) vals_mat[cell_idx[ok], k] <- v[ok]
    }
    
    rm(dt); gc()
  }
  
  stack <- terra::rast(replicate(length(RP_LABELS), r_ref, simplify = FALSE))
  names(stack) <- colnames(vals_mat)
  stack <- terra::setValues(stack, vals_mat)
  
  terra::writeCDF(
    stack,
    filename  = OUT_NC_BBOX,
    overwrite = OVERWRITE_STAGE_C,
    varname   = "precip_idf_pr24",
    unit      = "mm"
  )
  
  message("[WRITE] ", OUT_NC_BBOX)
  invisible(OUT_NC_BBOX)
}

# ----------------------------
# RUN ALL
# ----------------------------
run_all <- function() {
  stopifnot(fs::file_exists(REF_NC_CONUS))
  
  # Build lookup on bbox-cropped ref (mask enforced)
  grid <- build_lookup_from_ref_bbox()
  
  # QA: print bbox extent + non-NA cell counts
  r_ref <- get_ref_bbox_raster()
  message("[REF bbox extent] ", paste(as.vector(terra::ext(r_ref)), collapse = ", "))
  message(sprintf("[REF bbox ncell=%d | nonNA=%d]",
                  terra::ncell(r_ref),
                  sum(!is.na(terra::values(r_ref, mat = FALSE)))))
  
  if (RUN_STAGE_A) run_stageA(grid)
  if (RUN_STAGE_B) run_stageB()
  if (RUN_STAGE_C) run_stageC()
  invisible(TRUE)
}

# ---- RUN ----
run_all()

# Quick check plot (optional)
# rp <- terra::rast(OUT_NC_BBOX)
# plot(rp[[5]])
