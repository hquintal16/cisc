# # Calc max daily HI ----
# # STAGE A (HEAT INDEX): Daily max time series by year (standalone, resume-safe)
# #
# # - Builds reference grid lookup from REF_NC (cell index + lon/lat)
# # - Reads daily hourly .qs (YYYY-MM-DD.qs) with heat_index (numeric or NA)
# # - Computes DAILY MAX heat index per location
# # - Writes one .qs per year with:
# #     grid: cell (terra cell index), lon, lat
# #     days: Date vector
# #     hi_max: float matrix [n_cells x n_days]
# #
# # Robustness:
# # - Resume-safe based on actual files on disk
# # - Atomic writes: write temp -> rename
# # - Cleans leftover temp files
# # - Parallel optional; sequential recommended on flaky network I/O
# 
# # test<- terra::rast('V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation/_mask_target.tif')
# # test
# # plot(test)
# suppressPackageStartupMessages({
#   library(qs)
#   library(data.table)
#   library(fs)
#   library(terra)
#   library(float)
#   library(stringr)
#   library(progressr)
#   library(future)
#   library(future.apply)
# })
# 
# # USER SETTINGS
# IN_DIR <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/heatindex"
# REF_NC <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation/_mask_target.tif"
# 
# OUT_DIR    <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/04_rp"
# OUT_TS_DIR <- file.path(OUT_DIR, "heatindex_dailymax_timeseries_by_year")
# dir_create(OUT_TS_DIR)
# 
# LOOKUP_QS <- file.path(OUT_DIR, "heatindex_refgrid_lonlat_to_cell.qs")
# 
# # Resume/overwrite controls
# OVERWRITE_LOOKUP <- TRUE
# OVERWRITE_YEARS  <- TRUE
# 
# # Execution mode
# PARALLEL <- TRUE          # set FALSE for max stability
# WORKERS  <- 26L            # keep conservative on V: drive
# options(future.globals.maxSize = 8 * 1024^3)  # bump if needed (8GB)
# 
# QS_PRESET  <- "high"
# KEY_DIGITS <- 6  # lon/lat join key formatting
# 
# # Helpers
# kfmt <- function(x, digits = KEY_DIGITS) sprintf(paste0("%.", digits, "f"), x)
# 
# is_daily_qs <- function(p) {
#   p <- as.character(p)
#   grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}\\.qs$", fs::path_file(p))
# }
# day_from_file_chr <- function(p) {
#   p <- as.character(p)
#   sub("\\.qs$", "", fs::path_file(p))
# }
# out_year_file <- function(y) file.path(OUT_TS_DIR, sprintf("hi_dailymax_timeseries_%d.qs", y))
# 
# # Lookup from REF grid (cell index + lon/lat keys)
# build_lookup_from_ref <- function() {
#   if (!OVERWRITE_LOOKUP && fs::file_exists(LOOKUP_QS)) {
#     message("[LOOKUP] Exists: ", LOOKUP_QS)
#     grid <- as.data.table(qread(LOOKUP_QS))
#     if (!all(c("cell","lon","lat","lon_key","lat_key") %in% names(grid))) {
#       stop("[LOOKUP] Existing lookup missing required columns.")
#     }
#     setorder(grid, cell)
#     setkey(grid, lon_key, lat_key)
#     return(grid)
#   }
#   
#   message("[LOOKUP] Building lookup from REF_NC grid: ", REF_NC)
#   r <- terra::rast(REF_NC)[[1]]
#   
#   v <- terra::values(r, mat = FALSE)
#   cells <- which(!is.na(v))
#   if (length(cells) == 0) stop("REF_NC has 0 non-NA cells: ", REF_NC)
#   
#   xy <- terra::xyFromCell(r, cells)
#   grid <- data.table(
#     cell = as.integer(cells),
#     lon  = as.numeric(xy[,1]),
#     lat  = as.numeric(xy[,2])
#   )
#   grid[, `:=`(lon_key = kfmt(lon), lat_key = kfmt(lat))]
#   setorder(grid, cell)
#   setkey(grid, lon_key, lat_key)
#   
#   qsave(grid, LOOKUP_QS, preset = QS_PRESET)
#   message("[LOOKUP] Wrote: ", LOOKUP_QS, " | n_cells=", nrow(grid))
#   grid
# }
# 
# # Daily max vector aligned to GRID row order (length = nrow(grid))
# daily_max_vec <- function(day_file, grid, n_cells_grid) {
#   dt <- as.data.table(qread(day_file))
#   if (!all(c("lon","lat","heat_index") %in% names(dt))) {
#     stop("Missing columns in: ", day_file)
#   }
#   
#   dt[, `:=`(
#     lon = as.numeric(lon),
#     lat = as.numeric(lat),
#     heat_index = as.numeric(heat_index)
#   )]
#   
#   dt <- dt[is.finite(heat_index)]
#   if (nrow(dt) == 0) return(rep(NA_real_, n_cells_grid))  # all-NA day is normal
#   
#   dm <- dt[, .(hi_max = max(heat_index, na.rm = TRUE)), by = .(lon, lat)]
#   dm[, `:=`(lon_key = kfmt(lon), lat_key = kfmt(lat))]
#   setkey(dm, lon_key, lat_key)
#   
#   # IMPORTANT: dm[grid] preserves grid row count/order
#   m <- dm[grid, on = .(lon_key, lat_key)]
#   as.numeric(m$hi_max)
# }
# 
# # Atomic write utilities
# cleanup_tmp_files <- function(out_qs) {
#   base <- fs::path_ext_remove(out_qs)
#   patt <- paste0(fs::path_file(base), ".tmp_*")
#   junk <- fs::dir_ls(dirname(out_qs), glob = patt, fail = FALSE, type="file")
#   if (length(junk) > 0) fs::file_delete(junk)
# }
# 
# atomic_qsave <- function(obj, out_qs, preset = "high") {
#   tmp_qs <- paste0(out_qs, ".tmp_", Sys.getpid())
#   qs::qsave(obj, tmp_qs, preset = preset)
#   fs::file_move(tmp_qs, out_qs)
#   invisible(TRUE)
# }
# 
# # Process one year (writes one file)
# process_one_year <- function(year, files_year, days_year, grid) {
#   out_qs <- out_year_file(year)
#   
#   if (!OVERWRITE_YEARS && fs::file_exists(out_qs)) {
#     return(list(year=year, status="skip_exists", file=out_qs))
#   }
#   
#   cleanup_tmp_files(out_qs)
#   
#   n_cells <- nrow(grid)
#   n_days  <- length(files_year)
#   if (n_days == 0) return(list(year=year, status="skip_empty", file=out_qs))
#   
#   mat <- matrix(NA_real_, nrow = n_cells, ncol = n_days)
#   for (j in seq_len(n_days)) {
#     v <- daily_max_vec(files_year[j], grid, n_cells_grid = n_cells)
#     if (length(v) != n_cells) stop("Length mismatch day vec: got ", length(v), " expected ", n_cells)
#     mat[, j] <- v
#     if (j %% 25 == 0) gc()
#   }
#   
#   obj <- list(
#     meta = list(
#       variable  = "heat_index",
#       statistic = "daily_max",
#       units     = "C",
#       year      = year,
#       created_utc = format(Sys.time(), tz="UTC", usetz=TRUE),
#       ref_nc    = REF_NC,
#       key_digits = KEY_DIGITS
#     ),
#     grid = grid[, .(cell, lon, lat)],  # terra cell index into REF_NC
#     days = as.Date(days_year),
#     hi_max = float::fl(mat)
#   )
#   
#   atomic_qsave(obj, out_qs, preset = QS_PRESET)
#   
#   rm(obj, mat); gc()
#   list(year=year, status="ok", file=out_qs)
# }
# 
# # MAIN RUN (resume-safe via disk scan)
# run_stageA_heatindex_dailymax <- function() {
#   message("\n==================== STAGE A: Daily max HI time series by year ====================")
#   
#   grid <- build_lookup_from_ref()
#   
#   files <- fs::dir_ls(IN_DIR, type="file", glob="*.qs")
#   files <- sort(as.character(files))
#   files <- files[is_daily_qs(files)]
#   if (length(files) == 0) stop("No daily .qs files found in: ", IN_DIR)
#   
#   day_chr <- vapply(files, day_from_file_chr, character(1))
#   days <- as.Date(day_chr, format = "%Y-%m-%d")
#   o <- order(days)
#   files <- files[o]
#   days  <- days[o]
#   
#   years <- as.integer(format(days, "%Y"))
#   years_u <- sort(unique(years))
#   
#   message(sprintf("[FILES] days=%d | years=%d (%d–%d) | n_cells(ref)=%d",
#                   length(files), length(years_u), min(years_u), max(years_u), nrow(grid)))
#   
#   by_year_idx <- split(seq_along(files), years)
#   
#   # --- disk-truth resume scan ---
#   existing_files <- fs::dir_ls(OUT_TS_DIR, type="file", glob="hi_dailymax_timeseries_*.qs", fail = FALSE)
#   existing_files <- as.character(existing_files)
#   
#   existing_years <- integer(0)
#   if (length(existing_files) > 0) {
#     existing_years <- as.integer(str_match(fs::path_file(existing_files),
#                                            "hi_dailymax_timeseries_([0-9]{4})\\.qs")[,2])
#     existing_years <- existing_years[is.finite(existing_years)]
#   }
#   existing_years <- sort(unique(existing_years))
#   
#   if (!OVERWRITE_YEARS) {
#     years_u <- setdiff(years_u, existing_years)
#     message(sprintf("[SCAN] existing years=%d | remaining=%d",
#                     length(existing_years), length(years_u)))
#   }
#   if (length(years_u) == 0) {
#     message("[DONE] Nothing to do: all yearly outputs exist.")
#     return(invisible(TRUE))
#   }
#   
#   # Plan
#   if (isTRUE(PARALLEL)) {
#     message(sprintf("[MODE] parallel multisession | workers=%d", WORKERS))
#     future::plan(multisession, workers = WORKERS)
#   } else {
#     message("[MODE] sequential")
#     future::plan(sequential)
#   }
#   
#   progressr::handlers(global = TRUE)
#   progressr::handlers("txtprogressbar")
#   
#   res <- progressr::with_progress({
#     p <- progressor(along = years_u)
#     future_lapply(years_u, function(y) {
#       on.exit(p(sprintf("year=%d", y)), add = TRUE)
#       
#       idx <- by_year_idx[[as.character(y)]]
#       
#       # per-year error isolation (won't kill whole run)
#       tryCatch(
#         process_one_year(
#           year = y,
#           files_year = files[idx],
#           days_year  = days[idx],
#           grid = grid
#         ),
#         error = function(e) {
#           list(year = y, status = "error", file = out_year_file(y), error = conditionMessage(e))
#         }
#       )
#     })
#   })
#   
#   res_dt <- rbindlist(res, fill=TRUE)
#   message("\n[STAGE A] Summary:")
#   print(res_dt[, .N, by=status][order(-N)])
#   
#   if (res_dt[status == "error", .N] > 0) {
#     message("\n[ERROR YEARS]")
#     print(res_dt[status == "error", .(year, error)][order(year)])
#   }
#   
#   message("[OUT] ", OUT_TS_DIR)
#   invisible(res_dt)
# }
# 
# # ---- RUN ----
# # run_stageA_heatindex_dailymax()
# # PARALLEL <- FALSE
# # run_stageA_heatindex_dailymax()
# 
# PARALLEL <- TRUE
# WORKERS <- 20
# run_stageA_heatindex_dailymax()
# 
# # Calc rps ----
# # Heat Index IDF (3-day duration) — CONUS-correct, end-to-end
# # Reference grid MUST be CONUS: 1940_01_era5_temperature_conus.nc
# #
# # Stage A: build per-year DAILY MAX heat index matrices on CONUS grid
# # Stage B: per-cell IDF on annual maxima of 3-day trailing mean of daily max (hi3)
# # Stage C: compile NetCDF stack on CONUS grid
# #
# # INPUT daily HI files:
# #   IN_HI_DIR/YYYY-MM-DD.qs with lon, lat, datetime_utc, heat_index
# #
# # OUTPUT:
# #   Stage A yearly: OUT_TS_DIR/hi_dailymax_timeseries_YYYY.qs
# #   Stage B chunks: OUT_RP_CHUNKS_DIR/rp_hi3_chunk_####.qs
# #   Stage C NetCDF: OUT_NC_CONUS
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
# # USER SETTINGS
# ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
# IN_HI_DIR <- file.path(ROOT, "data/out/01_vars/heatindex")
# 
# # CONUS reference (this fixes your extent issue)
# # REF_NC_CONUS <- "V:/users/hquintal/phd1_cluster_southeast/data/input/temperature/1940_01_era5_temperature_conus.nc"
# REF_NC_CONUS <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation/_mask_target.tif"
# test <- terra::rast(REF_NC_CONUS)
# test
# plot(test)
# OUT_DIR <- file.path(ROOT, "data/out/04_rp")
# OUT_TS_DIR <- file.path(OUT_DIR, "heatindex_dailymax_timeseries_by_year")
# OUT_RP_CHUNKS_DIR <- file.path(OUT_DIR, "heatindex_idf_hi3_return_periods_chunks")
# OUT_NC_CONUS <- file.path(OUT_DIR, "heatindex_idf_hi3_return_periods.nc")
# dir_create(OUT_TS_DIR); dir_create(OUT_RP_CHUNKS_DIR)
# 
# LOOKUP_QS <- file.path(OUT_DIR, "heatindex_refgrid_lonlat_to_cell.qs")
# 
# # Controls
# RUN_STAGE_A <- TRUE
# RUN_STAGE_B <- TRUE
# RUN_STAGE_C <- TRUE
# 
# OVERWRITE_LOOKUP <- FALSE
# OVERWRITE_STAGE_A <- FALSE
# OVERWRITE_STAGE_B <- TRUE
# OVERWRITE_STAGE_C <- TRUE
# 
# WORKERS_A <- 4L     # year-level IO heavy; keep conservative
# WORKERS_B <- 6L     # chunk fits
# CHUNK_SIZE <- 300L
# 
# QS_PRESET <- "high"
# KEY_DIGITS <- 6
# kfmt <- function(x, digits = KEY_DIGITS) sprintf(paste0("%.", digits, "f"), x)
# 
# # Return periods
# RP_LABELS  <- c(1, 2, 5, 10, 25, 50, 100, 200, 500, 1000)
# RP_COMPUTE <- c(1.01, 2, 5, 10, 25, 50, 100, 200, 500, 1000)
# 
# # Fit settings
# MIN_YEARS_REQUIRED <- 20L
# PVAL_THRESHOLD     <- 0.4
# CANDIDATE_DISTS    <- c("cosine", "hypsecant", "logistic", "norm", "uniform")
# 
# # Small utilities
# is_daily_qs <- function(p) grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}\\.qs$", fs::path_file(p))
# day_from_file_chr <- function(p) sub("\\.qs$", "", fs::path_file(p))
# year_from_day <- function(d) as.integer(format(as.Date(d), "%Y"))
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
# # Lookup: CONUS ref grid (cell index + lon/lat + join keys)
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
# # Stage A: Yearly DAILY MAX heat index matrices on CONUS grid
# daily_max_vec_one_day <- function(day_file, grid) {
#   dt <- as.data.table(qread(day_file))
#   dt[, `:=`(lon=as.numeric(lon), lat=as.numeric(lat), heat_index=as.numeric(heat_index))]
#   dt <- dt[is.finite(heat_index)]
#   if (nrow(dt) == 0) return(rep(NA_real_, nrow(grid)))
#   dm <- dt[, .(hi_max = max(heat_index, na.rm = TRUE)), by=.(lon,lat)]
#   dm[, `:=`(lon_key=kfmt(lon), lat_key=kfmt(lat))]
#   setkey(dm, lon_key, lat_key)
#   m <- dm[grid, on=.(lon_key, lat_key)]
#   as.numeric(m$hi_max)
# }
# 
# stageA_process_one_year <- function(y, files_year, days_year, grid) {
#   out_qs <- file.path(OUT_TS_DIR, sprintf("hi_dailymax_timeseries_%d.qs", y))
#   if (!OVERWRITE_STAGE_A && fs::file_exists(out_qs)) return(list(year=y, status="skip_exists", file=out_qs))
#   cleanup_tmp_files(out_qs)
#   
#   n_cells <- nrow(grid); n_days <- length(files_year)
#   mat <- matrix(NA_real_, nrow=n_cells, ncol=n_days)
#   for (j in seq_len(n_days)) {
#     v <- daily_max_vec_one_day(files_year[j], grid)
#     if (length(v) != n_cells) stop("Day vec length mismatch for ", files_year[j])
#     mat[, j] <- v
#     if (j %% 25 == 0) gc()
#   }
#   
#   obj <- list(
#     meta = list(variable="heat_index", statistic="daily_max", units="C", year=y, ref_nc=REF_NC_CONUS),
#     grid = grid[, .(cell, lon, lat)],
#     days = as.Date(days_year),
#     hi_max = float::fl(mat)
#   )
#   atomic_qsave(obj, out_qs, preset=QS_PRESET)
#   rm(obj, mat); gc()
#   list(year=y, status="ok", file=out_qs)
# }
# 
# run_stageA <- function(grid) {
#   files <- sort(as.character(fs::dir_ls(IN_HI_DIR, type="file", regexp="\\.qs$", fail=FALSE)))
#   files <- files[is_daily_qs(files)]
#   if (length(files) == 0) stop("No daily HI .qs files in: ", IN_HI_DIR)
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
#   existing <- fs::dir_ls(OUT_TS_DIR, type="file", regexp="^hi_dailymax_timeseries_[0-9]{4}\\.qs$", fail=FALSE)
#   done_years <- integer(0)
#   if (length(existing) > 0) done_years <- as.integer(str_match(fs::path_file(existing), "([0-9]{4})")[,1])
#   done_years <- done_years[is.finite(done_years)]
#   todo <- years_u
#   if (!OVERWRITE_STAGE_A) todo <- setdiff(todo, done_years)
#   
#   message(sprintf("[Stage A] years total=%d | done=%d | todo=%d | n_cells(CONUS)=%d",
#                   length(years_u), length(done_years), length(todo), nrow(grid)))
#   if (length(todo) == 0) return(invisible(TRUE))
#   
#   future::plan(multisession, workers=WORKERS_A)
#   progressr::handlers(global=TRUE); progressr::handlers("txtprogressbar")
#   
#   res <- progressr::with_progress({
#     p <- progressr::progressor(along=todo)
#     future_lapply(todo, function(y) {
#       on.exit(p(sprintf("year=%d", y)), add=TRUE)
#       idx <- by_year[[as.character(y)]]
#       tryCatch(stageA_process_one_year(y, files[idx], days[idx], grid),
#                error=function(e) list(year=y, status="error", error=conditionMessage(e)))
#     })
#   })
#   print(rbindlist(res, fill=TRUE)[, .N, by=status][order(-N)])
#   invisible(TRUE)
# }
# 
# # Stage B: IDF on annual maxima of 3-day trailing mean of daily max (hi3)
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
# hi3_trailing_mean_strict <- function(xmat) {
#   nr <- nrow(xmat); nt <- ncol(xmat)
#   out <- matrix(NA_real_, nr, nt)
#   if (nt < 3) return(out)
#   a <- xmat[,1:(nt-2),drop=FALSE]
#   b <- xmat[,2:(nt-1),drop=FALSE]
#   c <- xmat[,3:nt,drop=FALSE]
#   ok <- is.finite(a)&is.finite(b)&is.finite(c)
#   out[,3:nt] <- ifelse(ok, (a+b+c)/3, NA_real_)
#   out
# }
# row_max_na <- function(m) apply(m, 1, function(v) if (all(!is.finite(v))) NA_real_ else max(v, na.rm=TRUE))
# 
# run_stageB <- function() {
#   
#   # Robust Stage A discovery (no regexp in dir_ls)
#   all_qs <- fs::dir_ls(OUT_TS_DIR, type = "file", regexp = "\\.qs$", fail = FALSE)
#   all_qs <- sort(as.character(all_qs))
#   year_files <- all_qs[grepl("^hi_dailymax_timeseries_[0-9]{4}\\.qs$", fs::path_file(all_qs))]
#   
#   if (length(year_files) == 0) {
#     message("[DEBUG] OUT_TS_DIR = ", OUT_TS_DIR)
#     message("[DEBUG] files (first 25):\n", paste(head(fs::path_file(all_qs), 25), collapse = "\n"))
#     stop("No Stage A files in: ", OUT_TS_DIR)
#   }
#   
#   years <- as.integer(stringr::str_match(fs::path_file(year_files),
#                                          "^hi_dailymax_timeseries_([0-9]{4})\\.qs$")[,2])
#   ord <- order(years)
#   year_files <- year_files[ord]
#   years <- years[ord]
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
#   # Resume scan (robust)
#   all_chunks <- fs::dir_ls(OUT_RP_CHUNKS_DIR, type="file", regexp="\\.qs$", fail=FALSE)
#   all_chunks <- sort(as.character(all_chunks))
#   done_files <- all_chunks[grepl("^rp_hi3_chunk_[0-9]{4}\\.qs$", fs::path_file(all_chunks))]
#   
#   done_ids <- integer(0)
#   if (length(done_files) > 0) {
#     done_ids <- as.integer(stringr::str_match(fs::path_file(done_files),
#                                               "^rp_hi3_chunk_([0-9]{4})\\.qs$")[,2])
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
#     out_qs <- file.path(OUT_RP_CHUNKS_DIR, sprintf("rp_hi3_chunk_%04d.qs", cid))
#     if (!OVERWRITE_STAGE_B && fs::file_exists(out_qs)) return(list(chunk=cid,status="skip_exists"))
#     
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
#       # Map by cell id (safe)
#       g <- as.data.table(obj$grid); g[, cell := as.integer(cell)]
#       pos <- match(grid_chunk$cell, g$cell)
#       if (anyNA(pos)) stop("Grid mismatch in year ", years[j])
#       
#       xmat <- as.matrix(obj$hi_max)[pos, , drop=FALSE]
#       hi3  <- hi3_trailing_mean_strict(xmat)
#       ann_max[, j] <- row_max_na(hi3)
#       
#       rm(obj, g, pos, xmat, hi3); gc()
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
#   # Parallel (no progressr to avoid the handler warning)
#   future::plan(multisession, workers = WORKERS_B)
#   res <- future.apply::future_lapply(todo, process_chunk, future.seed = TRUE)
#   
#   res_dt <- rbindlist(res, fill=TRUE)
#   print(res_dt[, .N, by=status][order(-N)])
#   invisible(TRUE)
# }
# 
# # Stage C: Compile NetCDF on CONUS ref grid (uses cell indices on CONUS grid)
# run_stageC <- function() {
#   
#   # robust scan (no regexp in dir_ls)
#   all_files <- fs::dir_ls(OUT_RP_CHUNKS_DIR, type = "file", fail = FALSE)
#   all_files <- sort(as.character(all_files))
#   
#   rp_files <- all_files[grepl("^rp_hi3_chunk_[0-9]{4}\\.qs$", fs::path_file(all_files))]
#   
#   message(sprintf("[Stage C] found %d rp chunk files in: %s", length(rp_files), OUT_RP_CHUNKS_DIR))
#   if (length(rp_files) == 0) {
#     message("[DEBUG] first 25 files actually in dir:\n",
#             paste(head(fs::path_file(all_files), 25), collapse = "\n"))
#     stop("No rp chunk files in: ", OUT_RP_CHUNKS_DIR)
#   }
#   
#   # CONUS reference grid
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
#     
#     miss <- setdiff(required_cols, names(dt))
#     if (length(miss) > 0) stop("Missing columns in ", fs::path_file(f), ": ", paste(miss, collapse = ", "))
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
#     varname   = "heatindex_idf_hi3",
#     unit      = "C"
#   )
#   
#   message("[WRITE] ", OUT_NC_CONUS)
#   invisible(OUT_NC_CONUS)
# }
# 
# # RUN ALL
# run_all <- function() {
#   stopifnot(fs::file_exists(REF_NC_CONUS))
#   grid <- build_lookup_from_conus_ref()
#   
#   # quick extent print (should be CONUS)
#   message("[REF CONUS] ", paste(as.vector(terra::ext(terra::rast(REF_NC_CONUS)[[1]])), collapse=", "))
#   
#   if (RUN_STAGE_A) run_stageA(grid)
#   if (RUN_STAGE_B) run_stageB()
#   if (RUN_STAGE_C) run_stageC()
#   invisible(TRUE)
# }
# 
# # ---- RUN ----
# run_all()
# 
# # run_stageB()
# # run_stageC()
# 
# test <- terra::rast('V:/users/hquintal/phd1_cluster_southeast_v2/data/out/04_rp/heatindex_idf_hi3_return_periods.nc')
# plot(test[[2]])

# Southeast ----
# =============================================================================
# Heat Index IDF (1-day duration) — BBOX-cropped, end-to-end
#
# Updates you requested:
#  1) Crop + mask HEAT INDEX workflow to:
#       BBOX_LONLAT <- c(xmin = -95, xmax = -75, ymin = 22, ymax = 40)
#     using REF_NC as the mask/reference grid (crop to bbox; preserve NA mask).
#  2) Compute MAX DAILY heat index per location (same as before, but bbox-filtered).
#  3) Compute RETURN PERIODS for 1-day heat index (annual maxima of daily max),
#     replacing the current 3-day (hi3) logic.
#
# Outputs:
#  - Stage A yearly:  OUT_TS_DIR/hi_dailymax_timeseries_YYYY.qs
#  - Stage B chunks:  OUT_RP_CHUNKS_DIR/rp_hi1_chunk_####.qs
#  - Stage C NetCDF:  OUT_NC_BBOX  (stack of rp layers on bbox ref grid)
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

# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

IN_HI_DIR <- file.path(ROOT, "data/out/01_vars/heatindex")

# Reference mask/grid (lon/lat)
REF_NC <- file.path(ROOT, "data/out/01_vars/precipitation/_mask_target.tif")
# test <- terra::rast(REF_NC)
# plot(test)
# BBOX in lon/lat
BBOX_LONLAT <- c(xmin = -95, xmax = -75, ymin = 22, ymax = 40)

OUT_DIR <- file.path(ROOT, "data/out/04_rp")

# NOTE: rename outputs to reflect 1-day IDF (hi1)
OUT_TS_DIR        <- file.path(OUT_DIR, "heatindex_dailymax_timeseries_by_year_bbox_hi1")
OUT_RP_CHUNKS_DIR <- file.path(OUT_DIR, "heatindex_idf_hi1_return_periods_chunks_bbox")
OUT_NC_BBOX       <- file.path(OUT_DIR, "heatindex_idf_hi1_return_periods_bbox.nc")

dir_create(OUT_TS_DIR)
dir_create(OUT_RP_CHUNKS_DIR)

LOOKUP_QS <- file.path(OUT_DIR, "heatindex_refgrid_lonlat_to_cell_bbox_hi1.qs")

# Controls
RUN_STAGE_A <- TRUE
RUN_STAGE_B <- TRUE
RUN_STAGE_C <- TRUE

OVERWRITE_LOOKUP   <- TRUE
OVERWRITE_STAGE_A  <- TRUE
OVERWRITE_STAGE_B  <- TRUE
OVERWRITE_STAGE_C  <- TRUE

# Parallelism
WORKERS_A  <- 10L     # IO-heavy stage
WORKERS_B  <- 10L     # compute-heavy stage
CHUNK_SIZE <- 300L

QS_PRESET  <- "high"
KEY_DIGITS <- 6
options(future.globals.maxSize = 8 * 1024^3)

# Return periods
RP_LABELS  <- c(1, 2, 5, 10, 25, 50, 100, 200, 500, 1000)
RP_COMPUTE <- c(1.01, 2, 5, 10, 25, 50, 100, 200, 500, 1000)

# Fit settings
MIN_YEARS_REQUIRED <- 20L
PVAL_THRESHOLD     <- 0.4
CANDIDATE_DISTS    <- c("cosine", "hypsecant", "logistic", "norm", "uniform")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
kfmt <- function(x, digits = KEY_DIGITS) sprintf(paste0("%.", digits, "f"), x)

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

# -----------------------------------------------------------------------------
# Stage 0: Build bbox-cropped reference raster + lookup
# -----------------------------------------------------------------------------
get_ref_bbox_raster <- function() {
  stopifnot(fs::file_exists(REF_NC))
  r0 <- terra::rast(REF_NC)[[1]]
  
  # Crop to bbox; keep NA mask (mask is already in r0 values)
  e <- terra::ext(BBOX_LONLAT["xmin"], BBOX_LONLAT["xmax"], BBOX_LONLAT["ymin"], BBOX_LONLAT["ymax"])
  r1 <- terra::crop(r0, e, snap = "out")  # snap out for safety on boundaries
  
  # Optional hard-trim to exact bbox (crop already does), but keep as-is
  # r1 <- terra::trim(r1) # DON'T trim; would shrink to non-NA which you don't want
  
  r1
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
  
  message("[LOOKUP] Building lookup from REF_NC bbox grid: ", REF_NC)
  r <- get_ref_bbox_raster()
  
  v <- terra::values(r, mat = FALSE)
  cells <- which(!is.na(v))
  if (length(cells) == 0) stop("REF bbox raster has 0 non-NA cells (mask empty in bbox).")
  
  xy <- terra::xyFromCell(r, cells)
  grid <- data.table(
    cell = as.integer(cells),         # cell index into the *bbox raster* used later in Stage C
    lon  = as.numeric(xy[,1]),
    lat  = as.numeric(xy[,2])
  )
  grid[, `:=`(lon_key = kfmt(lon), lat_key = kfmt(lat))]
  setorder(grid, cell)
  setkey(grid, lon_key, lat_key)
  
  qsave(grid, LOOKUP_QS, preset = QS_PRESET)
  message("[LOOKUP] Wrote: ", LOOKUP_QS, " | n_cells=", nrow(grid))
  grid
}

# -----------------------------------------------------------------------------
# Stage A: Daily max time series by year (bbox-filtered)
# -----------------------------------------------------------------------------
daily_max_vec_one_day_bbox <- function(day_file, grid) {
  dt <- as.data.table(qread(day_file))
  if (!all(c("lon","lat","heat_index") %in% names(dt))) {
    stop("Missing columns in: ", day_file, " (need lon, lat, heat_index)")
  }
  
  dt[, `:=`(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    heat_index = as.numeric(heat_index)
  )]
  
  # BBOX filter first (fast), then finite
  dt <- dt[
    lon >= BBOX_LONLAT["xmin"] & lon <= BBOX_LONLAT["xmax"] &
      lat >= BBOX_LONLAT["ymin"] & lat <= BBOX_LONLAT["ymax"]
  ]
  dt <- dt[is.finite(heat_index)]
  
  if (nrow(dt) == 0) return(rep(NA_real_, nrow(grid)))
  
  dm <- dt[, .(hi_max = max(heat_index, na.rm = TRUE)), by = .(lon, lat)]
  dm[, `:=`(lon_key = kfmt(lon), lat_key = kfmt(lat))]
  setkey(dm, lon_key, lat_key)
  
  # dm[grid] preserves grid order/length
  m <- dm[grid, on = .(lon_key, lat_key)]
  as.numeric(m$hi_max)
}

stageA_process_one_year <- function(y, files_year, days_year, grid) {
  out_qs <- file.path(OUT_TS_DIR, sprintf("hi_dailymax_timeseries_%d.qs", y))
  if (!OVERWRITE_STAGE_A && fs::file_exists(out_qs)) return(list(year = y, status = "skip_exists", file = out_qs))
  cleanup_tmp_files(out_qs)
  
  n_cells <- nrow(grid)
  n_days  <- length(files_year)
  if (n_days == 0) return(list(year = y, status = "skip_empty", file = out_qs))
  
  mat <- matrix(NA_real_, nrow = n_cells, ncol = n_days)
  for (j in seq_len(n_days)) {
    v <- daily_max_vec_one_day_bbox(files_year[j], grid)
    if (length(v) != n_cells) stop("Day vec length mismatch for ", files_year[j])
    mat[, j] <- v
    if (j %% 25 == 0) gc()
  }
  
  obj <- list(
    meta = list(
      variable    = "heat_index",
      statistic   = "daily_max",
      duration    = "1-day",
      units       = "C",
      year        = y,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      ref_mask    = REF_NC,
      bbox_lonlat = BBOX_LONLAT,
      key_digits  = KEY_DIGITS
    ),
    grid  = grid[, .(cell, lon, lat)], # cell index into bbox ref raster
    days  = as.Date(days_year),
    hi_max = float::fl(mat)
  )
  
  atomic_qsave(obj, out_qs, preset = QS_PRESET)
  rm(obj, mat); gc()
  list(year = y, status = "ok", file = out_qs)
}

run_stageA <- function(grid) {
  message("\n==================== STAGE A: Daily max HI by year (BBOX) ====================")
  
  files <- sort(as.character(fs::dir_ls(IN_HI_DIR, type = "file", glob = "*.qs", fail = FALSE)))
  files <- files[is_daily_qs(files)]
  if (length(files) == 0) stop("No daily HI .qs files in: ", IN_HI_DIR)
  
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
    done_years <- as.integer(str_match(fs::path_file(existing), "^hi_dailymax_timeseries_([0-9]{4})\\.qs$")[,2])
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
    future_lapply(todo, function(y) {
      on.exit(p(sprintf("year=%d", y)), add = TRUE)
      idx <- by_year[[as.character(y)]]
      tryCatch(
        stageA_process_one_year(y, files[idx], days[idx], grid),
        error = function(e) list(year = y, status = "error", file = file.path(OUT_TS_DIR, sprintf("hi_dailymax_timeseries_%d.qs", y)),
                                 error = conditionMessage(e))
      )
    })
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

# -----------------------------------------------------------------------------
# Stage B: IDF on ANNUAL MAXIMA of DAILY MAX (hi1)
# -----------------------------------------------------------------------------
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
  message("\n==================== STAGE B: Return periods (hi1 annual maxima) ====================")
  
  # Stage A discovery
  all_qs <- fs::dir_ls(OUT_TS_DIR, type = "file", regexp = "\\.qs$", fail = FALSE)
  all_qs <- sort(as.character(all_qs))
  year_files <- all_qs[grepl("^hi_dailymax_timeseries_[0-9]{4}\\.qs$", fs::path_file(all_qs))]
  if (length(year_files) == 0) stop("No Stage A files in: ", OUT_TS_DIR)
  
  years <- as.integer(str_match(fs::path_file(year_files), "^hi_dailymax_timeseries_([0-9]{4})\\.qs$")[,2])
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
  
  # Resume scan
  all_chunks <- fs::dir_ls(OUT_RP_CHUNKS_DIR, type = "file", regexp = "\\.qs$", fail = FALSE)
  all_chunks <- sort(as.character(all_chunks))
  done_files <- all_chunks[grepl("^rp_hi1_chunk_[0-9]{4}\\.qs$", fs::path_file(all_chunks))]
  done_ids <- integer(0)
  if (length(done_files) > 0) {
    done_ids <- as.integer(str_match(fs::path_file(done_files), "^rp_hi1_chunk_([0-9]{4})\\.qs$")[,2])
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
    out_qs <- file.path(OUT_RP_CHUNKS_DIR, sprintf("rp_hi1_chunk_%04d.qs", cid))
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
      
      # Convert float32 -> base double matrix FIRST (robust), then subset rows
      xmat_all <- float::dbl(obj$hi_max)              # [n_cells x n_days] base matrix
      xmat     <- xmat_all[pos, , drop = FALSE]       # [n_chunk x n_days]
      
      # HI1: annual maxima of DAILY MAX directly (no 3-day trailing mean)
      ann_max[, j] <- row_max_na(xmat)
      
      rm(xmat_all, xmat)
      
      rm(obj, g, pos, xmat); gc()
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

# -----------------------------------------------------------------------------
# Stage C: Compile NetCDF on bbox ref grid (cell indices refer to bbox raster)
# -----------------------------------------------------------------------------
run_stageC <- function() {
  message("\n==================== STAGE C: NetCDF compile (bbox) ====================")
  
  all_files <- fs::dir_ls(OUT_RP_CHUNKS_DIR, type = "file", fail = FALSE)
  all_files <- sort(as.character(all_files))
  rp_files <- all_files[grepl("^rp_hi1_chunk_[0-9]{4}\\.qs$", fs::path_file(all_files))]
  
  message(sprintf("[Stage C] found %d rp chunk files in: %s", length(rp_files), OUT_RP_CHUNKS_DIR))
  if (length(rp_files) == 0) stop("No rp chunk files in: ", OUT_RP_CHUNKS_DIR)
  
  # bbox reference raster (cropped to bbox, preserves mask)
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
    varname   = "heatindex_idf_hi1",
    unit      = "C"
  )
  
  message("[WRITE] ", OUT_NC_BBOX)
  invisible(OUT_NC_BBOX)
}

# -----------------------------------------------------------------------------
# RUN ALL
# -----------------------------------------------------------------------------
run_all <- function() {
  stopifnot(fs::file_exists(REF_NC))
  
  # Build lookup on bbox-cropped ref
  grid <- build_lookup_from_ref_bbox()
  
  # Print bbox ref extent
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

