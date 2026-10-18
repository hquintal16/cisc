# # =============================================================================
# # Cluster exceedances BY YEAR for 6 datasets (heatindex + precipitation)
# #
# # KEY DESIGN:
# #   1) Keep ONLY the needed columns for clustering:
# #        x, y, time_utc, intensity, rp_max_years (+ optional metadata)
# #   2) Build a NORMALIZED integer cube:
# #        row_ix = index of y in sorted unique(y)
# #        col_ix = index of x in sorted unique(x)
# #        t_ix   = integer timestep index at native cadence (e.g., 24 hr)
# #      -> Adjacent timesteps are distance 1
# #      -> Gaps in time produce large distances (e.g. day 1 -> day 10 => 9)
# #   3) Apply RP floor BEFORE clustering.
# #
# # Output:
# #   V:/users/hquintal/phd1_cluster_southeast_v2/data/out/06_clusters/hi{HI_FLOOR}_pr{PR_FLOOR}/{dataset_tag}
# # =============================================================================
# 
# suppressPackageStartupMessages({
#   library(data.table)
#   library(qs)
#   library(fs)
#   library(dbscan)
#   library(lubridate)
#   library(progressr)
#   library(future)
#   library(future.apply)
#   library(stringr)
# })
# 
# # ── CONFIG ───────────────────────────────────────────────────────────────────
# MU  <- 4L
# EPS <- 4.0
# 
# BASE <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out"
# 
# HI_FLOOR_YR <- 25
# PR_FLOOR_YR <- 25
# 
# OUT_BASE <- fs::path(BASE, sprintf("06_clusters/hi%d_pr%d", HI_FLOOR_YR, PR_FLOOR_YR))
# dir_create(OUT_BASE, recurse = TRUE)
# 
# DATASETS <- list(
#   # -------------------------
#   # PRECIP (PR24)
#   # -------------------------
#   list(tag="precipitation_0.25",
#        in_dir=fs::path(BASE,"05_exceedances/precipitation/pr24/0.25"),
#        out_dir=fs::path(OUT_BASE,"precipitation_0.25"),
#        time_step_hours=24L, floor_yr=PR_FLOOR_YR),
# 
#   list(tag="precipitation_0.17_bilinear",
#        in_dir=fs::path(BASE,"05_exceedances/precipitation/pr24/0.17_bilinear"),
#        out_dir=fs::path(OUT_BASE,"precipitation_0.17_bilinear"),
#        time_step_hours=24L, floor_yr=PR_FLOOR_YR),
# 
#   list(tag="precipitation_0.17_max",
#        in_dir=fs::path(BASE,"05_exceedances/precipitation/pr24/0.17_max"),
#        out_dir=fs::path(OUT_BASE,"precipitation_0.17_max"),
#        time_step_hours=24L, floor_yr=PR_FLOOR_YR),
# 
#   # NEW: precip 0.12°
#   list(tag="precipitation_0.12_bilinear",
#        in_dir=fs::path(BASE,"05_exceedances/precipitation/pr24/0.12_bilinear"),
#        out_dir=fs::path(OUT_BASE,"precipitation_0.12_bilinear"),
#        time_step_hours=24L, floor_yr=PR_FLOOR_YR),
#   
#   list(tag="precipitation_0.12_max",
#        in_dir=fs::path(BASE,"05_exceedances/precipitation/pr24/0.12_max"),
#        out_dir=fs::path(OUT_BASE,"precipitation_0.12_max"),
#        time_step_hours=24L, floor_yr=PR_FLOOR_YR),
#   
#   
#   # -------------------------
#   # HEAT INDEX
#   # -------------------------
#   list(tag="heatindex_0.25",
#        in_dir=fs::path(BASE,"05_exceedances/heatindex/0.25"),
#        out_dir=fs::path(OUT_BASE,"heatindex_0.25"),
#        time_step_hours=24L, floor_yr=HI_FLOOR_YR),
# 
#   list(tag="heatindex_0.44_bilinear",
#        in_dir=fs::path(BASE,"05_exceedances/heatindex/0.44_bilinear"),
#        out_dir=fs::path(OUT_BASE,"heatindex_0.44_bilinear"),
#        time_step_hours=24L, floor_yr=HI_FLOOR_YR),
# 
#   list(tag="heatindex_0.44_max",
#        in_dir=fs::path(BASE,"05_exceedances/heatindex/0.44_max"),
#        out_dir=fs::path(OUT_BASE,"heatindex_0.44_max"),
#        time_step_hours=24L, floor_yr=HI_FLOOR_YR),
#   
#   # NEW: HI 0.34°
#   list(tag="heatindex_0.34_bilinear",
#        in_dir=fs::path(BASE,"05_exceedances/heatindex/0.34_bilinear"),
#        out_dir=fs::path(OUT_BASE,"heatindex_0.34_bilinear"),
#        time_step_hours=24L, floor_yr=HI_FLOOR_YR),
#   
#   list(tag="heatindex_0.34_max",
#        in_dir=fs::path(BASE,"05_exceedances/heatindex/0.34_max"),
#        out_dir=fs::path(OUT_BASE,"heatindex_0.34_max"),
#        time_step_hours=24L, floor_yr=HI_FLOOR_YR)
# )
# 
# WORKERS <- 20L
# plan(multisession, workers = WORKERS)
# handlers(global = TRUE)
# 
# # ── HELPERS ──────────────────────────────────────────────────────────────────
# .log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
#                           sprintf(...), "\n", sep="")
# 
# .ensure_posix_utc <- function(tt) {
#   if (inherits(tt, "POSIXct")) return(as.POSIXct(tt, tz="UTC"))
#   if (inherits(tt, "Date"))   return(as.POSIXct(tt, tz="UTC"))
#   as.POSIXct(tt, tz="UTC", origin="1970-01-01")
# }
# 
# .get_time_utc <- function(dt) {
#   if ("time_utc" %in% names(dt)) return(.ensure_posix_utc(dt$time_utc))
#   if ("time" %in% names(dt)) return(.ensure_posix_utc(dt$time))
#   if ("datetime_utc" %in% names(dt)) return(.ensure_posix_utc(dt$datetime_utc))
#   if ("datetime" %in% names(dt)) return(.ensure_posix_utc(dt$datetime))
#   if ("day" %in% names(dt)) return(as.POSIXct(as.Date(dt$day), tz="UTC"))
#   if ("date" %in% names(dt)) return(as.POSIXct(as.Date(dt$date), tz="UTC"))
#   stop("No time column found (need time_utc/time/datetime/day/date).")
# }
# 
# .standardize_xy <- function(dt) {
#   if (all(c("x","y") %in% names(dt))) return(invisible(NULL))
#   if (all(c("lon","lat") %in% names(dt))) { setnames(dt,c("lon","lat"),c("x","y")); return(invisible(NULL)) }
#   if (all(c("longitude","latitude") %in% names(dt))) { setnames(dt,c("longitude","latitude"),c("x","y")); return(invisible(NULL)) }
#   stop("Missing x/y or lon/lat.")
# }
# 
# # intensity column: keep existing; we don't need to rename here
# .guess_intensity_col <- function(dt, tag) {
#   if (grepl("^precipitation", tag, ignore.case=TRUE)) {
#     for (cc in c("pr24_intensity_mm","intensity_mm","mm","precip_mm","prcp_mm","pr","ppt","value","intensity")) {
#       if (cc %in% names(dt)) return(cc)
#     }
#   } else {
#     for (cc in c("hi3_intensity","hi_intensity","heatindex","HeatIndex","value","intensity")) {
#       if (cc %in% names(dt)) return(cc)
#     }
#   }
#   # last resort: a single remaining numeric column
#   ignore <- c("x","y","lon","lat","time","time_utc","datetime","datetime_utc","day","date",
#               "rp","rp_max_years","rp_max_layer","resolution_deg")
#   num_cols <- names(dt)[vapply(dt, is.numeric, TRUE)]
#   num_cols <- setdiff(num_cols, ignore)
#   if (length(num_cols)==1) return(num_cols[1])
#   stop("Cannot infer intensity column for ", tag)
# }
# 
# .get_rp_col <- function(dt) {
#   if ("rp_max_years" %in% names(dt)) return("rp_max_years")
#   if ("rp" %in% names(dt)) return("rp")
#   stop("Missing return period column (need rp_max_years or rp).")
# }
# 
# .apply_floor <- function(dt, rp_col, floor_yr) {
#   rp <- suppressWarnings(as.numeric(dt[[rp_col]]))
#   dt[is.finite(rp) & rp >= floor_yr]
# }
# 
# # NORMALIZED integer indices: adjacent grid/time steps => distance 1
# .make_normalized_cube <- function(dt, time_step_hours) {
#   # spatial index by rank in sorted unique coords
#   ux <- sort(unique(dt$x))
#   uy <- sort(unique(dt$y))
#   x_map <- data.table(x = ux, col_ix = seq_along(ux))
#   y_map <- data.table(y = uy, row_ix = seq_along(uy))
#   dt <- x_map[dt, on="x"]
#   dt <- y_map[dt, on="y"]
#   
#   # time index: anchored to file-min time, in units of native timestep hours
#   t0 <- min(dt$time_utc, na.rm=TRUE)
#   dh <- as.numeric(difftime(dt$time_utc, t0, units="hours"))
#   dt[, t_ix := as.integer(round(dh / time_step_hours)) + 1L]
#   
#   cube <- as.matrix(dt[, .(row_ix, col_ix, t_ix)])
#   list(dt = dt, cube = cube)
# }
# 
# .add_cluster_id <- function(dt) {
#   dt[, cluster_id := NA_character_]
#   lab <- dt[cluster > 0]
#   if (!nrow(lab)) return(dt)
#   
#   meta <- lab[, {
#     ft <- min(time_utc, na.rm=TRUE)
#     sub <- .SD[time_utc == ft]
#     data.table(first_time=ft,
#                first_y=min(y, na.rm=TRUE),
#                first_x=min(x, na.rm=TRUE))
#   }, by=cluster]
#   
#   meta[, date := format(as.Date(first_time, tz="UTC"), "%Y-%m-%d")]
#   setorder(meta, date, first_time, first_y, first_x, cluster)
#   meta[, seq_id := ave(seq_len(.N), date, FUN=seq_along)]
#   meta[, cluster_id := sprintf("%s-%03d", date, seq_id)]
#   
#   dt[meta, on="cluster", cluster_id := i.cluster_id]
#   dt[]
# }
# 
# .get_year_from_filename <- function(path) {
#   bn <- fs::path_file(path)
#   yy <- suppressWarnings(as.integer(str_match(bn, "(19|20)\\d{2}")[,1]))
#   if (is.na(yy)) stop("Could not parse year from filename: ", bn)
#   yy
# }
# 
# .discover_year_files <- function(in_dir) {
#   files <- fs::dir_ls(in_dir, type="file", glob="*.qs", fail=FALSE)
#   files <- sort(as.character(files))
#   files[grepl("(19|20)\\d{2}", fs::path_file(files))]
# }
# 
# .cluster_one_year_file <- function(in_path, out_path, tag, time_step_hours, floor_yr) {
#   
#   dt0 <- as.data.table(qread(in_path))
#   if (!nrow(dt0)) {
#     qsave(data.table(), out_path)
#     return(list(n=0L, n_clusters=0L, status="empty"))
#   }
#   
#   .standardize_xy(dt0)
#   dt0[, time_utc := .get_time_utc(dt0)]
#   
#   rp_col <- .get_rp_col(dt0)
#   inten_col <- .guess_intensity_col(dt0, tag)
#   
#   # Keep only the columns we actually need (+ helpful provenance)
#   keep <- intersect(c("x","y","time_utc", inten_col, rp_col, "rp_max_layer", "resolution_deg"), names(dt0))
#   dt <- dt0[, ..keep]
#   setnames(dt, inten_col, "intensity")
#   if (rp_col != "rp_max_years") setnames(dt, rp_col, "rp_max_years")
#   
#   dt <- dt[is.finite(x) & is.finite(y) & is.finite(intensity) & is.finite(rp_max_years) & !is.na(time_utc)]
#   if (!nrow(dt)) {
#     qsave(data.table(), out_path)
#     return(list(n=0L, n_clusters=0L, status="all_filtered_xytime"))
#   }
#   
#   n_before <- nrow(dt)
#   dt <- .apply_floor(dt, rp_col="rp_max_years", floor_yr=floor_yr)
#   n_after <- nrow(dt)
#   if (!n_after) {
#     qsave(data.table(), out_path)
#     return(list(n=0L, n_clusters=0L, status=sprintf("floor_filtered(%d->%d)", n_before, n_after)))
#   }
#   
#   # Normalized cube (this is the critical change)
#   tmp <- .make_normalized_cube(dt, time_step_hours = time_step_hours)
#   dt  <- tmp$dt
#   cube <- tmp$cube
#   
#   ds <- dbscan::dbscan(cube, eps=EPS, minPts=MU)
#   dt[, `:=`(
#     cluster = as.integer(ds$cluster),
#     eps = EPS, mu = MU,
#     time_step_hours = as.integer(time_step_hours),
#     floor_yr = as.numeric(floor_yr),
#     dataset_tag = tag
#   )]
#   
#   dt <- .add_cluster_id(dt)
#   
#   ncl <- uniqueN(dt$cluster[dt$cluster > 0])
#   qsave(dt, out_path)
#   list(n=nrow(dt), n_clusters=ncl, status="ok")
# }
# 
# # ── RUN ──────────────────────────────────────────────────────────────────────
# all_summaries <- list()
# 
# for (cfg in DATASETS) {
#   
#   tag <- cfg$tag
#   dir_create(cfg$out_dir, recurse=TRUE)
#   
#   .log("============================================================")
#   .log("DATASET: %s | floor=%dyr | dt=%dhr", tag, cfg$floor_yr, cfg$time_step_hours)
#   .log("IN : %s", cfg$in_dir)
#   .log("OUT: %s", cfg$out_dir)
#   
#   files <- .discover_year_files(cfg$in_dir)
#   if (!length(files)) {
#     warning("No files: ", cfg$in_dir)
#     next
#   }
#   
#   years <- vapply(files, .get_year_from_filename, integer(1))
#   ord <- order(years); files <- files[ord]; years <- years[ord]
#   
#   jobs <- data.table(
#     in_path = files,
#     out_path = fs::path(cfg$out_dir, fs::path_file(files)),
#     year = years
#   )
#   
#   ds_summary <- with_progress({
#     p <- progressor(steps=nrow(jobs))
#     parts <- future_lapply(seq_len(nrow(jobs)), function(i) {
#       res <- tryCatch(
#         .cluster_one_year_file(
#           in_path = jobs$in_path[i],
#           out_path = jobs$out_path[i],
#           tag = tag,
#           time_step_hours = cfg$time_step_hours,
#           floor_yr = cfg$floor_yr
#         ),
#         error=function(e) list(n=NA_integer_, n_clusters=NA_integer_, status=paste0("error: ", e$message))
#       )
#       p(sprintf("%s %d", tag, jobs$year[i]))
#       data.table(dataset=tag, file=fs::path_file(jobs$in_path[i]), year=jobs$year[i],
#                  n=res$n, n_clusters=res$n_clusters, status=res$status)
#     })
#     rbindlist(parts, use.names=TRUE, fill=TRUE)
#   })
#   
#   sum_path <- fs::path(cfg$out_dir, sprintf("cluster_summary_%s.csv", tag))
#   fwrite(ds_summary, sum_path)
#   .log("Wrote: %s", sum_path)
#   
#   all_summaries[[tag]] <- ds_summary
# }
# 
# if (length(all_summaries)) {
#   combo <- rbindlist(all_summaries, use.names=TRUE, fill=TRUE)
#   combo_path <- fs::path(OUT_BASE, "cluster_summary_ALL_DATASETS.csv")
#   fwrite(combo, combo_path)
#   .log("Wrote combined summary: %s", combo_path)
# }
# 
# .log("DONE. eps=%.2f | mu=%d | HI_FLOOR=%d | PR_FLOOR=%d", EPS, MU, HI_FLOOR_YR, PR_FLOOR_YR)

# Southeast ----
# =============================================================================
# Cluster exceedances BY YEAR for heatindex + precipitation
#
# UPDATED FOR NEW EXCEEDANCE WORKFLOWS:
#   1) Heat index now uses HI1 exceedances (single-day daily max), not HI3.
#   2) Precipitation uses PR24 bbox return periods / updated exceedance outputs.
#   3) Dataset tags / target resolutions are NO LONGER hard-coded.
#      They are discovered from the STM-driven ref-grid maps:
#        - data/out/03_cov_model/heatindex/dataset_to_refgrid_heatindex.qs
#        - data/out/03_cov_model/precipitation/dataset_to_refgrid_precipitation.qs
#   4) Exceedance input folders are inferred from dataset_tag:
#        heatindex_<res>[_bilinear|_max]      -> data/out/05_exceedances/heatindex/<...>
#        precipitation_<res>[_bilinear|_max]  -> data/out/05_exceedances/precipitation/pr24/<...>
#
# KEY DESIGN:
#   1) Keep ONLY the needed columns for clustering:
#        x, y, time_utc, intensity, rp_max_years (+ optional metadata)
#   2) Build a NORMALIZED integer cube:
#        row_ix = index of y in sorted unique(y)
#        col_ix = index of x in sorted unique(x)
#        t_ix   = integer timestep index at native cadence (e.g., 24 hr)
#      -> Adjacent timesteps are distance 1
#      -> Gaps in time produce large distances
#   3) Apply RP floor BEFORE clustering.
#
# Output:
#   V:/users/hquintal/phd1_cluster_southeast_v2/data/out/06_clusters/hi{HI_FLOOR}_pr{PR_FLOOR}/{dataset_tag}
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(fs)
  library(dbscan)
  library(lubridate)
  library(progressr)
  library(future)
  library(future.apply)
  library(stringr)
})

# ── CONFIG ───────────────────────────────────────────────────────────────────
MU  <- 4L
EPS <- 4.0

BASE <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out"

HI_FLOOR_YR <- 50
PR_FLOOR_YR <- 50

OUT_BASE <- fs::path(BASE, sprintf("06_clusters/hi%d_pr%d", HI_FLOOR_YR, PR_FLOOR_YR))
dir_create(OUT_BASE, recurse = TRUE)

# STM-driven dataset maps from ref-grid builder scripts
HI_MAP_QS <- fs::path(BASE, "03_cov_model/heatindex/dataset_to_refgrid_heatindex.qs")
PR_MAP_QS <- fs::path(BASE, "03_cov_model/precipitation/dataset_to_refgrid_precipitation.qs")

WORKERS <- 20L
plan(multisession, workers = WORKERS)
handlers(global = TRUE)

# ── HELPERS ──────────────────────────────────────────────────────────────────
.log <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
}

.ensure_posix_utc <- function(tt) {
  if (inherits(tt, "POSIXct")) return(as.POSIXct(tt, tz = "UTC"))
  if (inherits(tt, "Date"))   return(as.POSIXct(tt, tz = "UTC"))
  as.POSIXct(tt, tz = "UTC", origin = "1970-01-01")
}

.get_time_utc <- function(dt) {
  if ("time_utc" %in% names(dt)) return(.ensure_posix_utc(dt$time_utc))
  if ("time" %in% names(dt)) return(.ensure_posix_utc(dt$time))
  if ("datetime_utc" %in% names(dt)) return(.ensure_posix_utc(dt$datetime_utc))
  if ("datetime" %in% names(dt)) return(.ensure_posix_utc(dt$datetime))
  if ("day" %in% names(dt)) return(as.POSIXct(as.Date(dt$day), tz = "UTC"))
  if ("date" %in% names(dt)) return(as.POSIXct(as.Date(dt$date), tz = "UTC"))
  stop("No time column found (need time_utc/time/datetime/day/date).")
}

.standardize_xy <- function(dt) {
  if (all(c("x","y") %in% names(dt))) return(invisible(NULL))
  if (all(c("lon","lat") %in% names(dt))) {
    setnames(dt, c("lon","lat"), c("x","y"))
    return(invisible(NULL))
  }
  if (all(c("longitude","latitude") %in% names(dt))) {
    setnames(dt, c("longitude","latitude"), c("x","y"))
    return(invisible(NULL))
  }
  stop("Missing x/y or lon/lat.")
}

.guess_intensity_col <- function(dt, tag) {
  if (grepl("^precipitation", tag, ignore.case = TRUE)) {
    for (cc in c("pr24_intensity_mm", "intensity_mm", "mm", "precip_mm", "prcp_mm",
                 "precipitation", "precip", "pr", "ppt", "value", "intensity")) {
      if (cc %in% names(dt)) return(cc)
    }
  } else {
    for (cc in c("hi1_intensity", "hi3_intensity", "hi_intensity", "heatindex",
                 "heat_index", "HeatIndex", "value", "intensity")) {
      if (cc %in% names(dt)) return(cc)
    }
  }
  
  ignore <- c("x","y","lon","lat","time","time_utc","datetime","datetime_utc",
              "day","date","rp","rp_max_years","rp_max_layer","resolution_deg")
  num_cols <- names(dt)[vapply(dt, is.numeric, TRUE)]
  num_cols <- setdiff(num_cols, ignore)
  if (length(num_cols) == 1) return(num_cols[1])
  
  stop("Cannot infer intensity column for ", tag)
}

.get_rp_col <- function(dt) {
  if ("rp_max_years" %in% names(dt)) return("rp_max_years")
  if ("rp" %in% names(dt)) return("rp")
  stop("Missing return period column (need rp_max_years or rp).")
}

.apply_floor <- function(dt, rp_col, floor_yr) {
  rp <- suppressWarnings(as.numeric(dt[[rp_col]]))
  dt[is.finite(rp) & rp >= floor_yr]
}

.make_normalized_cube <- function(dt, time_step_hours) {
  ux <- sort(unique(dt$x))
  uy <- sort(unique(dt$y))
  
  x_map <- data.table(x = ux, col_ix = seq_along(ux))
  y_map <- data.table(y = uy, row_ix = seq_along(uy))
  
  dt <- x_map[dt, on = "x"]
  dt <- y_map[dt, on = "y"]
  
  t0 <- min(dt$time_utc, na.rm = TRUE)
  dh <- as.numeric(difftime(dt$time_utc, t0, units = "hours"))
  dt[, t_ix := as.integer(round(dh / time_step_hours)) + 1L]
  
  cube <- as.matrix(dt[, .(row_ix, col_ix, t_ix)])
  list(dt = dt, cube = cube)
}

.add_cluster_id <- function(dt) {
  dt[, cluster_id := NA_character_]
  lab <- dt[cluster > 0]
  if (!nrow(lab)) return(dt)
  
  meta <- lab[, {
    ft <- min(time_utc, na.rm = TRUE)
    sub <- .SD[time_utc == ft]
    data.table(
      first_time = ft,
      first_y = min(sub$y, na.rm = TRUE),
      first_x = min(sub$x, na.rm = TRUE)
    )
  }, by = cluster]
  
  meta[, date := format(as.Date(first_time, tz = "UTC"), "%Y-%m-%d")]
  setorder(meta, date, first_time, first_y, first_x, cluster)
  meta[, seq_id := ave(seq_len(.N), date, FUN = seq_along)]
  meta[, cluster_id := sprintf("%s-%03d", date, seq_id)]
  
  dt[meta, on = "cluster", cluster_id := i.cluster_id]
  dt[]
}

.get_year_from_filename <- function(path) {
  bn <- fs::path_file(path)
  yy <- suppressWarnings(as.integer(str_match(bn, "(19|20)\\d{2}")[,1]))
  if (is.na(yy)) stop("Could not parse year from filename: ", bn)
  yy
}

.discover_year_files <- function(in_dir, prefix = NULL) {
  files <- fs::dir_ls(in_dir, type = "file", glob = "*.qs", fail = FALSE)
  files <- sort(as.character(files))
  
  if (!is.null(prefix)) {
    files <- files[grepl(paste0("^", prefix, ".*(19|20)\\d{2}\\.qs$"), fs::path_file(files))]
  } else {
    files <- files[grepl("(19|20)\\d{2}", fs::path_file(files))]
  }
  files
}

# Build dataset configs from STM-driven dataset maps
.build_dataset_configs <- function() {
  stopifnot(file.exists(HI_MAP_QS), file.exists(PR_MAP_QS))
  
  hi_map <- as.data.table(qread(HI_MAP_QS))
  pr_map <- as.data.table(qread(PR_MAP_QS))
  
  req_cols <- c("dataset_tag", "res_deg")
  if (!all(req_cols %in% names(hi_map))) stop("HI_MAP_QS missing required columns: ", paste(req_cols, collapse = ", "))
  if (!all(req_cols %in% names(pr_map))) stop("PR_MAP_QS missing required columns: ", paste(req_cols, collapse = ", "))
  
  hi_tags <- sort(unique(hi_map$dataset_tag))
  pr_tags <- sort(unique(pr_map$dataset_tag))
  
  hi_cfg <- rbindlist(lapply(hi_tags, function(tag) {
    data.table(
      tag = tag,
      in_dir = fs::path(BASE, "05_exceedances/heatindex", sub("^heatindex_", "", tag)),
      out_dir = fs::path(OUT_BASE, tag),
      time_step_hours = 24L,
      floor_yr = HI_FLOOR_YR,
      prefix = "hi1_exceedances_maxrp_"
    )
  }))
  
  pr_cfg <- rbindlist(lapply(pr_tags, function(tag) {
    data.table(
      tag = tag,
      in_dir = fs::path(BASE, "05_exceedances/precipitation/pr24", sub("^precipitation_", "", tag)),
      out_dir = fs::path(OUT_BASE, tag),
      time_step_hours = 24L,
      floor_yr = PR_FLOOR_YR,
      prefix = "pr24_exceedances_maxrp_"
    )
  }))
  
  cfg <- rbindlist(list(pr_cfg, hi_cfg), use.names = TRUE, fill = TRUE)
  cfg[]
}

.cluster_one_year_file <- function(in_path, out_path, tag, time_step_hours, floor_yr) {
  dt0 <- as.data.table(qread(in_path))
  
  if (!nrow(dt0)) {
    qsave(data.table(), out_path)
    return(list(n = 0L, n_clusters = 0L, status = "empty"))
  }
  
  .standardize_xy(dt0)
  dt0[, time_utc := .get_time_utc(dt0)]
  
  rp_col <- .get_rp_col(dt0)
  inten_col <- .guess_intensity_col(dt0, tag)
  
  keep <- intersect(c("x","y","time_utc", inten_col, rp_col, "rp_max_layer", "resolution_deg", "resample_method"), names(dt0))
  dt <- dt0[, ..keep]
  setnames(dt, inten_col, "intensity")
  if (rp_col != "rp_max_years") setnames(dt, rp_col, "rp_max_years")
  
  dt <- dt[is.finite(x) & is.finite(y) & is.finite(intensity) &
             is.finite(rp_max_years) & !is.na(time_utc)]
  
  if (!nrow(dt)) {
    qsave(data.table(), out_path)
    return(list(n = 0L, n_clusters = 0L, status = "all_filtered_xytime"))
  }
  
  n_before <- nrow(dt)
  dt <- .apply_floor(dt, rp_col = "rp_max_years", floor_yr = floor_yr)
  n_after <- nrow(dt)
  
  if (!n_after) {
    qsave(data.table(), out_path)
    return(list(n = 0L, n_clusters = 0L, status = sprintf("floor_filtered(%d->%d)", n_before, n_after)))
  }
  
  tmp <- .make_normalized_cube(dt, time_step_hours = time_step_hours)
  dt   <- tmp$dt
  cube <- tmp$cube
  
  ds <- dbscan::dbscan(cube, eps = EPS, minPts = MU)
  
  dt[, `:=`(
    cluster = as.integer(ds$cluster),
    eps = EPS,
    mu = MU,
    time_step_hours = as.integer(time_step_hours),
    floor_yr = as.numeric(floor_yr),
    dataset_tag = tag
  )]
  
  dt <- .add_cluster_id(dt)
  
  ncl <- uniqueN(dt$cluster[dt$cluster > 0])
  qsave(dt, out_path)
  list(n = nrow(dt), n_clusters = ncl, status = "ok")
}

# ── BUILD DATASET CONFIGS ────────────────────────────────────────────────────
DATASETS <- .build_dataset_configs()

.log("Discovered %d datasets from STM-driven ref-grid maps.", nrow(DATASETS))
print(DATASETS)

# ── RUN ──────────────────────────────────────────────────────────────────────
all_summaries <- list()

for (ii in seq_len(nrow(DATASETS))) {
  cfg <- as.list(DATASETS[ii])
  
  tag <- cfg$tag
  dir_create(cfg$out_dir, recurse = TRUE)
  
  .log("============================================================")
  .log("DATASET: %s | floor=%dyr | dt=%dhr", tag, cfg$floor_yr, cfg$time_step_hours)
  .log("IN : %s", cfg$in_dir)
  .log("OUT: %s", cfg$out_dir)
  
  if (!fs::dir_exists(cfg$in_dir)) {
    warning("Input directory does not exist: ", cfg$in_dir)
    next
  }
  
  files <- .discover_year_files(cfg$in_dir, prefix = cfg$prefix)
  if (!length(files)) {
    warning("No files: ", cfg$in_dir)
    next
  }
  
  years <- vapply(files, .get_year_from_filename, integer(1))
  ord <- order(years)
  files <- files[ord]
  years <- years[ord]
  
  jobs <- data.table(
    in_path = files,
    out_path = fs::path(cfg$out_dir, fs::path_file(files)),
    year = years
  )
  
  ds_summary <- with_progress({
    p <- progressor(steps = nrow(jobs))
    parts <- future_lapply(seq_len(nrow(jobs)), function(i) {
      res <- tryCatch(
        .cluster_one_year_file(
          in_path = jobs$in_path[i],
          out_path = jobs$out_path[i],
          tag = tag,
          time_step_hours = cfg$time_step_hours,
          floor_yr = cfg$floor_yr
        ),
        error = function(e) list(
          n = NA_integer_,
          n_clusters = NA_integer_,
          status = paste0("error: ", e$message)
        )
      )
      
      p(sprintf("%s %d", tag, jobs$year[i]))
      data.table(
        dataset = tag,
        file = fs::path_file(jobs$in_path[i]),
        year = jobs$year[i],
        n = res$n,
        n_clusters = res$n_clusters,
        status = res$status
      )
    })
    rbindlist(parts, use.names = TRUE, fill = TRUE)
  })
  
  sum_path <- fs::path(cfg$out_dir, sprintf("cluster_summary_%s.csv", tag))
  fwrite(ds_summary, sum_path)
  .log("Wrote: %s", sum_path)
  
  all_summaries[[tag]] <- ds_summary
}

if (length(all_summaries)) {
  combo <- rbindlist(all_summaries, use.names = TRUE, fill = TRUE)
  combo_path <- fs::path(OUT_BASE, "cluster_summary_ALL_DATASETS.csv")
  fwrite(combo, combo_path)
  .log("Wrote combined summary: %s", combo_path)
}

.log("DONE. eps=%.2f | mu=%d | HI_FLOOR=%d | PR_FLOOR=%d", EPS, MU, HI_FLOOR_YR, PR_FLOOR_YR)