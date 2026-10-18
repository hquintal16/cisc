# suppressPackageStartupMessages({
#   library(data.table)
#   library(qs)
#   library(fs)
#   library(lubridate)
#   library(terra)
#   
#   # parallel + progress
#   library(future)
#   library(future.apply)
#   library(progressr)
# })
# 
# source("V:/users/hquintal/phd1_cluster_southeast_v2/scripts/14_eval_config.R")
# 
# ref_map <- as.data.table(qread(file.path(REF_DIR, "ref_map.qs")))
# setkey(ref_map, dataset_tag)
# 
# .discover_hi_pr_folders <- function() {
#   all_dirs <- dir_ls(CLUST_BASE, type = "directory", fail = FALSE)
#   hi_pr <- all_dirs[grepl("^hi\\d+_pr\\d+$", path_file(all_dirs))]
#   if (!length(hi_pr)) stop("No hi*_pr* folders under: ", CLUST_BASE)
#   sort(as.character(hi_pr))
# }
# 
# .discover_dataset_dirs <- function(hi_pr_dir) {
#   ds <- dir_ls(hi_pr_dir, type = "directory", fail = FALSE)
#   ds <- ds[path_file(ds) %in% DATASET_TAGS]
#   sort(as.character(ds))
# }
# 
# # Your QS schema is stable: x, y, time_utc, cluster, cluster_id
# .read_clustered_points <- function(qs_path) {
#   dt <- as.data.table(qread(qs_path))
#   if (!nrow(dt)) return(dt[0])
#   
#   need <- c("x","y","time_utc","cluster","cluster_id")
#   miss <- setdiff(need, names(dt))
#   if (length(miss)) stop("Missing columns in ", qs_path, ": ", paste(miss, collapse=", "))
#   
#   dt <- dt[is.finite(x) & is.finite(y)]
#   dt <- dt[is.finite(cluster) & cluster > 0]
#   dt <- dt[!is.na(cluster_id) & nzchar(as.character(cluster_id))]
#   if (!nrow(dt)) return(dt[0])
#   
#   # date from time_utc
#   dt[, date := as.Date(time_utc, tz = "UTC")]
#   dt <- dt[!is.na(date)]
#   
#   dt[, .(x, y, date)]
# }
# 
# .points_to_modeled_daily_one_file <- function(qs_path, r_ref, lookup_dt) {
#   pts <- .read_clustered_points(qs_path)
#   if (!nrow(pts)) return(NULL)
#   
#   # cell id on canonical reference raster
#   pts[, cell_id := terra::cellFromXY(r_ref, cbind(x, y))]
#   pts <- pts[!is.na(cell_id)]
#   if (!nrow(pts)) return(NULL)
#   
#   # join lookup (cell_id -> CWA)
#   pts <- lookup_dt[pts, on = "cell_id", nomatch = 0L]
#   pts <- pts[!is.na(CWA) & nzchar(CWA)]
#   if (!nrow(pts)) return(NULL)
#   
#   unique(pts[, .(CWA, date)][, modeled := 1L], by = c("CWA","date"))
# }
# 
# .build_modeled_one <- function(ds_dir, scenario, dataset_tag, out_qs) {
#   if (file_exists(out_qs) && !OVERWRITE_MODEL) {
#     .log("[MODEL] exists: %s", out_qs)
#     return(out_qs)
#   }
#   
#   # canonical artifacts for this dataset_tag
#   ref_row <- ref_map[J(dataset_tag)]
#   if (!nrow(ref_row)) stop("No ref_map entry for dataset_tag=", dataset_tag)
#   
#   ref_tif   <- ref_row$ref_tif[1]
#   lookup_qs <- file.path(LOOKUP_DIR, sprintf("lookup_xy_to_cwa_%s.qs", dataset_tag))
#   if (!file_exists(lookup_qs)) stop("Missing lookup: ", lookup_qs)
#   
#   r_ref <- terra::rast(ref_tif)
#   
#   lookup_dt <- as.data.table(qread(lookup_qs))
#   stopifnot(all(c("cell_id","CWA") %in% names(lookup_dt)))
#   lookup_dt[, cell_id := as.integer(cell_id)]
#   lookup_dt[, CWA := toupper(trimws(as.character(CWA)))]
#   setkey(lookup_dt, cell_id)
#   
#   files <- dir_ls(ds_dir, type = "file", glob = "*.qs", fail = FALSE)
#   files <- sort(as.character(files))
#   if (!length(files)) {
#     qsave(data.table(), out_qs)
#     .log("[MODEL] no input qs files -> wrote empty: %s", out_qs)
#     return(out_qs)
#   }
#   
#   # NOTE: keep this sequential to avoid nested parallelism (we parallelize at scenario x dataset level)
#   parts <- lapply(files, function(f) {
#     tryCatch(
#       .points_to_modeled_daily_one_file(f, r_ref, lookup_dt),
#       error = function(e) {
#         .log("[MODEL][ERR] %s | %s", path_file(f), e$message)
#         NULL
#       }
#     )
#   })
#   
#   out <- rbindlist(Filter(Negate(is.null), parts), use.names = TRUE, fill = TRUE)
#   if (!nrow(out)) {
#     qsave(data.table(), out_qs)
#     .log("[MODEL] wrote empty: %s", out_qs)
#     return(out_qs)
#   }
#   
#   out <- out[, .(modeled = 1L), by = .(CWA, date)]
#   setorder(out, date, CWA)
#   
#   dir_create(path_dir(out_qs), recurse = TRUE)
#   tmp <- paste0(out_qs, ".tmp")
#   if (file_exists(tmp)) file_delete(tmp)
#   qsave(out, tmp); file_move(tmp, out_qs)
#   
#   .log("[MODEL] wrote: %s (rows=%d)", out_qs, nrow(out))
#   out_qs
# }
# 
# # ----------------------------
# # MAIN (PARALLEL over scenario x dataset)
# # ----------------------------
# scenarios <- .discover_hi_pr_folders()
# .log("Scenarios found: %d", length(scenarios))
# 
# # Build task table
# tasks <- rbindlist(lapply(scenarios, function(sc_dir) {
#   scenario <- path_file(sc_dir)
#   ds_dirs <- .discover_dataset_dirs(sc_dir)
#   if (!length(ds_dirs)) return(NULL)
#   
#   rbindlist(lapply(ds_dirs, function(ds_dir) {
#     dataset_tag <- path_file(ds_dir)
#     out_qs <- file.path(
#       MODEL_WFO_DAILY_DIR,
#       scenario,
#       sprintf("modeled_wfo_daily_%s.qs", dataset_tag)
#     )
#     data.table(
#       sc_dir = sc_dir,
#       scenario = scenario,
#       ds_dir = ds_dir,
#       dataset_tag = dataset_tag,
#       out_qs = out_qs
#     )
#   }))
# }), use.names = TRUE, fill = TRUE)
# 
# if (is.null(tasks) || !nrow(tasks)) {
#   .log("[MODEL] no tasks discovered. Nothing to do.")
#   quit(save = "no")
# }
# 
# # Pre-skip outputs that already exist (resume-safe)
# if (!OVERWRITE_MODEL) {
#   tasks[, exists := file_exists(out_qs)]
#   .log("[MODEL] tasks total=%d | already_done=%d | todo=%d",
#        nrow(tasks), sum(tasks$exists), sum(!tasks$exists))
#   
#   # FIX: data.table requires explicit column reference (avoid DT[!exists] scoping trap)
#   tasks <- tasks[exists == FALSE]
#   tasks[, exists := NULL]
#   
#   if (!nrow(tasks)) {
#     .log("[MODEL] all outputs already exist; nothing to do.")
#     quit(save = "no")
#   }
# } else {
#   .log("[MODEL] OVERWRITE_MODEL=TRUE; will rebuild all %d task(s).", nrow(tasks))
# }
# 
# # Workers (override by options if you want)
# WORKERS_MODEL <- getOption("cisc.workers_model", 8L)
# 
# plan(multisession, workers = WORKERS_MODEL)
# handlers(global = TRUE)
# progressr::handlers("txtprogressbar")
# 
# dir_create(MODEL_WFO_DAILY_DIR, recurse = TRUE)
# 
# with_progress({
#   p <- progressor(along = 1:nrow(tasks))
#   
#   out_files <- future_lapply(1:nrow(tasks), function(i) {
#     row <- tasks[i]
#     
#     p(sprintf("%s | %s", row$scenario, row$dataset_tag))
#     
#     dir_create(path_dir(row$out_qs), recurse = TRUE)
#     
#     .build_modeled_one(
#       ds_dir      = row$ds_dir,
#       scenario    = row$scenario,
#       dataset_tag = row$dataset_tag,
#       out_qs      = row$out_qs
#     )
#   })
#   
#   .log("[MODEL] finished %d task(s).", length(out_files))
# })
# 
# plan(sequential)

# 20260224 ----
suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(fs)
  library(lubridate)
  library(terra)
  
  library(future)
  library(future.apply)
  library(progressr)
})

source("V:/users/hquintal/phd1_cluster_southeast_v2/scripts/14_eval_config.R")

dataset_map <- as.data.table(qread(file.path(REF_DIR, "dataset_to_refgrid.qs")))
setkey(dataset_map, dataset_tag)

.discover_hi_pr_folders <- function() {
  all_dirs <- dir_ls(CLUST_BASE, type = "directory", fail = FALSE)
  hi_pr <- all_dirs[grepl("^hi\\d+_pr\\d+$", path_file(all_dirs))]
  if (!length(hi_pr)) stop("No hi*_pr* folders under: ", CLUST_BASE)
  sort(as.character(hi_pr))
}

.discover_dataset_dirs <- function(hi_pr_dir) {
  ds <- dir_ls(hi_pr_dir, type = "directory", fail = FALSE)
  ds <- ds[path_file(ds) %in% DATASET_TAGS]
  sort(as.character(ds))
}

.read_clustered_points <- function(qs_path) {
  dt <- as.data.table(qread(qs_path))
  if (!nrow(dt)) return(dt[0])
  
  need <- c("x","y","time_utc","cluster","cluster_id")
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("Missing columns in ", qs_path, ": ", paste(miss, collapse=", "))
  
  dt <- dt[is.finite(x) & is.finite(y)]
  dt <- dt[is.finite(cluster) & cluster > 0]
  dt <- dt[!is.na(cluster_id) & nzchar(as.character(cluster_id))]
  if (!nrow(dt)) return(dt[0])
  
  dt[, date := as.Date(time_utc, tz = "UTC")]
  dt <- dt[!is.na(date)]
  dt[, .(x, y, date)]
}

# Cache keyed by ref_grid_tag (resolution) for max reuse
.worker_cache <- new.env(parent = emptyenv())

.get_assets_for_refgrid <- function(ref_grid_tag, ref_tif) {
  key <- paste0("assets__", ref_grid_tag)
  if (exists(key, envir=.worker_cache, inherits=FALSE)) {
    return(get(key, envir=.worker_cache))
  }
  
  lookup_qs <- file.path(LOOKUP_DIR_COUNTY, sprintf("lookup_xy_to_county_%s.qs", ref_grid_tag))
  if (!file_exists(lookup_qs)) stop("Missing lookup: ", lookup_qs)
  
  r_ref <- terra::rast(ref_tif)
  
  lookup_dt <- as.data.table(qread(lookup_qs))
  stopifnot(all(c("cell_id","GEOID") %in% names(lookup_dt)))
  lookup_dt[, cell_id := as.integer(cell_id)]
  lookup_dt[, GEOID := as.character(GEOID)]
  lookup_dt <- lookup_dt[!is.na(GEOID) & nzchar(GEOID)]
  setkey(lookup_dt, cell_id)
  
  assets <- list(r_ref = r_ref, lookup_dt = lookup_dt)
  assign(key, assets, envir=.worker_cache)
  assets
}

.points_to_modeled_daily_one_file <- function(qs_path, assets) {
  pts <- .read_clustered_points(qs_path)
  if (!nrow(pts)) return(NULL)
  
  pts[, cell_id := terra::cellFromXY(assets$r_ref, cbind(x, y))]
  pts <- pts[!is.na(cell_id)]
  if (!nrow(pts)) return(NULL)
  
  pts <- assets$lookup_dt[pts, on="cell_id", nomatch=0L]
  pts <- pts[!is.na(GEOID) & nzchar(GEOID)]
  if (!nrow(pts)) return(NULL)
  
  unique(pts[, .(GEOID, date)][, modeled := 1L], by=c("GEOID","date"))
}

.write_atomic_qs <- function(obj, out_path) {
  dir_create(path_dir(out_path), recurse = TRUE)
  tmp <- paste0(out_path, ".tmp")
  if (file_exists(tmp)) file_delete(tmp)
  qsave(obj, tmp)
  file_move(tmp, out_path)
  invisible(out_path)
}

.build_modeled_one <- function(ds_dir, scenario, dataset_tag, out_qs) {
  if (file_exists(out_qs) && !OVERWRITE_MODEL_COUNTY) {
    .log("[MODEL COUNTY] exists: %s", out_qs)
    return(out_qs)
  }
  
  row <- dataset_map[J(dataset_tag)]
  if (!nrow(row)) stop("No dataset_map row for dataset_tag=", dataset_tag)
  
  ref_grid_tag <- row$ref_grid_tag[1]
  ref_tif      <- row$ref_tif[1]
  
  assets <- .get_assets_for_refgrid(ref_grid_tag, ref_tif)
  
  files <- dir_ls(ds_dir, type="file", glob="*.qs", fail=FALSE)
  files <- sort(as.character(files))
  if (!length(files)) {
    .write_atomic_qs(data.table(), out_qs)
    .log("[MODEL COUNTY] no input files -> wrote empty: %s", out_qs)
    return(out_qs)
  }
  
  parts <- lapply(files, function(f) {
    tryCatch(.points_to_modeled_daily_one_file(f, assets),
             error=function(e) { .log("[MODEL COUNTY][ERR] %s | %s", path_file(f), e$message); NULL })
  })
  
  out <- rbindlist(Filter(Negate(is.null), parts), use.names=TRUE, fill=TRUE)
  if (!nrow(out)) {
    .write_atomic_qs(data.table(), out_qs)
    .log("[MODEL COUNTY] wrote empty: %s", out_qs)
    return(out_qs)
  }
  
  out <- out[, .(modeled=1L), by=.(GEOID, date)]
  setorder(out, date, GEOID)
  
  .write_atomic_qs(out, out_qs)
  .log("[MODEL COUNTY] wrote: %s (rows=%d)", out_qs, nrow(out))
  out_qs
}

# ----------------------------
# MAIN (parallel over scenario × dataset)
# ----------------------------
scenarios <- .discover_hi_pr_folders()
.log("Scenarios found: %d", length(scenarios))

tasks <- rbindlist(lapply(scenarios, function(sc_dir) {
  scenario <- path_file(sc_dir)
  ds_dirs <- .discover_dataset_dirs(sc_dir)
  if (!length(ds_dirs)) return(NULL)
  
  rbindlist(lapply(ds_dirs, function(ds_dir) {
    dataset_tag <- path_file(ds_dir)
    out_qs <- file.path(MODEL_COUNTY_DAILY_DIR, scenario,
                        sprintf("modeled_county_daily_%s.qs", dataset_tag))
    data.table(sc_dir=sc_dir, scenario=scenario, ds_dir=ds_dir, dataset_tag=dataset_tag, out_qs=out_qs)
  }))
}), use.names=TRUE, fill=TRUE)

if (is.null(tasks) || !nrow(tasks)) {
  .log("[MODEL COUNTY] no tasks discovered.")
  quit(save="no")
}

if (!OVERWRITE_MODEL_COUNTY) {
  tasks[, exists := file_exists(out_qs)]
  .log("[MODEL COUNTY] tasks total=%d | already_done=%d | todo=%d",
       nrow(tasks), sum(tasks$exists), sum(!tasks$exists))
  tasks <- tasks[exists == FALSE]
  tasks[, exists := NULL]
  if (!nrow(tasks)) {
    .log("[MODEL COUNTY] nothing to do.")
    quit(save="no")
  }
}

plan(multisession, workers = min(WORKERS_MODEL, nrow(tasks)))
handlers(global = TRUE)
handlers("txtprogressbar")

with_progress({
  p <- progressor(along = 1:nrow(tasks))
  future_lapply(seq_len(nrow(tasks)), function(i) {
    row <- tasks[i]
    p(sprintf("%s | %s", row$scenario, row$dataset_tag))
    dir_create(path_dir(row$out_qs), recurse=TRUE)
    .build_modeled_one(row$ds_dir, row$scenario, row$dataset_tag, row$out_qs)
  }, future.seed=TRUE,
  future.packages=c("data.table","qs","terra","fs","lubridate"))
})

plan(sequential)
.log("✅ modeled county cache build complete")