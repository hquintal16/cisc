# scripts/08_write_event_tifs_from_clusters_hi_pr_allfolders_v1.R
# =============================================================================
# Write per-event GeoTIFFs (intensity cube + max RP) from clustered exceedance .qs
#
# Discovers ALL clusters across:
#   V:/users/hquintal/phd1_cluster_southeast_v2/data/out/06_clusters/hi{X}_pr{Y}/
# where X,Y are floors (e.g., 1,2,5,10).
#
# For each dataset folder inside hiX_prY (6 sets):
#   - precipitation_0.25
#   - precipitation_0.17_bilinear
#   - precipitation_0.17_max
#   - heatindex_0.25
#   - heatindex_0.44_bilinear
#   - heatindex_0.44_max
#
# Outputs (cropped to event extent) to:
#   data/out/08_rasters/{heatindex|precipitation}/rp{FLOOR}/res{RES_TAG}/{intensity|maxrp}/<cluster_id>.tif
#
# Notes:
# - Uses MASK_TIF as the authoritative CONUS template (EPSG:4326).
# - Builds target grids by resampling MASK_TIF to the dataset resolution:
#     0.25° native, 0.441501° (HI resampled), 0.167403858° (PR resampled)
# - Temporal resolution is preserved; layers are unique event times.
# - RP is taken from rp_max_years if present, else rp.
# - Intensity column:
#     precipitation: pr24_intensity_mm (preferred)
#     heatindex:     hi3_intensity (preferred)
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(data.table)
  library(qs)
  library(fs)
  library(stringr)
  library(lubridate)
  library(future)
  library(future.apply)
  library(progressr)
})

# ── CONFIG ────────────────────────────────────────────────────────────────────
BASE         <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out"
CLUST_BASE   <- fs::path(BASE, "06_clusters")
OUT_BASE     <- fs::path(BASE, "08_rasters")

MASK_TIF     <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/temperature/_mask_target.tif"

# target spatial resolutions (degrees)
RES_NATIVE_025 <- c(0.25, 0.25)
RES_HI_044     <- c(0.441501, 0.441501)               # <-- updated as requested
RES_PR_017     <- c(0.167403858, 0.167403858)         # keep exact

# cropped output buffer (in cells of the target grid)
BUFFER_CELLS <- 1L

# parallel
WORKERS   <- 30L
OVERWRITE <- FALSE

# DB: output datatypes
DTYPE_INT <- "FLT4S"
DTYPE_RP  <- "FLT4S"

terraOptions(progress = 0, todisk = TRUE)
Sys.setenv(GDAL_CACHEMAX = "1024")

plan(multisession, workers = WORKERS)
handlers(global = TRUE)

.log <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

# ── HELPERS: template grids from MASK_TIF ─────────────────────────────────────
.make_target_grid <- function(mask_r, target_res) {
  stopifnot(length(target_res) == 2)
  r0 <- terra::rast(ext(mask_r), resolution = target_res, crs = crs(mask_r))
  terra::resample(mask_r, r0, method = "near")
}

.get_ref_template <- local({
  cache <- new.env(parent = emptyenv())
  function(dataset_tag) {
    if (exists(dataset_tag, envir = cache, inherits = FALSE)) {
      return(get(dataset_tag, envir = cache, inherits = FALSE))
    }
    
    mask_r <- terra::rast(MASK_TIF)
    if (terra::nlyr(mask_r) > 1) mask_r <- mask_r[[1]]
    
    target_res <- switch(
      dataset_tag,
      "heatindex_0.25"               = RES_NATIVE_025,
      "precipitation_0.25"           = RES_NATIVE_025,
      "heatindex_0.44_bilinear"      = RES_HI_044,
      "heatindex_0.44_max"           = RES_HI_044,
      "precipitation_0.17_bilinear"  = RES_PR_017,
      "precipitation_0.17_max"       = RES_PR_017,
      stop("Unknown dataset_tag for template resolution: ", dataset_tag)
    )
    
    ref <- .make_target_grid(mask_r, target_res)
    assign(dataset_tag, ref, envir = cache)
    ref
  }
})

# ── HELPERS: robust column parsing ────────────────────────────────────────────
.standardize_xy <- function(dt) {
  nms <- names(dt)
  if (all(c("x","y") %in% nms)) return(invisible(NULL))
  if (all(c("lon","lat") %in% nms)) { setnames(dt, c("lon","lat"), c("x","y")); return(invisible(NULL)) }
  if (all(c("longitude","latitude") %in% nms)) { setnames(dt, c("longitude","latitude"), c("x","y")); return(invisible(NULL)) }
  stop("Missing x/y (or lon/lat). Have: ", paste(head(nms, 40), collapse = ", "))
}

.get_time <- function(dt) {
  if ("time_utc" %in% names(dt)) return(dt$time_utc)
  if ("time" %in% names(dt)) return(dt$time)
  if ("datetime_utc" %in% names(dt)) return(dt$datetime_utc)
  if ("datetime" %in% names(dt)) return(dt$datetime)
  if ("day" %in% names(dt)) return(dt$day)
  if ("date" %in% names(dt)) return(dt$date)
  stop("No time column found (expected time_utc/time/datetime/day/date).")
}

.as_posix_utc <- function(x) {
  if (inherits(x, "POSIXct")) return(as.POSIXct(x, tz = "UTC"))
  if (inherits(x, "Date"))   return(as.POSIXct(x, tz = "UTC"))
  # try parse
  xx <- suppressWarnings(ymd_hms(x, tz = "UTC"))
  if (all(is.na(xx))) xx <- suppressWarnings(ymd(x, tz = "UTC"))
  if (all(is.na(xx))) xx <- suppressWarnings(as.POSIXct(x, tz = "UTC"))
  xx
}

.get_rp <- function(dt) {
  if ("rp_max_years" %in% names(dt)) return(suppressWarnings(as.numeric(dt$rp_max_years)))
  if ("rp" %in% names(dt)) return(suppressWarnings(as.numeric(dt$rp)))
  return(rep(NA_real_, nrow(dt)))
}

.get_intensity_col <- function(dt, dataset_tag) {
  if (grepl("^precipitation", dataset_tag, ignore.case = TRUE)) {
    # user-specified
    if ("pr24_intensity_mm" %in% names(dt)) return("pr24_intensity_mm")
    # fallbacks
    for (cc in c("intensity_mm","mm","precip_mm","pr_mm","pr","ppt","value","intensity")) {
      if (cc %in% names(dt)) return(cc)
    }
  } else {
    if ("hi3_intensity" %in% names(dt)) return("hi3_intensity")
    for (cc in c("hi_intensity","heatindex","HeatIndex","value","intensity")) {
      if (cc %in% names(dt)) return(cc)
    }
  }
  stop("Could not find intensity column for dataset_tag=", dataset_tag,
       ". Available cols: ", paste(names(dt), collapse = ", "))
}

# ── HELPERS: cropping extents in cell units ───────────────────────────────────
.extent_from_cells_buffered <- function(ref_r, cells, buffer_cells = 1L) {
  if (!length(cells)) return(terra::ext(ref_r))
  xy  <- terra::xyFromCell(ref_r, cells)
  xrg <- range(xy[,1], na.rm = TRUE); yrg <- range(xy[,2], na.rm = TRUE)
  rs  <- terra::res(ref_r)
  ex  <- terra::ext(xrg[1]-buffer_cells*rs[1], xrg[2]+buffer_cells*rs[1],
                    yrg[1]-buffer_cells*rs[2], yrg[2]+buffer_cells*rs[2])
  terra::intersect(ex, terra::ext(ref_r))
}

.write_gtiff <- function(r, out_tif, overwrite = TRUE, datatype = "FLT4S") {
  terra::writeRaster(
    r, out_tif,
    overwrite = overwrite,
    datatype = datatype,
    gdal = c("COMPRESS=LZW", "PREDICTOR=2", "TILED=YES")
  )
  out_tif
}

# ── WRITERS ───────────────────────────────────────────────────────────────────
.write_intensity_cube <- function(dt_cl, dataset_tag, out_tif, overwrite = TRUE) {
  # dt_cl must have: x,y,time,intensity
  tt <- sort(unique(dt_cl$time))
  tt <- tt[!is.na(tt)]
  if (!length(tt)) stop("No valid times for intensity cube.")
  
  ref <- .get_ref_template(dataset_tag)
  
  # crop template to event footprint
  base_cells <- terra::cellFromXY(ref, as.matrix(dt_cl[, .(x,y)]))
  base_cells <- base_cells[!is.na(base_cells)]
  if (!length(base_cells)) stop("No valid cells from XY mapping (intensity).")
  
  ex <- .extent_from_cells_buffered(ref, base_cells, buffer_cells = BUFFER_CELLS)
  tmpl <- terra::crop(ref, ex, snap = "out")
  
  cube <- terra::rast(tmpl, nlyrs = length(tt))
  names(cube) <- format(tt, "t%Y%m%d_%H%M%S")
  
  dt <- copy(dt_cl)
  dt[, lyr := match(time, tt)]
  dt[, cell := terra::cellFromXY(tmpl, as.matrix(dt[, .(x,y)]))]
  
  dt <- dt[!is.na(cell) & !is.na(lyr) & is.finite(intensity),
           .(intensity = max(intensity, na.rm = TRUE)),
           by = .(cell, lyr)]
  if (!nrow(dt)) stop("No data after (cell,layer) aggregation for intensity.")
  
  ncell <- terra::ncell(tmpl)
  V <- matrix(NA_real_, nrow = ncell, ncol = length(tt))
  V[(dt$lyr - 1L) * ncell + dt$cell] <- dt$intensity
  terra::values(cube) <- V
  
  .write_gtiff(cube, out_tif, overwrite = overwrite, datatype = DTYPE_INT)
}

.write_maxrp <- function(dt_cl, dataset_tag, out_tif, overwrite = TRUE) {
  # dt_cl must have: x,y,rp
  ref <- .get_ref_template(dataset_tag)
  
  base_cells <- terra::cellFromXY(ref, as.matrix(dt_cl[, .(x,y)]))
  base_cells <- base_cells[!is.na(base_cells)]
  if (!length(base_cells)) stop("No valid cells from XY mapping (RP).")
  
  ex <- .extent_from_cells_buffered(ref, base_cells, buffer_cells = BUFFER_CELLS)
  tmpl <- terra::crop(ref, ex, snap = "out")
  
  dt <- copy(dt_cl)
  dt[, cell := terra::cellFromXY(tmpl, as.matrix(dt[, .(x,y)]))]
  dt <- dt[!is.na(cell) & is.finite(rp),
           .(rp = max(rp, na.rm = TRUE)),
           by = cell]
  if (!nrow(dt)) stop("No data after cell aggregation (RP).")
  
  r <- tmpl
  vals <- rep(NA_real_, terra::ncell(tmpl))
  vals[dt$cell] <- dt$rp
  terra::values(r) <- vals
  names(r) <- "max_rp"
  
  .write_gtiff(r, out_tif, overwrite = overwrite, datatype = DTYPE_RP)
}

# ── DISCOVERY ─────────────────────────────────────────────────────────────────
# Find hiX_prY folders (filter by basename, not full path)
all_dirs <- fs::dir_ls(CLUST_BASE, type = "directory", fail = FALSE)
all_dirs <- sort(as.character(all_dirs))

hi_pr_folders <- all_dirs[grepl("^hi\\d+_pr\\d+$", fs::path_file(all_dirs))]
hi_pr_folders <- sort(hi_pr_folders)

if (!length(hi_pr_folders)) {
  stop(
    "No hi*_pr* folders found under: ", CLUST_BASE, "\n",
    "Found these dirs:\n- ", paste(fs::path_file(all_dirs), collapse = "\n- ")
  )
}

# For each hiX_prY folder, discover dataset subfolders
.discover_dataset_folders <- function(hi_pr_dir) {
  ds <- fs::dir_ls(hi_pr_dir, type = "directory", fail = FALSE)
  ds <- sort(as.character(ds))
  # keep only known patterns
  ds[fs::path_file(ds) %in% c(
    "precipitation_0.25","precipitation_0.17_bilinear","precipitation_0.17_max",
    "heatindex_0.25","heatindex_0.44_bilinear","heatindex_0.44_max"
  )]
}

# parse floors from folder name
.parse_floors <- function(hi_pr_dir) {
  bn <- fs::path_file(hi_pr_dir)
  m <- stringr::str_match(bn, "^hi(\\d+)_pr(\\d+)$")
  list(hi = as.integer(m[,2]), pr = as.integer(m[,3]))
}

# resolution tag for output folders
.res_tag <- function(dataset_tag) {
  if (dataset_tag %in% c("heatindex_0.25","precipitation_0.25")) return("res0.25")
  if (dataset_tag %in% c("heatindex_0.44_bilinear","heatindex_0.44_max")) return("res0.44")
  if (dataset_tag %in% c("precipitation_0.17_bilinear","precipitation_0.17_max")) return("res0.17")
  "resUNK"
}

.var_group <- function(dataset_tag) {
  if (grepl("^heatindex", dataset_tag)) return("heatindex")
  "precipitation"
}

.method_tag <- function(dataset_tag) {
  if (grepl("_bilinear$", dataset_tag)) return("bilinear")
  if (grepl("_max$", dataset_tag)) return("max")
  "native"
}

# output directory builder
.out_dirs_for <- function(dataset_tag, floor_yr) {
  varg <- .var_group(dataset_tag)
  rest <- .res_tag(dataset_tag)
  meth <- .method_tag(dataset_tag)
  
  out_int <- fs::path(OUT_BASE, varg, sprintf("rp%d", floor_yr), sprintf("%s_%s", rest, meth), "intensity")
  out_rp  <- fs::path(OUT_BASE, varg, sprintf("rp%d", floor_yr), sprintf("%s_%s", rest, meth), "maxrp")
  fs::dir_create(c(out_int, out_rp), recurse = TRUE)
  list(intensity = out_int, maxrp = out_rp)
}

# build task table: each cluster_id within each qs year file
.tasks <- list()

for (hi_pr_dir in hi_pr_folders) {
  floors <- .parse_floors(hi_pr_dir)
  
  ds_dirs <- .discover_dataset_folders(hi_pr_dir)
  if (!length(ds_dirs)) next
  
  for (ds_dir in ds_dirs) {
    dataset_tag <- fs::path_file(ds_dir)
    
    # choose floor to label outputs by variable type
    floor_yr <- if (grepl("^heatindex", dataset_tag)) floors$hi else floors$pr
    
    qs_files <- fs::dir_ls(ds_dir, type = "file", glob = "*.qs", fail = FALSE)
    qs_files <- qs_files[!grepl("summary|characteristics", fs::path_file(qs_files), ignore.case = TRUE)]
    qs_files <- sort(as.character(qs_files))
    if (!length(qs_files)) next
    
    for (f in qs_files) {
      dt <- try(as.data.table(qread(f)), silent = TRUE)
      if (inherits(dt, "try-error") || !nrow(dt)) next
      if (!("cluster_id" %in% names(dt))) next
      
      cids <- unique(dt$cluster_id)
      cids <- cids[!is.na(cids) & nzchar(as.character(cids))]
      if (!length(cids)) next
      
      key <- paste(hi_pr_dir, dataset_tag, f, sep = "||")
      .tasks[[key]] <- data.table(
        hi_pr_dir   = hi_pr_dir,
        dataset_tag = dataset_tag,
        floor_yr    = floor_yr,
        qs_path     = f,
        cluster_id  = as.character(cids)
      )
    }
  }
}

if (!length(.tasks)) stop("No cluster tasks discovered under: ", CLUST_BASE)
task_dt <- rbindlist(.tasks, use.names = TRUE, fill = TRUE)

# Resume-safe skipping if OVERWRITE=FALSE
if (!OVERWRITE) {
  task_dt[, out_rp_tif := {
    out_dirs <- .out_dirs_for(dataset_tag, floor_yr)
    fs::path(out_dirs$maxrp, sprintf("%s.tif", cluster_id))
  }, by = .(dataset_tag, floor_yr)]
  task_dt <- task_dt[!fs::file_exists(out_rp_tif)]
  task_dt[, out_rp_tif := NULL]
}

.log("Queued clusters: %d | workers=%d | overwrite=%s", nrow(task_dt), WORKERS, OVERWRITE)

# ── PROCESS ONE TASK ──────────────────────────────────────────────────────────
.process_one <- function(qs_path, target_cluster_id, dataset_tag, floor_yr) {
  
  dt <- as.data.table(qread(qs_path))
  if (!nrow(dt)) return(NULL)
  
  .standardize_xy(dt)
  
  # subset to target cluster_id (robust + unambiguous)
  dt <- dt[!is.na(cluster_id) & as.character(cluster_id) == as.character(target_cluster_id)]
  if (!nrow(dt)) return(NULL)
  
  # time + rp + intensity
  time_raw <- .get_time(dt)
  tt <- .as_posix_utc(time_raw)
  
  rp <- .get_rp(dt)
  
  inten_col <- .get_intensity_col(dt, dataset_tag)
  intensity <- suppressWarnings(as.numeric(dt[[inten_col]]))
  
  # enforce valid rows
  keep <- is.finite(dt$x) & is.finite(dt$y) & !is.na(tt) & is.finite(intensity)
  
  dt2 <- data.table(
    x = dt$x[keep],
    y = dt$y[keep],
    time = tt[keep],
    intensity = intensity[keep],
    rp = rp[keep]
  )
  if (!nrow(dt2)) return(NULL)
  
  # output dirs
  out_dirs <- .out_dirs_for(dataset_tag, floor_yr)
  out_int_tif <- fs::path(out_dirs$intensity, sprintf("%s.tif", target_cluster_id))
  out_rp_tif  <- fs::path(out_dirs$maxrp,    sprintf("%s.tif", target_cluster_id))
  
  # write
  if (OVERWRITE || !fs::file_exists(out_int_tif)) {
    .write_intensity_cube(dt2, dataset_tag, out_int_tif, overwrite = OVERWRITE)
  }
  if (OVERWRITE || !fs::file_exists(out_rp_tif)) {
    .write_maxrp(dt2, dataset_tag, out_rp_tif, overwrite = OVERWRITE)
  }
  
  target_cluster_id
}

# ── PARALLEL RUN WITH PROGRESS ────────────────────────────────────────────────
with_progress({
  p <- progressor(steps = nrow(task_dt))
  res <- future_lapply(seq_len(nrow(task_dt)), function(i) {
    row <- task_dt[i]
    out <- tryCatch(
      .process_one(
        qs_path            = row$qs_path,
        target_cluster_id  = row$cluster_id,
        dataset_tag        = row$dataset_tag,
        floor_yr           = row$floor_yr
      ),
      error = function(e) {
        .log("ERROR | %s | %s | %s", row$dataset_tag, row$cluster_id, e$message)
        NULL
      }
    )
    p(sprintf("%s | %s", row$dataset_tag, row$cluster_id))
    out
  })
})

.log("✅ DONE. Wrote GeoTIFFs to: %s", OUT_BASE)