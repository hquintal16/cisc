# Extract monthly .nc to daily .qs ----
# ERA5 monthly NetCDF cubes (.nc) -> DAILY .qs tables (one file per UTC day)
#
# CHANGE:
# - Mask is now built by rasterizing an ecoregion polygon shapefile onto the ERA5 template grid.
# - Writes OUT_DIR/_mask_target.tif once, then workers load+apply it.

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(data.table)
  library(qs)
  library(fs)
  library(glue)
  library(progressr)
  library(future)
  library(future.apply)
})

# INPUTS ----
IN_DIR  <- "V:/users/hquintal/phd1_cluster_southeast/data/input/precipitation"
OUT_DIR <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation"
FILE_GLOB <- "*.nc"

TARGET_UNITS <- "mm"     # precip: m -> mm
KEEP_ZERO <- TRUE
DROP_NA   <- TRUE

# NEW: polygon mask source
ECOREGION_SHP <- "V:/users/hquintal/phd1_cluster_southeast/data/input/regional_aggregation/ecoregion/Aggr_Ecoregions_2015/Aggr_Ecoregions_2015.shp"
# Optional: filter to a subset (leave NULL to keep all 9)
ECOREGION_KEEP <- NULL  # e.g., c("CPL","SAP") or NULL

MASK_DILATE_CELLS <- 1L
N_WORKERS <- 36L

terra::terraOptions(progress = 0)

# Helpers ----
convert_units <- function(values, from_unit, target_unit) {
  
  if (is.na(from_unit) || !nzchar(from_unit) ||
      is.na(target_unit) || !nzchar(target_unit)) {
    return(values)
  }
  
  from_unit   <- trimws(from_unit)
  target_unit <- trimws(target_unit)
  
  if (identical(from_unit, target_unit)) return(values)
  
  # Precipitation
  if (from_unit == "m" && target_unit == "mm") {
    return(values * 1000)
  }
  
  warning(sprintf(
    "No conversion rule from '%s' to '%s'. Leaving unchanged.",
    from_unit, target_unit
  ))
  
  values
}

# UTC calendar date label: "YYYY-MM-DD"
date_key_utc <- function(x_utc) {
  x_utc  <- as.POSIXct(x_utc, tz = "UTC")
  x_date <- as.Date(x_utc, tz = "UTC")
  format(x_date, "%Y-%m-%d")
}

ym_from_filename <- function(nc_file) {
  base <- fs::path_file(nc_file)
  m <- regexec("([0-9]{4})[_-]([0-9]{2})", base)
  g <- regmatches(base, m)[[1]]
  if (length(g) >= 3) return(sprintf("%s-%s", g[2], g[3]))
  base
}

# ----------------------------------------------------------------------
# NEW: Build mask by rasterizing polygon onto template grid
# ----------------------------------------------------------------------
build_polygon_mask <- function(template_layer, shp_path,
                               keep_codes = NULL,
                               field = "WSA9",
                               dilate_cells = 0L) {
  
  stopifnot(inherits(template_layer, "SpatRaster"))
  
  # Read polygons (sf), optionally filter
  eco <- sf::read_sf(shp_path)
  
  if (!is.null(keep_codes)) {
    eco <- eco[eco[[field]] %in% keep_codes, , drop = FALSE]
    if (nrow(eco) == 0) stop("No polygons after filtering keep_codes on field='", field, "'.")
  }
  
  # Ensure valid geometries (important for rasterize)
  eco <- sf::st_make_valid(eco)
  
  # Transform polygons to template CRS
  # template CRS is typically lon/lat; eco is Albers in your example
  eco <- sf::st_transform(eco, crs = terra::crs(template_layer))
  
  # Convert to terra vect
  v <- terra::vect(eco)
  
  # Rasterize onto template grid:
  # - background = NA
  # - field not strictly needed; we'll burn 1s
  r0 <- terra::rast(template_layer)
  m  <- terra::rasterize(v, r0, field = 1, background = NA)
  
  # Optional dilation in target grid space
  if (isTRUE(dilate_cells > 0L)) {
    k <- 2L * as.integer(dilate_cells) + 1L
    w <- matrix(1, k, k)
    m <- terra::focal(m, w = w, fun = max, na.rm = TRUE)
  }
  
  # Force strict 0/1 with NA outside
  m <- terra::mask(m, m, maskvalues = 0)
  m[m > 0] <- 1
  m
}

# Core: one monthly file -> daily .qs outputs ----
process_one_nc_to_daily_qs <- function(nc_file, out_dir,
                                       target_units = NULL,
                                       keep_zero = TRUE,
                                       drop_na = TRUE,
                                       mask_file = NULL) {
  
  r <- terra::rast(nc_file)
  
  tt <- terra::time(r)
  if (is.null(tt) || length(tt) != terra::nlyr(r)) {
    stop("Could not read a valid time vector from NetCDF via terra::time(). File: ", nc_file)
  }
  tt <- as.POSIXct(tt, tz = "UTC")
  attr(tt, "tzone") <- "UTC"
  
  # Units (best-effort)
  from_unit <- NA_character_
  u <- tryCatch(terra::units(r), error = function(e) NULL)
  if (!is.null(u) && length(u) >= 1) from_unit <- u[[1]]
  
  # Apply mask from disk (robust across workers)
  if (!is.null(mask_file)) {
    m <- terra::rast(mask_file)
    # Guard: if mask is not aligned for any reason, force alignment
    if (!terra::compareGeom(r[[1]], m, stopOnError = FALSE)) {
      m <- terra::resample(m, r[[1]], method = "near")
      m[m > 0] <- 1
    }
    r <- terra::mask(r, m)
  }
  
  # Split layer indices by UTC date
  dk <- date_key_utc(tt)
  idx_by_day <- split(seq_along(tt), dk)
  
  for (d in names(idx_by_day)) {
    idx <- idx_by_day[[d]]
    sub <- r[[idx]]
    
    df <- terra::as.data.frame(sub, xy = TRUE, na.rm = TRUE)
    if (nrow(df) == 0) next
    
    dt <- data.table::as.data.table(df)
    data.table::setnames(dt, c("x", "y"), c("lon", "lat"))
    val_cols <- setdiff(names(dt), c("lon", "lat"))
    
    dt_long <- data.table::melt(
      dt,
      id.vars = c("lon", "lat"),
      measure.vars = val_cols,
      variable.name = "layer",
      value.name = "intensity",
      variable.factor = FALSE
    )
    
    # Map layer order -> datetime
    layer_map <- data.table::data.table(layer = val_cols, t_index = seq_along(val_cols))
    dt_long <- layer_map[dt_long, on = "layer"]
    dt_long[, datetime := tt[idx][t_index]]
    attr(dt_long$datetime, "tzone") <- "UTC"
    
    dt_long[, datetime_utc := format(datetime, "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")]
    
    # Convert units
    if (!is.null(target_units)) {
      dt_long[, intensity := convert_units(intensity, from_unit = from_unit, target_unit = target_units)]
      dt_long[, units := target_units]
    } else {
      dt_long[, units := from_unit]
    }
    
    # Clean
    dt_long[, c("layer", "t_index") := NULL]
    if (drop_na)     dt_long <- dt_long[!is.na(intensity)]
    if (!keep_zero)  dt_long <- dt_long[intensity != 0]
    if (nrow(dt_long) == 0) next
    
    out_file <- fs::path(out_dir, glue::glue("{d}.qs"))
    
    if (fs::file_exists(out_file)) {
      old <- qs::qread(out_file)
      old <- data.table::as.data.table(old)
      dt_long <- data.table::rbindlist(list(old, dt_long), use.names = TRUE, fill = TRUE)
    }
    
    data.table::setkey(dt_long, lon, lat, datetime_utc)
    dt_long <- unique(dt_long)
    
    qs::qsave(dt_long, out_file, preset = "high")
  }
  
  invisible(TRUE)
}

# DISCOVER FILES + PRECOMPUTE MASK FILE ONCE ----
fs::dir_create(OUT_DIR)

nc_files <- fs::dir_ls(IN_DIR, glob = FILE_GLOB, type = "file")
stopifnot(length(nc_files) >= 1)
nc_files <- sort(nc_files)

message(glue::glue("[NFILES] {length(nc_files)}"))

MASK_FILE <- fs::path(OUT_DIR, "_mask_target.tif")

# If you previously wrote a bad mask and want to force rebuild, uncomment:
if (fs::file_exists(MASK_FILE)) fs::file_delete(MASK_FILE)

if (!fs::file_exists(MASK_FILE)) {
  tmpl <- terra::rast(nc_files[[1]])[[1]]
  
  mask_target <- build_polygon_mask(
    template_layer = tmpl,
    shp_path = ECOREGION_SHP,
    keep_codes = ECOREGION_KEEP,
    field = "WSA9",
    dilate_cells = MASK_DILATE_CELLS
  )
  
  terra::writeRaster(mask_target, MASK_FILE, overwrite = TRUE)
  message(glue::glue("[MASK] wrote {MASK_FILE}"))
  
  # quick sanity
  s <- terra::global(mask_target, fun = "sum", na.rm = TRUE)
  message(glue::glue("[MASK] sum(1s) = {as.numeric(s[1,1])}"))
} else {
  message(glue::glue("[MASK] using existing {MASK_FILE}"))
}

# Optional quick plot check
# plot(terra::rast(MASK_FILE))

# PARALLEL RUN ACROSS MONTHS ----
progressr::handlers(global = TRUE)
future::plan(future::multisession, workers = N_WORKERS)

results <- progressr::with_progress({
  p <- progressr::progressor(along = nc_files)
  
  future.apply::future_lapply(
    nc_files,
    function(f) {
      ym <- ym_from_filename(f)
      base <- fs::path_file(f)
      p(message = glue::glue("{ym}  |  {base}"))
      
      process_one_nc_to_daily_qs(
        nc_file = f,
        out_dir = OUT_DIR,
        target_units = TARGET_UNITS,
        keep_zero = KEEP_ZERO,
        drop_na = DROP_NA,
        mask_file = MASK_FILE
      )
      
      list(file = f, ym = ym, ok = TRUE, error = NA_character_)
    },
    future.seed = TRUE
  )
})

# Summarize
res_dt <- data.table::rbindlist(results, fill = TRUE)
bad <- res_dt[ok == FALSE | !is.na(error)]
if (nrow(bad) > 0) {
  message("\n[FAILURES]")
  print(bad)
  stop("Some months failed. See failures above.")
}

message("[DONE] All months processed.")

# # Parallel "bad qs" check (zstd / corrupt files) -------------------------------
# suppressPackageStartupMessages({
#   library(fs)
#   library(qs)
#   library(future)
#   library(future.apply)
#   library(progressr)
#   library(data.table)
# })
# 
# OUT_DIR    <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation"
# N_WORKERS  <- 36L
# 
# qs_files <- fs::dir_ls(OUT_DIR, glob = "*.qs", type = "file")
# message("[CHECK] n_qs = ", length(qs_files))
# 
# progressr::handlers(global = TRUE)
# future::plan(future::multisession, workers = N_WORKERS)
# 
# res <- progressr::with_progress({
#   p <- progressr::progressor(along = qs_files)
#   
#   future.apply::future_lapply(
#     qs_files,
#     function(f) {
#       p(message = fs::path_file(f))
#       t0 <- Sys.time()
#       
#       ok <- TRUE
#       msg <- NA_character_
#       
#       tryCatch(
#         {
#           qs::qread(f)  # just test decompress/read
#         },
#         error = function(e) {
#           ok  <<- FALSE
#           msg <<- conditionMessage(e)
#         }
#       )
#       
#       list(
#         file = as.character(f),
#         ok = ok,
#         msg = msg,
#         elapsed_s = as.numeric(difftime(Sys.time(), t0, units = "secs"))
#       )
#     },
#     future.seed = TRUE
#   )
# })
# 
# future::plan(future::sequential)
# 
# res_dt <- data.table::rbindlist(res, fill = TRUE)
# bad_dt <- res_dt[ok == FALSE]
# 
# message("[CHECK DONE] bad = ", nrow(bad_dt))
# 
# # View the bad ones (and their error messages)
# if (nrow(bad_dt) > 0) {
#   print(bad_dt[order(-elapsed_s)][, .(file, msg, elapsed_s)])
# }
# 
# # Optional: delete bad files (so they get regenerated cleanly)
# # if (nrow(bad_dt) > 0) fs::file_delete(bad_dt$file)
# 
# # Optional: quarantine bad files instead of deleting
# # if (nrow(bad_dt) > 0) {
# #   qdir <- fs::path(OUT_DIR, "_corrupt_qs_quarantine")
# #   fs::dir_create(qdir)
# #   for (f in bad_dt$file) {
# #     fs::file_move(f, fs::path(qdir, fs::path_file(f)))
# #   }
# # }
# 
# library(fs)
# 
# bad_files <- c(
#   "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation/1982-02-05.qs",
#   "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation/1943-04-05.qs",
#   "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation/1940-02-05.qs",
#   "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation/2010-02-01.qs",
#   "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation/1970-06-05.qs"
# )
# 
# fs::file_delete(bad_files)
# 
# # qdir <- fs::path("V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/precipitation", "_corrupt_qs_quarantine")
# # fs::dir_create(qdir)
# # 
# # stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
# # for (f in bad_files) {
# #   if (!fs::file_exists(f)) next
# #   fs::file_move(f, fs::path(qdir, paste0(fs::path_ext_remove(fs::path_file(f)), "_", stamp, ".qs")))
# # }