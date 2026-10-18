# Extract monthly .nc to daily .qs ----
#
# scripts/nc_months_to_daily_qs_temperature_parallel_full.R
# ERA5 monthly NetCDF cubes (.nc) -> DAILY .qs tables (one file per UTC day)
#
# Outputs:
#   OUT_DIR/1940-01-01.qs ... OUT_DIR/2023-12-31.qs
#
# Requirements:
# - Aggregate by UTC calendar date (00:00–23:00 UTC)
# - Parallelize across monthly files (N_WORKERS)
# - Use on-disk mask file (_mask_target.tif) for robustness across workers
# - Rerun-safe (append+dedup), though monthlies should not overlap dates
# - Add datetime_utc string column for explicit timestamps
#
# CHANGE (requested):
# - Rasterize Aggr_Ecoregions_2015 polygons to the ERA5 template grid
#   and use that raster as the mask (instead of Atlas14/Köppen rasters).

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
IN_DIR  <- "V:/users/hquintal/phd1_cluster_southeast/data/input/temperature"
OUT_DIR <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars/temperature"
FILE_GLOB <- "*.nc"

TARGET_UNITS <- "C"     # temperature: K -> C
KEEP_ZERO <- TRUE
DROP_NA   <- TRUE

# NEW: polygon mask source
ECOREGION_SHP <- "V:/users/hquintal/phd1_cluster_southeast/data/input/regional_aggregation/ecoregion/Aggr_Ecoregions_2015/Aggr_Ecoregions_2015.shp"
# Optional: filter to subset by WSA9 codes (NULL keeps all 9)
ECOREGION_KEEP <- NULL  # e.g., c("SAP","CPL") or NULL
ECOREGION_FIELD <- "WSA9"

MASK_DILATE_CELLS <- 1L
N_WORKERS <- 36L

terra::terraOptions(progress = 0)

# Helpers ----
convert_units <- function(values, from_unit, target_unit) {
  if (is.na(from_unit) || !nzchar(from_unit) || is.na(target_unit) || !nzchar(target_unit)) return(values)
  from_unit   <- trimws(from_unit)
  target_unit <- trimws(target_unit)
  if (identical(from_unit, target_unit)) return(values)
  
  if (from_unit == "K" && target_unit == "C") return(values - 273.15)
  
  warning(sprintf("No conversion rule from '%s' to '%s'. Leaving unchanged.", from_unit, target_unit))
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

# NEW: polygon -> raster mask on ERA5 template grid
build_polygon_mask <- function(template_layer, shp_path,
                               keep_codes = NULL,
                               field = "WSA9",
                               dilate_cells = 0L) {
  
  stopifnot(inherits(template_layer, "SpatRaster"))
  eco <- sf::read_sf(shp_path)
  
  if (!is.null(keep_codes)) {
    eco <- eco[eco[[field]] %in% keep_codes, , drop = FALSE]
    if (nrow(eco) == 0) stop("No polygons after filtering keep_codes on field='", field, "'.")
  }
  
  eco <- sf::st_make_valid(eco)
  eco <- sf::st_transform(eco, crs = terra::crs(template_layer))
  v   <- terra::vect(eco)
  
  # Burn 1s wherever polygon exists; NA outside
  r0 <- terra::rast(template_layer)
  m  <- terra::rasterize(v, r0, field = 1, background = NA)
  
  if (isTRUE(dilate_cells > 0L)) {
    k <- 2L * as.integer(dilate_cells) + 1L
    w <- matrix(1, k, k)
    m <- terra::focal(m, w = w, fun = max, na.rm = TRUE)
  }
  
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
    
    # Guard: if mask grid isn't aligned (cheap insurance)
    if (!terra::compareGeom(r[[1]], m, stopOnError = FALSE)) {
      m <- terra::resample(m, r[[1]], method = "near")
      m[m > 0] <- 1
    }
    
    r <- terra::mask(r, m)
  }
  
  # Split layer indices by UTC date
  dk <- date_key_utc(tt)
  idx_by_day <- split(seq_along(tt), dk)
  
  # For each day: extract those layers, melt, write
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
    
    # Explicit string timestamp
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
    
    # Rerun-safe append+dedup
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

# If you previously wrote a bad/old mask and want to force rebuild, uncomment:
# if (fs::file_exists(MASK_FILE)) fs::file_delete(MASK_FILE)

if (!fs::file_exists(MASK_FILE)) {
  tmpl <- terra::rast(nc_files[[1]])[[1]]
  
  mask_target <- build_polygon_mask(
    template_layer = tmpl,
    shp_path = ECOREGION_SHP,
    keep_codes = ECOREGION_KEEP,
    field = ECOREGION_FIELD,
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