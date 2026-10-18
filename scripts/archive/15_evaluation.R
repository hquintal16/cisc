# =============================================================================
# Evaluate NOAA Storm Events rasters vs clustered exceedances using WFO extents
#
# Outputs (resume-safe):
#   data/out/09_performance/noaa_wfo_daily/*.csv
#   data/out/09_performance/lookup_xy_to_cwa/lookup_xy_to_cwa_<dataset_tag>.qs
#   data/out/09_performance/modeled_wfo_daily/<scenario>/modeled_wfo_daily_<dataset_tag>.qs
#   data/out/09_performance/confusion_daily/confusion_daily_<var>_<event>.qs
#   data/out/09_performance/metrics/metrics_<var>_<event>.csv (+ .qs)
#
# Observed (NOAA): 1 if any Storm Events record overlaps that WFO on that day
# Modeled (clusters): 1 if any clustered point (cluster>0 AND non-empty cluster_id)
#                     occurs inside that WFO on that day
#
# NOTE on AUC: with binary 0/1 "scores", ROC-AUC is degenerate; we report
#   AUC ~= (Sensitivity + Specificity) / 2  (balanced accuracy)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(fs)
  library(stringr)
  library(lubridate)
  library(sf)
  library(terra)
  library(future)
  library(future.apply)
  library(progressr)
})

# ----------------------------
# USER CONFIG
# ----------------------------
BASE_V2    <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out"
CLUST_BASE <- fs::path(BASE_V2, "06_clusters")   # contains hi1_pr1/, hi2_pr2/, ...
OUT_BASE   <- fs::path(BASE_V2, "09_performance")
dir_create(OUT_BASE, recurse = TRUE)

# WFO shapefile (NAD83)
WFO_SHP <- "V:/users/hquintal/phd1_cluster_southeast/data/input/regional_aggregation/nws_wfo_cws/w_05mr24/w_05mr24.shp"

# Cluster dataset folders expected under each scenario hiX_prY
DATASET_TAGS <- c(
  "precipitation_0.25","precipitation_0.17_bilinear","precipitation_0.17_max",
  "heatindex_0.25","heatindex_0.44_bilinear","heatindex_0.44_max"
)

.lookup_has_cell_id <- function(lookup_qs) {
  if (!fs::file_exists(lookup_qs)) return(FALSE)
  x <- tryCatch(qread(lookup_qs), error = function(e) NULL)
  if (is.null(x)) return(FALSE)
  is.data.frame(x) && ("cell_id" %in% names(x)) && ("CWA" %in% names(x))
}

# ----------------------------
# TEMPLATE RASTERS (one per resolution family)
# Use FIRST tif in each directory you provided.
# ----------------------------

TEMPLATE_ROOT_PR <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/08_rasters/precipitation/rp1"
TEMPLATE_ROOT_HI <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/08_rasters/heatindex/rp1"

# Map dataset_tag -> directory that contains representative rasters at that resolution
.template_dir_for_dataset <- function(dataset_tag) {
  if (dataset_tag == "precipitation_0.17_bilinear") return(fs::path(TEMPLATE_ROOT_PR, "res0.17_bilinear"))
  if (dataset_tag == "precipitation_0.17_max")      return(fs::path(TEMPLATE_ROOT_PR, "res0.17_max"))
  if (dataset_tag == "precipitation_0.25")          return(fs::path(TEMPLATE_ROOT_PR, "res0.25_native"))
  
  if (dataset_tag == "heatindex_0.44_bilinear")     return(fs::path(TEMPLATE_ROOT_HI, "res0.44_bilinear"))
  if (dataset_tag == "heatindex_0.44_max")          return(fs::path(TEMPLATE_ROOT_HI, "res0.44_max"))
  if (dataset_tag == "heatindex_0.25")              return(fs::path(TEMPLATE_ROOT_HI, "res0.25_native"))
  
  stop("No template directory rule for dataset_tag='", dataset_tag, "'.")
}

.get_first_tif_in_dir <- function(dir_path) {
  if (!fs::dir_exists(dir_path)) stop("Template dir does not exist: ", dir_path)
  tif <- fs::dir_ls(fs::path(dir_path, "maxrp"), type = "file", glob = "*.tif", recurse = TRUE, fail = FALSE)
  tif <- sort(as.character(tif))
  if (!length(tif)) stop("No .tif found under: ", dir_path)
  tif[1]
}

.get_template_rast <- function(dataset_tag) {
  d <- .template_dir_for_dataset(dataset_tag)
  f <- .get_first_tif_in_dir(d)
  r <- terra::rast(f)
  if (is.na(terra::crs(r, proj = TRUE))) stop("Template raster CRS missing: ", f)
  r
}

# Time window (strings can be "YYYY", "YYYY-MM", or "YYYY-MM-DD"; leave NA for auto)
START_DATE <- NA_character_  # e.g. "1995" or "2023-01" or "2023-01-01"
END_DATE   <- NA_character_  # e.g. "2023" or "2023-04" or "2023-04-15"

# Parallelism (default for stages that use future_lapply)
WORKERS <- 10L
plan(multisession, workers = WORKERS)

# progress (elapsed + ETA)
handlers(global = TRUE)
handlers("cli")
options(progressr.format = "[:bar] :percent | elapsed: :elapsed | eta: :eta | :message")
options(progressr.show_after = 0)

# Performance / resume controls
OVERWRITE_NOAA       <- FALSE
OVERWRITE_MODELED    <- FALSE
OVERWRITE_CONFUSION  <- TRUE
OVERWRITE_METRICS    <- TRUE

# NEW: lookup cache for (x,y)->CWA per dataset_tag (huge speedup)
LOOKUP_CACHE_DIR <- fs::path(OUT_BASE, "lookup_xy_to_cwa")
fs::dir_create(LOOKUP_CACHE_DIR, recurse = TRUE)

OVERWRITE_LOOKUP <- FALSE

terraOptions(progress = 0)

# ----------------------------
# LOGGING
# ----------------------------
.t0 <- Sys.time()
.log <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      sprintf(...), "\n", sep = "")
}

# ----------------------------
# DATE RANGE PARSING
# ----------------------------
.parse_date_any <- function(x, end = FALSE) {
  if (is.na(x) || !nzchar(x)) return(NA_Date_)
  x <- trimws(x)
  
  # YYYY
  if (grepl("^\\d{4}$", x)) {
    y <- as.integer(x)
    return(if (!end) as.Date(sprintf("%04d-01-01", y)) else as.Date(sprintf("%04d-12-31", y)))
  }
  
  # YYYY-MM
  if (grepl("^\\d{4}-\\d{2}$", x)) {
    y <- as.integer(substr(x, 1, 4))
    m <- as.integer(substr(x, 6, 7))
    if (!end) return(as.Date(sprintf("%04d-%02d-01", y, m)))
    return(as.Date(ceiling_date(as.Date(sprintf("%04d-%02d-01", y, m)), "month") - days(1)))
  }
  
  # YYYY-MM-DD
  if (grepl("^\\d{4}-\\d{2}-\\d{2}$", x)) return(as.Date(x))
  
  stop("Could not parse date: ", x, " (expected YYYY, YYYY-MM, or YYYY-MM-DD)")
}

.apply_window <- function(dt, date_col = "date") {
  if (!nrow(dt)) return(dt)
  d0 <- .parse_date_any(START_DATE, end = FALSE)
  d1 <- .parse_date_any(END_DATE,   end = TRUE)
  
  if (is.na(d0)) d0 <- min(dt[[date_col]], na.rm = TRUE)
  if (is.na(d1)) d1 <- max(dt[[date_col]], na.rm = TRUE)
  
  dt[dt[[date_col]] >= d0 & dt[[date_col]] <= d1]
}

# ----------------------------
# WFO PREP
# ----------------------------
.wfo_load <- function() {
  wfo <- sf::st_read(WFO_SHP, quiet = TRUE)
  wfo_4326 <- sf::st_transform(wfo, 4326)
  wfo_4326 <- wfo_4326[, c("CWA","WFO","STATE","REGION"), drop = FALSE]
  wfo_4326
}

WFO_SF   <- .wfo_load()
WFO_V    <- terra::vect(WFO_SF)  # SpatVector for fast terra::extract()
WFO_KEYS <- data.table(CWA = toupper(trimws(as.character(WFO_SF$CWA))))

.log("[WFO] features=%d | unique CWA=%d", nrow(WFO_SF), length(unique(WFO_KEYS$CWA)))

# ----------------------------
# DISCOVERY HELPERS
# ----------------------------
.discover_hi_pr_folders <- function() {
  all_dirs <- fs::dir_ls(CLUST_BASE, type = "directory", fail = FALSE)
  all_dirs <- sort(as.character(all_dirs))
  hi_pr <- all_dirs[grepl("^hi\\d+_pr\\d+$", fs::path_file(all_dirs))]
  hi_pr <- sort(hi_pr)
  if (!length(hi_pr)) stop("No hi*_pr* folders found under: ", CLUST_BASE)
  hi_pr
}

.parse_floors <- function(hi_pr_dir) {
  bn <- fs::path_file(hi_pr_dir)
  m <- stringr::str_match(bn, "^hi(\\d+)_pr(\\d+)$")
  list(hi = as.integer(m[,2]), pr = as.integer(m[,3]))
}

.discover_dataset_dirs <- function(hi_pr_dir) {
  ds <- fs::dir_ls(hi_pr_dir, type = "directory", fail = FALSE)
  ds <- sort(as.character(ds))
  ds <- ds[fs::path_file(ds) %in% DATASET_TAGS]
  ds
}

.var_group <- function(dataset_tag) {
  if (grepl("^heatindex", dataset_tag, ignore.case = TRUE)) "heatindex" else "precipitation"
}

# =============================================================================
# NOAA STAGE A: WFO x day observed presence (from CSVs)
# =============================================================================

NOAA_EVENTS_DIR <- "V:/users/hquintal/phd1_cluster_southeast/data/output/04_noaa"

NOAA_FILES <- c(
  fs::path(NOAA_EVENTS_DIR, "NOAA_events.csv"),
  fs::path(NOAA_EVENTS_DIR, "NOAA_events_2009_2015.csv"),
  fs::path(NOAA_EVENTS_DIR, "NOAA_events_2009-2015.csv")
)
NOAA_FILES <- NOAA_FILES[fs::file_exists(NOAA_FILES)]
if (!length(NOAA_FILES)) stop("No NOAA CSV files found in: ", NOAA_EVENTS_DIR)

# Map evaluation categories -> Storm Events EVENT_TYPE labels
NOAA_EVENT_MAP <- list(
  heatindex = list(
    excess_heat = c("Excessive Heat"),
    heat        = c("Heat")
  ),
  precipitation = list(
    flash_flood    = c("Flash Flood"),
    flood          = c("Flood"),
    heavy_rain     = c("Heavy Rain"),
    tropical_storm = c("Hurricane", "Tropical Depression", "Tropical Storm", "Typhoon")
  )
)

.pick_col <- function(nms, candidates_regex) {
  hit <- NA_character_
  for (pat in candidates_regex) {
    m <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
    if (length(m)) { hit <- m[1]; break }
  }
  hit
}

.parse_noaa_time <- function(x) {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = "UTC"))
  if (is.numeric(x)) return(as.POSIXct(as.Date(x, origin = "1899-12-30"), tz = "UTC"))
  
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  if (all(is.na(x))) return(as.POSIXct(rep(NA, length(x)), tz = "UTC"))
  
  suppressWarnings(lubridate::parse_date_time(
    x,
    orders = c(
      "mdy HMS", "mdy HM", "mdy",
      "ymd HMS", "ymd HM", "ymd",
      "mdyHMS", "mdyHM",
      "ymdHMS", "ymdHM"
    ),
    tz = "UTC"
  ))
}

.read_one_noaa_csv_min <- function(f) {
  hdr <- names(data.table::fread(f, nrows = 0, showProgress = FALSE))
  if (!length(hdr)) stop("Could not read header for: ", f)
  
  col_wfo   <- .pick_col(hdr, c("^WFO$", "WFO", "CWA"))
  col_type  <- .pick_col(hdr, c("^EVENT_TYPE$", "EVENT_TYPE", "EVENT TYPE"))
  col_start <- .pick_col(hdr, c("^START\\.TIMESTAMP\\.UTC$", "START\\.TIMESTAMP", "^BEGIN", "START", "BEGIN_DATE", "BEGIN_DATE_TIME"))
  col_end   <- .pick_col(hdr, c("^END\\.TIMESTAMP\\.UTC$", "END\\.TIMESTAMP", "^END", "END_DATE", "END_DATE_TIME"))
  
  if (anyNA(c(col_wfo, col_type, col_start))) {
    stop(sprintf(
      "Required columns not found in %s.\nHeader:\n%s",
      fs::path_file(f), paste(hdr, collapse = ", ")
    ))
  }
  
  sel <- unique(na.omit(c(col_wfo, col_type, col_start, col_end)))
  
  dt <- data.table::fread(f, select = sel, showProgress = FALSE) |> as.data.table()
  setnames(dt, col_wfo,   "WFO")
  setnames(dt, col_type,  "EVENT_TYPE")
  setnames(dt, col_start, "start_raw")
  if (!is.na(col_end)) setnames(dt, col_end, "end_raw") else dt[, end_raw := NA_character_]
  
  dt[, WFO := toupper(trimws(as.character(WFO)))]
  dt[, EVENT_TYPE := trimws(as.character(EVENT_TYPE))]
  dt[, start_time := .parse_noaa_time(start_raw)]
  dt[, end_time   := .parse_noaa_time(end_raw)]
  dt[is.na(end_time) & !is.na(start_time), end_time := start_time]
  
  dt <- dt[!is.na(WFO) & nzchar(WFO) & !is.na(EVENT_TYPE) & nzchar(EVENT_TYPE) & !is.na(start_time)]
  if (!nrow(dt)) return(dt[0])
  
  dt[, `:=`(
    start_date = as.Date(start_time, tz = "UTC"),
    end_date   = as.Date(end_time,   tz = "UTC")
  )]
  dt[is.na(end_date), end_date := start_date]
  dt[end_date < start_date, `:=`(start_date = end_date, end_date = start_date)]
  
  dt[, .(WFO, EVENT_TYPE, start_date, end_date)]
}

.read_noaa_events_min <- function() {
  .log("[NOAA-CSV] reading %d file(s)...", length(NOAA_FILES))
  dt_list <- lapply(NOAA_FILES, function(f) {
    .log("[NOAA-CSV] read: %s", fs::path_file(f))
    tryCatch(.read_one_noaa_csv_min(f), error = function(e) {
      .log("[NOAA-CSV][ERR] %s | %s", fs::path_file(f), e$message)
      NULL
    })
  })
  dt <- rbindlist(Filter(Negate(is.null), dt_list), use.names = TRUE, fill = TRUE)
  if (!nrow(dt)) stop("NOAA CSV read produced 0 rows after filtering.")
  dt[, WFO := toupper(trimws(as.character(WFO)))]
  dt[]
}

.noaa_events_to_wfo_daily_csv <- function(noaa_dt, event_types, out_csv, label) {
  
  et <- tolower(event_types)
  dt <- noaa_dt[tolower(EVENT_TYPE) %in% et]
  if (!nrow(dt)) {
    warning("No NOAA rows for: ", label, " event_types={", paste(event_types, collapse=", "), "}")
    return(invisible(NULL))
  }
  
  d0 <- .parse_date_any(START_DATE, end = FALSE)
  d1 <- .parse_date_any(END_DATE,   end = TRUE)
  if (!is.na(d0)) dt <- dt[end_date >= d0]
  if (!is.na(d1)) dt <- dt[start_date <= d1]
  if (!nrow(dt)) return(invisible(NULL))
  
  minD <- min(dt$start_date, na.rm = TRUE)
  maxD <- max(dt$end_date,   na.rm = TRUE)
  all_days <- seq(minD, maxD, by = "day")
  
  days <- CJ(WFO = sort(unique(dt$WFO)), date = all_days, unique = TRUE)
  days[, `:=`(dstart = date, dend = date)]
  setkey(days, WFO, dstart, dend)
  
  ev <- dt[, .(WFO, estart = start_date, eend = end_date)]
  setkey(ev, WFO, estart, eend)
  
  ov <- foverlaps(days, ev,
                  by.x = c("WFO","dstart","dend"),
                  by.y = c("WFO","estart","eend"),
                  nomatch = 0L)
  
  obs <- ov[, .(observed = 1L), by = .(CWA = WFO, date)]
  
  full <- CJ(CWA = WFO_KEYS$CWA, date = all_days, unique = TRUE)
  out  <- merge(full, obs, by = c("CWA","date"), all.x = TRUE)
  out[, observed := fifelse(is.na(observed), 0L, observed)]
  
  out <- .apply_window(out, "date")
  setorder(out, date, CWA)
  
  fs::dir_create(fs::path_dir(out_csv), recurse = TRUE)
  fwrite(out, out_csv)
  .log("[NOAA-CSV] wrote: %s  (rows=%d; days=%d)", out_csv, nrow(out), length(unique(out$date)))
  
  invisible(out_csv)
}

# ---- Read NOAA once
noaa_all <- .read_noaa_events_min()

# Remap legacy airport IDs -> CWAs (expand if needed)
WFO_REMAP <- c("EYW" = "KEY")
noaa_all[WFO %in% names(WFO_REMAP), WFO := unname(WFO_REMAP[WFO])]

.log("[NOAA-CSV] unique WFO in events: %d", length(unique(noaa_all$WFO)))
.log("[WFO-SHP] unique CWA in shapefile: %d", length(unique(WFO_KEYS$CWA)))
.log("[NOAA-CSV] WFO not in shapefile: %d",
     length(setdiff(unique(noaa_all$WFO), unique(WFO_KEYS$CWA))))

# =============================================================================
# MODELED STAGE B: clusters -> WFO x day modeled presence (FAST via lookup)
# =============================================================================

.standardize_xy <- function(dt) {
  nms <- names(dt)
  if (all(c("x","y") %in% nms)) return(dt)
  if (all(c("lon","lat") %in% nms)) { setnames(dt, c("lon","lat"), c("x","y")); return(dt) }
  if (all(c("longitude","latitude") %in% nms)) { setnames(dt, c("longitude","latitude"), c("x","y")); return(dt) }
  stop("No x/y or lon/lat columns in clustered qs.")
}

.get_date_from_cluster_dt <- function(dt) {
  if ("day"  %in% names(dt)) return(as.Date(dt$day))
  if ("date" %in% names(dt)) return(as.Date(dt$date))
  if ("time_utc" %in% names(dt)) return(as.Date(as.POSIXct(dt$time_utc, tz="UTC")))
  if ("time"     %in% names(dt)) return(as.Date(as.POSIXct(dt$time,     tz="UTC")))
  stop("No recognizable day/date/time column in clustered qs.")
}

.build_xy_to_cwa_lookup <- function(dataset_tag, out_lookup_qs) {
  
  if (!OVERWRITE_LOOKUP && .lookup_has_cell_id(out_lookup_qs)) {
    .log("[LOOKUP] exists (skip): %s", out_lookup_qs)
    return(invisible(out_lookup_qs))
  }
  
  if (fs::file_exists(out_lookup_qs) && !OVERWRITE_LOOKUP) {
    .log("[LOOKUP] stale lookup (no cell_id) -> rebuilding: %s", out_lookup_qs)
  }
  
  # Template raster defines the grid for this dataset_tag (fixed grid)
  r <- .get_template_rast(dataset_tag)
  n_cells <- terra::ncell(r)
  
  .log("[LOOKUP] template grid: %s | ncell=%d", dataset_tag, n_cells)
  
  # Cell centers -> points
  xy <- terra::xyFromCell(r, seq_len(n_cells))
  pts <- terra::vect(as.data.frame(xy), geom = c("x","y"), crs = terra::crs(r, proj = TRUE))
  
  # Extract WFO polygon attributes at each cell-center point
  ex <- terra::extract(WFO_V, pts)
  if (!nrow(ex) || !("CWA" %in% names(ex))) {
    stop("[LOOKUP] terra::extract did not return CWA. Check WFO_V fields.")
  }
  
  lookup <- data.table(
    cell_id = seq_len(n_cells),
    x = xy[, 1],
    y = xy[, 2],
    CWA = toupper(trimws(as.character(ex$CWA)))
  )
  
  # Drop non-matches (ocean/outside WFO coverage)
  lookup <- lookup[!is.na(CWA) & nzchar(CWA)]
  
  tmp <- paste0(out_lookup_qs, ".tmp")
  if (fs::file_exists(tmp)) fs::file_delete(tmp)
  qsave(lookup, tmp); fs::file_move(tmp, out_lookup_qs)
  
  .log("[LOOKUP] wrote: %s  (rows=%d; matched=%d/%d)",
       out_lookup_qs, nrow(lookup), nrow(lookup), n_cells)
  
  invisible(out_lookup_qs)
}

.clustered_presence_wfo_daily_one_file <- function(qs_path, lookup_dt, r_template) {
  
  dt <- as.data.table(qread(qs_path))
  if (!nrow(dt)) return(NULL)
  
  dt <- .standardize_xy(dt)
  
  # keep only clustered points
  if ("cluster" %in% names(dt))    dt <- dt[is.finite(cluster) & cluster > 0]
  if ("cluster_id" %in% names(dt)) dt <- dt[!is.na(cluster_id) & nzchar(as.character(cluster_id))]
  if (!nrow(dt)) return(NULL)
  
  dt[, date := .get_date_from_cluster_dt(dt)]
  dt <- dt[is.finite(x) & is.finite(y) & !is.na(date)]
  if (!nrow(dt)) return(NULL)
  
  # compute raster cell index for each point (lon/lat in EPSG:4326)
  pts <- terra::vect(dt[, .(x, y)], geom = c("x","y"), crs = terra::crs(r_template, proj = TRUE))
  dt[, cell_id := terra::cellFromXY(r_template, terra::crds(pts))]
  
  dt <- dt[!is.na(cell_id)]
  if (!nrow(dt)) return(NULL)
  
  # join cell_id -> CWA
  lookup_dt[, cell_id := as.integer(cell_id)]
  setkey(lookup_dt, cell_id)
  dt <- lookup_dt[dt, on = "cell_id"]     # brings in CWA from lookup_dt
  
  # after join: CWA is in dt$CWA (not i.CWA)
  dt <- dt[!is.na(CWA) & nzchar(CWA)]
  if (!nrow(dt)) return(NULL)
  
  unique(dt[, .(CWA, date)][, modeled := 1L], by = c("CWA","date"))
}

.build_modeled_wfo_daily <- function(out_qs, ds_dir, dataset_tag) {
  
  if (fs::file_exists(out_qs) && !OVERWRITE_MODELED) {
    .log("[MODEL] exists (skip): %s", out_qs)
    return(invisible(out_qs))
  }
  
  qs_files <- fs::dir_ls(ds_dir, type = "file", glob = "*.qs", fail = FALSE)
  qs_files <- sort(as.character(qs_files))
  if (!length(qs_files)) return(invisible(NULL))
  
  r_template <- .get_template_rast(dataset_tag)
  
  # Template raster + lookup for this dataset_tag
  r_template <- .get_template_rast(dataset_tag)
  
  lookup_qs <- fs::path(LOOKUP_CACHE_DIR, sprintf("lookup_xy_to_cwa_%s.qs", dataset_tag))
  .build_xy_to_cwa_lookup(dataset_tag, lookup_qs)
  
  lookup_dt <- as.data.table(qread(lookup_qs))
  if (!nrow(lookup_dt)) stop("[MODEL] lookup empty for dataset_tag=", dataset_tag)
  lookup_dt[, CWA := toupper(trimws(as.character(CWA)))]
  lookup_dt[, cell_id := as.integer(cell_id)]
  setkey(lookup_dt, cell_id)
  
  .log("[MODEL] %s | files=%d | lookup_cells=%d",
       ds_dir, length(qs_files), nrow(lookup_dt))
  
  parts <- lapply(qs_files, function(f) {
    tryCatch(.clustered_presence_wfo_daily_one_file(f, lookup_dt, r_template),
             error = function(e) {
               .log("[MODEL][ERR] %s | %s", fs::path_file(f), e$message)
               NULL
             })
  })
  
  dt <- rbindlist(Filter(Negate(is.null), parts), use.names = TRUE, fill = TRUE)
  if (!nrow(dt)) return(invisible(NULL))
  
  dt <- dt[, .(modeled = 1L), by = .(CWA, date)]
  dt[, CWA := toupper(trimws(as.character(CWA)))]
  setorder(dt, date, CWA)
  
  tmp <- paste0(out_qs, ".tmp")
  if (fs::file_exists(tmp)) fs::file_delete(tmp)
  qsave(dt, tmp); fs::file_move(tmp, out_qs)
  
  .log("[MODEL] wrote: %s  (rows=%d)", out_qs, nrow(dt))
  invisible(out_qs)
}

# =============================================================================
# CONFUSION + METRICS
# =============================================================================

.confusion_daily <- function(obs_dt, mod_dt) {
  all_dates <- sort(unique(c(obs_dt$date, mod_dt$date)))
  grid <- CJ(date = all_dates, CWA = WFO_KEYS$CWA, unique = TRUE)
  
  obs <- merge(grid, obs_dt, by = c("date","CWA"), all.x = TRUE)
  mod <- merge(grid, mod_dt, by = c("date","CWA"), all.x = TRUE)
  
  obs[, observed := fifelse(is.na(observed), 0L, as.integer(observed))]
  mod[, modeled  := fifelse(is.na(modeled),  0L, as.integer(modeled))]
  
  x <- merge(obs[, .(date,CWA,observed)], mod[, .(date,CWA,modeled)], by = c("date","CWA"))
  x[, `:=`(
    TP = as.integer(observed == 1L & modeled == 1L),
    FP = as.integer(observed == 0L & modeled == 1L),
    FN = as.integer(observed == 1L & modeled == 0L),
    TN = as.integer(observed == 0L & modeled == 0L)
  )]
  x
}

.metrics_from_counts <- function(TP, FP, FN, TN) {
  N <- TP + FP + FN + TN
  if (!is.finite(N) || N <= 0) {
    return(list(
      TP=TP,FP=FP,FN=FN,TN=TN,N=N,
      hit_rate=NA_real_, fpr=NA_real_, fnr=NA_real_, precision=NA_real_,
      f1=NA_real_, f2=NA_real_, kappa=NA_real_,
      sensitivity=NA_real_, specificity=NA_real_, accuracy=NA_real_, auc=NA_real_
    ))
  }
  
  sensitivity <- if ((TP+FN) > 0) TP/(TP+FN) else NA_real_
  specificity <- if ((TN+FP) > 0) TN/(TN+FP) else NA_real_
  precision   <- if ((TP+FP) > 0) TP/(TP+FP) else NA_real_
  
  fpr <- if ((FP+TN) > 0) FP/(FP+TN) else NA_real_
  fnr <- if ((FN+TP) > 0) FN/(FN+TP) else NA_real_
  accuracy <- (TP+TN)/N
  
  f1 <- if (is.finite(precision) && is.finite(sensitivity) && (precision+sensitivity) > 0)
    2*precision*sensitivity/(precision+sensitivity) else NA_real_
  
  beta <- 2
  f2 <- if (is.finite(precision) && is.finite(sensitivity) && ((beta^2*precision)+sensitivity) > 0)
    (1+beta^2)*precision*sensitivity/((beta^2*precision)+sensitivity) else NA_real_
  
  p_obs1 <- (TP+FN)/N
  p_obs0 <- (TN+FP)/N
  p_mod1 <- (TP+FP)/N
  p_mod0 <- (TN+FN)/N
  pe <- p_obs1*p_mod1 + p_obs0*p_mod0
  po <- accuracy
  kappa <- if (is.finite(pe) && (1-pe) > 0) (po - pe)/(1 - pe) else NA_real_
  
  auc <- if (is.finite(sensitivity) && is.finite(specificity)) (sensitivity + specificity)/2 else NA_real_
  
  list(
    TP=TP,FP=FP,FN=FN,TN=TN,N=N,
    hit_rate=sensitivity, fpr=fpr, fnr=fnr, precision=precision,
    f1=f1, f2=f2, kappa=kappa,
    sensitivity=sensitivity, specificity=specificity, accuracy=accuracy, auc=auc
  )
}

# =============================================================================
# RUN PIPELINE
# =============================================================================

.log("============================================================")
.log("WORKERS: %d", WORKERS)

hi_pr_folders <- .discover_hi_pr_folders()
.log("Scenarios found: %d", length(hi_pr_folders))

# ---- Build NOAA observed caches (CSV outputs) ----
NOAA_CACHE_DIR <- fs::path(OUT_BASE, "noaa_wfo_daily")
dir_create(NOAA_CACHE_DIR, recurse = TRUE)

noaa_cache_paths <- list()

.write_noaa_if_needed <- function(var_group, evt, out_csv) {
  if (fs::file_exists(out_csv) && !OVERWRITE_NOAA) {
    .log("[NOAA-CSV] exists (skip): %s", out_csv)
    return(invisible(out_csv))
  }
  .noaa_events_to_wfo_daily_csv(
    noaa_dt     = noaa_all,
    event_types = NOAA_EVENT_MAP[[var_group]][[evt]],
    out_csv     = out_csv,
    label       = paste(var_group, evt)
  )
  invisible(out_csv)
}

for (evt in names(NOAA_EVENT_MAP$heatindex)) {
  out_csv <- fs::path(NOAA_CACHE_DIR, sprintf("noaa_wfo_daily_heatindex_%s.csv", evt))
  .write_noaa_if_needed("heatindex", evt, out_csv)
  noaa_cache_paths[[paste0("heatindex__", evt)]] <- out_csv
}

for (evt in names(NOAA_EVENT_MAP$precipitation)) {
  out_csv <- fs::path(NOAA_CACHE_DIR, sprintf("noaa_wfo_daily_precipitation_%s.csv", evt))
  .write_noaa_if_needed("precipitation", evt, out_csv)
  noaa_cache_paths[[paste0("precipitation__", evt)]] <- out_csv
}

# ---- Build modeled cache jobs list (scenario x dataset_tag) ----
MODEL_CACHE_DIR <- fs::path(OUT_BASE, "modeled_wfo_daily")
dir_create(MODEL_CACHE_DIR, recurse = TRUE)

scenario_jobs <- list()

for (sc_dir in hi_pr_folders) {
  floors <- .parse_floors(sc_dir)
  sc_name <- fs::path_file(sc_dir)
  
  ds_dirs <- .discover_dataset_dirs(sc_dir)
  if (!length(ds_dirs)) next
  
  for (ds_dir in ds_dirs) {
    dataset_tag <- fs::path_file(ds_dir)
    varg <- .var_group(dataset_tag)
    floor_yr <- if (varg == "heatindex") floors$hi else floors$pr
    
    out_qs <- fs::path(MODEL_CACHE_DIR, sc_name, sprintf("modeled_wfo_daily_%s.qs", dataset_tag))
    dir_create(fs::path_dir(out_qs), recurse = TRUE)
    
    scenario_jobs[[length(scenario_jobs)+1L]] <- list(
      scenario    = sc_name,
      dataset_tag = dataset_tag,
      var_group   = varg,
      floor_yr    = floor_yr,
      ds_dir      = ds_dir,
      out_qs      = out_qs
    )
  }
}

scenario_jobs_dt <- rbindlist(scenario_jobs, use.names = TRUE, fill = TRUE)
if (!nrow(scenario_jobs_dt)) stop("No dataset folders found under scenarios in: ", CLUST_BASE)
.log("Modeled cache jobs discovered: %d", nrow(scenario_jobs_dt))

# ---- Build lookup tables once per dataset_tag (fast, reused across scenarios) ----
unique_tags <- sort(unique(scenario_jobs_dt$dataset_tag))
.log("[LOOKUP] unique dataset_tags=%d", length(unique_tags))

for (tag in unique_tags) {
  out_lookup_qs <- fs::path(LOOKUP_CACHE_DIR, sprintf("lookup_xy_to_cwa_%s.qs", tag))
  .build_xy_to_cwa_lookup(tag, out_lookup_qs)
}

# ---- Build modeled caches (RESUME SAFE + SEQUENTIAL to avoid FutureInterruptError) ----
scenario_jobs_dt[, needs_run := OVERWRITE_MODELED | !fs::file_exists(out_qs)]
todo <- scenario_jobs_dt[needs_run == TRUE]

.log("Modeled cache jobs total=%d | todo=%d | done=%d",
     nrow(scenario_jobs_dt), nrow(todo), nrow(scenario_jobs_dt) - nrow(todo))

if (nrow(todo) > 0) {
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(sequential)
  
  with_progress({
    p <- progressor(steps = nrow(todo))
    ok <- logical(nrow(todo))
    
    for (i in seq_len(nrow(todo))) {
      row <- todo[i]
      p(sprintf("[%d/%d] %s | %s", i, nrow(todo), row$scenario, row$dataset_tag))
      
      ok[i] <- tryCatch({
        .build_modeled_wfo_daily(
          out_qs      = row$out_qs,
          ds_dir      = row$ds_dir,
          dataset_tag = row$dataset_tag
        )
        TRUE
      }, error = function(e) {
        .log("[MODEL][JOB-ERR] %s | %s | %s", row$scenario, row$dataset_tag, e$message)
        FALSE
      })
    }
    
    .log("[MODEL] finished: ok=%d | failed=%d", sum(ok), sum(!ok))
  })
  
} else {
  .log("[MODEL] nothing to do (all modeled caches exist).")
}

# Restore multisession for eval stage (optional, but explicit)
plan(multisession, workers = WORKERS)

# ---- Confusion + metrics per hazard event ----
CONF_DIR <- fs::path(OUT_BASE, "confusion_daily")
MET_DIR  <- fs::path(OUT_BASE, "metrics")
dir_create(c(CONF_DIR, MET_DIR), recurse = TRUE)

.eval_one_event <- function(var_group_in, event_name, obs_csv_path) {
  
  if (!fs::file_exists(obs_csv_path)) {
    warning("Missing observed cache: ", obs_csv_path)
    return(NULL)
  }
  
  obs_dt <- fread(obs_csv_path)
  if (!nrow(obs_dt)) return(NULL)
  obs_dt[, date := as.Date(date)]
  obs_dt[, CWA  := toupper(trimws(as.character(CWA)))]
  obs_dt[, observed := as.integer(observed)]
  obs_dt <- .apply_window(obs_dt, "date")
  
  jobs <- scenario_jobs_dt[var_group == var_group_in]
  if (!nrow(jobs)) return(NULL)
  
  .log("[EVAL] %s | %s | jobs=%d", var_group_in, event_name, nrow(jobs))
  
  parts <- with_progress({
    p <- progressor(steps = nrow(jobs))
    
    future_lapply(seq_len(nrow(jobs)), function(i) {
      row <- jobs[i]
      p(sprintf("%s | %s", row$scenario, row$dataset_tag))
      
      if (!fs::file_exists(row$out_qs)) return(NULL)
      
      mod_dt <- as.data.table(qread(row$out_qs))
      if (!nrow(mod_dt)) return(NULL)
      
      mod_dt[, date := as.Date(date)]
      mod_dt[, CWA  := toupper(trimws(as.character(CWA)))]
      mod_dt <- mod_dt[, .(modeled = 1L), by = .(CWA, date)]
      mod_dt <- .apply_window(mod_dt, "date")
      
      conf <- .confusion_daily(obs_dt, mod_dt)
      conf[, `:=`(
        scenario    = row$scenario,
        dataset_tag = row$dataset_tag,
        floor_yr    = row$floor_yr,
        var_group   = row$var_group,
        event       = event_name
      )]
      conf <- .apply_window(conf, "date")
      if (!nrow(conf)) return(NULL)
      
      agg_wfo <- conf[, .(TP=sum(TP), FP=sum(FP), FN=sum(FN), TN=sum(TN)),
                      by = .(scenario, dataset_tag, floor_yr, var_group, event, CWA)]
      
      agg_all <- conf[, .(TP=sum(TP), FP=sum(FP), FN=sum(FN), TN=sum(TN)),
                      by = .(scenario, dataset_tag, floor_yr, var_group, event)][, CWA := "ALL_WFO"]
      
      agg <- rbindlist(list(agg_wfo, agg_all), use.names = TRUE, fill = TRUE)
      
      met <- agg[, {
        m <- .metrics_from_counts(TP, FP, FN, TN)
        as.list(m)
      }, by = .(scenario, dataset_tag, floor_yr, var_group, event, CWA)]
      
      list(conf = conf, met = met)
    })
  })
  
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(NULL)
  
  conf_all <- rbindlist(lapply(parts, `[[`, "conf"), use.names = TRUE, fill = TRUE)
  met_all  <- rbindlist(lapply(parts, `[[`, "met"),  use.names = TRUE, fill = TRUE)
  
  conf_out   <- fs::path(CONF_DIR, sprintf("confusion_daily_%s_%s.qs", var_group, event_name))
  met_out_csv <- fs::path(MET_DIR, sprintf("metrics_%s_%s.csv", var_group, event_name))
  met_out_qs  <- fs::path(MET_DIR, sprintf("metrics_%s_%s.qs",  var_group, event_name))
  
  if (OVERWRITE_CONFUSION || !fs::file_exists(conf_out)) {
    tmp <- paste0(conf_out, ".tmp")
    if (fs::file_exists(tmp)) fs::file_delete(tmp)
    qsave(conf_all, tmp); fs::file_move(tmp, conf_out)
    .log("[EVAL] wrote: %s  (rows=%d)", conf_out, nrow(conf_all))
  }
  
  if (OVERWRITE_METRICS || !fs::file_exists(met_out_qs)) {
    tmp <- paste0(met_out_qs, ".tmp")
    if (fs::file_exists(tmp)) fs::file_delete(tmp)
    qsave(met_all, tmp); fs::file_move(tmp, met_out_qs)
    fwrite(met_all, met_out_csv)
    .log("[EVAL] wrote: %s  (rows=%d)", met_out_csv, nrow(met_all))
  }
  
  invisible(list(conf = conf_out, metrics = met_out_csv))
}

# ---- Run all events ----
for (evt in names(NOAA_EVENT_MAP$heatindex)) {
  key <- paste0("heatindex__", evt)
  .eval_one_event("heatindex", evt, noaa_cache_paths[[key]])
}

for (evt in names(NOAA_EVENT_MAP$precipitation)) {
  key <- paste0("precipitation__", evt)
  .eval_one_event("precipitation", evt, noaa_cache_paths[[key]])
}

.dt <- difftime(Sys.time(), .t0, units = "mins")
.log("============================================================")
.log("DONE. Elapsed: %.2f minutes", as.numeric(.dt))
.log("Outputs: %s", OUT_BASE)

test <- terra::rast('V:/users/hquintal/phd1_cluster_southeast_v2/data/out/08_rasters/heatindex/rp1/res0.25_native/maxrp/1940-06-08-001.tif')
plot(test)
