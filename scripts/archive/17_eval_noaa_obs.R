# scripts/09perf_01_build_noaa_county_days_definitive.R
# =============================================================================
# NOAA county-day builder (DEFINITIVE; restart-safe)
#
# Produces:
#   1) Per-tag definitive county-day presence:
#        data/out/09_performance/noaa_county_days_definitive/noaa_county_days_<tag>.{qs,csv}
#   2) Per-group definitive county-day presence (dedup across tags within group):
#        .../noaa_county_days__<group>.{qs,csv}   e.g., __all_heat, __all_precip
#   3) Authoritative lookup table (tag + group keys):
#        .../noaa_county_day_lookup.{qs,csv}
#
# Key behaviors:
#   - Multi-stage location resolution: COUNTY -> WFO -> STATE
#   - Treat invalid county GEOIDs as missing so fallbacks still happen
#   - Supports `date` OR (`date_start`,`date_end`) (expanded to daily rows)
#   - Deduplicates by (key, GEOID, date) for group keys (prevents double counting)
#   - Optional failed-row audit output
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(lubridate)
  library(stringr)
  library(sf)
  library(tigris)
  library(qs)
})

options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)

# ----------------------------
# SHARED HELPERS (identical across scripts 1–5)
# ----------------------------
.norm_digits <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_extract(x, "\\d+")
  x[!nzchar(x)] <- NA_character_
  x
}

.norm_geoid5 <- function(x) {
  x <- .norm_digits(x)
  ok <- !is.na(x) & nchar(x) <= 5
  out <- rep(NA_character_, length(x))
  out[ok] <- stringr::str_pad(x[ok], 5, pad = "0")
  out
}

.norm_statefp2 <- function(x) {
  x <- .norm_digits(x)
  ok <- !is.na(x) & nchar(x) <= 2
  out <- rep(NA_character_, length(x))
  out[ok] <- stringr::str_pad(x[ok], 2, pad = "0")
  out
}

.norm_wfo <- function(x) {
  x <- toupper(as.character(x))
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "[^A-Z0-9]+", "")
  x[!nzchar(x)] <- NA_character_
  x
}

event_tag_norm <- function(x) {
  x <- tolower(as.character(x))
  x <- stringr::str_replace_all(x, "[^a-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  x
}

.log <- function(fmt, ...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))

.write_atomic_qs <- function(obj, out_path, preset="high") {
  dir_create(path_dir(out_path), recurse = TRUE)
  tmp <- paste0(out_path, ".tmp")
  if (file_exists(tmp)) file_delete(tmp)
  qsave(obj, tmp, preset = preset)
  file_move(tmp, out_path)
  invisible(out_path)
}

.write_atomic_csv <- function(dt, out_path) {
  dir_create(path_dir(out_path), recurse = TRUE)
  tmp <- paste0(out_path, ".tmp")
  if (file_exists(tmp)) file_delete(tmp)
  fwrite(dt, tmp)
  file_move(tmp, out_path)
  invisible(out_path)
}

# ----------------------------
# USER SETTINGS
# ----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

NOAA_EVENT_SUMMARY_DIR <- file.path(ROOT, "data/out/08_noaa/event_summaries")

# WFO polygons (for county-centroid -> WFO)
WFO_SHP <- "V:/users/hquintal/phd1_cluster_southeast/data/input/regional_aggregation/nws_wfo_cws/w_05mr24/w_05mr24.shp"

# OUTPUTS (canonical)
OUT_DIR <- file.path(ROOT, "data/out/09_performance/noaa_county_days_definitive")
dir_create(OUT_DIR, recurse = TRUE)

# Lookup outputs
LOOKUP_DIR <- file.path(ROOT, "data/out/09_performance/lookup_county_to_wfo")
dir_create(LOOKUP_DIR, recurse = TRUE)
LOOKUP_QS <- file.path(LOOKUP_DIR, "lookup_county_to_wfo_2022_epsg5070.qs")

CRS_ALBERS <- 5070

# Clip window
EVAL_START <- as.Date("1996-01-01")
EVAL_END   <- as.Date("2023-12-31")

# GROUP DEFINITIONS (TOGGLES LIVE HERE)
EVENT_GROUPS <- list(
  all_heat = c("heat", "excessive_heat"
               #, "drought"
               ),
  all_precip = c(
    "flash_flood","flood","heavy_rain","tropical_storm",
    "hurricane_typhoon","tropical_depression","coastal_flood"
    # ,"heavy_snow"
  )
)

# Optional additional tags not in groups
EXTRA_TAGS <- character()

HEAT_TAGS   <- event_tag_norm(EVENT_GROUPS$all_heat)
PRECIP_TAGS <- event_tag_norm(EVENT_GROUPS$all_precip)
TAGS_ALL    <- sort(unique(c(unlist(EVENT_GROUPS, use.names = FALSE), EXTRA_TAGS)))
TAGS_ALL    <- event_tag_norm(TAGS_ALL)

OVERWRITE <- TRUE
WRITE_FAILED_ROWS <- TRUE  # writes per-tag FAILED rows for audit (optional)

# ----------------------------
# Build county->WFO centroid lookup (ONE WFO per county)
# ----------------------------
.build_lookup_county_to_wfo <- function() {
  drop_statefp <- c("02","15","72","60","66","69","78")
  
  COUNTIES <- tigris::counties(cb = TRUE, year = 2022, class = "sf")
  STATES   <- tigris::states(cb = TRUE, year = 2022, class = "sf")
  COUNTIES <- COUNTIES[!(COUNTIES$STATEFP %in% drop_statefp), ]
  STATES   <- STATES[!(STATES$STATEFP %in% drop_statefp), ]
  
  COUNTIES <- st_transform(COUNTIES, CRS_ALBERS)
  STATES   <- st_transform(STATES, CRS_ALBERS)
  
  if (!file_exists(WFO_SHP)) stop("Missing WFO shapefile: ", WFO_SHP)
  WFO <- st_read(WFO_SHP, quiet = TRUE)
  WFO <- st_make_valid(WFO)
  WFO <- st_transform(WFO, CRS_ALBERS)
  
  # Crop to CONUS
  conus_union <- st_union(st_geometry(STATES))
  WFO <- WFO[st_intersects(WFO, conus_union, sparse = FALSE), ]
  
  # Prefer explicit WFO field; otherwise fall back to common candidates
  .pick_wfo_id_field <- function(wfo_sf) {
    candidates <- c("WFO","CWA","wfo","cwa","WFO_ID","SITE","ID","NAME")
    hit <- intersect(candidates, names(wfo_sf))
    if (length(hit)) return(hit[[1]])
    stop("Could not find WFO/CWA id field in WFO shapefile.")
  }
  WFO_ID_FIELD <- if ("WFO" %in% names(WFO)) "WFO" else .pick_wfo_id_field(WFO)
  WFO[[WFO_ID_FIELD]] <- .norm_wfo(WFO[[WFO_ID_FIELD]])
  
  cent <- st_point_on_surface(st_make_valid(COUNTIES))
  cent_join <- st_join(cent, WFO[, WFO_ID_FIELD, drop = FALSE], join = st_within, left = TRUE)
  
  lut <- data.table(
    GEOID   = as.character(COUNTIES$GEOID),
    STATEFP = as.character(COUNTIES$STATEFP),
    wfo     = .norm_wfo(cent_join[[WFO_ID_FIELD]])
  )
  
  n_missing <- sum(is.na(lut$wfo))
  .log("[LOOKUP] county->wfo built: counties=%d | missing_wfo=%d", nrow(lut), n_missing)
  
  list(
    county_to_wfo     = lut,
    wfo_to_counties   = lut[!is.na(wfo), .(GEOID = unique(GEOID)), by = wfo],
    state_to_counties = lut[, .(GEOID = unique(GEOID)), by = STATEFP]
  )
}

LOOKUP <- if (file_exists(LOOKUP_QS) && !OVERWRITE) {
  .log("[LOOKUP] using existing: %s", LOOKUP_QS)
  qread(LOOKUP_QS)
} else {
  .log("[LOOKUP] building county->wfo lookup via centroids...")
  x <- .build_lookup_county_to_wfo()
  .write_atomic_qs(x, LOOKUP_QS)
  .log("[LOOKUP] wrote: %s", LOOKUP_QS)
  x
}

setkey(LOOKUP$wfo_to_counties, wfo)
setkey(LOOKUP$state_to_counties, STATEFP)

# ----------------------------
# Read + expand one tag to UNIQUE county-days
# Returns: (event_tag, GEOID, date, observed=1)
# ----------------------------
.read_tag_to_county_days <- function(tag) {
  tag0 <- event_tag_norm(tag)
  f <- file.path(NOAA_EVENT_SUMMARY_DIR, paste0(tag0, ".csv"))
  if (!file_exists(f)) {
    warning("Missing event summary CSV: ", f)
    return(data.table(event_tag=character(), GEOID=character(), date=as.Date(character()), observed=integer()))
  }
  
  dt <- fread(f, showProgress = FALSE)
  dt[, event_tag := tag0]
  
  # ---- dates -> dt_days with explicit date column ----
  if ("date" %in% names(dt) && any(nzchar(as.character(dt$date)))) {
    dt[, date := as.Date(as.character(date))]
    dt <- dt[!is.na(date)]
    dt <- dt[date >= EVAL_START & date <= EVAL_END]
    dt_days <- dt[, .(
      event_tag,
      date,
      county_geoid = if ("county_geoid" %in% names(dt)) county_geoid else NA_character_,
      wfo          = if ("wfo" %in% names(dt)) wfo else if ("WFO" %in% names(dt)) get("WFO") else if ("cwa" %in% names(dt)) cwa else if ("CWA" %in% names(dt)) get("CWA") else NA_character_,
      statefp      = if ("statefp" %in% names(dt)) statefp else if ("STATEFP" %in% names(dt)) get("STATEFP") else if ("state_fips" %in% names(dt)) state_fips else NA_character_
    )]
  } else {
    if (!all(c("date_start","date_end") %in% names(dt))) {
      stop("Tag ", tag0, " missing both `date` and (`date_start`,`date_end`).")
    }
    dt[, date_start := as.Date(date_start)]
    dt[, date_end   := as.Date(date_end)]
    dt <- dt[!is.na(date_start) & !is.na(date_end)]
    
    dt[, `:=`(lo = pmin(date_start, date_end), hi = pmax(date_start, date_end))]
    dt[, `:=`(date_start = lo, date_end = hi)][, c("lo","hi") := NULL]
    
    dt <- dt[date_end >= EVAL_START & date_start <= EVAL_END]
    if (nrow(dt)) {
      dt[, date_start := pmax(date_start, EVAL_START)]
      dt[, date_end   := pmin(date_end,   EVAL_END)]
    }
    
    dt[, rid := .I]
    dt_days <- dt[, .(
      date = unlist(Map(seq.Date, date_start, date_end, MoreArgs=list(by="day")))
    ), by = .(
      rid,
      event_tag,
      county_geoid = if ("county_geoid" %in% names(dt)) county_geoid else NA_character_,
      wfo          = if ("wfo" %in% names(dt)) wfo else if ("WFO" %in% names(dt)) get("WFO") else if ("cwa" %in% names(dt)) cwa else if ("CWA" %in% names(dt)) get("CWA") else NA_character_,
      statefp      = if ("statefp" %in% names(dt)) statefp else if ("STATEFP" %in% names(dt)) get("STATEFP") else if ("state_fips" %in% names(dt)) state_fips else NA_character_
    )]
  }
  
  # Normalize location fields ONCE (length-preserving)
  dt_days[, county_geoid_norm := .norm_geoid5(county_geoid)]
  dt_days[, wfo_norm         := .norm_wfo(wfo)]
  dt_days[, statefp_norm     := .norm_statefp2(statefp)]
  
  # Validity flags: invalid county GEOID counts as missing
  dt_days[, has_county := !is.na(county_geoid_norm)]
  dt_days[, has_wfo    := !is.na(wfo_norm)]
  dt_days[, has_state  := !is.na(statefp_norm)]
  
  out_list <- list()
  
  # 1) direct county
  x <- dt_days[has_county == TRUE, .(event_tag, GEOID = county_geoid_norm, date)]
  x <- x[!is.na(GEOID) & !is.na(date)]
  if (nrow(x)) out_list[["county"]] <- x
  
  # 2) WFO -> counties (fallback only when county missing/invalid)
  x0 <- dt_days[has_county == FALSE & has_wfo == TRUE, .(event_tag, wfo = wfo_norm, date)]
  x0 <- x0[!is.na(wfo) & !is.na(date)]
  if (nrow(x0)) {
    x <- merge(x0, LOOKUP$wfo_to_counties, by="wfo", all.x=TRUE, allow.cartesian=TRUE)
    x <- x[!is.na(GEOID), .(event_tag, GEOID, date)]
    if (nrow(x)) out_list[["wfo"]] <- x
  }
  
  # 3) State -> counties (fallback only when county and wfo missing/invalid)
  x0 <- dt_days[has_county == FALSE & has_wfo == FALSE & has_state == TRUE,
                .(event_tag, STATEFP = statefp_norm, date)]
  x0 <- x0[!is.na(STATEFP) & !is.na(date)]
  if (nrow(x0)) {
    x <- merge(x0, LOOKUP$state_to_counties, by="STATEFP", all.x=TRUE, allow.cartesian=TRUE)
    x <- x[!is.na(GEOID), .(event_tag, GEOID, date)]
    if (nrow(x)) out_list[["state"]] <- x
  }
  
  # Diagnostics for rows that could not be resolved
  n_in   <- nrow(dt_days)
  n_dead <- dt_days[has_county==FALSE & has_wfo==FALSE & has_state==FALSE, .N]
  .log("[DEF:%s] dt_days=%d | dead(no county/wfo/state)=%d", tag0, n_in, n_dead)
  
  if (WRITE_FAILED_ROWS && n_dead > 0) {
    bad <- dt_days[has_county==FALSE & has_wfo==FALSE & has_state==FALSE]
    fwrite(bad, file.path(OUT_DIR, sprintf("noaa_county_days_%s__FAILED.csv", tag0)))
  }
  
  out <- rbindlist(out_list, use.names=TRUE, fill=TRUE)
  if (!nrow(out)) {
    return(data.table(event_tag=character(), GEOID=character(), date=as.Date(character()), observed=integer()))
  }
  
  out <- unique(out, by = c("event_tag","GEOID","date"))
  out[, observed := 1L]
  setorder(out, GEOID, date)
  out[]
}

# ----------------------------
# Build + save per-tag definitive files
# ----------------------------
.log("[TAGS] heat:   %s", paste(HEAT_TAGS, collapse=", "))
.log("[TAGS] precip: %s", paste(PRECIP_TAGS, collapse=", "))
.log("[TAGS] total unique=%d", length(TAGS_ALL))

for (tag in TAGS_ALL) {
  tag0 <- event_tag_norm(tag)
  out_qs  <- file.path(OUT_DIR, sprintf("noaa_county_days_%s.qs", tag0))
  out_csv <- file.path(OUT_DIR, sprintf("noaa_county_days_%s.csv", tag0))
  
  if (file_exists(out_qs) && !OVERWRITE) {
    .log("[DEF] exists (skip): %s", path_file(out_qs))
    next
  }
  
  .log("[DEF] building tag=%s", tag0)
  x <- .read_tag_to_county_days(tag0)
  
  .write_atomic_qs(x, out_qs)
  .write_atomic_csv(x, out_csv)
  
  if (nrow(x)) {
    mx <- x[, .N, by = GEOID][, max(N)]
    .log("[DEF] tag=%s | rows=%d | unique_counties=%d | max_county_days=%d",
         tag0, nrow(x), uniqueN(x$GEOID), mx)
  } else {
    .log("[DEF] tag=%s | EMPTY", tag0)
  }
}

# ----------------------------
# Group-level outputs (dedup across tags within group at county-day)
# ----------------------------
.write_group_union <- function(tags, group_name) {
  tags <- event_tag_norm(tags)
  paths <- file.path(OUT_DIR, sprintf("noaa_county_days_%s.qs", tags))
  paths <- paths[file_exists(paths)]
  if (!length(paths)) {
    warning("[DEF] no per-tag qs files found for group: ", group_name)
    return(invisible(NULL))
  }
  
  x <- rbindlist(lapply(paths, qread), use.names=TRUE, fill=TRUE)
  u <- unique(x[, .(GEOID, date)], by=c("GEOID","date"))
  u[, `:=`(event_tag = event_tag_norm(group_name), observed = 1L)]
  out_qs  <- file.path(OUT_DIR, sprintf("noaa_county_days__%s.qs", event_tag_norm(group_name)))
  out_csv <- file.path(OUT_DIR, sprintf("noaa_county_days__%s.csv", event_tag_norm(group_name)))
  
  .write_atomic_qs(u, out_qs)
  .write_atomic_csv(u, out_csv)
  .log("[DEF] wrote group union: %s (rows=%d)", path_file(out_qs), nrow(u))
}

.write_group_union(HEAT_TAGS, "all_heat")
.write_group_union(PRECIP_TAGS, "all_precip")

# ----------------------------
# Authoritative lookup table (tag + group keys)
# ----------------------------
.build_lookup_table <- function() {
  tag_paths <- file.path(OUT_DIR, sprintf("noaa_county_days_%s.qs", TAGS_ALL))
  tag_paths <- tag_paths[file_exists(tag_paths)]
  
  tag_dt <- rbindlist(lapply(tag_paths, qread), use.names=TRUE, fill=TRUE)
  tag_dt[, event_tag := event_tag_norm(event_tag)]
  tag_dt[, `:=`(event_group = NA_character_, key = event_tag)]
  tag_dt <- tag_dt[, .(key, event_tag, event_group, GEOID, date, observed)]
  
  grp_rows <- rbindlist(lapply(names(EVENT_GROUPS), function(g) {
    tags <- event_tag_norm(EVENT_GROUPS[[g]])
    # union + dedup across tags at county-day
    y <- rbindlist(lapply(tags, function(tg) {
      p <- file.path(OUT_DIR, sprintf("noaa_county_days_%s.qs", tg))
      if (file_exists(p)) qread(p) else data.table()
    }), use.names=TRUE, fill=TRUE)
    if (!nrow(y)) return(data.table())
    u <- unique(as.data.table(y)[, .(GEOID, date)], by=c("GEOID","date"))
    u[, `:=`(
      key = event_tag_norm(g),
      event_tag = NA_character_,
      event_group = event_tag_norm(g),
      observed = 1L
    )]
    u[, .(key, event_tag, event_group, GEOID, date, observed)]
  }), use.names=TRUE, fill=TRUE)
  
  out <- rbindlist(list(tag_dt, grp_rows), use.names=TRUE, fill=TRUE)
  out <- unique(out, by=c("key","GEOID","date","event_tag","event_group"))
  setorder(out, key, GEOID, date)
  out[]
}

lookup <- .build_lookup_table()
LOOKUP_OUT_QS  <- file.path(OUT_DIR, "noaa_county_day_lookup.qs")
LOOKUP_OUT_CSV <- file.path(OUT_DIR, "noaa_county_day_lookup.csv")
.write_atomic_qs(lookup, LOOKUP_OUT_QS)
.write_atomic_csv(lookup, LOOKUP_OUT_CSV)
.log("[DEF] wrote lookup: %s (rows=%d)", path_file(LOOKUP_OUT_QS), nrow(lookup))

.log("✅ definitive county-day build complete: %s", OUT_DIR)