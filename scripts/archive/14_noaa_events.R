# 1. Extract events from GZ files ----
# =============================================================================
# Read NOAA Storm Events "details" yearly CSV.GZ files into R
#
# Input:
#   V:/users/hquintal/phd1_cluster_southeast/data/input/noaa_events/
#
# Output:
#   V:/users/hquintal/phd1_cluster_southeast_v2/data/out/08_noaa/
#       stormevents_details_all.qs
#       stormevents_details_all.csv.gz  (optional)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(stringr)
  library(qs)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
ROOT_IN  <- "V:/users/hquintal/phd1_cluster_southeast"
ROOT_OUT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

IN_DIR  <- file.path(ROOT_IN,  "data/input/noaa_events")
OUT_DIR <- file.path(ROOT_OUT, "data/out/08_noaa")
dir_create(OUT_DIR)

# Optional year filter (based on filename year after _dYYYY_)
YEAR_MIN <- NA_integer_   # e.g. 1980L
YEAR_MAX <- NA_integer_   # e.g. 2023L

# Write outputs?
WRITE_QS     <- TRUE
WRITE_CSV_GZ <- FALSE

OUT_QS     <- file.path(OUT_DIR, "stormevents_details_all.qs")
OUT_CSV_GZ <- file.path(OUT_DIR, "stormevents_details_all.csv.gz")

# -----------------------------
# Logging helper
# -----------------------------
.log <- function(fmt, ...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")))
  cat(sprintf(fmt, ...))
  cat("\n")
}

# -----------------------------
# Extract year from filename
# -----------------------------
extract_year_from_filename <- function(path) {
  # Example:
  # StormEvents_details-ftp_v1.0_d1950_c20210803.csv.gz
  m <- str_match(basename(path), "_d(\\d{4})_c\\d+\\.csv\\.gz$")
  if (is.na(m[1,2])) return(NA_integer_)
  as.integer(m[1,2])
}

# -----------------------------
# Reader
# -----------------------------
read_one_details_gz <- function(path) {
  dt <- fread(
    path,
    showProgress = FALSE,
    fill = TRUE,
    quote = "",
    colClasses = "character"
  )
  setnames(dt, tolower(names(dt)))
  dt[, source_file := basename(path)]
  dt
}

# -----------------------------
# Discover files
# -----------------------------
files <- dir_ls(
  IN_DIR,
  regexp = "StormEvents_details-ftp_v1\\.0_.*\\.csv\\.gz$",
  type = "file"
)

if (length(files) == 0) {
  stop("No StormEvents_details-*.csv.gz files found in: ", IN_DIR)
}

years <- vapply(files, extract_year_from_filename, integer(1))
ord <- order(years, files)
files <- files[ord]
years <- years[ord]

.log("Found %d file(s) spanning %s–%s.",
     length(files),
     suppressWarnings(min(years, na.rm = TRUE)),
     suppressWarnings(max(years, na.rm = TRUE)))

# Optional year filter
if (!is.na(YEAR_MIN)) {
  keep <- is.na(years) | years >= YEAR_MIN
  files <- files[keep]; years <- years[keep]
}
if (!is.na(YEAR_MAX)) {
  keep <- is.na(years) | years <= YEAR_MAX
  files <- files[keep]; years <- years[keep]
}

.log("After filtering: %d file(s).", length(files))

# -----------------------------
# Read + Combine
# -----------------------------
dt_list <- vector("list", length(files))

for (i in seq_along(files)) {
  .log("Reading [%d/%d]: %s", i, length(files), basename(files[i]))
  dt_list[[i]] <- read_one_details_gz(files[i])
}

dt_all <- rbindlist(dt_list, fill = TRUE, use.names = TRUE)

.log("Combined rows: %s | columns: %d",
     format(nrow(dt_all), big.mark = ","),
     ncol(dt_all))

# -----------------------------
# Optional type coercion (safe guards)
# -----------------------------
if ("begin_yearmonth" %in% names(dt_all))
  dt_all[, begin_yearmonth := as.integer(begin_yearmonth)]

if ("cz_fips" %in% names(dt_all))
  dt_all[, cz_fips := as.integer(cz_fips)]

if ("injuries_direct" %in% names(dt_all))
  dt_all[, injuries_direct := as.integer(injuries_direct)]

if ("deaths_direct" %in% names(dt_all))
  dt_all[, deaths_direct := as.integer(deaths_direct)]

# -----------------------------
# Write outputs
# -----------------------------
if (WRITE_QS) {
  .log("Writing QS → %s", OUT_QS)
  qsave(dt_all, OUT_QS, preset = "high")
}

if (WRITE_CSV_GZ) {
  .log("Writing CSV.GZ → %s", OUT_CSV_GZ)
  fwrite(
    dt_all,
    cmd = paste("gzip -c >", shQuote(OUT_CSV_GZ)),
    row.names = FALSE
  )
}

.log("Done.")

# 2. Per-event files ----
# =============================================================================
# Write one CSV per NOAA StormEvents event_type, expanding cz_type == "Z"
# rows to counties using a deterministic county->WFO crosswalk:
#   county point-on-surface ∈ WFO polygon
#
# Updated:
#   - Robust multi-pass datetime parsing + future-year (-100y) correction
#   - Date window clamping + duration cap + optional daily expansion
#   - Robust county assignment with tiered fallbacks
#
# County assignment rules:
#   cz_type == "C" (County/Parish):
#     1) county GEOID from state_fips + cz_fips (validated against county layer)
#     2) cz_name + state_fips text match to county NAME in that state
#     3) fallback to WFO expansion (state_scoped if state_fips known; else wfo_scoped)
#     4) fallback to all counties in state
#
#   cz_type == "Z" (Forecast Zone):
#     - Expand using WFO crosswalk (state_scoped or wfo_scoped)
#     - If missing wfo OR expansion yields 0, fallback to all counties in state
#
# Outputs:
#   .../data/out/08_noaa/event_summaries/<event_type>.csv
# QA:
#   .../data/out/08_noaa/event_summaries/_qa/*
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(fs)
  library(sf)
  library(lubridate)
  library(stringi)
  library(tigris)
})

options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)

# -----------------------------
# USER SETTINGS
# -----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
IN_QS <- file.path(ROOT, "data/out/08_noaa/stormevents_details_all.qs")

WFO_SHP <- "V:/users/hquintal/phd1_cluster_southeast/data/input/regional_aggregation/nws_wfo_cws/w_05mr24/w_05mr24.shp"

COUNTY_SHP  <- NA_character_  # local county shapefile w/ GEOID+STATEFP+NAME; else uses tigris
TIGRIS_YEAR <- 2023L
TIGRIS_CB   <- TRUE

OUTDIR <- file.path(ROOT, "data/out/08_noaa/event_summaries")
dir_create(OUTDIR)
QA_DIR <- file.path(OUTDIR, "_qa")
dir_create(QA_DIR)

JOIN_EPSG <- 5070L  # NAD83 / CONUS Albers

# Z expansion mode
EXPAND_Z_MODE <- "wfo_scoped"  # "state_scoped" or "wfo_scoped"
stopifnot(EXPAND_Z_MODE %in% c("state_scoped", "wfo_scoped"))

# Filter out rows whose state_fips not present in county layer (drops territories, ocean basins, etc.)
DROP_STATEFIPS_NOT_IN_COUNTIES <- TRUE

# Evaluation window (optional). Events are clamped to overlap this window.
EVAL_START <- as.IDate(NA_character_)  # e.g., as.IDate("1996-01-01")
EVAL_END   <- as.IDate(NA_character_)  # e.g., as.IDate("2023-12-31")

# Cap event durations after clamping (avoid multi-month expansions)
MAX_EVENT_DAYS <- 31L

# Expand to one row per day (recommended for daily evaluation workflows)
EXPAND_TO_DAILY <- TRUE

# QA toggles
QA_CHECK_NAME_MATCH <- TRUE
QA_HARD_STOP_ON_MISSING_GEOID <- FALSE

# Future-year correction ceiling
FUTURE_YEAR_MAX <- 2026L  # or lubridate::year(Sys.Date())

# -----------------------------
# Helpers
# -----------------------------
.log <- function(fmt, ...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")))
  cat(sprintf(fmt, ...))
  cat("\n")
}

safe_filename <- function(x) {
  x <- tolower(x)
  x <- gsub('"', "", x, fixed = TRUE)
  x <- stri_trans_general(x, "Latin-ASCII")
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

strip_quotes <- function(x) {
  x <- as.character(x)
  x <- trimws(gsub('"', "", x, fixed = TRUE))
  x[x == ""] <- NA_character_
  x
}

.std_name <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub('"', "", x, fixed = TRUE)
  x <- gsub("[^A-Z0-9 ]", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- gsub("\\bSAINT\\b", "ST", x)
  x <- gsub("\\bSTE\\b", "ST", x)
  x
}

# Multi-pass parse + future-year correction (no cutoff_2000 dependency)
parse_dt <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  
  dt1 <- suppressWarnings(lubridate::parse_date_time(
    x, orders = c("d-b-y H:M:S", "d-b-Y H:M:S"), tz = "UTC"
  ))
  dt2 <- suppressWarnings(lubridate::parse_date_time(
    x, orders = c("m/d/y H:M:S", "m/d/Y H:M:S", "Y-m-d H:M:S", "Y-m-d H:M"), tz = "UTC"
  ))
  
  dt <- dt1
  dt[is.na(dt)] <- dt2[is.na(dt)]
  
  yy <- suppressWarnings(lubridate::year(dt))
  bad <- !is.na(dt) & !is.na(yy) & (yy > FUTURE_YEAR_MAX)
  if (any(bad)) dt[bad] <- dt[bad] %m-% lubridate::years(100)
  
  dt
}

# state text -> STATEFP (2-digit) using tigris fips_codes
.state_to_fips <- function(state_chr) {
  s <- toupper(trimws(gsub('"', "", as.character(state_chr), fixed = TRUE)))
  fc <- tigris::fips_codes
  idx <- match(s, toupper(fc$state))
  out <- fc$state_code[idx]
  out[is.na(out) | out == ""] <- NA_character_
  sprintf("%02d", suppressWarnings(as.integer(out)))
}

clamp_event_dates <- function(DT) {
  DT <- copy(DT)
  
  DT[, `:=`(
    flag_dt_missing_start = is.na(date_start),
    flag_dt_missing_end   = is.na(date_end),
    flag_dt_end_before_start = (!is.na(date_start) & !is.na(date_end) & date_end < date_start)
  )]
  DT[flag_dt_end_before_start == TRUE, date_end := date_start]
  
  DT[, `:=`(date_start2 = date_start, date_end2 = date_end)]
  
  if (!is.na(EVAL_START)) {
    DT[!is.na(date_start2) & date_start2 < EVAL_START, date_start2 := EVAL_START]
    DT[!is.na(date_end2)   & date_end2   < EVAL_START, date_end2   := EVAL_START]
  }
  if (!is.na(EVAL_END)) {
    DT[!is.na(date_end2)   & date_end2   > EVAL_END, date_end2   := EVAL_END]
    DT[!is.na(date_start2) & date_start2 > EVAL_END, date_start2 := EVAL_END]
  }
  
  DT[, flag_no_overlap_after_clamp := (!is.na(date_start2) & !is.na(date_end2) & date_end2 < date_start2)]
  DT[, event_days := as.integer(date_end2 - date_start2) + 1L]
  
  DT[, flag_truncated_max_days := FALSE]
  DT[!is.na(event_days) & event_days > MAX_EVENT_DAYS, `:=`(
    date_end2 = date_start2 + (MAX_EVENT_DAYS - 1L),
    event_days = MAX_EVENT_DAYS,
    flag_truncated_max_days = TRUE
  )]
  
  DT
}

expand_daily <- function(DT, date_start_col = "date_start2", date_end_col = "date_end2") {
  DT <- DT[!is.na(get(date_start_col)) & !is.na(get(date_end_col))]
  if (nrow(DT) == 0) return(DT[0])
  
  out_list <- vector("list", nrow(DT))
  for (i in seq_len(nrow(DT))) {
    s <- DT[[date_start_col]][i]
    e <- DT[[date_end_col]][i]
    if (is.na(s) || is.na(e)) next
    dd <- seq(s, e, by = "1 day")
    tmp <- DT[i]
    tmp <- tmp[rep(1L, length(dd))]
    tmp[, date := dd]
    out_list[[i]] <- tmp
  }
  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}

# =============================================================================
# 1) Build county -> WFO crosswalk via point-on-surface containment
# =============================================================================
.log("EXPAND_Z_MODE = %s", EXPAND_Z_MODE)
.log("Reading WFO polygons: %s", WFO_SHP)
wfo_sf <- st_read(WFO_SHP, quiet = TRUE)
if (!("WFO" %in% names(wfo_sf)) && ("CWA" %in% names(wfo_sf))) wfo_sf$WFO <- wfo_sf$CWA
stopifnot("WFO" %in% names(wfo_sf))

if (!is.na(COUNTY_SHP)) {
  .log("Reading counties: %s", COUNTY_SHP)
  cty_sf <- st_read(COUNTY_SHP, quiet = TRUE)
} else {
  .log("Reading counties from tigris::counties(year=%d, cb=%s)", TIGRIS_YEAR, TIGRIS_CB)
  cty_sf <- tigris::counties(year = TIGRIS_YEAR, cb = TIGRIS_CB, class = "sf")
}
stopifnot(all(c("GEOID", "STATEFP", "NAME") %in% names(cty_sf)))

# County lookup
cty_lu <- data.table(
  county_geoid      = as.character(cty_sf$GEOID),
  county_state_fips = as.character(cty_sf$STATEFP),
  county_name_raw   = as.character(cty_sf$NAME)
)
cty_lu[, county_name_std := .std_name(county_name_raw)]
setkey(cty_lu, county_geoid)
fwrite(cty_lu, file.path(QA_DIR, "county_lookup_geoid_state_name.csv"))

# Project + valid geometries
wfo_sf <- st_make_valid(wfo_sf)
cty_sf <- st_make_valid(cty_sf)
wfo_j <- st_transform(wfo_sf, JOIN_EPSG) |> st_buffer(0)
cty_j <- st_transform(cty_sf, JOIN_EPSG) |> st_buffer(0)
cty_pt <- st_point_on_surface(cty_j)

.log("Assigning each county to exactly one WFO via st_within(county_point, wfo_polygon)")
idx <- st_within(cty_pt, wfo_j, sparse = TRUE)
pairs <- data.table(
  cty_i = rep(seq_along(idx), lengths(idx)),
  wfo_i = unlist(idx, use.names = FALSE)
)

pairs[, wfo := as.character(wfo_j$WFO[wfo_i])]
pairs[, county_geoid := as.character(cty_j$GEOID[cty_i])]
pairs[, county_state_fips := as.character(cty_j$STATEFP[cty_i])]
setorder(pairs, county_geoid, wfo)
pairs <- pairs[, .SD[1], by = county_geoid]

wfo_county_all <- unique(pairs[, .(wfo, county_geoid)])
setkey(wfo_county_all, wfo)

wfo_county_state <- unique(pairs[, .(wfo, state_fips = county_state_fips, county_geoid)])
setkey(wfo_county_state, wfo, state_fips)

valid_statefips <- sort(unique(as.character(cty_j$STATEFP)))

fwrite(wfo_county_all,   file.path(QA_DIR, "wfo_county_crosswalk_all.csv"))
fwrite(wfo_county_state, file.path(QA_DIR, "wfo_state_county_crosswalk.csv"))

# State -> all counties (for fallback)
state_counties <- cty_lu[, .(county_geoid), by = .(state_fips = county_state_fips)]
setkey(state_counties, state_fips)

# Name+state lookup for county matching
cty_name_state <- cty_lu[, .(
  state_fips = county_state_fips,
  county_geoid,
  cz_name_std = county_name_std
)]
setkey(cty_name_state, state_fips, cz_name_std)

# =============================================================================
# 2) Read StormEvents and prep base summary table
# =============================================================================
.log("Reading StormEvents .qs: %s", IN_QS)
dt <- qread(IN_QS)
setDT(dt)
setnames(dt, tolower(names(dt)))

for (cc in c("cz_type", "wfo", "event_type", "state", "cz_name",
             "state_fips", "cz_fips", "begin_date_time", "end_date_time")) {
  if (cc %in% names(dt)) dt[, (cc) := strip_quotes(get(cc))]
}

stopifnot(all(c("begin_date_time", "end_date_time") %in% names(dt)))

# Parse datetimes
dt[, datetime_start := parse_dt(begin_date_time)]
dt[, datetime_end   := parse_dt(end_date_time)]
dt[is.na(datetime_end) & !is.na(datetime_start), datetime_end := datetime_start]

dt[, date_start := as.IDate(datetime_start)]
dt[, date_end   := as.IDate(datetime_end)]
dt[is.na(date_end) & !is.na(date_start), date_end := date_start]

# Normalize state_fips
dt[, state_fips_int := suppressWarnings(as.integer(state_fips))]
dt[, state_fips := fifelse(!is.na(state_fips_int), sprintf("%02d", state_fips_int), NA_character_)]

# Fill missing state_fips from state name if possible
dt[is.na(state_fips) & !is.na(state) & nzchar(state), state_fips := .state_to_fips(state)]

dt[, cz_fips_int := suppressWarnings(as.integer(cz_fips))]

# Future-year QA (after correction this should be empty)
dt[, `:=`(year_start = lubridate::year(datetime_start), year_end = lubridate::year(datetime_end))]
qa_future <- dt[
  (!is.na(year_start) & year_start > FUTURE_YEAR_MAX) |
    (!is.na(year_end)   & year_end   > FUTURE_YEAR_MAX)
]
fwrite(
  qa_future[, .(begin_date_time, end_date_time, year_start, year_end, event_id, episode_id, state, event_type, source_file)],
  file.path(QA_DIR, "qa_future_year_rows_after_parse.csv")
)

# Quick counts directly from qs
fwrite(dt[, .N, by = .(state, event_type)][order(state, -N)],
       file.path(QA_DIR, "event_type_counts_by_state_from_qs.csv"))
fwrite(dt[, .(n_rows = .N, n_event_types = uniqueN(event_type)), by = state][order(-n_rows)],
       file.path(QA_DIR, "event_type_richness_by_state_from_qs.csv"))

# County GEOID for cz_type == "C" (single county when cz_fips>0)
dt[, county_geoid := NA_character_]
dt[cz_type == "C" & !is.na(state_fips) & !is.na(cz_fips_int) & cz_fips_int > 0,
   county_geoid := paste0(state_fips, sprintf("%03d", cz_fips_int))]

keep_cols <- c(
  "event_type", "event_id", "episode_id", "wfo", "cz_type",
  "state", "state_fips", "cz_fips", "county_geoid", "cz_name",
  "date_start", "date_end", "datetime_start", "datetime_end", "source_file"
)
keep_cols <- keep_cols[keep_cols %in% names(dt)]
dt_sum <- dt[, ..keep_cols]

# Optional drop state_fips not in county layer
if (DROP_STATEFIPS_NOT_IN_COUNTIES) {
  before <- nrow(dt_sum)
  dt_sum <- dt_sum[state_fips %in% valid_statefips]
  .log("DROP_STATEFIPS_NOT_IN_COUNTIES=TRUE: kept %s / %s rows after state filter",
       format(nrow(dt_sum), big.mark=","), format(before, big.mark=","))
}

# Date clamp + QA
dt_sum <- clamp_event_dates(dt_sum)
fwrite(dt_sum[, .N, by = cz_type][order(-N)], file.path(QA_DIR, "cz_type_counts.csv"))
fwrite(
  dt_sum[, .(
    n = .N,
    missing_start = sum(flag_dt_missing_start),
    missing_end = sum(flag_dt_missing_end),
    end_before_start = sum(flag_dt_end_before_start),
    truncated_max_days = sum(flag_truncated_max_days),
    no_overlap_after_clamp = sum(flag_no_overlap_after_clamp)
  )],
  file.path(QA_DIR, "qa_date_ranges_summary.csv")
)

# QA: Validate cz_type=="C" county GEOIDs (exist + state match + optional name check)
C_raw <- dt_sum[cz_type == "C" & !is.na(county_geoid) & nzchar(county_geoid)]
C_qa <- merge(
  C_raw,
  cty_lu[, .(county_geoid, county_state_fips, county_name_std)],
  by = "county_geoid",
  all.x = TRUE
)
C_qa[, cz_name_std := .std_name(cz_name)]
C_qa[, flag_missing_geoid := is.na(county_state_fips)]
C_qa[, flag_state_mismatch := !flag_missing_geoid & !is.na(state_fips) & (state_fips != county_state_fips)]
if (QA_CHECK_NAME_MATCH) {
  C_qa[, flag_name_mismatch := !flag_missing_geoid &
         !is.na(cz_name_std) & nzchar(cz_name_std) &
         !is.na(county_name_std) & nzchar(county_name_std) &
         (cz_name_std != county_name_std)
  ]
} else {
  C_qa[, flag_name_mismatch := FALSE]
}
qa_county_summary <- C_qa[, .(
  n_rows = .N,
  n_unique_geoid = uniqueN(county_geoid),
  missing_geoid = sum(flag_missing_geoid),
  state_mismatch = sum(flag_state_mismatch),
  name_mismatch = sum(flag_name_mismatch)
)]
fwrite(qa_county_summary, file.path(QA_DIR, "qa_county_geoid_summary.csv"))
fwrite(C_qa[flag_missing_geoid | flag_state_mismatch | flag_name_mismatch],
       file.path(QA_DIR, "qa_county_geoid_problem_rows.csv"))
.log("QA C-rows validated GEOID: rows=%s | missing_geoid=%s | state_mismatch=%s | name_mismatch=%s",
     format(qa_county_summary$n_rows, big.mark=","),
     format(qa_county_summary$missing_geoid, big.mark=","),
     format(qa_county_summary$state_mismatch, big.mark=","),
     format(qa_county_summary$name_mismatch, big.mark=","))

if (QA_HARD_STOP_ON_MISSING_GEOID && qa_county_summary$missing_geoid > 0) {
  stop("Some cz_type=='C' county GEOIDs do not exist in county layer. See QA outputs.")
}

# =============================================================================
# 3) Write per-event_type CSVs with Z expansion + tiered C fallbacks
# =============================================================================
evs <- sort(unique(dt_sum$event_type))
.log("Unique event_type count: %d", length(evs))

qa_event <- vector("list", length(evs))
names(qa_event) <- evs

for (ev in evs) {
  out_name <- paste0(safe_filename(ev), ".csv")
  out_path <- file.path(OUTDIR, out_name)
  
  x <- dt_sum[event_type == ev]
  x <- x[flag_no_overlap_after_clamp == FALSE]
  
  # -----------------------------
  # C rows: GEOID -> name+state -> WFO -> state
  # -----------------------------
  xC <- x[cz_type == "C"]
  
  # (1) GEOID validated
  xC_geoid <- xC[!is.na(county_geoid) & nzchar(county_geoid)]
  if (nrow(xC_geoid) > 0) {
    xC_geoid <- merge(
      xC_geoid,
      cty_lu[, .(county_geoid)],
      by = "county_geoid",
      all.x = FALSE
    )
    xC_geoid[, county_assign_method := "C_native_geoid_validated"]
  }
  
  needC <- xC
  if (nrow(xC_geoid) > 0) needC <- needC[!(event_id %in% unique(xC_geoid$event_id))]
  
  # (2) name+state match (cz_name within state)
  xC_name <- needC[0]
  if (nrow(needC) > 0) {
    tmp <- copy(needC)
    tmp[, cz_name_std := .std_name(cz_name)]
    tmp <- tmp[!is.na(state_fips) & nzchar(state_fips) &
                 !is.na(cz_name_std) & nzchar(cz_name_std)]
    if (nrow(tmp) > 0) {
      xC_name <- merge(
        tmp,
        cty_name_state,
        by = c("state_fips", "cz_name_std"),
        all.x = FALSE,
        allow.cartesian = TRUE
      )
      
      # --- FIX: reconcile county_geoid after merge suffixing (.x/.y) ---
      if ("county_geoid.y" %in% names(xC_name)) {
        xC_name[, county_geoid := county_geoid.y]
        xC_name[, `:=`(county_geoid.x = NULL, county_geoid.y = NULL)]
      }
      # ---------------------------------------------------------------
      
      setorder(xC_name, event_id, county_geoid)
      xC_name <- xC_name[, .SD[1], by = .(event_id, episode_id, datetime_start, datetime_end, date_start2, date_end2, state_fips, cz_fips, cz_type, wfo, state, cz_name, event_type, source_file)]
      xC_name[, county_assign_method := "C_fallback_name_state"]
    }
  }
  
  if (nrow(xC_name) > 0) needC <- needC[!(event_id %in% unique(xC_name$event_id))]
  if (nrow(xC_geoid) > 0) needC <- needC[!(event_id %in% unique(xC_geoid$event_id))]
  
  # (3) fallback to WFO
  xC_wfo <- needC[0]
  if (nrow(needC) > 0) {
    tmp <- needC[!is.na(wfo) & nzchar(wfo)]
    if (nrow(tmp) > 0) {
      tmp_state <- tmp[!is.na(state_fips) & nzchar(state_fips)]
      tmp_nost  <- tmp[is.na(state_fips) | !nzchar(state_fips)]
      
      xC_wfo1 <- tmp_state[0]
      if (nrow(tmp_state) > 0) {
        xC_wfo1 <- merge(
          tmp_state,
          wfo_county_state,
          by.x = c("wfo", "state_fips"),
          by.y = c("wfo", "state_fips"),
          all.x = FALSE,
          allow.cartesian = TRUE
        )
        xC_wfo1[, county_assign_method := "C_fallback_wfo_state_scoped"]
      }
      
      xC_wfo2 <- tmp_nost[0]
      if (nrow(tmp_nost) > 0) {
        xC_wfo2 <- merge(
          tmp_nost,
          wfo_county_all,
          by = "wfo",
          all.x = FALSE,
          allow.cartesian = TRUE
        )
        xC_wfo2[, county_assign_method := "C_fallback_wfo_wfo_scoped"]
      }
      
      xC_wfo <- rbindlist(list(xC_wfo1, xC_wfo2), use.names = TRUE, fill = TRUE)
      if (nrow(xC_wfo) > 0 && "county_geoid.y" %in% names(xC_wfo)) {
        xC_wfo[, county_geoid := county_geoid.y]
        xC_wfo[, `:=`(county_geoid.x = NULL, county_geoid.y = NULL)]
      }
    }
  }
  
  if (nrow(xC_wfo) > 0) needC <- needC[!(event_id %in% unique(xC_wfo$event_id))]
  
  # (4) fallback to all counties in state
  xC_state <- needC[0]
  if (nrow(needC) > 0) {
    tmp <- needC[!is.na(state_fips) & nzchar(state_fips)]
    if (nrow(tmp) > 0) {
      xC_state <- merge(tmp, state_counties, by = "state_fips", allow.cartesian = TRUE)
      xC_state[, county_geoid := county_geoid.y]
      xC_state[, `:=`(county_geoid.x = NULL, county_geoid.y = NULL)]
      xC_state[, county_assign_method := "C_fallback_all_state_counties"]
    }
  }
  
  # -----------------------------
  # Z rows: expand + fallback to state if needed
  # -----------------------------
  xZ <- x[cz_type == "Z"]
  
  xZ_wfo <- xZ[!is.na(wfo) & nzchar(wfo)]
  xZ_exp <- xZ_wfo[0]
  if (nrow(xZ_wfo) > 0) {
    if (EXPAND_Z_MODE == "state_scoped") {
      tmp <- xZ_wfo[!is.na(state_fips) & nzchar(state_fips)]
      if (nrow(tmp) > 0) {
        xZ_exp <- merge(
          tmp,
          wfo_county_state,
          by.x = c("wfo", "state_fips"),
          by.y = c("wfo", "state_fips"),
          all.x = FALSE,
          allow.cartesian = TRUE
        )
        xZ_exp[, county_assign_method := "Z_wfo_expand_state_scoped"]
      }
    } else {
      xZ_exp <- merge(
        xZ_wfo,
        wfo_county_all,
        by = "wfo",
        all.x = FALSE,
        allow.cartesian = TRUE
      )
      xZ_exp[, county_assign_method := "Z_wfo_expand_wfo_scoped"]
    }
    if (nrow(xZ_exp) > 0 && "county_geoid.y" %in% names(xZ_exp)) {
      xZ_exp[, county_geoid := county_geoid.y]
      xZ_exp[, `:=`(county_geoid.x = NULL, county_geoid.y = NULL)]
    }
  }
  
  # Z fallback: missing wfo OR wfo present but produced no expansions
  xZ_fb <- xZ[0]
  if (nrow(xZ) > 0) {
    have_exp_ids <- if (nrow(xZ_exp) > 0) unique(xZ_exp$event_id) else character()
    needZ <- xZ[ (is.na(wfo) | !nzchar(wfo)) | !(event_id %in% have_exp_ids) ]
    needZ <- needZ[!is.na(state_fips) & nzchar(state_fips)]
    if (nrow(needZ) > 0) {
      xZ_fb <- merge(needZ, state_counties, by = "state_fips", allow.cartesian = TRUE)
      xZ_fb[, county_geoid := county_geoid.y]
      xZ_fb[, `:=`(county_geoid.x = NULL, county_geoid.y = NULL)]
      xZ_fb[, county_assign_method := "Z_fallback_all_state_counties"]
    }
  }
  
  # -----------------------------
  # Combine outputs
  # -----------------------------
  out <- rbindlist(
    list(xC_geoid, xC_name, xC_wfo, xC_state, xZ_exp, xZ_fb),
    use.names = TRUE, fill = TRUE
  )
  
  # Optional: daily expansion
  if (EXPAND_TO_DAILY) {
    out <- out[!is.na(date_start2) & !is.na(date_end2)]
    out <- expand_daily(out, "date_start2", "date_end2")
  }
  
  # Write per-event file
  keep_out <- intersect(
    c("event_type","event_id","episode_id","cz_type","wfo","state","state_fips","cz_fips","cz_name",
      "county_geoid","county_assign_method",
      "datetime_start","datetime_end","date_start","date_end","date_start2","date_end2","date",
      "event_days","flag_truncated_max_days","source_file"),
    names(out)
  )
  out <- out[, ..keep_out]
  fwrite(out, out_path)
  
  # QA per event_type
  qa_methods <- out[, .N, by = county_assign_method][order(-N)]
  fwrite(qa_methods, file.path(QA_DIR, paste0("qa_methods_", safe_filename(ev), ".csv")))
  
  qa_event[[ev]] <- data.table(
    event_type = ev,
    rows_in = nrow(x),
    rows_out = nrow(out),
    C_geoid = sum(out$county_assign_method == "C_native_geoid_validated", na.rm=TRUE),
    C_name  = sum(out$county_assign_method == "C_fallback_name_state", na.rm=TRUE),
    C_wfo   = sum(grepl("^C_fallback_wfo", out$county_assign_method), na.rm=TRUE),
    C_state = sum(out$county_assign_method == "C_fallback_all_state_counties", na.rm=TRUE),
    Z_exp   = sum(grepl("^Z_wfo_expand", out$county_assign_method), na.rm=TRUE),
    Z_state = sum(out$county_assign_method == "Z_fallback_all_state_counties", na.rm=TRUE),
    truncated = sum(out$flag_truncated_max_days %in% TRUE, na.rm=TRUE)
  )
  
  .log("Wrote %s | out_rows=%s | C_geoid=%s | C_name=%s | C_wfo=%s | C_state=%s | Z_exp=%s | Z_state=%s",
       out_name,
       format(nrow(out), big.mark=","),
       format(qa_event[[ev]]$C_geoid, big.mark=","),
       format(qa_event[[ev]]$C_name, big.mark=","),
       format(qa_event[[ev]]$C_wfo, big.mark=","),
       format(qa_event[[ev]]$C_state, big.mark=","),
       format(qa_event[[ev]]$Z_exp, big.mark=","),
       format(qa_event[[ev]]$Z_state, big.mark=","))
}

qa_event_dt <- rbindlist(qa_event, use.names = TRUE, fill = TRUE)
fwrite(qa_event_dt[order(-rows_out)], file.path(QA_DIR, "qa_per_event_type_output_counts.csv"))

.log("Done. Output folder: %s", OUTDIR)
.log("QA outputs: %s", QA_DIR)

