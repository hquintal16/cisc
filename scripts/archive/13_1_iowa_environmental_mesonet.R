# VTEC + SIG ~ Mesonet ----
# scripts/fig_01_wwa_pheno_sig_frequency_all_years.R
# =============================================================================
# WWA (IEM/Mesonet) daily frequency maps + per-location time series
# Output per (PHENOM, SIG): PNG map + QS sparse time series
#
# Naming:
#   wwa_freq_<PHENOM>_<SIG>_<phenom_name>_<sig_name>.png
#   wwa_freq_<PHENOM>_<SIG>_<phenom_name>_<sig_name>.qs
# Example:
#   wwa_freq_FA_W_areal_flood_warning.png / .qs
#
# "time series per location" output:
#   sparse table of (STATE_ZONE, date, present=1) for that (PHENOM,SIG).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(lubridate)
  library(ggplot2)
  library(ragg)
  library(tigris)
  library(fs)
  library(qs)
  library(future)
  library(future.apply)
  library(rvest)
  library(xml2)
})

options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)

# ----------------------------
# USER SETTINGS
# ----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

MESONET_BASE <- file.path(ROOT, "data/in/mesonet")  # contains 2002_all ... 2023_all
ZONE_SHP     <- file.path(ROOT, "data/in/shapefile/nws_public_forecast_zone/z_18mr25.shp")
WFO_SHP      <- "V:/users/hquintal/phd1_cluster_southeast/data/input/regional_aggregation/nws_wfo_cws/w_05mr24/w_05mr24.shp"

OUT_DIR <- file.path(ROOT, "figures/fig_01_wwa")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

IEM_VTEC_COUNTS_URL <- "https://mesonet.agron.iastate.edu/vtec/yearly_counts.php"
IEM_LOOKUP_QS  <- file.path(OUT_DIR, "iem_lookup_phenom_sig.qs")
IEM_LOOKUP_CSV <- file.path(OUT_DIR, "iem_lookup_phenom_sig.csv")
REFRESH_IEM_LOOKUP <- FALSE

# bbox in degrees (lon/lat)
BBOX_LONLAT <- c(xmin = -95, xmax = -75, ymin = 22, ymax = 40)

CRS_ALBERS <- 5070

# plotting
RES_MAIN  <- 200
WIDTH_PX  <- 2600
HEIGHT_PX <- 1700

# which SIG to output (you can add Y, S later)
SIG_KEEP <- c("W")

# performance / caching
USE_PARALLEL    <- TRUE
WORKERS         <- 10  # adjust
OVERWRITE_CACHE <- FALSE

CACHE_QS <- file.path(OUT_DIR, "wwa_cache_ugcday_bbox_2002_2023.qs")
LOOKUP_CSV <- file.path(OUT_DIR, "wwa_lookup_pheno_sig.csv")
LOOKUP_QS  <- file.path(OUT_DIR, "wwa_lookup_pheno_sig.qs")

# ----------------------------
# Helpers
# ----------------------------
.log <- function(fmt, ...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))

.parse_dt <- function(x) {
  # your sample is "YYYY-MM-DD HH:MM"
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  suppressWarnings(lubridate::ymd_hm(x, tz = "UTC"))
}

.expand_days <- function(d0, d1) {
  if (is.na(d0) || is.na(d1)) return(as.Date(character()))
  if (d1 < d0) return(as.Date(character()))
  if (d0 == d1) return(d0)
  seq.Date(d0, d1, by = "day")
}

# Drop 3rd character: ALZ019 -> AL019, OHC027 -> OH027, etc.
.drop_third_char <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x == ""] <- NA_character_
  out <- x
  ok <- !is.na(x) & nchar(x) >= 3
  out[ok] <- paste0(substr(x[ok], 1, 2), substr(x[ok], 4, nchar(x[ok])))
  out
}

.slug <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# ----------------------------
# IEM-driven PHENOM/SIG lookup (cached)
# ----------------------------
.build_iem_lookup <- function(url, out_qs, out_csv, refresh = FALSE) {
  if (!refresh && file.exists(out_qs)) {
    .log("Loading IEM lookup cache: %s", out_qs)
    return(as.data.table(qread(out_qs)))
  }
  
  .log("Fetching IEM lookup table: %s", url)
  # yearly_counts.php includes a big table where first columns are:
  # P (PHENOM), S (SIG), Phenomena (name), Significance (name)
  pg <- tryCatch(read_html(url), error = function(e) NULL)
  if (is.null(pg)) stop("Failed to read IEM yearly_counts page: ", url)
  
  tabs <- html_elements(pg, "table")
  if (!length(tabs)) stop("No tables found on IEM yearly_counts page.")
  
  # Find the table that contains header row starting with "P" and "S"
  tab_dt <- NULL
  for (tb in tabs) {
    x <- tryCatch(html_table(tb, fill = TRUE), error = function(e) NULL)
    if (is.null(x) || !ncol(x)) next
    nms <- names(x)
    if (length(nms) >= 4 && toupper(nms[1]) == "P" && toupper(nms[2]) == "S") {
      tab_dt <- as.data.table(x)
      break
    }
  }
  if (is.null(tab_dt)) stop("Could not locate the IEM VTEC counts table with columns P,S,...")
  
  setnames(tab_dt, names(tab_dt), toupper(names(tab_dt)))
  
  need <- c("P","S","PHENOMENA","SIGNIFICANCE")
  if (!all(need %in% names(tab_dt))) {
    stop("IEM table missing expected columns: ", paste(setdiff(need, names(tab_dt)), collapse=", "))
  }
  
  lk <- tab_dt[, .(
    PHENOM = toupper(trimws(P)),
    SIG    = toupper(trimws(S)),
    PHENOMENA = trimws(PHENOMENA),
    SIGNIFICANCE = trimws(SIGNIFICANCE)
  )]
  
  lk <- lk[nchar(PHENOM) == 2 & nchar(SIG) == 1]
  lk <- unique(lk)
  
  # Standardize text to slugs for filenames
  lk[, phenom_slug := .slug(PHENOMENA)]
  lk[, sig_slug    := .slug(SIGNIFICANCE)]
  
  fwrite(lk, out_csv)
  qsave(lk, out_qs)
  .log("Wrote IEM lookup: %s", out_csv)
  .log("Wrote IEM lookup: %s", out_qs)
  lk
}

iem_lk <- .build_iem_lookup(
  url = IEM_VTEC_COUNTS_URL,
  out_qs = IEM_LOOKUP_QS,
  out_csv = IEM_LOOKUP_CSV,
  refresh = REFRESH_IEM_LOOKUP
)

# Build a PHENOM -> slug map from IEM
iem_phenom_dict <- setNames(iem_lk$phenom_slug, iem_lk$PHENOM)
iem_sig_dict    <- setNames(iem_lk$sig_slug,    iem_lk$SIG)

# ---- VTEC PHENOM dictionary (common set; unknowns fall back to code) ----
# You can extend this list as you encounter new codes.
PHENOM_DICT <- c(
  AF="ashfall",
  AS="air_stagnation",
  AV="avalanche",
  BH="beach_hazards",
  BZ="blizzard",
  CF="coastal_flood",
  DF="debris_flow",
  DS="dust_storm",
  DU="blowing_dust",
  EC="extreme_cold",
  EH="excessive_heat",
  EW="extreme_wind",
  FA="areal_flood",
  FF="flash_flood",
  FG="dense_fog",
  FL="flood",
  FR="frost",
  FW="fire_weather",
  GL="gale",
  HF="hurricane_force_wind",
  HI="inland_hurricane",
  HS="heavy_snow",
  HT="heat",
  HU="hurricane",
  HW="high_wind",
  HY="hydrologic_outlook",
  HZ="hard_freeze",
  IP="sleet",
  IS="ice_storm",
  LE="lake_effect_snow",
  LO="low_water",
  LS="lakeshore_flood",
  MA="marine",
  RB="small_craft_for_rough_bar",
  RP="rip_current",
  SC="small_craft",
  SE="hazardous_seas",
  SI="small_craft_for_winds",
  SM="dense_smoke",
  SR="storm",
  SS="storm_surge",
  SU="high_surf",
  SV="severe_thunderstorm",
  TO="tornado",
  TR="tropical_storm",
  TS="tsunami",
  TY="typhoon",
  UP="heavy_freezing_spray",
  WC="wind_chill",
  WI="wind",
  WS="winter_storm",
  WW="winter_weather",
  ZF="freezing_fog",
  ZR="freezing_rain"
)

SIG_DICT <- c(
  W = "warning",
  A = "advisory",
  Y = "watch",
  S = "statement",
  O = "outlook",
  N = "synopsis"
)

.pheno_name <- function(code) {
  code <- toupper(code)
  out <- unname(PHENOM_DICT[code])
  if (!is.na(out)) return(out)
  
  # IEM fallback (covers TI, BW, BS, BH, etc.)
  out2 <- unname(iem_phenom_dict[code])
  if (!is.na(out2)) return(out2)
  
  paste0("phenom_", .slug(code))
}

.sig_name <- function(sig) {
  sig <- toupper(sig)
  out <- unname(SIG_DICT[sig])
  if (!is.na(out)) return(out)
  
  out2 <- unname(iem_sig_dict[sig])
  if (!is.na(out2)) return(out2)
  
  paste0("sig_", .slug(sig))
}

.theme_map <- function() {
  theme_minimal(base_size = 20) +
    theme(
      panel.grid.major = element_line(linewidth = 0.25),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 20),
      legend.text  = element_text(size = 20),
      legend.key.width  = grid::unit(1.6, "in"),
      legend.key.height = grid::unit(0.18, "in"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

.make_map <- function(zones_sf, states_sf, wfo_sf, out_png) {
  p <- ggplot() +
    geom_sf(data = zones_sf, aes(fill = freq_days), color = NA) +
    geom_sf(data = wfo_sf, fill = NA, color = "grey55", linewidth = 0.25) +
    geom_sf(data = states_sf, fill = NA, color = "black", linewidth = 0.35) +
    scale_fill_gradient(
      low = "grey90",
      high = "blue",
      na.value = "white",
      name = "Frequency"
    ) +
    coord_sf(
      crs   = sf::st_crs(CRS_ALBERS),
      datum = sf::st_crs(4326),  # lon/lat degree labels
      xlim  = c(st_bbox(zones_sf)$xmin, st_bbox(zones_sf)$xmax),
      ylim  = c(st_bbox(zones_sf)$ymin, st_bbox(zones_sf)$ymax),
      expand = FALSE
    ) +
    .theme_map()
  
  ragg::agg_png(out_png, width = WIDTH_PX, height = HEIGHT_PX, res = RES_MAIN)
  print(p)
  dev.off()
  .log("Wrote: %s", out_png)
}

# ----------------------------
# Read + prep geometries (States + WFO) in Albers
# ----------------------------
STATES <- tigris::states(cb = TRUE, year = 2022, class = "sf")
drop_statefp <- c("02","15","72","60","66","69","78")
STATES <- STATES[!(STATES$STATEFP %in% drop_statefp), ]
STATES <- st_make_valid(STATES)
STATES <- st_transform(STATES, CRS_ALBERS)
conus_union <- st_union(st_geometry(STATES))

if (!file.exists(WFO_SHP)) stop("Missing WFO shapefile: ", WFO_SHP)
WFO <- st_read(WFO_SHP, quiet = TRUE)
WFO <- st_make_valid(WFO)
WFO <- st_transform(WFO, CRS_ALBERS)
WFO <- WFO[st_intersects(WFO, conus_union, sparse = FALSE), ]

# Zones in Albers and crop/mask in Albers (avoids lon/lat planar intersection warning)
if (!file.exists(ZONE_SHP)) stop("Missing zone shapefile: ", ZONE_SHP)
zones <- st_read(ZONE_SHP, quiet = TRUE)
zones <- st_make_valid(zones)
if (!("STATE_ZONE" %in% names(zones))) stop("ZONE_SHP missing STATE_ZONE field.")
zones$STATE_ZONE <- toupper(trimws(as.character(zones$STATE_ZONE)))

bbox_ll  <- st_as_sfc(st_bbox(BBOX_LONLAT, crs = st_crs(4326)))
bbox_alb <- st_transform(bbox_ll, CRS_ALBERS)

zones_alb  <- st_transform(zones, CRS_ALBERS)
zones_crop <- st_intersection(zones_alb, bbox_alb)

STATES_crop <- STATES[st_intersects(STATES, bbox_alb, sparse = FALSE), ]
WFO_crop    <- WFO[st_intersects(WFO, bbox_alb, sparse = FALSE), ]

zones_keys <- unique(zones_crop$STATE_ZONE)
.log("Zones after crop/mask: %d", nrow(zones_crop))

# ----------------------------
# Build / load cached ugc-day table within bbox
# ----------------------------
if (file.exists(CACHE_QS) && !OVERWRITE_CACHE) {
  .log("Loading cache: %s", CACHE_QS)
  ugcdays_bbox <- as.data.table(qread(CACHE_QS))
} else {
  # Find all CSVs
  if (!dir_exists(MESONET_BASE)) stop("Missing MESONET_BASE: ", MESONET_BASE)
  csv_files <- dir_ls(MESONET_BASE, recurse = TRUE, type = "file", glob = "*.csv")
  if (!length(csv_files)) stop("No CSV files found under: ", MESONET_BASE)
  .log("Found %d CSV files", length(csv_files))
  
  if (USE_PARALLEL) {
    plan(multisession, workers = WORKERS)
    .log("Parallel enabled: workers=%d", WORKERS)
  } else {
    plan(sequential)
    .log("Parallel disabled.")
  }
  
  # Per-file: read, normalize, expand to ugc-day, transform to STATE_ZONE, filter to bbox keys
  .read_expand_one <- function(f, zones_keys, sig_keep) {
    x <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
    if (is.null(x)) return(NULL)
    
    need <- c("PHENOM", "SIG", "ISSUED", "EXPIRED", "NWS_UGC")
    if (!all(need %in% names(x))) return(NULL)
    
    x <- x[, ..need]
    x[, `:=`(
      PHENOM  = toupper(trimws(PHENOM)),
      SIG     = toupper(trimws(SIG)),
      NWS_UGC = toupper(trimws(NWS_UGC))
    )]
    x <- x[SIG %in% sig_keep]
    x[NWS_UGC == "", NWS_UGC := NA_character_]
    x <- x[!is.na(NWS_UGC)]
    
    x[, issued_dt  := .parse_dt(ISSUED)]
    x[, expired_dt := .parse_dt(EXPIRED)]
    x <- x[!is.na(issued_dt) & !is.na(expired_dt)]
    x <- x[expired_dt >= issued_dt]
    if (!nrow(x)) return(NULL)
    
    x[, `:=`(issued_date = as.Date(issued_dt), expired_date = as.Date(expired_dt))]
    
    day_list <- Map(.expand_days, x$issued_date, x$expired_date)
    out <- data.table(
      PHENOM  = rep(x$PHENOM, lengths(day_list)),
      SIG     = rep(x$SIG, lengths(day_list)),
      NWS_UGC = rep(x$NWS_UGC, lengths(day_list)),
      date    = as.Date(unlist(day_list, use.names = FALSE))
    )
    
    if (!nrow(out)) return(NULL)
    setkey(out, PHENOM, SIG, NWS_UGC, date)
    out <- unique(out)
    
    out[, STATE_ZONE := .drop_third_char(NWS_UGC)]
    out <- out[!is.na(STATE_ZONE)]
    out[, NWS_UGC := NULL]
    
    # Filter to bbox zones immediately (major speed/memory win)
    out <- out[STATE_ZONE %in% zones_keys]
    if (!nrow(out)) return(NULL)
    
    setkey(out, PHENOM, SIG, STATE_ZONE, date)
    unique(out)
  }
  
  ugc_list <- future_lapply(
    csv_files,
    .read_expand_one,
    zones_keys = zones_keys,
    sig_keep = SIG_KEEP
  )
  
  ugc_list <- Filter(Negate(is.null), ugc_list)
  if (!length(ugc_list)) stop("No ugc-day rows produced after processing all CSVs.")
  
  ugcdays_bbox <- rbindlist(ugc_list, use.names = TRUE, fill = TRUE)
  setkey(ugcdays_bbox, PHENOM, SIG, STATE_ZONE, date)
  ugcdays_bbox <- unique(ugcdays_bbox)
  
  .log("ugc-day rows in bbox (unique): %d", nrow(ugcdays_bbox))
  qsave(ugcdays_bbox, CACHE_QS)
  .log("Wrote cache: %s", CACHE_QS)
}

# ----------------------------
# Build lookup table (PHENOM + SIG expansions) + write
# ----------------------------
ph <- sort(unique(ugcdays_bbox$PHENOM))
sg <- sort(unique(ugcdays_bbox$SIG))
lookup <- CJ(PHENOM = ph, SIG = sg, unique = TRUE)
lookup[, phenom_name := vapply(PHENOM, .pheno_name, character(1))]
lookup[, sig_name    := vapply(SIG, .sig_name, character(1))]
lookup[, label := sprintf("%s.%s  (%s %s)", PHENOM, SIG, phenom_name, sig_name)]

fwrite(lookup, LOOKUP_CSV)
qsave(lookup, LOOKUP_QS)
.log("Wrote lookup: %s", LOOKUP_CSV)
.log("Wrote lookup: %s", LOOKUP_QS)

# ----------------------------
# Aggregate frequency + write per (PHENOM, SIG): PNG + QS time series
# ----------------------------
# Frequency: count unique days per zone
counts <- ugcdays_bbox[, .(freq_days = .N), by = .(PHENOM, SIG, STATE_ZONE)]
setkey(counts, PHENOM, SIG, STATE_ZONE)

pairs <- unique(counts[, .(PHENOM, SIG)])
setorder(pairs, PHENOM, SIG)

for (i in seq_len(nrow(pairs))) {
  phc <- pairs$PHENOM[i]
  sig <- pairs$SIG[i]
  
  phenom_name <- .pheno_name(phc)
  sig_name    <- .sig_name(sig)
  
  base <- sprintf("wwa_freq_%s_%s_%s_%s", phc, sig, phenom_name, sig_name)
  out_png <- file.path(OUT_DIR, paste0(base, ".png"))
  out_qs  <- file.path(OUT_DIR, paste0(base, ".qs"))
  
  # ---- time series per location (sparse) ----
  # store (STATE_ZONE, date, present=1) for this phenom/sig
  ts_dt <- ugcdays_bbox[PHENOM == phc & SIG == sig, .(STATE_ZONE, date)]
  ts_dt[, present := 1L]
  setkey(ts_dt, STATE_ZONE, date)
  ts_dt <- unique(ts_dt)
  
  qsave(ts_dt, out_qs)
  .log("Wrote TS: %s  (%d rows)", out_qs, nrow(ts_dt))
  
  # ---- map ----
  cc <- counts[PHENOM == phc & SIG == sig, .(STATE_ZONE, freq_days)]
  zones_join <- merge(zones_crop, cc, by = "STATE_ZONE", all.x = TRUE, sort = FALSE)
  
  .make_map(zones_join, STATES_crop, WFO_crop, out_png)
}

.log("Done.")