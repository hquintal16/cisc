# suppressPackageStartupMessages({
#   library(fs)
#   library(data.table)
#   library(lubridate)
# })
# 
# ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
# 
# # Canonical mask defining full CONUS extent + CRS
# MASK_TIF <- file.path(ROOT, "data/out/01_vars/precipitation/_mask_target.tif")
# 
# # STM cov model CSVs
# HI_COV_CSV <- file.path(ROOT, "data/out/03_cov_model/heatindex/heatindex_daily_optimal_cov_model.csv")
# PR_COV_CSV <- file.path(ROOT, "data/out/03_cov_model/precipitation/precipitation_daily_optimal_cov_model.csv")
# 
# # Output base
# OUT_BASE <- file.path(ROOT, "data/out/09_performance")
# REF_DIR  <- file.path(OUT_BASE, "ref_grids")
# LOOKUP_DIR <- file.path(OUT_BASE, "lookup_xy_to_cwa")
# NOAA_WFO_DAILY_DIR <- file.path(OUT_BASE, "noaa_wfo_daily")
# MODEL_WFO_DAILY_DIR <- file.path(OUT_BASE, "modeled_wfo_daily")
# CONF_DIR <- file.path(OUT_BASE, "confusion_daily")
# MET_DIR  <- file.path(OUT_BASE, "metrics")
# 
# dir_create(c(REF_DIR, LOOKUP_DIR, NOAA_WFO_DAILY_DIR, MODEL_WFO_DAILY_DIR, CONF_DIR, MET_DIR), recurse = TRUE)
# 
# # CWA polygon shapefile used for lookup (points -> CWA)
# CWA_SHP <- "V:/users/hquintal/phd1_cluster_southeast/data/input/regional_aggregation/nws_wfo_cws/w_05mr24/w_05mr24.shp"
# CWA_FIELD_CANDIDATES <- c("CWA")  # it's there; keep it strict
# 
# # Cluster scenarios
# CLUST_BASE <- file.path(ROOT, "data/out/06_clusters")
# 
# DATASET_TAGS <- c(
#   "precipitation_0.25","precipitation_0.17_bilinear","precipitation_0.17_max","precipitation_0.12_bilinear","precipitation_0.12_max",
#   "heatindex_0.25","heatindex_0.44_bilinear","heatindex_0.44_max","heatindex_0.34_bilinear","heatindex_0.34_max"
# )
# 
# # STM months used to compute resampling resolution
# HI_STM_MONTHS <- 6:8
# PR_STM_MONTHS <- 1:12
# BASE_RES <- 0.25
# 
# # Join rounding (only used if you later merge by x/y; cell_id join is preferred)
# ROUND_DIGITS <- 6
# 
# # Parallelism (use sequential for the modeled-cache stage if you want max stability)
# WORKERS <- 8L
# 
# # Overwrite flags
# OVERWRITE_REF    <- FALSE
# OVERWRITE_LOOKUP <- FALSE
# OVERWRITE_NOAA   <- FALSE
# OVERWRITE_MODEL  <- FALSE
# OVERWRITE_EVAL   <- TRUE
# 
# # NOAA Events CSVs (your existing exports)
# NOAA_EVENTS_DIR <- "V:/users/hquintal/phd1_cluster_southeast/data/output/04_noaa"
# NOAA_FILES <- c(
#   file.path(NOAA_EVENTS_DIR, "NOAA_events.csv"),
#   file.path(NOAA_EVENTS_DIR, "NOAA_events_2009_2015.csv"),
#   file.path(NOAA_EVENTS_DIR, "NOAA_events_2009-2015.csv")
# )
# NOAA_FILES <- NOAA_FILES[file_exists(NOAA_FILES)]
# if (!length(NOAA_FILES)) stop("No NOAA CSV files found in: ", NOAA_EVENTS_DIR)
# 
# NOAA_EVENT_MAP <- list(
#   heatindex = list(
#     excess_heat = c("Excessive Heat"),
#     heat        = c("Heat")
#   ),
#   precipitation = list(
#     flash_flood    = c("Flash Flood"),
#     flood          = c("Flood"),
#     heavy_rain     = c("Heavy Rain"),
#     tropical_storm = c("Hurricane", "Tropical Depression", "Tropical Storm", "Typhoon")
#   )
# )
# 
# .log <- function(fmt, ...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))

# 20260224 ----
suppressPackageStartupMessages({
  library(fs)
  library(data.table)
  library(lubridate)
})

ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

# Canonical mask defining full CONUS extent + CRS
MASK_TIF <- file.path(ROOT, "data/out/01_vars/precipitation/_mask_target.tif")

# STM cov model CSVs
HI_COV_CSV <- file.path(ROOT, "data/out/03_cov_model/heatindex/heatindex_daily_optimal_cov_model.csv")
PR_COV_CSV <- file.path(ROOT, "data/out/03_cov_model/precipitation/precipitation_daily_optimal_cov_model.csv")

# Output base
OUT_BASE <- file.path(ROOT, "data/out/09_performance")

REF_DIR          <- file.path(OUT_BASE, "ref_grids")
LOOKUP_DIR_COUNTY <- file.path(OUT_BASE, "lookup_xy_to_county")
NOAA_COUNTY_DAILY_DIR  <- file.path(OUT_BASE, "noaa_county_daily")
MODEL_COUNTY_DAILY_DIR <- file.path(OUT_BASE, "modeled_county_daily")
CONF_DIR_COUNTY  <- file.path(OUT_BASE, "confusion_daily_county")
MET_DIR_COUNTY   <- file.path(OUT_BASE, "metrics_county")

dir_create(
  c(REF_DIR, LOOKUP_DIR_COUNTY, NOAA_COUNTY_DAILY_DIR, MODEL_COUNTY_DAILY_DIR, CONF_DIR_COUNTY, MET_DIR_COUNTY),
  recurse = TRUE
)

# Cluster scenarios
CLUST_BASE <- file.path(ROOT, "data/out/06_clusters")

DATASET_TAGS <- c(
  "precipitation_0.25","precipitation_0.17_bilinear","precipitation_0.17_max","precipitation_0.12_bilinear","precipitation_0.12_max",
  "heatindex_0.25","heatindex_0.44_bilinear","heatindex_0.44_max","heatindex_0.34_bilinear","heatindex_0.34_max"
)

# STM months used to compute resampling resolution
HI_STM_MONTHS <- 6:8
PR_STM_MONTHS <- 1:12
BASE_RES <- 0.25

# Parallelism defaults (override with options(cisc.workers_model=..., cisc.workers_eval=...))
WORKERS_MODEL <- getOption("cisc.workers_model", 40L)
WORKERS_EVAL  <- getOption("cisc.workers_eval",  40L)

# Overwrite flags
OVERWRITE_REF         <- FALSE
OVERWRITE_LOOKUP      <- FALSE
OVERWRITE_NOAA_COUNTY <- FALSE
OVERWRITE_MODEL_COUNTY<- FALSE
OVERWRITE_EVAL_COUNTY <- TRUE

# Optional evaluation window (set NA to disable)
EVAL_START_DATE <- as.Date("2014-01-01")
EVAL_END_DATE   <- as.Date("2023-12-31")

# NOAA Events CSVs (your existing exports)
NOAA_EVENTS_DIR <- "V:/users/hquintal/phd1_cluster_southeast/data/output/04_noaa"
NOAA_FILES <- c(
  file.path(NOAA_EVENTS_DIR, "NOAA_events.csv"),
  file.path(NOAA_EVENTS_DIR, "NOAA_events_2009_2015.csv"),
  file.path(NOAA_EVENTS_DIR, "NOAA_events_2009-2015.csv")
)
NOAA_FILES <- NOAA_FILES[file_exists(NOAA_FILES)]
if (!length(NOAA_FILES)) stop("No NOAA CSV files found in: ", NOAA_EVENTS_DIR)

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

# Counties: drop AK/HI + territories by STATEFP
DROP_STATEFP <- c("02","15","60","66","69","72","78")

.log <- function(fmt, ...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))