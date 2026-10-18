# scripts/01_make_heatindex_qs_parallel_validonly.R
# Build daily heat index .qs files from paired dewpoint + temperature .qs inputs
#
# Inputs (per day):
#   .../dewpoint/YYYY-MM-DD.qs      (dewpoint in Celsius; column: intensity)
#   .../temperature/YYYY-MM-DD.qs   (temperature in Celsius; column: intensity)
#
# Output (per day):
#   .../heatindex/YYYY-MM-DD.qs     (heat index in Celsius; NA outside valid regime)
#
# Heat Index method:
#   - Compute RH (%) from T and Td via Magnus approximation
#   - Compute NOAA Rothfusz regression (heat index) in Fahrenheit
#   - Valid only when:
#       * Td <= T  (otherwise NA)
#       * RH in [0, 100]
#       * T >= 80F and RH >= 40% (otherwise NA)
#   - Output converted back to Celsius
#
# Parallelization:
#   - Parallel over days using 18 workers

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(fs)
  library(stringr)
  library(future)
  library(future.apply)
})

# USER INPUTS ----
BASE_DIR <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/01_vars"
DP_DIR   <- file.path(BASE_DIR, "dewpoint")
T_DIR    <- file.path(BASE_DIR, "temperature")
OUT_DIR  <- file.path(BASE_DIR, "heatindex")

N_WORKERS <- 26L
OVERWRITE <- TRUE

# if TRUE, also keep temperature/dewpoint/RH columns in output for QA
KEEP_QA_COLS <- FALSE

# Helpers ----
.msg <- function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))

# RH from T and Td (Celsius) using Magnus approximation
rh_from_t_td_c <- function(t_c, td_c) {
  a <- 17.625
  b <- 243.04
  es <- 6.1094 * exp(a * t_c  / (b + t_c))
  e  <- 6.1094 * exp(a * td_c / (b + td_c))
  rh <- 100 * (e / es)
  rh[rh < 0 | rh > 100] <- NA_real_
  rh
}

# NOAA Rothfusz regression (Fahrenheit). Returns NA outside valid regime.
heat_index_noaa_f_validonly <- function(t_f, rh) {
  hi <- rep(NA_real_, length(t_f))
  
  # valid regime (plus non-missing)
  valid <- !is.na(t_f) & !is.na(rh) & t_f >= 80 & rh >= 40
  if (!any(valid)) return(hi)
  
  T <- t_f[valid]
  R <- rh[valid]
  
  hi_val <- -42.379 +
    2.04901523*T +
    10.14333127*R -
    0.22475541*T*R -
    0.00683783*T^2 -
    0.05481717*R^2 +
    0.00122874*T^2*R +
    0.00085282*T*R^2 -
    0.00000199*T^2*R^2
  
  # Adjustment for very low humidity
  adj1 <- (R < 13 & T >= 80 & T <= 112)
  if (any(adj1)) {
    hi_val[adj1] <- hi_val[adj1] -
      ((13 - R[adj1]) / 4) *
      sqrt((17 - abs(T[adj1] - 95)) / 17)
  }
  
  # Adjustment for very high humidity
  adj2 <- (R > 85 & T >= 80 & T <= 87)
  if (any(adj2)) {
    hi_val[adj2] <- hi_val[adj2] +
      ((R[adj2] - 85) / 10) *
      ((87 - T[adj2]) / 5)
  }
  
  hi[valid] <- hi_val
  hi
}

# Discover paired files ----
dp_files <- dir_ls(DP_DIR, glob = "*.qs")
t_files  <- dir_ls(T_DIR,  glob = "*.qs")

dp_days <- str_replace(path_file(dp_files), "\\.qs$", "")
t_days  <- str_replace(path_file(t_files),  "\\.qs$", "")

days <- sort(intersect(dp_days, t_days))
if (length(days) == 0) stop("No matching YYYY-MM-DD.qs filenames found between dewpoint/ and temperature/.")

dir_create(OUT_DIR)
.msg("Found %d paired days.", length(days))

# Parallel setup ----
plan(multisession, workers = N_WORKERS)
options(future.globals.maxSize = 12 * 1024^3)  # 12GB; adjust if needed

# Per-day worker ----
process_one_day <- function(day) {
  
  out_path <- file.path(OUT_DIR, paste0(day, ".qs"))
  if (!OVERWRITE && file_exists(out_path)) {
    return(list(day = day, status = "skip_exists", n = NA_integer_, na_hi = NA_integer_))
  }
  
  dp_path <- file.path(DP_DIR, paste0(day, ".qs"))
  t_path  <- file.path(T_DIR,  paste0(day, ".qs"))
  
  dp <- as.data.table(qread(dp_path))
  tt <- as.data.table(qread(t_path))
  
  # Expect your schema exactly (from your example):
  # keys: lon, lat, datetime_utc ; value: intensity ; units: "C"
  needed <- c("lon", "lat", "datetime_utc", "intensity")
  miss_dp <- setdiff(needed, names(dp))
  miss_tt <- setdiff(needed, names(tt))
  if (length(miss_dp)) stop("Missing columns in dewpoint: ", paste(miss_dp, collapse = ", "))
  if (length(miss_tt)) stop("Missing columns in temperature: ", paste(miss_tt, collapse = ", "))
  
  setkey(dp, lon, lat, datetime_utc)
  setkey(tt, lon, lat, datetime_utc)
  
  j <- tt[dp, nomatch = 0L]
  if (nrow(j) == 0) {
    return(list(day = day, status = "no_overlap", n = 0L, na_hi = 0L))
  }
  
  # Extract Celsius values
  temp_c <- j$intensity
  dew_c  <- j$i.intensity
  
  # Implausible Td > T => NA
  dew_c[dew_c > temp_c] <- NA_real_
  
  # RH from T and Td (C)
  rh <- rh_from_t_td_c(temp_c, dew_c)
  
  # Convert T to F
  temp_f <- temp_c * 9/5 + 32
  
  # Heat index in F (valid-only)
  hi_f <- heat_index_noaa_f_validonly(temp_f, rh)
  
  # Convert HI to C
  hi_c <- (hi_f - 32) * 5/9
  
  # Output: keep a single heat_index column; allow many NAs
  out <- data.table(
    lon = j$lon,
    lat = j$lat,
    datetime_utc = j$datetime_utc,
    datetime = as.POSIXct(j$datetime_utc, tz = "UTC"),
    heat_index = hi_c,
    units = "C"
  )
  
  if (KEEP_QA_COLS) {
    out[, `:=`(temp_c = temp_c, dewpoint_c = dew_c, rh = rh)]
  }
  
  setkey(out, lon, lat, datetime_utc)
  
  qsave(out, out_path, preset = "high")
  
  list(
    day = day,
    status = "ok",
    n = nrow(out),
    na_hi = sum(is.na(out$heat_index))
  )
}

# Run ----
results <- future_lapply(days, process_one_day, future.seed = TRUE)

res_dt <- rbindlist(lapply(results, as.data.table), fill = TRUE)
.msg("Done. ok=%d | skipped=%d | no_overlap=%d",
     res_dt[status == "ok", .N],
     res_dt[status == "skip_exists", .N],
     res_dt[status == "no_overlap", .N])

# quick NA diagnostics
if (res_dt[status == "ok", .N] > 0) {
  res_ok <- res_dt[status == "ok"]
  res_ok[, na_frac := na_hi / n]
  .msg("NA fraction (ok days): median=%.3f | p95=%.3f | max=%.3f",
       median(res_ok$na_frac, na.rm = TRUE),
       quantile(res_ok$na_frac, 0.95, na.rm = TRUE),
       max(res_ok$na_frac, na.rm = TRUE))
}

print(res_dt[, .N, by = status][order(status)])
