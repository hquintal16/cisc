# scripts/09perf_04_eval_confusion_and_metrics_county_multiwindow.R
# =============================================================================
# County evaluation (restart-safe):
# - Reads observed NOAA definitive county-days from lookup table (Script 1)
# - Reads modeled county daily caches (*.qs)
# - Writes one part file per modeled job then combines (restart-safe)
# - Computes metrics per GEOID × window_tag across yearly + fixed windows
#
# Outputs:
#   data/out/09_performance/metrics_county/metrics_<var_group>_<event>_multiwindow.{qs,csv}
#   data/out/09_performance/metrics_county/parts/part_<var_group>_<event>__<scenario>__<dataset_tag>.qs
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(fs)
  library(lubridate)
  library(future)
  library(future.apply)
  library(progressr)
  library(stringr)
})

source("V:/users/hquintal/phd1_cluster_southeast_v2/scripts/14_eval_config.R")

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
  qsave(obj, tmp, preset=preset)
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
WORKERS_EVAL <- getOption("cisc.workers_eval", 40L)

OVERWRITE_EVAL <- if (exists("OVERWRITE_EVAL")) OVERWRITE_EVAL else FALSE
RUN_SEQUENTIAL_DEBUG <- FALSE

YEAR_START <- 1996L
YEAR_END   <- 2023L

FIXED_WINDOWS <- list(
  list(tag="19960101_20231231", start=as.Date("1996-01-01"), end=as.Date("2023-12-31")),
  list(tag="20040101_20231231", start=as.Date("2004-01-01"), end=as.Date("2023-12-31")),
  list(tag="20140101_20231231", start=as.Date("2014-01-01"), end=as.Date("2023-12-31")),
  list(tag="20190101_20231231", start=as.Date("2019-01-01"), end=as.Date("2023-12-31"))
)

GLOBAL_START <- as.Date(sprintf("%d-01-01", YEAR_START))
GLOBAL_END   <- as.Date(sprintf("%d-12-31", YEAR_END))

# Observed NOAA definitive lookup (from Script 1)
NOAA_COUNTY_DAYS_DEFINITIVE_DIR <- file.path(ROOT, "data/out/09_performance/noaa_county_days_definitive")
NOAA_LOOKUP_QS <- file.path(NOAA_COUNTY_DAYS_DEFINITIVE_DIR, "noaa_county_day_lookup.qs")

# County universe (for true TN and complete coverage)
COUNTY_LUT_QS <- file.path(ROOT, "data/out/09_performance/lookup_county_to_wfo/lookup_county_to_wfo_2022_epsg5070.qs")

if (!exists("MET_DIR")) MET_DIR <- file.path(ROOT, "data/out/09_performance/metrics_county")
dir_create(MET_DIR, recurse = TRUE)

PART_DIR <- file.path(MET_DIR, "parts")
dir_create(PART_DIR, recurse = TRUE)

# ----------------------------
# County universe
# ----------------------------
if (!file_exists(COUNTY_LUT_QS)) stop("Missing COUNTY_LUT_QS: ", COUNTY_LUT_QS)
COUNTY_UNIVERSE <- local({
  x <- qread(COUNTY_LUT_QS)
  sort(unique(as.character(x$county_to_wfo$GEOID)))
})
.log("County universe: %d GEOIDs", length(COUNTY_UNIVERSE))

# ----------------------------
# Discover modeled caches (county)
# ----------------------------
.discover_modeled_jobs <- function() {
  qs_files <- dir_ls(
    MODEL_COUNTY_DAILY_DIR,
    type = "file",
    recurse = TRUE,
    glob = "*.qs",
    fail = FALSE
  )
  qs_files <- qs_files[grepl("^modeled_county_daily_.*\\.qs$", path_file(qs_files))]
  qs_files <- sort(as.character(qs_files))
  if (!length(qs_files)) return(data.table())
  
  rbindlist(lapply(qs_files, function(f) {
    scenario <- path_file(path_dir(f))
    dataset_tag <- sub("^modeled_county_daily_", "", sub("\\.qs$", "", path_file(f)))
    list(
      scenario = scenario,
      dataset_tag = dataset_tag,
      var_group = if (grepl("^heatindex", dataset_tag, ignore.case = TRUE)) "heatindex" else "precipitation",
      mod_path = f
    )
  }), use.names = TRUE, fill = TRUE)
}

# ----------------------------
# Windows table
# ----------------------------
.build_windows <- function() {
  yearly <- rbindlist(lapply(YEAR_START:YEAR_END, function(y) {
    data.table(
      window_tag = sprintf("%d", y),
      window_type = "year",
      year = y,
      start = as.Date(sprintf("%d-01-01", y)),
      end   = as.Date(sprintf("%d-12-31", y))
    )
  }))
  
  fixed <- rbindlist(lapply(FIXED_WINDOWS, function(w) {
    data.table(
      window_tag = w$tag,
      window_type = "fixed",
      year = NA_integer_,
      start = w$start,
      end   = w$end
    )
  }))
  
  windows <- rbindlist(list(yearly, fixed), use.names = TRUE, fill = TRUE)
  windows[, total_days := as.integer(end - start) + 1L]
  windows[]
}
WINDOWS <- .build_windows()

# ----------------------------
# Observed lookup: discover keys + read observed set for a key
# ----------------------------
.discover_noaa_keys_from_lookup <- function(lookup_qs) {
  if (!file_exists(lookup_qs)) stop("Missing NOAA lookup qs: ", lookup_qs)
  lk <- as.data.table(qread(lookup_qs))
  stopifnot(all(c("key","GEOID","date","observed") %in% names(lk)))
  sort(unique(as.character(lk$key)))
}

# FIXED: remove `!!` and avoid `key` name collision with column `key`
.read_obs_from_lookup <- function(key_value, lookup_qs) {
  lk <- as.data.table(qread(lookup_qs))
  stopifnot(all(c("key","GEOID","date","observed") %in% names(lk)))
  
  lk <- lk[key == key_value, .(GEOID, date, observed)]
  
  lk[, GEOID := .norm_geoid5(GEOID)]
  lk[, date := as.Date(date)]
  lk[, observed := as.integer(observed)]
  lk <- lk[!is.na(GEOID) & !is.na(date)]
  lk <- lk[date >= GLOBAL_START & date <= GLOBAL_END]
  lk <- lk[observed > 0L]
  unique(lk[, .(GEOID, date)])
}

.key_to_var_group <- function(k) {
  k <- tolower(as.character(k))
  if (k %in% c("all_heat","heat","excessive_heat")) return("heatindex")
  "precipitation"
}

# ----------------------------
# Readers (modeled sparse)
# ----------------------------
.read_mod_sparse <- function(path) {
  if (!file_exists(path)) return(data.table(GEOID=character(), date=as.Date(character())))
  dt <- as.data.table(qread(path))
  if (!nrow(dt)) return(data.table(GEOID=character(), date=as.Date(character())))
  stopifnot(all(c("GEOID","date","modeled") %in% names(dt)))
  dt[, GEOID := .norm_geoid5(GEOID)]
  dt[, date := as.Date(date)]
  dt[, modeled := as.integer(modeled)]
  dt <- dt[!is.na(GEOID) & !is.na(date)]
  dt <- dt[date >= GLOBAL_START & date <= GLOBAL_END]
  dt <- dt[modeled > 0L]
  unique(dt[, .(GEOID, date)])
}

# ----------------------------
# Metrics
# ----------------------------
.metrics_from_counts <- function(TP, FP, FN, TN) {
  N <- TP + FP + FN + TN
  precision <- if ((TP+FP)>0) TP/(TP+FP) else NA_real_
  recall    <- if ((TP+FN)>0) TP/(TP+FN) else NA_real_
  specificity <- if ((TN+FP)>0) TN/(TN+FP) else NA_real_
  accuracy  <- if (N>0) (TP+TN)/N else NA_real_
  fpr <- if ((FP+TN)>0) FP/(FP+TN) else NA_real_
  fnr <- if ((FN+TP)>0) FN/(FN+TP) else NA_real_
  
  beta <- 2
  f2 <- if (is.finite(precision) && is.finite(recall) && (beta^2*precision + recall) > 0) {
    (1 + beta^2) * precision * recall / (beta^2 * precision + recall)
  } else NA_real_
  f1 <- if (is.finite(precision) && is.finite(recall) && (precision+recall)>0) {
    2*precision*recall/(precision+recall)
  } else NA_real_
  auc_balanced <- if (is.finite(recall) && is.finite(specificity)) (recall + specificity)/2 else NA_real_
  
  list(
    TP=TP,FP=FP,FN=FN,TN=TN,N=N,
    precision=precision, recall=recall,
    specificity=specificity, accuracy=accuracy,
    f1=f1, f2=f2, fpr=fpr, fnr=fnr,
    auc=auc_balanced
  )
}

.compute_counts_all_windows <- function(obs_set, mod_set) {
  
  # full key grid: ALL counties × ALL windows
  keys <- CJ(GEOID = COUNTY_UNIVERSE, window_tag = WINDOWS$window_tag, unique = TRUE)
  setkey(keys, GEOID, window_tag)
  
  # Build date->window mapping once from all relevant dates in global range
  all_dates <- data.table(date = seq.Date(GLOBAL_START, GLOBAL_END, by="day"))
  all_dates[, year := year(date)]
  date_windows <- all_dates[year >= YEAR_START & year <= YEAR_END, .(
    date,
    window_tag_year = as.character(year)
  )]
  for (w in FIXED_WINDOWS) {
    tag <- w$tag
    date_windows[, (tag) := as.integer(date >= w$start & date <= w$end)]
  }
  
  .attach_windows <- function(set_dt) {
    if (!nrow(set_dt)) return(data.table(GEOID=character(), date=as.Date(character()), window_tag=character()))
    x <- merge(set_dt, date_windows, by="date", all.x=TRUE)
    
    out <- list()
    out[[1]] <- x[!is.na(window_tag_year), .(GEOID, date, window_tag = window_tag_year)]
    for (w in FIXED_WINDOWS) {
      tag <- w$tag
      out[[length(out)+1L]] <- x[get(tag) == 1L, .(GEOID, date, window_tag = tag)]
    }
    rbindlist(out, use.names=TRUE, fill=TRUE)
  }
  
  obs_w <- .attach_windows(obs_set)
  mod_w <- .attach_windows(mod_set)
  
  obs_n <- if (nrow(obs_w)) obs_w[, .(obs_n = .N), by=.(GEOID, window_tag)] else data.table(GEOID=character(), window_tag=character(), obs_n=integer())
  mod_n <- if (nrow(mod_w)) mod_w[, .(mod_n = .N), by=.(GEOID, window_tag)] else data.table(GEOID=character(), window_tag=character(), mod_n=integer())
  
  tp_set <- if (nrow(obs_set) && nrow(mod_set)) {
    setkey(obs_set, GEOID, date)
    setkey(mod_set, GEOID, date)
    inter <- obs_set[mod_set, nomatch=0L]
    .attach_windows(inter)
  } else data.table(GEOID=character(), date=as.Date(character()), window_tag=character())
  tp_n <- if (nrow(tp_set)) tp_set[, .(tp_n = .N), by=.(GEOID, window_tag)] else data.table(GEOID=character(), window_tag=character(), tp_n=integer())
  
  setkey(obs_n, GEOID, window_tag)
  setkey(mod_n, GEOID, window_tag)
  setkey(tp_n,  GEOID, window_tag)
  
  out <- copy(keys)
  out <- out[obs_n]
  out <- out[mod_n]
  out <- out[tp_n]
  
  out[is.na(obs_n), obs_n := 0L]
  out[is.na(mod_n), mod_n := 0L]
  out[is.na(tp_n),  tp_n  := 0L]
  
  win_days <- WINDOWS[, .(window_tag, total_days, window_type, year, start, end)]
  out <- merge(out, win_days, by="window_tag", all.x=TRUE)
  
  out[, FP := pmax.int(0L, mod_n - tp_n)]
  out[, FN := pmax.int(0L, obs_n - tp_n)]
  out[, TP := tp_n]
  out[, TN := pmax.int(0L, total_days - TP - FP - FN)]
  out[]
}

# ----------------------------
# One modeled job -> one part file (restart-safe)
# ----------------------------
.part_path <- function(var_group_value, event_name, scenario, dataset_tag) {
  safe_scn <- gsub("[^A-Za-z0-9_\\-]+", "_", scenario)
  safe_dsn <- gsub("[^A-Za-z0-9_\\-]+", "_", dataset_tag)
  file.path(PART_DIR, sprintf("part_%s_%s__%s__%s.qs", var_group_value, event_name, safe_scn, safe_dsn))
}

.eval_one_job_part <- function(obs_set, var_group_value, event_name, job_row) {
  
  out_part <- .part_path(var_group_value, event_name, job_row$scenario, job_row$dataset_tag)
  if (file_exists(out_part) && !OVERWRITE_EVAL) return(out_part)
  
  mod_set <- .read_mod_sparse(job_row$mod_path)
  cnt <- .compute_counts_all_windows(obs_set, mod_set)
  if (!nrow(cnt) || !("window_tag" %in% names(cnt))) return(NA_character_)
  
  met <- cnt[, as.list(.metrics_from_counts(TP, FP, FN, TN)),
             by=.(GEOID, window_tag, window_type, year, start, end)]
  met[, `:=`(
    scenario    = job_row$scenario,
    dataset_tag = job_row$dataset_tag,
    var_group   = var_group_value,
    event       = event_name
  )]
  
  .write_atomic_qs(met, out_part)
  out_part
}

# ----------------------------
# Evaluate one observed set against all modeled jobs
# ----------------------------
.eval_one_event_from_set <- function(var_group_value, event_name, obs_set, jobs_dt) {
  
  out_qs  <- file.path(MET_DIR, sprintf("metrics_%s_%s_multiwindow.qs",  var_group_value, event_name))
  out_csv <- file.path(MET_DIR, sprintf("metrics_%s_%s_multiwindow.csv", var_group_value, event_name))
  
  jobs <- jobs_dt[var_group == var_group_value]
  if (!nrow(jobs)) return(invisible(FALSE))
  
  .log("[EVAL] %s | %s | modeled_jobs=%d | obs_rows=%d", var_group_value, event_name, nrow(jobs), nrow(obs_set))
  
  if (RUN_SEQUENTIAL_DEBUG) {
    part_files <- character(nrow(jobs))
    for (i in seq_len(nrow(jobs))) {
      part_files[i] <- .eval_one_job_part(obs_set, var_group_value, event_name, jobs[i])
    }
  } else {
    part_files <- with_progress({
      p <- progressor(steps = nrow(jobs))
      future_lapply(seq_len(nrow(jobs)), function(i) {
        row <- jobs[i]
        p(sprintf("%s | %s", row$scenario, row$dataset_tag))
        .eval_one_job_part(obs_set, var_group_value, event_name, row)
      }, future.seed = TRUE, future.packages = c("data.table","qs","lubridate","fs","stringr"))
    })
    part_files <- unlist(part_files, use.names = FALSE)
  }
  
  part_files <- part_files[is.character(part_files) & nzchar(part_files) & file_exists(part_files)]
  if (!length(part_files)) return(invisible(FALSE))
  
  parts <- lapply(part_files, function(pth) as.data.table(qread(pth)))
  met_all <- rbindlist(parts, use.names = TRUE, fill = TRUE)
  
  setcolorder(met_all, c("scenario","dataset_tag","var_group","event","window_type","window_tag","year","start","end","GEOID",
                         "TP","FP","FN","TN","N","precision","recall","specificity","accuracy","f1","f2","fpr","fnr","auc"))
  
  if (!OVERWRITE_EVAL && file_exists(out_qs)) {
    .log("[EVAL] final exists (skip write): %s", path_file(out_qs))
    return(invisible(TRUE))
  }
  
  .write_atomic_qs(met_all, out_qs)
  .write_atomic_csv(met_all, out_csv)
  .log("[EVAL] wrote: %s (rows=%d)", path_file(out_csv), nrow(met_all))
  
  invisible(TRUE)
}

# =============================================================================
# MAIN
# =============================================================================
handlers(global = TRUE)
handlers("txtprogressbar")

if (RUN_SEQUENTIAL_DEBUG) {
  plan(sequential)
} else {
  plan(multisession, workers = WORKERS_EVAL)
}

jobs_dt <- .discover_modeled_jobs()
.log("modeled jobs: %d", nrow(jobs_dt))
if (!nrow(jobs_dt)) stop("No modeled county caches found in: ", MODEL_COUNTY_DAILY_DIR)

if (!file_exists(NOAA_LOOKUP_QS)) stop("Missing NOAA lookup: ", NOAA_LOOKUP_QS)
keys <- .discover_noaa_keys_from_lookup(NOAA_LOOKUP_QS)
.log("NOAA keys found in lookup: %d", length(keys))

for (k in keys) {
  vg  <- .key_to_var_group(k)
  evt <- k
  obs_set <- .read_obs_from_lookup(k, NOAA_LOOKUP_QS)
  if (!nrow(obs_set)) {
    warning("[EVAL] empty observed set for key=", k)
    next
  }
  .eval_one_event_from_set(vg, evt, obs_set, jobs_dt)
}

plan(sequential)
.log("✅ county evaluation complete")