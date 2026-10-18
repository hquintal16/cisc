# # ~~~~~ ----
# # Either single or nested temporal ----
# # scripts/03_covariance_models_daily.R
# # Fit optimal covariance models to DAILY experimental covariance (space + time)
# # - UNWEIGHTED SSE (np ignored)
# # - Forces negative covariances to 0 BEFORE fitting
# # - Fills missing lags with 0 covariances on fixed lag grids BEFORE fitting
# # - Time model candidates:
# #     (1) nested temporal: alpha * Gaussian(at1) + (1-alpha) * Exponential(at2)
# #     (2) single Exponential(at)
# #     (3) single Gaussian(at)
# # - Selects optimal model by minimized total SSE (space + time)
# # NOTE: plotting can later ignore the filled zeros; this script uses them ONLY for fitting.
# 
# suppressPackageStartupMessages({
#   library(dplyr)
#   library(stringr)
#   library(tibble)
#   library(readr)
#   library(progressr)
#   library(future)
#   library(future.apply)
#   library(tidyr)
# })
# 
# # USER SETTINGS ----
# WORKERS <- 26
# future::plan(multisession, workers = WORKERS)
# 
# IN_BASE  <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/02_exp_cov"
# OUT_BASE <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/03_cov_model"
# 
# VARS <- c("heatindex", "precipitation")
# DOMAINS_KEEP <- c("CONUS")
# 
# # ---- LAG GRIDS (used for zero-filling missing lags) ----
# base_dist <- 0.25
# rLag_grid <- base_dist * c(0:4, 6, 8, 12)   # 0,0.25,0.5,0.75,1,1.5,2,3
# tLag_grid <- c(0:4, 6, 8, 12)               # 0,1,2,3,4,6,8,12
# 
# # range of plausible values
# # (used for nested time model + space range; single-time uses subsets)
# init3  <- c(at1 = 20, at2 = 10, Ar = 1)
# lower3 <- c(at1 = 1,  at2 = 1,  Ar = 0.25)
# upper3 <- c(at1 = 200, at2 = 200, Ar = 100)
# 
# # single-time model bounds
# init2  <- c(at = 20, Ar = 1)
# lower2 <- c(at = 1,  Ar = 0.25)
# upper2 <- c(at = 200, Ar = 100)
# 
# alpha_grid <- seq(0.01, 0.99, by = 0.01)
# 
# # Helpers ----
# safe_read_cov <- function(path, lag_col_guess = c("rLag", "tLag", "lag")) {
#   x <- tryCatch(read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE),
#                 error = function(e) NULL)
#   if (is.null(x) || nrow(x) == 0) return(NULL)
#   
#   nm <- names(x)
#   
#   lag_col <- intersect(lag_col_guess, nm)
#   lag_col <- if (length(lag_col) == 0) nm[1] else lag_col[1]
#   
#   cov_col <- if ("Cov" %in% nm) "Cov" else if ("cov" %in% nm) "cov" else nm[2]
#   np_col  <- if ("np" %in% nm) "np" else NA_character_
#   
#   out <- tibble(
#     lag = as.numeric(x[[lag_col]]),
#     cov = as.numeric(x[[cov_col]]),
#     np  = if (!is.na(np_col)) as.numeric(x[[np_col]]) else NA_real_
#   ) %>%
#     filter(is.finite(lag), is.finite(cov)) %>%
#     arrange(lag)
#   
#   if (nrow(out) == 0) return(NULL)
#   out
# }
# 
# safe_fn <- function(fn) {
#   function(p) {
#     val <- tryCatch(fn(p), error = function(e) NA_real_)
#     if (!is.finite(val)) val <- .Machine$double.xmax / 2
#     val
#   }
# }
# 
# SSE <- function(obs, pred) {
#   sum((obs - pred)^2, na.rm = TRUE)
# }
# 
# # ---- Clamp negatives to 0 + fill missing lags with 0 on a fixed grid ----
# standardize_cov_curve <- function(df, lag_grid) {
#   # Clamp negative covariances to 0
#   df2 <- df %>%
#     transmute(lag = as.numeric(lag),
#               cov = pmax(0, as.numeric(cov))) %>%
#     filter(is.finite(lag), is.finite(cov)) %>%
#     group_by(lag) %>%
#     summarise(cov = mean(cov), .groups = "drop")  # collapse duplicates
#   
#   out <- tibble(lag = lag_grid) %>%
#     left_join(df2, by = "lag") %>%
#     mutate(cov = replace_na(cov, 0)) %>%
#     arrange(lag)
#   
#   out
# }
# 
# # Covariance model kernels (sill fixed)
# pred_space <- function(lag, sill, Ar, type = c("gaussian", "exponential")) {
#   type <- match.arg(type)
#   if (type == "gaussian") sill * exp(-3 * (lag^2) / (Ar^2)) else sill * exp(-3 * lag / Ar)
# }
# 
# pred_time_single <- function(lag, sill, at, type = c("gaussian", "exponential")) {
#   type <- match.arg(type)
#   if (type == "gaussian") sill * exp(-3 * (lag^2) / (at^2)) else sill * exp(-3 * lag / at)
# }
# 
# pred_time_nested <- function(lag, sill, alpha, at1, at2) {
#   # alpha * Gaussian(at1) + (1-alpha) * Exponential(at2)
#   sill * (alpha * exp(-3 * (lag^2) / (at1^2)) + (1 - alpha) * exp(-3 * lag / at2))
# }
# 
# discover_daily_pairs <- function(var_folder) {
#   cov_dir <- file.path(IN_BASE, var_folder, "daily")
#   if (!dir.exists(cov_dir)) return(tibble())
#   
#   sp_files <- list.files(cov_dir, pattern = "_cov_space\\.txt$", full.names = TRUE)
#   
#   rx <- "^([a-zA-Z0-9]+)_([A-Za-z0-9]+)_day_([0-9]{4}-[0-9]{2}-[0-9]{2})_cov_space\\.txt$"
#   
#   tibble(space_file = sp_files) %>%
#     mutate(
#       fname     = basename(space_file),
#       matched   = str_match(fname, rx),
#       prefix    = matched[, 2],
#       domain    = matched[, 3],
#       day       = matched[, 4],
#       time_file = file.path(cov_dir, str_replace(fname, "_cov_space\\.txt$", "_cov_time.txt"))
#     ) %>%
#     filter(!is.na(day), file.exists(time_file), domain %in% DOMAINS_KEEP)
# }
# 
# fit_one_pair <- function(space_file, time_file, var_name, domain, day) {
#   
#   sp_raw <- safe_read_cov(space_file, lag_col_guess = c("rLag", "lag"))
#   tm_raw <- safe_read_cov(time_file,  lag_col_guess = c("tLag", "lag"))
#   if (is.null(sp_raw) || is.null(tm_raw)) return(NULL)
#   
#   # Standardize curves: negative -> 0, then fill missing lags with 0 on fixed grids
#   space_df <- standardize_cov_curve(sp_raw, rLag_grid)
#   time_df  <- standardize_cov_curve(tm_raw, tLag_grid)
#   
#   # Need a valid sill at lag 0 (after clamping)
#   if (!any(space_df$lag == 0) || !any(time_df$lag == 0)) return(NULL)
#   
#   sill <- space_df$cov[which(space_df$lag == 0)[1]]
#   if (!is.finite(sill) || sill <= 0) return(NULL)
#   
#   # Track how much was original vs filled (for diagnostics)
#   n_space_lags_calc <- nrow(sp_raw)
#   n_time_lags_calc  <- nrow(tm_raw)
#   n_space_lags_used <- nrow(space_df)
#   n_time_lags_used  <- nrow(time_df)
#   
#   # Candidate model families
#   space_models <- c("exponential", "gaussian")
#   time_models  <- c("nested_gauss_exp", "exponential", "gaussian")  # temporal candidates
#   
#   best <- NULL
#   best_sse_total <- Inf
#   
#   for (space_model in space_models) {
#     
#     # Time model: NESTED
#     for (alpha in alpha_grid) {
#       
#       sse_total_fn <- function(p) {
#         at1 <- p[1]; at2 <- p[2]; Ar <- p[3]
#         
#         sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
#         tm_pred <- pred_time_nested(time_df$lag, sill, alpha, at1, at2)
#         
#         SSE(space_df$cov, sp_pred) + SSE(time_df$cov, tm_pred)
#       }
#       
#       opt <- optim(
#         par    = init3,
#         fn     = safe_fn(sse_total_fn),
#         method = "L-BFGS-B",
#         lower  = lower3,
#         upper  = upper3
#       )
#       
#       if (is.finite(opt$value) && all(is.finite(opt$par)) && opt$value < best_sse_total) {
#         at1 <- opt$par[1]; at2 <- opt$par[2]; Ar <- opt$par[3]
#         
#         sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
#         tm_pred <- pred_time_nested(time_df$lag, sill, alpha, at1, at2)
#         
#         sse_space <- SSE(space_df$cov, sp_pred)
#         sse_time  <- SSE(time_df$cov, tm_pred)
#         sse_total <- sse_space + sse_time
#         
#         stm1 <- Ar / (alpha * at1 + (1 - alpha) * at2)
#         stm2 <- alpha * (Ar / at1) + (1 - alpha) * (Ar / at2)
#         
#         best_sse_total <- sse_total
#         best <- tibble(
#           variable          = var_name,
#           domain            = domain,
#           day               = day,
#           sill              = sill,
#           
#           # params (keep consistent columns across model types)
#           alpha1            = alpha,
#           alpha2            = 1 - alpha,
#           at1               = at1,
#           at2               = at2,
#           at                = NA_real_,   # single-time param (unused here)
#           Ar                = Ar,
#           
#           stm1              = stm1,
#           stm2              = stm2,
#           
#           sse_space         = sse_space,
#           sse_time          = sse_time,
#           sse_total         = sse_total,
#           
#           converged         = (opt$convergence == 0),
#           model_space       = space_model,
#           model_time        = "nested_gauss_exp",
#           
#           n_space_lags_calc = n_space_lags_calc,
#           n_time_lags_calc  = n_time_lags_calc,
#           n_space_lags_used = n_space_lags_used,
#           n_time_lags_used  = n_time_lags_used,
#           
#           space_file        = space_file,
#           time_file         = time_file
#         )
#       }
#     }
#     
#     # Time model: SINGLE EXPONENTIAL
#     {
#       sse_total_fn <- function(p) {
#         at <- p[1]; Ar <- p[2]
#         sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
#         tm_pred <- pred_time_single(time_df$lag, sill, at, type = "exponential")
#         SSE(space_df$cov, sp_pred) + SSE(time_df$cov, tm_pred)
#       }
#       
#       opt <- optim(
#         par    = init2,
#         fn     = safe_fn(sse_total_fn),
#         method = "L-BFGS-B",
#         lower  = lower2,
#         upper  = upper2
#       )
#       
#       if (is.finite(opt$value) && all(is.finite(opt$par)) && opt$value < best_sse_total) {
#         at <- opt$par[1]; Ar <- opt$par[2]
#         
#         sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
#         tm_pred <- pred_time_single(time_df$lag, sill, at, type = "exponential")
#         
#         sse_space <- SSE(space_df$cov, sp_pred)
#         sse_time  <- SSE(time_df$cov, tm_pred)
#         sse_total <- sse_space + sse_time
#         
#         stm1 <- Ar / at
#         stm2 <- Ar / at
#         
#         best_sse_total <- sse_total
#         best <- tibble(
#           variable          = var_name,
#           domain            = domain,
#           day               = day,
#           sill              = sill,
#           
#           alpha1            = NA_real_,
#           alpha2            = NA_real_,
#           at1               = NA_real_,
#           at2               = NA_real_,
#           at                = at,
#           Ar                = Ar,
#           
#           stm1              = stm1,
#           stm2              = stm2,
#           
#           sse_space         = sse_space,
#           sse_time          = sse_time,
#           sse_total         = sse_total,
#           
#           converged         = (opt$convergence == 0),
#           model_space       = space_model,
#           model_time        = "exponential",
#           
#           n_space_lags_calc = n_space_lags_calc,
#           n_time_lags_calc  = n_time_lags_calc,
#           n_space_lags_used = n_space_lags_used,
#           n_time_lags_used  = n_time_lags_used,
#           
#           space_file        = space_file,
#           time_file         = time_file
#         )
#       }
#     }
#     
#     # Time model: SINGLE GAUSSIAN
#     {
#       sse_total_fn <- function(p) {
#         at <- p[1]; Ar <- p[2]
#         sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
#         tm_pred <- pred_time_single(time_df$lag, sill, at, type = "gaussian")
#         SSE(space_df$cov, sp_pred) + SSE(time_df$cov, tm_pred)
#       }
#       
#       opt <- optim(
#         par    = init2,
#         fn     = safe_fn(sse_total_fn),
#         method = "L-BFGS-B",
#         lower  = lower2,
#         upper  = upper2
#       )
#       
#       if (is.finite(opt$value) && all(is.finite(opt$par)) && opt$value < best_sse_total) {
#         at <- opt$par[1]; Ar <- opt$par[2]
#         
#         sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
#         tm_pred <- pred_time_single(time_df$lag, sill, at, type = "gaussian")
#         
#         sse_space <- SSE(space_df$cov, sp_pred)
#         sse_time  <- SSE(time_df$cov, tm_pred)
#         sse_total <- sse_space + sse_time
#         
#         stm1 <- Ar / at
#         stm2 <- Ar / at
#         
#         best_sse_total <- sse_total
#         best <- tibble(
#           variable          = var_name,
#           domain            = domain,
#           day               = day,
#           sill              = sill,
#           
#           alpha1            = NA_real_,
#           alpha2            = NA_real_,
#           at1               = NA_real_,
#           at2               = NA_real_,
#           at                = at,
#           Ar                = Ar,
#           
#           stm1              = stm1,
#           stm2              = stm2,
#           
#           sse_space         = sse_space,
#           sse_time          = sse_time,
#           sse_total         = sse_total,
#           
#           converged         = (opt$convergence == 0),
#           model_space       = space_model,
#           model_time        = "gaussian",
#           
#           n_space_lags_calc = n_space_lags_calc,
#           n_time_lags_calc  = n_time_lags_calc,
#           n_space_lags_used = n_space_lags_used,
#           n_time_lags_used  = n_time_lags_used,
#           
#           space_file        = space_file,
#           time_file         = time_file
#         )
#       }
#     }
#   }
#   
#   best
# }
# 
# # Run ----
# handler_txtprogressbar(format = ":current/:total [:bar] :elapsed | eta: :eta")
# 
# for (var_folder in VARS) {
#   
#   pairs <- discover_daily_pairs(var_folder)
#   if (nrow(pairs) == 0) {
#     message("[SKIP] No cov files found for: ", var_folder)
#     next
#   }
#   
#   message("[RUN] ", var_folder, " (n=", nrow(pairs), " days)")
#   
#   results_list <- NULL
#   with_progress({
#     p <- progressor(along = seq_len(nrow(pairs)))
#     results_list <- future_lapply(seq_len(nrow(pairs)), function(i) {
#       p()
#       row <- pairs[i, ]
#       fit_one_pair(
#         space_file = row$space_file,
#         time_file  = row$time_file,
#         var_name   = var_folder,
#         domain     = row$domain,
#         day        = row$day
#       )
#     })
#   })
#   
#   master <- bind_rows(results_list)
#   
#   out_dir <- file.path(OUT_BASE, var_folder)
#   dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#   
#   out_csv <- file.path(out_dir, paste0(var_folder, "_daily_optimal_cov_model.csv"))
#   
#   if (nrow(master) > 0) {
#     master <- master %>% arrange(domain, day, sse_total)
#     write_csv(master, out_csv)
#     message("→ Wrote: ", out_csv)
#   } else {
#     warning("[WARN] No successful fits recorded for: ", var_folder)
#   }
# }

# ~~~~~ ----
# Southeast ----
# scripts/03_covariance_models_daily_SOUTHEAST.R
# Fit optimal covariance models to DAILY experimental covariance (space + time)
# - UNWEIGHTED SSE (np ignored)
# - Forces negative covariances to 0 BEFORE fitting
# - Fills missing lags with 0 covariances on fixed lag grids BEFORE fitting
# - Time model candidates:
#     (1) nested temporal: alpha * Gaussian(at1) + (1-alpha) * Exponential(at2)
#     (2) single Exponential(at)
#     (3) single Gaussian(at)
# - Selects optimal model by minimized total SSE (space + time)
#
# NEW:
# - Fits from Southeast experimental covariance folders:
#     02_exp_cov/heatindex/southeast
#     02_exp_cov/precipitation/southeast
# - Writes:
#     03_cov_model/<var>/<var>_daily_optimal_cov_model_southeast.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(readr)
  library(progressr)
  library(future)
  library(future.apply)
  library(tidyr)
  library(fs)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
WORKERS <- 26
future::plan(multisession, workers = WORKERS)

IN_BASE  <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/02_exp_cov"
OUT_BASE <- "V:/users/hquintal/phd1_cluster_southeast_v2/data/out/03_cov_model"

VARS <- c("heatindex", "precipitation")

# NEW: domain folder(s) under 02_exp_cov/<var>/
DOMAINS_KEEP <- c("southeast")

# ---- LAG GRIDS (used for zero-filling missing lags) ----
base_dist <- 0.25
rLag_grid <- base_dist * c(0:4, 6, 8, 12)   # 0,0.25,0.5,0.75,1,1.5,2,3
tLag_grid <- c(0:4, 6, 8, 12)               # 0,1,2,3,4,6,8,12

# range of plausible values
init3  <- c(at1 = 20, at2 = 10, Ar = 1)
lower3 <- c(at1 = 1,  at2 = 1,  Ar = 0.25)
upper3 <- c(at1 = 200, at2 = 200, Ar = 100)

# single-time model bounds
init2  <- c(at = 20, Ar = 1)
lower2 <- c(at = 1,  Ar = 0.25)
upper2 <- c(at = 200, Ar = 100)

alpha_grid <- seq(0.01, 0.99, by = 0.01)

# -----------------------------
# Helpers
# -----------------------------
safe_read_cov <- function(path, lag_col_guess = c("rLag", "tLag", "lag")) {
  x <- tryCatch(
    read.table(path, header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE),
    error = function(e) NULL
  )
  if (is.null(x) || nrow(x) == 0) return(NULL)
  
  nm <- names(x)
  
  lag_col <- intersect(lag_col_guess, nm)
  lag_col <- if (length(lag_col) == 0) nm[1] else lag_col[1]
  
  cov_col <- if ("Cov" %in% nm) "Cov" else if ("cov" %in% nm) "cov" else nm[2]
  np_col  <- if ("np" %in% nm) "np" else NA_character_
  
  out <- tibble(
    lag = as.numeric(x[[lag_col]]),
    cov = as.numeric(x[[cov_col]]),
    np  = if (!is.na(np_col)) as.numeric(x[[np_col]]) else NA_real_
  ) %>%
    filter(is.finite(lag), is.finite(cov)) %>%
    arrange(lag)
  
  if (nrow(out) == 0) return(NULL)
  out
}

safe_fn <- function(fn) {
  function(p) {
    val <- tryCatch(fn(p), error = function(e) NA_real_)
    if (!is.finite(val)) val <- .Machine$double.xmax / 2
    val
  }
}

SSE <- function(obs, pred) sum((obs - pred)^2, na.rm = TRUE)

# ---- Clamp negatives to 0 + fill missing lags with 0 on a fixed grid ----
standardize_cov_curve <- function(df, lag_grid) {
  df2 <- df %>%
    transmute(lag = as.numeric(lag),
              cov = pmax(0, as.numeric(cov))) %>%
    filter(is.finite(lag), is.finite(cov)) %>%
    group_by(lag) %>%
    summarise(cov = mean(cov), .groups = "drop")
  
  tibble(lag = lag_grid) %>%
    left_join(df2, by = "lag") %>%
    mutate(cov = replace_na(cov, 0)) %>%
    arrange(lag)
}

# Covariance model kernels (sill fixed)
pred_space <- function(lag, sill, Ar, type = c("gaussian", "exponential")) {
  type <- match.arg(type)
  if (type == "gaussian") sill * exp(-3 * (lag^2) / (Ar^2)) else sill * exp(-3 * lag / Ar)
}

pred_time_single <- function(lag, sill, at, type = c("gaussian", "exponential")) {
  type <- match.arg(type)
  if (type == "gaussian") sill * exp(-3 * (lag^2) / (at^2)) else sill * exp(-3 * lag / at)
}

pred_time_nested <- function(lag, sill, alpha, at1, at2) {
  sill * (alpha * exp(-3 * (lag^2) / (at1^2)) + (1 - alpha) * exp(-3 * lag / at2))
}

# -----------------------------
# NEW: discover daily pairs from 02_exp_cov/<var>/<domain_folder>
# -----------------------------
discover_daily_pairs <- function(var_folder) {
  out <- list()
  
  for (domain in DOMAINS_KEEP) {
    cov_dir <- file.path(IN_BASE, var_folder, domain)
    if (!dir.exists(cov_dir)) next
    
    sp_files <- list.files(cov_dir, pattern = "_cov_space\\.txt$", full.names = TRUE)
    
    # Accept both:
    #   heatindex_southeast_day_YYYY-MM-DD_cov_space.txt
    #   precip_southeast_day_YYYY-MM-DD_cov_space.txt
    # and similar, as long as it matches *_<domain>_day_<date>_cov_space.txt
    rx <- "^([A-Za-z0-9]+)_([A-Za-z0-9]+)_day_([0-9]{4}-[0-9]{2}-[0-9]{2})_cov_space\\.txt$"
    
    out[[domain]] <- tibble(space_file = sp_files) %>%
      mutate(
        fname   = basename(space_file),
        matched = str_match(fname, rx),
        prefix  = matched[, 2],
        domain  = matched[, 3],
        day     = matched[, 4],
        time_file = file.path(cov_dir, str_replace(fname, "_cov_space\\.txt$", "_cov_time.txt"))
      ) %>%
      filter(!is.na(day), file.exists(time_file), domain %in% DOMAINS_KEEP)
  }
  
  bind_rows(out)
}

fit_one_pair <- function(space_file, time_file, var_name, domain, day) {
  
  sp_raw <- safe_read_cov(space_file, lag_col_guess = c("rLag", "lag"))
  tm_raw <- safe_read_cov(time_file,  lag_col_guess = c("tLag", "lag"))
  if (is.null(sp_raw) || is.null(tm_raw)) return(NULL)
  
  # Standardize curves: negative -> 0, then fill missing lags with 0 on fixed grids
  space_df <- standardize_cov_curve(sp_raw, rLag_grid)
  time_df  <- standardize_cov_curve(tm_raw, tLag_grid)
  
  # Need a valid sill at lag 0 (after clamping)
  if (!any(space_df$lag == 0) || !any(time_df$lag == 0)) return(NULL)
  
  sill <- space_df$cov[which(space_df$lag == 0)[1]]
  if (!is.finite(sill) || sill <= 0) return(NULL)
  
  # Diagnostics
  n_space_lags_calc <- nrow(sp_raw)
  n_time_lags_calc  <- nrow(tm_raw)
  n_space_lags_used <- nrow(space_df)
  n_time_lags_used  <- nrow(time_df)
  
  space_models <- c("exponential", "gaussian")
  best <- NULL
  best_sse_total <- Inf
  
  for (space_model in space_models) {
    
    # Time model: NESTED
    for (alpha in alpha_grid) {
      
      sse_total_fn <- function(p) {
        at1 <- p[1]; at2 <- p[2]; Ar <- p[3]
        sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
        tm_pred <- pred_time_nested(time_df$lag, sill, alpha, at1, at2)
        SSE(space_df$cov, sp_pred) + SSE(time_df$cov, tm_pred)
      }
      
      opt <- optim(
        par    = init3,
        fn     = safe_fn(sse_total_fn),
        method = "L-BFGS-B",
        lower  = lower3,
        upper  = upper3
      )
      
      if (is.finite(opt$value) && all(is.finite(opt$par)) && opt$value < best_sse_total) {
        at1 <- opt$par[1]; at2 <- opt$par[2]; Ar <- opt$par[3]
        
        sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
        tm_pred <- pred_time_nested(time_df$lag, sill, alpha, at1, at2)
        
        sse_space <- SSE(space_df$cov, sp_pred)
        sse_time  <- SSE(time_df$cov, tm_pred)
        sse_total <- sse_space + sse_time
        
        stm1 <- Ar / (alpha * at1 + (1 - alpha) * at2) # ratio of variance weighted average ranges
        stm2 <- alpha * (Ar / at1) + (1 - alpha) * (Ar / at2)# variance weighted average of the ratios
        
        best_sse_total <- sse_total
        best <- tibble(
          variable          = var_name,
          domain            = domain,
          day               = day,
          sill              = sill,
          
          alpha1            = alpha,
          alpha2            = 1 - alpha,
          at1               = at1,
          at2               = at2,
          at                = NA_real_,
          Ar                = Ar,
          
          stm1              = stm1,
          stm2              = stm2,
          
          sse_space         = sse_space,
          sse_time          = sse_time,
          sse_total         = sse_total,
          
          converged         = (opt$convergence == 0),
          model_space       = space_model,
          model_time        = "nested_gauss_exp",
          
          n_space_lags_calc = n_space_lags_calc,
          n_time_lags_calc  = n_time_lags_calc,
          n_space_lags_used = n_space_lags_used,
          n_time_lags_used  = n_time_lags_used,
          
          space_file        = space_file,
          time_file         = time_file
        )
      }
    }
    
    # Time model: SINGLE EXPONENTIAL
    {
      sse_total_fn <- function(p) {
        at <- p[1]; Ar <- p[2]
        sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
        tm_pred <- pred_time_single(time_df$lag, sill, at, type = "exponential")
        SSE(space_df$cov, sp_pred) + SSE(time_df$cov, tm_pred)
      }
      
      opt <- optim(
        par    = init2,
        fn     = safe_fn(sse_total_fn),
        method = "L-BFGS-B",
        lower  = lower2,
        upper  = upper2
      )
      
      if (is.finite(opt$value) && all(is.finite(opt$par)) && opt$value < best_sse_total) {
        at <- opt$par[1]; Ar <- opt$par[2]
        
        sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
        tm_pred <- pred_time_single(time_df$lag, sill, at, type = "exponential")
        
        sse_space <- SSE(space_df$cov, sp_pred)
        sse_time  <- SSE(time_df$cov, tm_pred)
        sse_total <- sse_space + sse_time
        
        stm1 <- Ar / at
        stm2 <- Ar / at
        
        best_sse_total <- sse_total
        best <- tibble(
          variable          = var_name,
          domain            = domain,
          day               = day,
          sill              = sill,
          
          alpha1            = NA_real_,
          alpha2            = NA_real_,
          at1               = NA_real_,
          at2               = NA_real_,
          at                = at,
          Ar                = Ar,
          
          stm1              = stm1,
          stm2              = stm2,
          
          sse_space         = sse_space,
          sse_time          = sse_time,
          sse_total         = sse_total,
          
          converged         = (opt$convergence == 0),
          model_space       = space_model,
          model_time        = "exponential",
          
          n_space_lags_calc = n_space_lags_calc,
          n_time_lags_calc  = n_time_lags_calc,
          n_space_lags_used = n_space_lags_used,
          n_time_lags_used  = n_time_lags_used,
          
          space_file        = space_file,
          time_file         = time_file
        )
      }
    }
    
    # Time model: SINGLE GAUSSIAN
    {
      sse_total_fn <- function(p) {
        at <- p[1]; Ar <- p[2]
        sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
        tm_pred <- pred_time_single(time_df$lag, sill, at, type = "gaussian")
        SSE(space_df$cov, sp_pred) + SSE(time_df$cov, tm_pred)
      }
      
      opt <- optim(
        par    = init2,
        fn     = safe_fn(sse_total_fn),
        method = "L-BFGS-B",
        lower  = lower2,
        upper  = upper2
      )
      
      if (is.finite(opt$value) && all(is.finite(opt$par)) && opt$value < best_sse_total) {
        at <- opt$par[1]; Ar <- opt$par[2]
        
        sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
        tm_pred <- pred_time_single(time_df$lag, sill, at, type = "gaussian")
        
        sse_space <- SSE(space_df$cov, sp_pred)
        sse_time  <- SSE(time_df$cov, tm_pred)
        sse_total <- sse_space + sse_time
        
        stm1 <- Ar / at
        stm2 <- Ar / at
        
        best_sse_total <- sse_total
        best <- tibble(
          variable          = var_name,
          domain            = domain,
          day               = day,
          sill              = sill,
          
          alpha1            = NA_real_,
          alpha2            = NA_real_,
          at1               = NA_real_,
          at2               = NA_real_,
          at                = at,
          Ar                = Ar,
          
          stm1              = stm1,
          stm2              = stm2,
          
          sse_space         = sse_space,
          sse_time          = sse_time,
          sse_total         = sse_total,
          
          converged         = (opt$convergence == 0),
          model_space       = space_model,
          model_time        = "gaussian",
          
          n_space_lags_calc = n_space_lags_calc,
          n_time_lags_calc  = n_time_lags_calc,
          n_space_lags_used = n_space_lags_used,
          n_time_lags_used  = n_time_lags_used,
          
          space_file        = space_file,
          time_file         = time_file
        )
      }
    }
  }
  
  best
}

# -----------------------------
# Run
# -----------------------------
handler_txtprogressbar(format = ":current/:total [:bar] :elapsed | eta: :eta")

for (var_folder in VARS) {
  
  pairs <- discover_daily_pairs(var_folder)
  if (nrow(pairs) == 0) {
    message("[SKIP] No cov files found for: ", var_folder, " (domains=", paste(DOMAINS_KEEP, collapse=", "), ")")
    next
  }
  
  message("[RUN] ", var_folder, " | domain(s)=", paste(DOMAINS_KEEP, collapse=", "),
          " | n_pairs=", nrow(pairs))
  
  results_list <- NULL
  with_progress({
    p <- progressor(along = seq_len(nrow(pairs)))
    results_list <- future_lapply(seq_len(nrow(pairs)), function(i) {
      p()
      row <- pairs[i, ]
      fit_one_pair(
        space_file = row$space_file,
        time_file  = row$time_file,
        var_name   = var_folder,
        domain     = row$domain,
        day        = row$day
      )
    })
  })
  
  master <- bind_rows(results_list)
  
  out_dir <- file.path(OUT_BASE, var_folder)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_csv <- file.path(out_dir, paste0(var_folder, "_daily_optimal_cov_model_southeast.csv"))
  
  if (nrow(master) > 0) {
    master <- master %>% arrange(domain, day, sse_total)
    write_csv(master, out_csv)
    message("→ Wrote: ", out_csv)
  } else {
    warning("[WARN] No successful fits recorded for: ", var_folder,
            " (domain(s)=", paste(DOMAINS_KEEP, collapse=", "), ")")
  }
}

# statistical testing ----
# =============================================================================
# Hypothesis testing for Space-Time Metrics (STM1, STM2)
#
# Purpose:
#   - Test whether STM distributions differ from 1
#   - HeatIndex: summertime only (June-August)
#   - Precipitation: year-round
#   - Run separately for STM1 and STM2
#
# Tests:
#   1) Shapiro-Wilk normality test
#      - Because n > 5000, use a random sample of 5000 observations for Shapiro
#   2) One-sample Wilcoxon signed-rank test against mu = 1
#
# Inputs:
#   V:/users/hquintal/phd1_cluster_southeast_v2/data/out/03_cov_model/heatindex/heatindex_daily_optimal_cov_model_southeast.csv
#   V:/users/hquintal/phd1_cluster_southeast_v2/data/out/03_cov_model/precipitation/precipitation_daily_optimal_cov_model_southeast.csv
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(tibble)
})

# -----------------------------
# Paths
# -----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

HEAT_FILE <- file.path(
  ROOT,
  "data/out/03_cov_model/heatindex/heatindex_daily_optimal_cov_model_southeast.csv"
)

PREC_FILE <- file.path(
  ROOT,
  "data/out/03_cov_model/precipitation/precipitation_daily_optimal_cov_model_southeast.csv"
)

# -----------------------------
# Helpers
# -----------------------------
parse_day_safe <- function(x) {
  suppressWarnings(parse_date_time(x, orders = c("Y-m-d", "m/d/Y", "m/d/y")) |> as.Date())
}

prepare_stm <- function(file_path, var_name, metric_col = c("stm1", "stm2")) {
  metric_col <- match.arg(metric_col)
  
  dt <- readr::read_csv(file_path, show_col_types = FALSE)
  
  req <- c("day", metric_col)
  if (!all(req %in% names(dt))) {
    stop(
      "Missing required columns (",
      paste(req, collapse = ", "),
      ") in: ",
      file_path
    )
  }
  
  dt %>%
    mutate(
      date = parse_day_safe(as.character(day)),
      year = year(date),
      month_num = month(date),
      variable = var_name,
      stm_name = metric_col,
      stm_val = as.numeric(.data[[metric_col]])
    ) %>%
    filter(!is.na(date), is.finite(stm_val))
}

run_stm_tests <- function(dat, variable_name, metric_name, summer_only = FALSE, shapiro_n = 5000) {
  sub <- dat %>%
    filter(variable == variable_name, stm_name == metric_name)
  
  if (summer_only) {
    sub <- sub %>% filter(month_num %in% 6:8)
  }
  
  x <- sub$stm_val
  x <- x[is.finite(x)]
  
  if (length(x) < 3) {
    return(tibble(
      variable = variable_name,
      stm_name = metric_name,
      subset = ifelse(summer_only, "summer", "all_months"),
      n = length(x),
      shapiro_n = NA_integer_,
      shapiro_p = NA_real_,
      wilcox_p = NA_real_,
      wilcox_statistic = NA_real_,
      median_stm = NA_real_,
      mean_stm = NA_real_,
      min_stm = NA_real_,
      max_stm = NA_real_,
      n_less_1 = NA_integer_,
      n_equal_1 = NA_integer_,
      n_greater_1 = NA_integer_
    ))
  }
  
  # Shapiro-Wilk cannot be run on n > 5000 in base R
  # so use a reproducible sample of up to 5000 observations
  set.seed(123)
  x_shapiro <- if (length(x) > shapiro_n) sample(x, shapiro_n) else x
  
  shapiro_res <- shapiro.test(x_shapiro)
  wilcox_res <- wilcox.test(
    x,
    mu = 1,
    exact = FALSE,
    conf.int = TRUE
  )
  
  tibble(
    variable = variable_name,
    stm_name = metric_name,
    subset = ifelse(summer_only, "summer", "all_months"),
    n = length(x),
    shapiro_n = length(x_shapiro),
    shapiro_p = shapiro_res$p.value,
    wilcox_p = wilcox_res$p.value,
    wilcox_statistic = unname(wilcox_res$statistic),
    median_stm = median(x, na.rm = TRUE),
    mean_stm = mean(x, na.rm = TRUE),
    min_stm = min(x, na.rm = TRUE),
    max_stm = max(x, na.rm = TRUE),
    n_less_1 = sum(x < 1, na.rm = TRUE),
    n_equal_1 = sum(x == 1, na.rm = TRUE),
    n_greater_1 = sum(x > 1, na.rm = TRUE)
  )
}

print_test_summary <- function(res_row) {
  cat("\n============================================================\n")
  cat("Variable:   ", res_row$variable, "\n", sep = "")
  cat("Metric:     ", res_row$stm_name, "\n", sep = "")
  cat("Subset:     ", res_row$subset, "\n", sep = "")
  cat("N:          ", res_row$n, "\n", sep = "")
  cat("Shapiro N:  ", res_row$shapiro_n, "\n", sep = "")
  cat("Median STM: ", signif(res_row$median_stm, 5), "\n", sep = "")
  cat("Mean STM:   ", signif(res_row$mean_stm, 5), "\n", sep = "")
  cat("Range STM:  [", signif(res_row$min_stm, 5), ", ", signif(res_row$max_stm, 5), "]\n", sep = "")
  cat("Count < 1:  ", res_row$n_less_1, "\n", sep = "")
  cat("Count = 1:  ", res_row$n_equal_1, "\n", sep = "")
  cat("Count > 1:  ", res_row$n_greater_1, "\n", sep = "")
  cat("Shapiro p:  ", signif(res_row$shapiro_p, 5), "\n", sep = "")
  cat("Wilcox p:   ", signif(res_row$wilcox_p, 5), "\n", sep = "")
  cat("W statistic:", signif(res_row$wilcox_statistic, 7), "\n", sep = "")
  cat("Interpretation:\n")
  cat("  - Shapiro-Wilk tests normality of the sampled STM values.\n")
  cat("  - Wilcoxon signed-rank tests whether STM is centered at 1.\n")
  cat("============================================================\n")
}

# -----------------------------
# Read data
# -----------------------------
stm_heat_1 <- prepare_stm(HEAT_FILE, "HeatIndex", "stm1")
stm_heat_2 <- prepare_stm(HEAT_FILE, "HeatIndex", "stm2")
stm_prec_1 <- prepare_stm(PREC_FILE, "Precipitation", "stm1")
stm_prec_2 <- prepare_stm(PREC_FILE, "Precipitation", "stm2")

raw_data <- bind_rows(
  stm_heat_1,
  stm_heat_2,
  stm_prec_1,
  stm_prec_2
)

# -----------------------------
# Run tests
# -----------------------------
stm_test_results <- bind_rows(
  run_stm_tests(raw_data, "HeatIndex",     "stm1", summer_only = TRUE),
  run_stm_tests(raw_data, "HeatIndex",     "stm2", summer_only = TRUE),
  run_stm_tests(raw_data, "Precipitation", "stm1", summer_only = FALSE),
  run_stm_tests(raw_data, "Precipitation", "stm2", summer_only = FALSE)
)

# -----------------------------
# Print concise results table
# -----------------------------
print(stm_test_results)

# -----------------------------
# Print readable summaries
# -----------------------------
invisible(lapply(seq_len(nrow(stm_test_results)), function(i) {
  print_test_summary(stm_test_results[i, ])
}))

# # ~~~~~ ----
# # NCA regions ----
# # scripts/03_covariance_models_daily_nca.R
# # ------------------------------------------------------------------------------
# # Fit optimal covariance models to DAILY experimental covariance (space + time)
# # for NATIONAL CLIMATE ASSESSMENT (NCA) REGIONS
# #
# # UPDATE (2026-02-23): BATCHED processing (default 1000 pairs at a time) to
# # avoid R aborts from excessive parallel memory/CPU pressure. Also:
# #   - restart-safe (skips already-fit region|day keys)
# #   - append-only output writing (no repeated full CSV read/merge per batch)
# #   - safe per-pair tryCatch wrapper so single bad files don’t kill the run
# #   - optional year window filtering remains
# #   - limits BLAS thread oversubscription (OMP/OPENBLAS/MKL) to 1 thread
# # ------------------------------------------------------------------------------
# 
# suppressPackageStartupMessages({
#   library(dplyr)
#   library(stringr)
#   library(tibble)
#   library(readr)
#   library(progressr)
#   library(future)
#   library(future.apply)
#   library(tidyr)
#   library(fs)
#   library(data.table)   # NEW (for fast append writes)
# })
# 
# # USER SETTINGS ----
# # IMPORTANT: This workload is heavy (hundreds of optim() calls per pair).
# # 35 workers commonly causes oversubscription + memory pressure -> abort.
# WORKERS <- 12L
# future::plan(multisession, workers = WORKERS)
# 
# # Prevent BLAS oversubscription inside each worker (common abort culprit)
# Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
# 
# ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
# 
# # NCA experimental covariance files live here (per your note)
# IN_BASE  <- file.path(ROOT, "data/out/03_cov_model/nca")
# 
# # outputs requested here:
# OUT_BASE <- file.path(ROOT, "data/out/03_cov_model/nca")
# 
# VARS <- c("heatindex", "precipitation")
# 
# # ---- BATCHING (NEW) ----
# BATCH_SIZE <- 1000L      # number of (space,time) pairs per batch
# 
# # ---- Optional year window ----
# EVAL_START_YEAR <- NA_integer_   # e.g., 1940L
# EVAL_END_YEAR   <- NA_integer_   # e.g., 2023L
# 
# # ---- LAG GRIDS (used for zero-filling missing lags) ----
# base_dist <- 0.25
# rLag_grid <- base_dist * c(0:4, 6, 8, 12)   # 0,0.25,0.5,0.75,1,1.5,2,3
# tLag_grid <- c(0:4, 6, 8, 12)               # 0,1,2,3,4,6,8,12
# 
# # range of plausible values
# # (used for nested time model + space range; single-time uses subsets)
# init3  <- c(at1 = 20, at2 = 10, Ar = 1)
# lower3 <- c(at1 = 1,  at2 = 1,  Ar = 0.25)
# upper3 <- c(at1 = 200, at2 = 200, Ar = 100)
# 
# # single-time model bounds
# init2  <- c(at = 20, Ar = 1)
# lower2 <- c(at = 1,  Ar = 0.25)
# upper2 <- c(at = 200, Ar = 100)
# 
# alpha_grid <- seq(0.01, 0.99, by = 0.01)
# 
# # Helpers ----
# safe_read_cov <- function(path, lag_col_guess = c("rLag", "tLag", "lag")) {
#   x <- tryCatch(read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE),
#                 error = function(e) NULL)
#   if (is.null(x) || nrow(x) == 0) return(NULL)
#   
#   nm <- names(x)
#   
#   lag_col <- intersect(lag_col_guess, nm)
#   lag_col <- if (length(lag_col) == 0) nm[1] else lag_col[1]
#   
#   cov_col <- if ("Cov" %in% nm) "Cov" else if ("cov" %in% nm) "cov" else nm[2]
#   np_col  <- if ("np" %in% nm) "np" else NA_character_
#   
#   out <- tibble(
#     lag = as.numeric(x[[lag_col]]),
#     cov = as.numeric(x[[cov_col]]),
#     np  = if (!is.na(np_col)) as.numeric(x[[np_col]]) else NA_real_
#   ) %>%
#     filter(is.finite(lag), is.finite(cov)) %>%
#     arrange(lag)
#   
#   if (nrow(out) == 0) return(NULL)
#   out
# }
# 
# safe_fn <- function(fn) {
#   function(p) {
#     val <- tryCatch(fn(p), error = function(e) NA_real_)
#     if (!is.finite(val)) val <- .Machine$double.xmax / 2
#     val
#   }
# }
# 
# SSE <- function(obs, pred) sum((obs - pred)^2, na.rm = TRUE)
# 
# # ---- Clamp negatives to 0 + fill missing lags with 0 on a fixed grid ----
# standardize_cov_curve <- function(df, lag_grid) {
#   df2 <- df %>%
#     transmute(lag = as.numeric(lag),
#               cov = pmax(0, as.numeric(cov))) %>%
#     filter(is.finite(lag), is.finite(cov)) %>%
#     group_by(lag) %>%
#     summarise(cov = mean(cov), .groups = "drop")
#   
#   tibble(lag = lag_grid) %>%
#     left_join(df2, by = "lag") %>%
#     mutate(cov = replace_na(cov, 0)) %>%
#     arrange(lag)
# }
# 
# # Covariance model kernels (sill fixed)
# pred_space <- function(lag, sill, Ar, type = c("gaussian", "exponential")) {
#   type <- match.arg(type)
#   if (type == "gaussian") sill * exp(-3 * (lag^2) / (Ar^2)) else sill * exp(-3 * lag / Ar)
# }
# 
# pred_time_single <- function(lag, sill, at, type = c("gaussian", "exponential")) {
#   type <- match.arg(type)
#   if (type == "gaussian") sill * exp(-3 * (lag^2) / (at^2)) else sill * exp(-3 * lag / at)
# }
# 
# pred_time_nested <- function(lag, sill, alpha, at1, at2) {
#   sill * (alpha * exp(-3 * (lag^2) / (at1^2)) + (1 - alpha) * exp(-3 * lag / at2))
# }
# 
# # ------------------------------------------------------------------------------
# # NCA discovery (Region may include underscores/hyphens/etc.)
# # Parse:
# #   <prefix>_<Region>_day_<YYYY-MM-DD>_cov_space.txt
# # Variable assignment comes from folder name, not prefix.
# # ------------------------------------------------------------------------------
# discover_daily_pairs_nca <- function(var_folder) {
#   cov_dir <- file.path(IN_BASE, var_folder, "daily")
#   if (!dir.exists(cov_dir)) return(tibble())
#   
#   sp_files <- dir_ls(cov_dir, type = "file", regexp = "_cov_space\\.txt$", fail = FALSE)
#   sp_files <- sort(as.character(sp_files))
#   
#   rx <- "^(.+?)_(.+?)_day_([0-9]{4}-[0-9]{2}-[0-9]{2})_cov_space\\.txt$"
#   
#   tibble(space_file = sp_files) %>%
#     mutate(
#       fname       = path_file(space_file),
#       matched     = str_match(fname, rx),
#       file_prefix = matched[, 2],
#       region      = matched[, 3],
#       day         = matched[, 4],
#       year        = suppressWarnings(as.integer(substr(day, 1, 4))),
#       time_file   = file.path(cov_dir, str_replace(fname, "_cov_space\\.txt$", "_cov_time.txt"))
#     ) %>%
#     filter(!is.na(region), !is.na(day), is.finite(year), file.exists(time_file)) %>%
#     select(space_file, time_file, region, day, year, file_prefix)
# }
# 
# fit_one_pair <- function(space_file, time_file, var_name, region, day) {
#   
#   sp_raw <- safe_read_cov(space_file, lag_col_guess = c("rLag", "lag"))
#   tm_raw <- safe_read_cov(time_file,  lag_col_guess = c("tLag", "lag"))
#   if (is.null(sp_raw) || is.null(tm_raw)) return(NULL)
#   
#   # Standardize curves: negative -> 0, fill missing lags with 0
#   space_df <- standardize_cov_curve(sp_raw, rLag_grid)
#   time_df  <- standardize_cov_curve(tm_raw, tLag_grid)
#   
#   # Need a valid sill at lag 0 (after clamping)
#   if (!any(space_df$lag == 0) || !any(time_df$lag == 0)) return(NULL)
#   
#   sill <- space_df$cov[which(space_df$lag == 0)[1]]
#   if (!is.finite(sill) || sill <= 0) return(NULL)
#   
#   # diagnostics
#   n_space_lags_calc <- nrow(sp_raw)
#   n_time_lags_calc  <- nrow(tm_raw)
#   n_space_lags_used <- nrow(space_df)
#   n_time_lags_used  <- nrow(time_df)
#   
#   space_models <- c("exponential", "gaussian")
#   
#   best <- NULL
#   best_sse_total <- Inf
#   
#   for (space_model in space_models) {
#     
#     # ---- Time model: NESTED (alpha * Gauss(at1) + (1-alpha) * Exp(at2)) ----
#     for (alpha in alpha_grid) {
#       
#       sse_total_fn <- function(p) {
#         at1 <- p[1]; at2 <- p[2]; Ar <- p[3]
#         sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
#         tm_pred <- pred_time_nested(time_df$lag, sill, alpha, at1, at2)
#         SSE(space_df$cov, sp_pred) + SSE(time_df$cov, tm_pred)
#       }
#       
#       opt <- optim(
#         par    = init3,
#         fn     = safe_fn(sse_total_fn),
#         method = "L-BFGS-B",
#         lower  = lower3,
#         upper  = upper3
#       )
#       
#       if (is.finite(opt$value) && all(is.finite(opt$par)) && opt$value < best_sse_total) {
#         at1 <- opt$par[1]; at2 <- opt$par[2]; Ar <- opt$par[3]
#         
#         sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
#         tm_pred <- pred_time_nested(time_df$lag, sill, alpha, at1, at2)
#         
#         sse_space <- SSE(space_df$cov, sp_pred)
#         sse_time  <- SSE(time_df$cov, tm_pred)
#         sse_total <- sse_space + sse_time
#         
#         stm1 <- Ar / (alpha * at1 + (1 - alpha) * at2)
#         stm2 <- alpha * (Ar / at1) + (1 - alpha) * (Ar / at2)
#         
#         best_sse_total <- sse_total
#         best <- tibble(
#           variable          = var_name,
#           region            = region,
#           day               = day,
#           sill              = sill,
#           
#           alpha1            = alpha,
#           alpha2            = 1 - alpha,
#           at1               = at1,
#           at2               = at2,
#           at                = NA_real_,
#           Ar                = Ar,
#           
#           stm1              = stm1,
#           stm2              = stm2,
#           
#           sse_space         = sse_space,
#           sse_time          = sse_time,
#           sse_total         = sse_total,
#           
#           converged         = (opt$convergence == 0),
#           model_space       = space_model,
#           model_time        = "nested_gauss_exp",
#           
#           n_space_lags_calc = n_space_lags_calc,
#           n_time_lags_calc  = n_time_lags_calc,
#           n_space_lags_used = n_space_lags_used,
#           n_time_lags_used  = n_time_lags_used,
#           
#           space_file        = space_file,
#           time_file         = time_file
#         )
#       }
#     }
#     
#     # ---- Time model: SINGLE EXPONENTIAL ----
#     {
#       sse_total_fn <- function(p) {
#         at <- p[1]; Ar <- p[2]
#         sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
#         tm_pred <- pred_time_single(time_df$lag, sill, at, type = "exponential")
#         SSE(space_df$cov, sp_pred) + SSE(time_df$cov, tm_pred)
#       }
#       
#       opt <- optim(
#         par    = init2,
#         fn     = safe_fn(sse_total_fn),
#         method = "L-BFGS-B",
#         lower  = lower2,
#         upper  = upper2
#       )
#       
#       if (is.finite(opt$value) && all(is.finite(opt$par)) && opt$value < best_sse_total) {
#         at <- opt$par[1]; Ar <- opt$par[2]
#         
#         sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
#         tm_pred <- pred_time_single(time_df$lag, sill, at, type = "exponential")
#         
#         sse_space <- SSE(space_df$cov, sp_pred)
#         sse_time  <- SSE(time_df$cov, tm_pred)
#         sse_total <- sse_space + sse_time
#         
#         best_sse_total <- sse_total
#         best <- tibble(
#           variable          = var_name,
#           region            = region,
#           day               = day,
#           sill              = sill,
#           
#           alpha1            = NA_real_,
#           alpha2            = NA_real_,
#           at1               = NA_real_,
#           at2               = NA_real_,
#           at                = at,
#           Ar                = Ar,
#           
#           stm1              = Ar / at,
#           stm2              = Ar / at,
#           
#           sse_space         = sse_space,
#           sse_time          = sse_time,
#           sse_total         = sse_total,
#           
#           converged         = (opt$convergence == 0),
#           model_space       = space_model,
#           model_time        = "exponential",
#           
#           n_space_lags_calc = n_space_lags_calc,
#           n_time_lags_calc  = n_time_lags_calc,
#           n_space_lags_used = n_space_lags_used,
#           n_time_lags_used  = n_time_lags_used,
#           
#           space_file        = space_file,
#           time_file         = time_file
#         )
#       }
#     }
#     
#     # ---- Time model: SINGLE GAUSSIAN ----
#     {
#       sse_total_fn <- function(p) {
#         at <- p[1]; Ar <- p[2]
#         sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
#         tm_pred <- pred_time_single(time_df$lag, sill, at, type = "gaussian")
#         SSE(space_df$cov, sp_pred) + SSE(time_df$cov, tm_pred)
#       }
#       
#       opt <- optim(
#         par    = init2,
#         fn     = safe_fn(sse_total_fn),
#         method = "L-BFGS-B",
#         lower  = lower2,
#         upper  = upper2
#       )
#       
#       if (is.finite(opt$value) && all(is.finite(opt$par)) && opt$value < best_sse_total) {
#         at <- opt$par[1]; Ar <- opt$par[2]
#         
#         sp_pred <- pred_space(space_df$lag, sill, Ar, type = space_model)
#         tm_pred <- pred_time_single(time_df$lag, sill, at, type = "gaussian")
#         
#         sse_space <- SSE(space_df$cov, sp_pred)
#         sse_time  <- SSE(time_df$cov, tm_pred)
#         sse_total <- sse_space + sse_time
#         
#         best_sse_total <- sse_total
#         best <- tibble(
#           variable          = var_name,
#           region            = region,
#           day               = day,
#           sill              = sill,
#           
#           alpha1            = NA_real_,
#           alpha2            = NA_real_,
#           at1               = NA_real_,
#           at2               = NA_real_,
#           at                = at,
#           Ar                = Ar,
#           
#           stm1              = Ar / at,
#           stm2              = Ar / at,
#           
#           sse_space         = sse_space,
#           sse_time          = sse_time,
#           sse_total         = sse_total,
#           
#           converged         = (opt$convergence == 0),
#           model_space       = space_model,
#           model_time        = "gaussian",
#           
#           n_space_lags_calc = n_space_lags_calc,
#           n_time_lags_calc  = n_time_lags_calc,
#           n_space_lags_used = n_space_lags_used,
#           n_time_lags_used  = n_time_lags_used,
#           
#           space_file        = space_file,
#           time_file         = time_file
#         )
#       }
#     }
#   }
#   
#   best
# }
# 
# # NEW: safe wrapper so one broken file doesn’t crash the whole run
# fit_one_pair_safe <- function(space_file, time_file, var_name, region, day) {
#   tryCatch(
#     fit_one_pair(space_file, time_file, var_name, region, day),
#     error = function(e) {
#       tibble(
#         variable = var_name, region = region, day = day,
#         sill = NA_real_,
#         alpha1 = NA_real_, alpha2 = NA_real_,
#         at1 = NA_real_, at2 = NA_real_, at = NA_real_, Ar = NA_real_,
#         stm1 = NA_real_, stm2 = NA_real_,
#         sse_space = NA_real_, sse_time = NA_real_, sse_total = NA_real_,
#         converged = FALSE,
#         model_space = NA_character_, model_time = NA_character_,
#         n_space_lags_calc = NA_integer_, n_time_lags_calc = NA_integer_,
#         n_space_lags_used = NA_integer_, n_time_lags_used = NA_integer_,
#         space_file = space_file, time_file = time_file,
#         error = paste0("fit_error: ", conditionMessage(e))
#       )
#     }
#   )
# }
# 
# # Run (BATCHED + RESTART SAFE) ----
# handler_txtprogressbar(format = ":current/:total [:bar] :elapsed | eta: :eta")
# 
# for (var_folder in VARS) {
#   
#   pairs_all <- discover_daily_pairs_nca(var_folder)
#   if (nrow(pairs_all) == 0) {
#     message("[SKIP] No cov files found for: ", var_folder)
#     next
#   }
#   
#   # optional year filtering
#   if (is.finite(EVAL_START_YEAR)) pairs_all <- pairs_all %>% filter(year >= EVAL_START_YEAR)
#   if (is.finite(EVAL_END_YEAR))   pairs_all <- pairs_all %>% filter(year <= EVAL_END_YEAR)
#   
#   if (nrow(pairs_all) == 0) {
#     message("[SKIP] No pairs after year filtering for: ", var_folder)
#     next
#   }
#   
#   out_dir <- file.path(OUT_BASE, var_folder, "model")
#   dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#   out_csv <- file.path(out_dir, paste0(var_folder, "_daily_optimal_cov_model_nca.csv"))
#   
#   # Restart-safe done keys (read ONCE)
#   done_key <- character(0)
#   if (file.exists(out_csv)) {
#     old <- tryCatch(readr::read_csv(out_csv, show_col_types = FALSE), error = function(e) NULL)
#     if (!is.null(old) && nrow(old) > 0 && all(c("region","day") %in% names(old))) {
#       done_key <- unique(paste0(old$region, "|", old$day))
#       message("[RESUME] existing fits: ", length(done_key))
#     }
#   }
#   
#   # Filter to TODO once
#   pairs_todo <- pairs_all %>%
#     mutate(key = paste0(region, "|", day)) %>%
#     filter(!(key %in% done_key)) %>%
#     select(-key)
#   
#   n_total <- nrow(pairs_all)
#   n_todo  <- nrow(pairs_todo)
#   message(sprintf("[DISCOVER] %s | total_pairs=%d | todo=%d | batch_size=%d | out=%s",
#                   var_folder, n_total, n_todo, BATCH_SIZE, out_csv))
#   
#   if (n_todo == 0) {
#     message("[DONE] nothing new to fit for: ", var_folder)
#     next
#   }
#   
#   # Batch indices
#   idx <- seq_len(n_todo)
#   batch_id <- ceiling(idx / BATCH_SIZE)
#   n_batches <- max(batch_id)
#   
#   for (b in seq_len(n_batches)) {
#     sel <- which(batch_id == b)
#     pairs_b <- pairs_todo[sel, , drop = FALSE]
#     
#     message(sprintf("[RUN] %s | batch %d/%d | n=%d",
#                     var_folder, b, n_batches, nrow(pairs_b)))
#     
#     results_list <- NULL
#     with_progress({
#       p <- progressor(along = seq_len(nrow(pairs_b)))
#       results_list <- future_lapply(seq_len(nrow(pairs_b)), function(i) {
#         p()
#         row <- pairs_b[i, ]
#         fit_one_pair_safe(
#           space_file = row$space_file,
#           time_file  = row$time_file,
#           var_name   = var_folder,
#           region     = row$region,
#           day        = row$day
#         )
#       }, future.seed = TRUE)
#     })
#     
#     master_new <- bind_rows(results_list)
#     
#     if (nrow(master_new) == 0) {
#       warning(sprintf("[WARN] %s | batch %d/%d | No rows returned", var_folder, b, n_batches))
#       rm(results_list); gc()
#       next
#     }
#     
#     # Append-only write (fast, avoids reading huge CSV each batch)
#     dt_new <- as.data.table(master_new)
#     
#     write_header <- !file.exists(out_csv)
#     fwrite(dt_new, out_csv, append = !write_header)
#     
#     # Update done_key in-memory so we don't double-do within a run
#     done_key <- unique(c(done_key, paste0(dt_new$region, "|", dt_new$day)))
#     
#     message(sprintf("→ Appended: %s | +%d rows | done_now=%d",
#                     out_csv, nrow(dt_new), length(done_key)))
#     
#     rm(results_list, master_new, dt_new)
#     gc()
#   }
# }