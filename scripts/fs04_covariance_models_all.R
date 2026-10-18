# ~~~~~~~ ----
# ~STACKED CURVES ACROSS ALL DAYS~ ----
# scripts/fig_01_cov_model_all_days_4panel.R
# Plot ALL daily covariance MODEL curves for ALL days (HI + PR) in a 4-panel grid
# - Reads *_daily_optimal_cov_model.csv for heatindex + precipitation
# - One curve per day (semi-transparent), no legends by default
# - Fixed x-axes: space 0–3 deg, time 0–12 hr
# - Clamps negative model covariances to 0
# - y-lims determined from max(model at lag=0) per variable (HI vs PR)
# - DOES NOT require overlapping HI/PR days

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# ---- USER SETTINGS ----
ROOT <- "/proj/sebastianlab/users/hquintal/cisc"

HI_MODEL_CSV <- file.path(ROOT, "data/out/03_cov_model/heatindex/heatindex_daily_optimal_cov_model_southeast.csv")
PR_MODEL_CSV <- file.path(ROOT, "data/out/03_cov_model/precipitation/precipitation_daily_optimal_cov_model_southeast.csv")

OUT_PNG <- file.path(ROOT, "figures","fs04_cov_model_all", "fs04_cov_model_all.png")
dir.create(dirname(OUT_PNG), recursive = TRUE, showWarnings = FALSE)

# Fixed x domains
X_SPACE <- c(0, 3)
X_TIME  <- c(0, 12)

# Curve resolution
N_LAG <- 400L

# Kernels (same conventions as your per-day script)
.space_kernel <- function(r, Ar, model_space) {
  ms <- tolower(trimws(model_space))
  if (ms == "gaussian") {
    exp(-3 * r^2 / (Ar^2))
  } else if (ms == "exponential") {
    exp(-3 * r / Ar)
  } else {
    stop("Unsupported model_space: ", model_space)
  }
}

.time_kernel <- function(t, model_time,
                         alpha1 = NA, alpha2 = NA,
                         at1 = NA, at2 = NA,
                         at = NA) {
  
  mt <- tolower(trimws(model_time))
  
  if (mt %in% c("nested_gauss_exp", "nested_gaussian_exponential", "nested_gauss_exponential")) {
    if (!is.finite(alpha1) || !is.finite(alpha2) || !is.finite(at1) || !is.finite(at2)) {
      stop("alpha1, alpha2, at1, at2 required for nested_gauss_exp time model.")
    }
    return(alpha1 * exp(-3 * (t^2) / (at1^2)) + alpha2 * exp(-3 * t / at2))
  }
  
  if (mt == "gaussian") {
    if (!is.finite(at)) stop("at required for gaussian time model.")
    return(exp(-3 * (t^2) / (at^2)))
  }
  
  if (mt == "exponential") {
    if (!is.finite(at)) stop("at required for exponential time model.")
    return(exp(-3 * t / at))
  }
  
  stop("Unsupported model_time: ", model_time)
}

.cov_curve_space <- function(r, sill, Ar, model_space) {
  sill * .space_kernel(r, Ar, model_space)
}

.cov_curve_time <- function(t, sill, model_time, alpha1, alpha2, at1, at2, at) {
  sill * .time_kernel(t, model_time, alpha1, alpha2, at1, at2, at)
}

# Read + sanitize daily model tables
.read_models <- function(path, var_name) {
  dt <- fread(path)
  
  # robust day parsing (handles "1/10/1940" etc.)
  dt[, day_date := as.IDate(as.Date(day))]
  
  # converged columns differ across outputs; accept either "converged" or converged_sp/tm
  if ("converged" %in% names(dt)) {
    dt[, converged_all := tolower(trimws(as.character(converged))) %in% c("true","t","1","yes")]
  } else if (all(c("converged_sp","converged_tm") %in% names(dt))) {
    dt[, converged_all :=
         tolower(trimws(as.character(converged_sp))) %in% c("true","t","1","yes") &
         tolower(trimws(as.character(converged_tm))) %in% c("true","t","1","yes")]
  } else {
    warning("No converged or converged_sp/tm columns found in: ", path, " (keeping all rows).")
    dt[, converged_all := TRUE]
  }
  
  dt <- dt[converged_all == TRUE & !is.na(day_date)]
  
  # one row per day: best sse_total if present
  if ("sse_total" %in% names(dt)) {
    setorder(dt, day_date, sse_total)
    dt <- dt[, .SD[1], by = day_date]
  } else {
    dt <- unique(dt, by = "day_date")
  }
  
  dt[, variable := var_name]
  dt[]
}

hi_mod <- .read_models(HI_MODEL_CSV, "heatindex")
pr_mod <- .read_models(PR_MODEL_CSV, "precipitation")

message(sprintf("[HI] days=%d", nrow(hi_mod)))
message(sprintf("[PR] days=%d", nrow(pr_mod)))

# Build model curves for ALL days (modeled only, no exp points)
lag_space <- seq(X_SPACE[1], X_SPACE[2], length.out = N_LAG)
lag_time  <- seq(X_TIME[1],  X_TIME[2],  length.out = N_LAG)

.make_curves_dt <- function(row) {
  sill <- row$sill
  Ar   <- row$Ar
  ms   <- row$model_space
  mt   <- row$model_time
  
  a1  <- if ("alpha1" %in% names(row)) row$alpha1 else NA_real_
  a2  <- if ("alpha2" %in% names(row)) row$alpha2 else NA_real_
  at1 <- if ("at1"    %in% names(row)) row$at1    else NA_real_
  at2 <- if ("at2"    %in% names(row)) row$at2    else NA_real_
  at  <- if ("at"     %in% names(row)) row$at     else NA_real_
  
  cs <- .cov_curve_space(lag_space, sill, Ar, ms)
  ct <- .cov_curve_time(lag_time,  sill, mt, a1, a2, at1, at2, at)
  
  cs <- pmax(cs, 0)
  ct <- pmax(ct, 0)
  
  id <- format(row$day_date, "%Y-%m-%d")
  
  rbind(
    data.table(variable=row$variable, day_date=row$day_date, day_id=id, panel="space", lag=lag_space, cov_mod=cs),
    data.table(variable=row$variable, day_date=row$day_date, day_id=id, panel="time",  lag=lag_time,  cov_mod=ct)
  )
}

.safe_apply <- function(dt) {
  out <- vector("list", nrow(dt))
  for (i in seq_len(nrow(dt))) {
    row <- dt[i]
    out[[i]] <- tryCatch(.make_curves_dt(row), error = function(e) NULL)
  }
  rbindlist(out, fill = TRUE)
}

curves <- rbindlist(list(.safe_apply(hi_mod), .safe_apply(pr_mod)), fill = TRUE)
if (nrow(curves) == 0) stop("No curves produced (all rows failed).")

# Y limits per variable based on max at lag==0 across space+time
ymax_by_var <- curves[lag == 0, .(ymax0 = max(cov_mod, na.rm = TRUE)), by = variable]
ymax_by_var[, ymax := ifelse(is.finite(ymax0) & ymax0 > 0, ymax0 * 1.05, 1e-12)]

hi_ymax <- ymax_by_var[variable == "heatindex", ymax]
pr_ymax <- ymax_by_var[variable == "precipitation", ymax]

# Panel helper
.make_panel_all <- function(dt, xlab, ylab, xlim_fixed, ylim_fixed) {
  ggplot(dt, aes(x = lag, y = cov_mod, group = day_id)) +
    geom_line(alpha = 0.03, linewidth = 0.5) +   # <-- change alpha here
    scale_x_continuous(limits = xlim_fixed, expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, ylim_fixed), clip = "off") +
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      plot.margin = margin(6, 10, 6, 10)
    )
}

# Split + plot
xlab_space <- "Spatial degrees"
xlab_time  <- "Hours"

p_hi_space <- .make_panel_all(curves[variable=="heatindex"     & panel=="space"], xlab_space, "Covariance (C²)",  X_SPACE, hi_ymax)
p_hi_time  <- .make_panel_all(curves[variable=="heatindex"     & panel=="time"],  xlab_time,  "Covariance (C²)",  X_TIME,  hi_ymax)
p_pr_space <- .make_panel_all(curves[variable=="precipitation" & panel=="space"], xlab_space, "Covariance (mm²)", X_SPACE, pr_ymax)
p_pr_time  <- .make_panel_all(curves[variable=="precipitation" & panel=="time"],  xlab_time,  "Covariance (mm²)", X_TIME,  pr_ymax)

fig <- wrap_plots(p_hi_space, p_hi_time, p_pr_space, p_pr_time, ncol = 2)

ggsave(OUT_PNG, fig, width = 12, height = 8, dpi = 300, device = "png")
message("[WROTE] ", OUT_PNG)
