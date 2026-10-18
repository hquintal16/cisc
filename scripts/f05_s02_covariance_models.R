# scripts/fig_03_cov_model_4panel_simple.R
# =============================================================================
# Batch 4-panel covariance-model plots for dates where both heat index and
# precipitation have:
#
#   1. a valid parsed date,
#   2. a converged model row, when a converged column is present,
#   3. existing spatial and temporal empirical-covariance files,
#   4. a date inside the requested inclusive date window.
#
# Output:
#   /proj/sebastianlab/users/hquintal/cisc/figures/f05_cov_model/
#     fig_03_cov_model_YYYY-MM-DD.png
#
# Plot rules:
#   - spatial x-axis fixed at 0-3 degrees
#   - temporal x-axis fixed at 0-12 hours
#   - y-axis starts at zero
#   - negative empirical and modeled covariance values clamped to zero
#   - SSE shown in upper-right corner
#   - panel tags (a)-(d)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(lubridate)
})

# =============================================================================
# USER SETTINGS
# =============================================================================

ROOT <- "/proj/sebastianlab/users/hquintal/cisc"

HI_MODEL_CSV <- file.path(
  ROOT,
  "data/out/03_cov_model/heatindex",
  "heatindex_daily_optimal_cov_model_southeast.csv"
)

PR_MODEL_CSV <- file.path(
  ROOT,
  "data/out/03_cov_model/precipitation",
  "precipitation_daily_optimal_cov_model_southeast.csv"
)

OUT_DIR <- file.path(
  ROOT,
  "figures/f05_cov_model"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Inclusive date window.
# Set either value to NA_character_ to disable that bound.
DATE_START <- "2023-06-01"
DATE_END <- "2023-08-31"

# Fixed x-axis domains.
X_SPACE <- c(0, 3)
X_TIME <- c(0, 12)

# Skip figures that already exist.
SKIP_EXISTING <- TRUE

# =============================================================================
# GLOBAL PLOT THEME
# =============================================================================

theme_set(
  theme_bw(
    base_size = 20,
    base_family = "Arial"
  ) +
    theme(
      text = element_text(
        family = "Arial",
        size = 20
      ),
      axis.title = element_text(
        family = "Arial",
        size = 20
      ),
      axis.text = element_text(
        family = "Arial",
        size = 12
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        color = "grey90"
      )
    )
)

# =============================================================================
# LOGGING
# =============================================================================

.log <- function(...) {
  message(
    sprintf(
      "[%s] %s",
      format(Sys.time(), "%H:%M:%S"),
      sprintf(...)
    )
  )
}

# =============================================================================
# GENERAL HELPERS
# =============================================================================

.parse_day <- function(x) {
  # Parse each value independently rather than selecting one parser for the
  # entire vector.
  
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  
  out <- rep(
    as.Date(NA),
    length(x)
  )
  
  # ISO date or ISO datetime:
  # YYYY-MM-DD
  # YYYY-MM-DD HH:MM:SS
  # YYYY-MM-DDTHH:MM:SS
  missing <- is.na(out) & !is.na(x)
  
  if (any(missing)) {
    out[missing] <- suppressWarnings(
      as.Date(
        substr(
          x[missing],
          1,
          10
        ),
        format = "%Y-%m-%d"
      )
    )
  }
  
  # Month/day/year.
  missing <- is.na(out) & !is.na(x)
  
  if (any(missing)) {
    out[missing] <- as.Date(
      suppressWarnings(
        mdy(
          x[missing],
          quiet = TRUE,
          tz = "UTC"
        )
      )
    )
  }
  
  # Day/month/year.
  missing <- is.na(out) & !is.na(x)
  
  if (any(missing)) {
    out[missing] <- as.Date(
      suppressWarnings(
        dmy(
          x[missing],
          quiet = TRUE,
          tz = "UTC"
        )
      )
    )
  }
  
  # Compact YYYYMMDD.
  missing <- is.na(out) & !is.na(x)
  
  if (any(missing)) {
    out[missing] <- suppressWarnings(
      as.Date(
        x[missing],
        format = "%Y%m%d"
      )
    )
  }
  
  out
}

.parse_bound <- function(x, label) {
  if (
    length(x) == 0L ||
    is.na(x) ||
    !nzchar(trimws(x))
  ) {
    return(as.Date(NA))
  }
  
  out <- .parse_day(x)
  
  if (is.na(out)) {
    stop(
      sprintf(
        "Could not parse %s: %s",
        label,
        x
      )
    )
  }
  
  out
}

.as_logical01 <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  
  if (is.numeric(x)) {
    out <- rep(NA, length(x))
    out[!is.na(x) & x == 1] <- TRUE
    out[!is.na(x) & x == 0] <- FALSE
    return(out)
  }
  
  x <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  out <- rep(
    NA,
    length(x)
  )
  
  out[
    x %in% c(
      "true",
      "t",
      "1",
      "yes",
      "y",
      "converged",
      "success"
    )
  ] <- TRUE
  
  out[
    x %in% c(
      "false",
      "f",
      "0",
      "no",
      "n",
      "failed",
      "failure",
      "not converged"
    )
  ] <- FALSE
  
  out
}

.normalize_cov_path <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  
  # Convert Windows separators to Unix separators.
  x <- gsub(
    "\\\\",
    "/",
    x
  )
  
  # Convert the original Windows project path to the current Longleaf root.
  x <- sub(
    "^V:/users/hquintal/phd1_cluster_southeast_v2/",
    paste0(ROOT, "/"),
    x,
    ignore.case = TRUE
  )
  
  # Also support paths previously written using the cisc project name.
  x <- sub(
    "^V:/users/hquintal/cisc/",
    paste0(ROOT, "/"),
    x,
    ignore.case = TRUE
  )
  
  # Convert relative paths to paths beneath ROOT.
  is_relative <- (
    !is.na(x) &
      nzchar(x) &
      !grepl("^/", x) &
      !grepl("^[A-Za-z]:/", x)
  )
  
  x[is_relative] <- file.path(
    ROOT,
    x[is_relative]
  )
  
  x
}

.fmt2 <- function(x) {
  if (
    length(x) == 0L ||
    is.na(x) ||
    !is.finite(as.numeric(x))
  ) {
    return("NA")
  }
  
  formatC(
    as.numeric(x),
    format = "f",
    digits = 2
  )
}

.require_columns <- function(dt, required, table_name) {
  missing <- setdiff(
    required,
    names(dt)
  )
  
  if (length(missing) > 0L) {
    stop(
      sprintf(
        "%s is missing required columns: %s",
        table_name,
        paste(
          missing,
          collapse = ", "
        )
      )
    )
  }
}

.print_date_summary <- function(dt, label) {
  n_parsed <- sum(
    !is.na(dt$day_date)
  )
  
  .log(
    "%s parsed dates: %d/%d",
    label,
    n_parsed,
    nrow(dt)
  )
  
  if (n_parsed > 0L) {
    .log(
      "%s date range: %s to %s",
      label,
      format(
        min(
          dt$day_date,
          na.rm = TRUE
        ),
        "%Y-%m-%d"
      ),
      format(
        max(
          dt$day_date,
          na.rm = TRUE
        ),
        "%Y-%m-%d"
      )
    )
  } else {
    .log(
      "%s date range: unavailable because no dates parsed",
      label
    )
  }
}

.print_unparsed_examples <- function(dt, label) {
  bad <- unique(
    dt[
      is.na(day_date) &
        !is.na(day),
      as.character(day)
    ]
  )
  
  if (length(bad) > 0L) {
    .log(
      "%s unparsed day examples:\n%s",
      label,
      paste(
        head(
          bad,
          10
        ),
        collapse = "\n"
      )
    )
  }
}

.print_missing_path_examples <- function(
    dt,
    exists_column,
    path_column,
    label,
    n = 5L
) {
  bad_paths <- unique(
    dt[
      get(exists_column) == FALSE |
        is.na(get(exists_column)),
      get(path_column)
    ]
  )
  
  bad_paths <- bad_paths[
    !is.na(bad_paths)
  ]
  
  if (length(bad_paths) > 0L) {
    .log(
      "%s missing path examples:\n%s",
      label,
      paste(
        head(
          bad_paths,
          n
        ),
        collapse = "\n"
      )
    )
  }
}

# =============================================================================
# COVARIANCE HELPERS
# =============================================================================

.read_exp_cov <- function(path) {
  if (
    length(path) == 0L ||
    is.na(path) ||
    !nzchar(path)
  ) {
    stop(
      "Empirical covariance path is empty."
    )
  }
  
  if (!file.exists(path)) {
    stop(
      "Empirical covariance file does not exist: ",
      path
    )
  }
  
  dt <- tryCatch(
    fread(path),
    error = function(e) {
      stop(
        "Failed to read empirical covariance file ",
        path,
        ": ",
        conditionMessage(e)
      )
    }
  )
  
  if (nrow(dt) == 0L) {
    stop(
      "Empirical covariance file contains zero rows: ",
      path
    )
  }
  
  if ("rLag" %in% names(dt)) {
    setnames(
      dt,
      "rLag",
      "lag"
    )
  }
  
  if ("tLag" %in% names(dt)) {
    setnames(
      dt,
      "tLag",
      "lag"
    )
  }
  
  if (!"lag" %in% names(dt)) {
    stop(
      "Expected a lag, rLag, or tLag column in: ",
      path
    )
  }
  
  if (!"Cov" %in% names(dt)) {
    stop(
      "Expected a Cov column in: ",
      path
    )
  }
  
  setnames(
    dt,
    "Cov",
    "cov_exp"
  )
  
  dt[
    ,
    `:=`(
      lag = as.numeric(lag),
      cov_exp = as.numeric(cov_exp)
    )
  ]
  
  dt[]
}

.space_kernel <- function(
    r,
    Ar,
    model_space
) {
  ms <- tolower(
    trimws(
      as.character(model_space)
    )
  )
  
  Ar <- as.numeric(Ar)
  
  if (
    !is.finite(Ar) ||
    Ar <= 0
  ) {
    stop(
      "Invalid spatial range parameter Ar: ",
      Ar
    )
  }
  
  if (ms == "gaussian") {
    return(
      exp(
        -3 * r^2 / Ar^2
      )
    )
  }
  
  if (ms == "exponential") {
    return(
      exp(
        -3 * r / Ar
      )
    )
  }
  
  stop(
    "Unsupported model_space: ",
    model_space,
    ". Expected gaussian or exponential."
  )
}

.time_kernel <- function(
    t,
    model_time,
    alpha1 = NA_real_,
    alpha2 = NA_real_,
    at1 = NA_real_,
    at2 = NA_real_,
    at = NA_real_
) {
  mt <- tolower(
    trimws(
      as.character(model_time)
    )
  )
  
  alpha1 <- as.numeric(alpha1)
  alpha2 <- as.numeric(alpha2)
  at1 <- as.numeric(at1)
  at2 <- as.numeric(at2)
  at <- as.numeric(at)
  
  nested_names <- c(
    "nested_gauss_exp",
    "nested_gaussian_exponential",
    "nested_gauss_exponential"
  )
  
  if (mt %in% nested_names) {
    required_values <- c(
      alpha1,
      alpha2,
      at1,
      at2
    )
    
    if (
      any(!is.finite(required_values)) ||
      at1 <= 0 ||
      at2 <= 0
    ) {
      stop(
        paste0(
          "Valid alpha1, alpha2, at1, and at2 values are required ",
          "for nested_gauss_exp."
        )
      )
    }
    
    return(
      alpha1 *
        exp(
          -3 * t^2 / at1^2
        ) +
        alpha2 *
        exp(
          -3 * t / at2
        )
    )
  }
  
  if (mt == "gaussian") {
    if (
      !is.finite(at) ||
      at <= 0
    ) {
      stop(
        "A valid positive at value is required for a Gaussian time model."
      )
    }
    
    return(
      exp(
        -3 * t^2 / at^2
      )
    )
  }
  
  if (mt == "exponential") {
    if (
      !is.finite(at) ||
      at <= 0
    ) {
      stop(
        "A valid positive at value is required for an exponential time model."
      )
    }
    
    return(
      exp(
        -3 * t / at
      )
    )
  }
  
  stop(
    paste0(
      "Unsupported model_time: ",
      model_time,
      ". Expected nested_gauss_exp, gaussian, or exponential."
    )
  )
}

.cov_curve_space <- function(
    r,
    sill,
    Ar,
    model_space
) {
  sill <- as.numeric(sill)
  
  if (!is.finite(sill)) {
    stop(
      "Invalid spatial sill: ",
      sill
    )
  }
  
  sill *
    .space_kernel(
      r = r,
      Ar = Ar,
      model_space = model_space
    )
}

.cov_curve_time <- function(
    t,
    sill,
    model_time,
    alpha1,
    alpha2,
    at1,
    at2,
    at
) {
  sill <- as.numeric(sill)
  
  if (!is.finite(sill)) {
    stop(
      "Invalid temporal sill: ",
      sill
    )
  }
  
  sill *
    .time_kernel(
      t = t,
      model_time = model_time,
      alpha1 = alpha1,
      alpha2 = alpha2,
      at1 = at1,
      at2 = at2,
      at = at
    )
}

.make_panel <- function(
    exp_dt,
    cov_fun,
    xlab,
    ylab,
    sse_value,
    xlim_fixed
) {
  exp_dt <- copy(exp_dt)
  
  exp_dt <- exp_dt[
    is.finite(lag) &
      is.finite(cov_exp)
  ]
  
  exp_dt[
    ,
    cov_exp := pmax(
      cov_exp,
      0
    )
  ]
  
  lag_seq <- seq(
    from = xlim_fixed[1],
    to = xlim_fixed[2],
    length.out = 400
  )
  
  modeled_values <- cov_fun(
    lag_seq
  )
  
  if (
    length(modeled_values) != length(lag_seq)
  ) {
    stop(
      "The covariance function returned an unexpected number of values."
    )
  }
  
  mod_dt <- data.table(
    lag = lag_seq,
    cov_mod = as.numeric(modeled_values)
  )
  
  mod_dt[
    ,
    cov_mod := pmax(
      cov_mod,
      0
    )
  ]
  
  ymax <- max(
    c(
      mod_dt$cov_mod,
      exp_dt$cov_exp,
      0
    ),
    na.rm = TRUE
  )
  
  if (
    !is.finite(ymax) ||
    ymax <= 0
  ) {
    ymax <- 1e-12
  } else {
    ymax <- ymax * 1.05
  }
  
  ggplot() +
    geom_line(
      data = mod_dt,
      mapping = aes(
        x = lag,
        y = cov_mod
      ),
      linewidth = 1.05
    ) +
    geom_point(
      data = exp_dt,
      mapping = aes(
        x = lag,
        y = cov_exp
      ),
      size = 2.4
    ) +
    annotate(
      geom = "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.4,
      label = paste0(
        "SSE = ",
        .fmt2(sse_value)
      ),
      size = 6,
      family = "Arial"
    ) +
    labs(
      x = xlab,
      y = ylab
    ) +
    scale_x_continuous(
      limits = xlim_fixed,
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    coord_cartesian(
      ylim = c(
        0,
        ymax
      ),
      clip = "off"
    ) +
    theme(
      plot.margin = margin(
        t = 6,
        r = 10,
        b = 6,
        l = 10
      )
    )
}

# =============================================================================
# READ MODEL TABLES
# =============================================================================

if (!file.exists(HI_MODEL_CSV)) {
  stop(
    "Heat-index model table not found: ",
    HI_MODEL_CSV
  )
}

if (!file.exists(PR_MODEL_CSV)) {
  stop(
    "Precipitation model table not found: ",
    PR_MODEL_CSV
  )
}

.log(
  "Reading heat-index model table: %s",
  HI_MODEL_CSV
)

hi_mod <- fread(
  HI_MODEL_CSV
)

.log(
  "Reading precipitation model table: %s",
  PR_MODEL_CSV
)

pr_mod <- fread(
  PR_MODEL_CSV
)

.log(
  "Raw rows: HI=%d | PR=%d",
  nrow(hi_mod),
  nrow(pr_mod)
)

# =============================================================================
# VALIDATE REQUIRED COLUMNS
# =============================================================================

base_required_columns <- c(
  "day",
  "space_file",
  "time_file",
  "sill",
  "Ar",
  "model_space",
  "model_time",
  "sse_space",
  "sse_time"
)

.require_columns(
  hi_mod,
  base_required_columns,
  "Heat-index model table"
)

.require_columns(
  pr_mod,
  base_required_columns,
  "Precipitation model table"
)

# Ensure optional time-model parameters exist so merged-column references are
# valid even when a table only uses one temporal model family.
optional_numeric_columns <- c(
  "alpha1",
  "alpha2",
  "at1",
  "at2",
  "at",
  "sse_total"
)

for (column_name in optional_numeric_columns) {
  if (!column_name %in% names(hi_mod)) {
    hi_mod[
      ,
      (column_name) := NA_real_
    ]
  }
  
  if (!column_name %in% names(pr_mod)) {
    pr_mod[
      ,
      (column_name) := NA_real_
    ]
  }
}

# =============================================================================
# CREATE day_date
# =============================================================================

# day_date is created here from the original CSV column named day.
hi_mod[
  ,
  day_date := .parse_day(day)
]

pr_mod[
  ,
  day_date := .parse_day(day)
]

.print_date_summary(
  hi_mod,
  "HI"
)

.print_date_summary(
  pr_mod,
  "PR"
)

.print_unparsed_examples(
  hi_mod,
  "HI"
)

.print_unparsed_examples(
  pr_mod,
  "PR"
)

# =============================================================================
# PARSE DATE WINDOW
# =============================================================================

.date_start <- .parse_bound(
  DATE_START,
  "DATE_START"
)

.date_end <- .parse_bound(
  DATE_END,
  "DATE_END"
)

if (
  !is.na(.date_start) &&
  !is.na(.date_end) &&
  .date_start > .date_end
) {
  stop(
    sprintf(
      "DATE_START (%s) occurs after DATE_END (%s).",
      format(
        .date_start,
        "%Y-%m-%d"
      ),
      format(
        .date_end,
        "%Y-%m-%d"
      )
    )
  )
}

window_start_label <- if (
  is.na(.date_start)
) {
  "-Inf"
} else {
  format(
    .date_start,
    "%Y-%m-%d"
  )
}

window_end_label <- if (
  is.na(.date_end)
) {
  "+Inf"
} else {
  format(
    .date_end,
    "%Y-%m-%d"
  )
}

.log(
  "Requested date window: [%s, %s]",
  window_start_label,
  window_end_label
)

# =============================================================================
# NORMALIZE FILE PATHS
# =============================================================================

hi_mod[
  ,
  space_file := .normalize_cov_path(space_file)
]

hi_mod[
  ,
  time_file := .normalize_cov_path(time_file)
]

pr_mod[
  ,
  space_file := .normalize_cov_path(space_file)
]

pr_mod[
  ,
  time_file := .normalize_cov_path(time_file)
]

# =============================================================================
# FILTER INVALID DATES
# =============================================================================

hi_before <- nrow(hi_mod)
pr_before <- nrow(pr_mod)

hi_mod <- hi_mod[
  !is.na(day_date)
]

pr_mod <- pr_mod[
  !is.na(day_date)
]

.log(
  paste0(
    "After parsed-date filter: ",
    "HI=%d/%d retained | PR=%d/%d retained"
  ),
  nrow(hi_mod),
  hi_before,
  nrow(pr_mod),
  pr_before
)

# =============================================================================
# FILTER CONVERGED MODEL ROWS
# =============================================================================

if ("converged" %in% names(hi_mod)) {
  hi_mod[
    ,
    converged_ok := .as_logical01(converged)
  ]
  
  .log(
    "HI convergence values: TRUE=%d | FALSE=%d | NA=%d",
    sum(
      hi_mod$converged_ok %in% TRUE
    ),
    sum(
      hi_mod$converged_ok %in% FALSE
    ),
    sum(
      is.na(
        hi_mod$converged_ok
      )
    )
  )
  
  hi_before <- nrow(hi_mod)
  
  hi_mod <- hi_mod[
    converged_ok %in% TRUE
  ]
  
  .log(
    "After HI convergence filter: %d/%d rows retained",
    nrow(hi_mod),
    hi_before
  )
} else {
  .log(
    "HI model table has no converged column; convergence filtering skipped"
  )
}

if ("converged" %in% names(pr_mod)) {
  pr_mod[
    ,
    converged_ok := .as_logical01(converged)
  ]
  
  .log(
    "PR convergence values: TRUE=%d | FALSE=%d | NA=%d",
    sum(
      pr_mod$converged_ok %in% TRUE
    ),
    sum(
      pr_mod$converged_ok %in% FALSE
    ),
    sum(
      is.na(
        pr_mod$converged_ok
      )
    )
  )
  
  pr_before <- nrow(pr_mod)
  
  pr_mod <- pr_mod[
    converged_ok %in% TRUE
  ]
  
  .log(
    "After PR convergence filter: %d/%d rows retained",
    nrow(pr_mod),
    pr_before
  )
} else {
  .log(
    "PR model table has no converged column; convergence filtering skipped"
  )
}

# =============================================================================
# CHECK EMPIRICAL-COVARIANCE FILES
# =============================================================================

hi_mod[
  ,
  space_exists := (
    !is.na(space_file) &
      file.exists(space_file)
  )
]

hi_mod[
  ,
  time_exists := (
    !is.na(time_file) &
      file.exists(time_file)
  )
]

pr_mod[
  ,
  space_exists := (
    !is.na(space_file) &
      file.exists(space_file)
  )
]

pr_mod[
  ,
  time_exists := (
    !is.na(time_file) &
      file.exists(time_file)
  )
]

.log(
  paste0(
    "HI covariance files: ",
    "space=%d/%d exist | time=%d/%d exist"
  ),
  sum(
    hi_mod$space_exists
  ),
  nrow(hi_mod),
  sum(
    hi_mod$time_exists
  ),
  nrow(hi_mod)
)

.log(
  paste0(
    "PR covariance files: ",
    "space=%d/%d exist | time=%d/%d exist"
  ),
  sum(
    pr_mod$space_exists
  ),
  nrow(pr_mod),
  sum(
    pr_mod$time_exists
  ),
  nrow(pr_mod)
)

.print_missing_path_examples(
  hi_mod,
  exists_column = "space_exists",
  path_column = "space_file",
  label = "HI missing spatial covariance paths"
)

.print_missing_path_examples(
  hi_mod,
  exists_column = "time_exists",
  path_column = "time_file",
  label = "HI missing temporal covariance paths"
)

.print_missing_path_examples(
  pr_mod,
  exists_column = "space_exists",
  path_column = "space_file",
  label = "PR missing spatial covariance paths"
)

.print_missing_path_examples(
  pr_mod,
  exists_column = "time_exists",
  path_column = "time_file",
  label = "PR missing temporal covariance paths"
)

hi_before <- nrow(hi_mod)
pr_before <- nrow(pr_mod)

hi_mod <- hi_mod[
  space_exists == TRUE &
    time_exists == TRUE
]

pr_mod <- pr_mod[
  space_exists == TRUE &
    time_exists == TRUE
]

.log(
  paste0(
    "After covariance-file filter: ",
    "HI=%d/%d retained | PR=%d/%d retained"
  ),
  nrow(hi_mod),
  hi_before,
  nrow(pr_mod),
  pr_before
)

# =============================================================================
# FILTER TO DATE WINDOW
# =============================================================================

if (!is.na(.date_start)) {
  hi_mod <- hi_mod[
    day_date >= .date_start
  ]
  
  pr_mod <- pr_mod[
    day_date >= .date_start
  ]
}

if (!is.na(.date_end)) {
  hi_mod <- hi_mod[
    day_date <= .date_end
  ]
  
  pr_mod <- pr_mod[
    day_date <= .date_end
  ]
}

.log(
  "Rows inside date window: HI=%d | PR=%d",
  nrow(hi_mod),
  nrow(pr_mod)
)

if (nrow(hi_mod) > 0L) {
  .log(
    "HI window range: %s to %s",
    format(
      min(
        hi_mod$day_date
      ),
      "%Y-%m-%d"
    ),
    format(
      max(
        hi_mod$day_date
      ),
      "%Y-%m-%d"
    )
  )
}

if (nrow(pr_mod) > 0L) {
  .log(
    "PR window range: %s to %s",
    format(
      min(
        pr_mod$day_date
      ),
      "%Y-%m-%d"
    ),
    format(
      max(
        pr_mod$day_date
      ),
      "%Y-%m-%d"
    )
  )
}

# =============================================================================
# CHOOSE ONE MODEL ROW PER DAY
# =============================================================================

if (
  "sse_total" %in% names(hi_mod) &&
  any(
    is.finite(
      hi_mod$sse_total
    )
  )
) {
  setorder(
    hi_mod,
    day_date,
    sse_total,
    na.last = TRUE
  )
} else {
  setorder(
    hi_mod,
    day_date
  )
}

if (
  "sse_total" %in% names(pr_mod) &&
  any(
    is.finite(
      pr_mod$sse_total
    )
  )
) {
  setorder(
    pr_mod,
    day_date,
    sse_total,
    na.last = TRUE
  )
} else {
  setorder(
    pr_mod,
    day_date
  )
}

hi_mod <- hi_mod[
  ,
  .SD[1L],
  by = day_date
]

pr_mod <- pr_mod[
  ,
  .SD[1L],
  by = day_date
]

.log(
  "Unique valid dates after deduplication: HI=%d | PR=%d",
  nrow(hi_mod),
  nrow(pr_mod)
)

# =============================================================================
# MERGE MATCHING HEAT-INDEX AND PRECIPITATION DAYS
# =============================================================================

pairs <- merge(
  x = hi_mod,
  y = pr_mod,
  by = "day_date",
  suffixes = c(
    "_hi",
    "_pr"
  ),
  all = FALSE,
  sort = TRUE
)

if (nrow(pairs) == 0L) {
  hi_dates <- sort(
    unique(
      hi_mod$day_date
    )
  )
  
  pr_dates <- sort(
    unique(
      pr_mod$day_date
    )
  )
  
  stop(
    sprintf(
      paste0(
        "No overlapping heat-index/precipitation dates found in ",
        "[%s, %s].\n",
        "HI valid dates in window: %d\n",
        "PR valid dates in window: %d\n",
        "HI first dates: %s\n",
        "PR first dates: %s"
      ),
      window_start_label,
      window_end_label,
      length(hi_dates),
      length(pr_dates),
      if (length(hi_dates) > 0L) {
        paste(
          format(
            head(
              hi_dates,
              5
            ),
            "%Y-%m-%d"
          ),
          collapse = ", "
        )
      } else {
        "none"
      },
      if (length(pr_dates) > 0L) {
        paste(
          format(
            head(
              pr_dates,
              5
            ),
            "%Y-%m-%d"
          ),
          collapse = ", "
        )
      } else {
        "none"
      }
    )
  )
}

.log(
  "Found %d overlapping dates in [%s, %s]",
  nrow(pairs),
  window_start_label,
  window_end_label
)

# =============================================================================
# RENDER FIGURES
# =============================================================================

xlab_space <- "Spatial degrees"
xlab_time <- "Hours"

n_written <- 0L
n_skipped <- 0L
n_failed <- 0L

for (i in seq_len(nrow(pairs))) {
  row <- pairs[i]
  
  day_date <- row$day_date
  
  day_label <- format(
    day_date,
    "%Y-%m-%d"
  )
  
  out_png <- file.path(
    OUT_DIR,
    paste0(
      "fig_03_cov_model_",
      day_label,
      ".png"
    )
  )
  
  if (
    SKIP_EXISTING &&
    file.exists(out_png)
  ) {
    n_skipped <- n_skipped + 1L
    
    .log(
      "[SKIP] %s already exists",
      day_label
    )
    
    next
  }
  
  ok <- tryCatch(
    {
      hi_space <- .read_exp_cov(
        row$space_file_hi
      )
      
      hi_time <- .read_exp_cov(
        row$time_file_hi
      )
      
      pr_space <- .read_exp_cov(
        row$space_file_pr
      )
      
      pr_time <- .read_exp_cov(
        row$time_file_pr
      )
      
      hi_space_fun <- function(r) {
        .cov_curve_space(
          r = r,
          sill = row$sill_hi,
          Ar = row$Ar_hi,
          model_space = row$model_space_hi
        )
      }
      
      hi_time_fun <- function(t) {
        .cov_curve_time(
          t = t,
          sill = row$sill_hi,
          model_time = row$model_time_hi,
          alpha1 = row$alpha1_hi,
          alpha2 = row$alpha2_hi,
          at1 = row$at1_hi,
          at2 = row$at2_hi,
          at = row$at_hi
        )
      }
      
      pr_space_fun <- function(r) {
        .cov_curve_space(
          r = r,
          sill = row$sill_pr,
          Ar = row$Ar_pr,
          model_space = row$model_space_pr
        )
      }
      
      pr_time_fun <- function(t) {
        .cov_curve_time(
          t = t,
          sill = row$sill_pr,
          model_time = row$model_time_pr,
          alpha1 = row$alpha1_pr,
          alpha2 = row$alpha2_pr,
          at1 = row$at1_pr,
          at2 = row$at2_pr,
          at = row$at_pr
        )
      }
      
      p_hi_space <- .make_panel(
        exp_dt = hi_space,
        cov_fun = hi_space_fun,
        xlab = xlab_space,
        ylab = expression(
          "Covariance (" * degree * C^2 * ")"
        ),
        sse_value = row$sse_space_hi,
        xlim_fixed = X_SPACE
      )
      
      p_hi_time <- .make_panel(
        exp_dt = hi_time,
        cov_fun = hi_time_fun,
        xlab = xlab_time,
        ylab = expression(
          "Covariance (" * degree * C^2 * ")"
        ),
        sse_value = row$sse_time_hi,
        xlim_fixed = X_TIME
      )
      
      p_pr_space <- .make_panel(
        exp_dt = pr_space,
        cov_fun = pr_space_fun,
        xlab = xlab_space,
        ylab = expression(
          "Covariance (mm"^2 * ")"
        ),
        sse_value = row$sse_space_pr,
        xlim_fixed = X_SPACE
      )
      
      p_pr_time <- .make_panel(
        exp_dt = pr_time,
        cov_fun = pr_time_fun,
        xlab = xlab_time,
        ylab = expression(
          "Covariance (mm"^2 * ")"
        ),
        sse_value = row$sse_time_pr,
        xlim_fixed = X_TIME
      )
      
      fig <- (
        p_hi_space |
          p_hi_time
      ) /
        (
          p_pr_space |
            p_pr_time
        ) +
        patchwork::plot_annotation(
          tag_levels = "a",
          tag_prefix = "(",
          tag_suffix = ")",
          theme = theme(
            plot.tag = element_text(
              face = "bold",
              size = 20,
              family = "Arial"
            ),
            plot.tag.position = c(
              0.02,
              0.98
            )
          )
        )
      
      ggsave(
        filename = out_png,
        plot = fig,
        width = 12,
        height = 8,
        units = "in",
        dpi = 300,
        bg = "white"
      )
      
      TRUE
    },
    error = function(e) {
      .log(
        "[FAIL] %s | %s",
        day_label,
        conditionMessage(e)
      )
      
      FALSE
    }
  )
  
  if (isTRUE(ok)) {
    n_written <- n_written + 1L
    
    .log(
      "[WROTE] %s",
      out_png
    )
  } else {
    n_failed <- n_failed + 1L
  }
}

.log(
  paste0(
    "[DONE] written=%d | skipped=%d | failed=%d | ",
    "output=%s"
  ),
  n_written,
  n_skipped,
  n_failed,
  OUT_DIR
)