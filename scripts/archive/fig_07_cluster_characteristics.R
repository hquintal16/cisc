# scripts/fig_07_plot_ifdc_curves.R
# ------------------------------------------------------------------------------
# Build a 2-panel Intensity–Frequency–Duration Curve figure:
#   - Left panel: heat index events
#   - Right panel: precipitation events
# Curves:
#   - Color = return period (1,2,5,10,25 yr)
#   - Linetype = scenario/resolution (e.g., 0.25, 0.44_max, 0.12_bilinear)
# IFDC construction (per variable × RP × scenario):
#   - Bin events by duration (quantile bins)
#   - For each bin: mean duration, mean intensity, and frequency (events/year)
#   - Plot intensity vs duration as a curve; frequency shown as point size along the curve
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(fs)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(ggplot2)
  library(scales)
})

# ----------------------------
# USER SETTINGS
# ----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

IN_DIR <- file.path(ROOT, "data/out/07_characteristics")
OUT_DIR <- file.path(ROOT, "figures/fig_07_cluster_characteristics")
OUT_PNG <- file.path(OUT_DIR, "fig_07_cluster_characteristics.png")

# Return periods to include (folders are assumed to contain these)
RPS_KEEP <- c(1, 2, 5, 10, 25)

# Which intensity column to use for curves
INTENSITY_COL <- "mean_intensity"

# Duration column (hours)
DURATION_COL <- "duration_hr"

# Number of duration bins per (variable × RP × scenario)
N_BINS <- 12

# ----------------------------
# Helpers
# ----------------------------
.parse_rp_from_folder <- function(path) {
  # expects .../hi1_pr1, hi2_pr2, hi5_pr5, hi10_pr10, hi25_pr25
  b <- basename(path)
  m <- str_match(b, "^hi(\\d+)_pr(\\d+)$")
  if (any(is.na(m))) return(NA_integer_)
  as.integer(m[,2])
}

.parse_scenario_from_cluster_set <- function(cluster_set) {
  # examples:
  #   heatindex_0.25
  #   heatindex_0.44_max
  #   precipitation_0.12_bilinear
  # keep everything after the variable prefix
  if (str_starts(cluster_set, "heatindex_")) {
    return(str_replace(cluster_set, "^heatindex_", ""))
  }
  if (str_starts(cluster_set, "precipitation_")) {
    return(str_replace(cluster_set, "^precipitation_", ""))
  }
  # fallback: after first underscore
  return(str_replace(cluster_set, "^[^_]+_", ""))
}

.safe_quantile_bins <- function(x, n_bins = 10) {
  x <- x[is.finite(x)]
  if (length(x) < 5) return(NULL)
  
  probs <- seq(0, 1, length.out = n_bins + 1)
  brks <- as.numeric(quantile(x, probs = probs, na.rm = TRUE, type = 7))
  brks <- unique(brks)
  
  # If quantiles collapse (too many ties), fall back to pretty breaks
  if (length(brks) < 4) {
    brks <- pretty(x, n = min(n_bins, 8))
    brks <- unique(brks)
  }
  if (length(brks) < 3) return(NULL)
  brks
}

.make_ifdc <- function(df, n_bins = 12) {
  # df already filtered to one variable × rp × scenario
  # returns binned table for plotting a curve + frequency encoding
  
  # parse time range for frequency denominator
  # (events/year over the span of years present in this group)
  if (!("start_time" %in% names(df))) df$start_time <- NA_character_
  st <- suppressWarnings(ymd_hms(df$start_time, quiet = TRUE))
  yrs <- year(st)
  yrs <- yrs[is.finite(yrs)]
  if (length(yrs) == 0) {
    # fallback: use 1 to avoid division by 0
    years_span <- 1
  } else {
    years_span <- max(yrs) - min(yrs) + 1
    years_span <- max(1, years_span)
  }
  
  # required columns
  dur <- df[[DURATION_COL]]
  inten <- df[[INTENSITY_COL]]
  
  ok <- is.finite(dur) & is.finite(inten) & dur > 0
  df <- df[ok, , drop = FALSE]
  if (nrow(df) < 5) return(NULL)
  
  brks <- .safe_quantile_bins(df[[DURATION_COL]], n_bins = n_bins)
  if (is.null(brks)) return(NULL)
  
  df$dur_bin <- cut(df[[DURATION_COL]], breaks = brks, include.lowest = TRUE, right = TRUE)
  
  out <- df %>%
    group_by(dur_bin) %>%
    summarise(
      duration_hr = mean(.data[[DURATION_COL]], na.rm = TRUE),
      intensity   = mean(.data[[INTENSITY_COL]], na.rm = TRUE),
      n_events    = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(
      years_span = years_span,
      freq_per_year = n_events / years_span
    ) %>%
    arrange(duration_hr)
  
  # drop bins that might be NA due to cuts
  out <- out %>% filter(is.finite(duration_hr), is.finite(intensity), n_events > 0)
  if (nrow(out) < 2) return(NULL)
  
  out
}

# ----------------------------
# Discover files
# ----------------------------
dir_create(OUT_DIR)

rp_dirs <- dir_ls(IN_DIR, type = "directory", recurse = FALSE)
rp_dirs <- rp_dirs[str_detect(basename(rp_dirs), "^hi\\d+_pr\\d+$")]
if (length(rp_dirs) == 0) stop("No RP folders found under: ", IN_DIR)

rp_tbl <- tibble(
  rp_dir = rp_dirs,
  rp_yr  = vapply(rp_dirs, .parse_rp_from_folder, integer(1))
) %>%
  filter(!is.na(rp_yr), rp_yr %in% RPS_KEEP)

if (nrow(rp_tbl) == 0) stop("No RP folders matched RPS_KEEP under: ", IN_DIR)

# read all matching CSVs from each RP folder
all_files <- unlist(lapply(rp_tbl$rp_dir, function(d) {
  dir_ls(d, recurse = FALSE, type = "file", glob = "*.csv")
}), use.names = FALSE)

# only cluster characteristics files (spelling varies in your note; keep broad)
all_files <- all_files[str_detect(basename(all_files), "^cluster_characterist|^cluster_characteristic")]

if (length(all_files) == 0) stop("No cluster characteristics CSVs found under: ", IN_DIR)

# ----------------------------
# Read + standardize
# ----------------------------
read_one <- function(f) {
  rp_yr <- .parse_rp_from_folder(path_dir(f))
  dt <- suppressWarnings(fread(f))
  if (nrow(dt) == 0) return(NULL)
  
  # Standardize expected columns
  if (!("cluster_set" %in% names(dt))) {
    # try to infer from filename if needed
    # e.g., cluster_characteristics_heatindex_0.25.csv
    bn <- basename(f)
    if (str_detect(bn, "heatindex")) dt$cluster_set <- str_replace(bn, "^.*heatindex_", "heatindex_")
    if (str_detect(bn, "precipitation")) dt$cluster_set <- str_replace(bn, "^.*precipitation_", "precipitation_")
  }
  
  # Determine variable from cluster_set
  dt <- as.data.frame(dt)
  dt$rp_yr_from_folder <- rp_yr
  dt$source_path <- f
  
  dt
}

raw_list <- lapply(all_files, read_one)
raw_list <- raw_list[!vapply(raw_list, is.null, logical(1))]
raw <- bind_rows(raw_list)

# Basic checks
need <- c("cluster_set", DURATION_COL, INTENSITY_COL)
miss <- setdiff(need, names(raw))
if (length(miss) > 0) stop("Missing required column(s): ", paste(miss, collapse = ", "))

# Derive variable + scenario + rp
dat <- raw %>%
  mutate(
    variable = case_when(
      str_starts(cluster_set, "heatindex") ~ "Heat index",
      str_starts(cluster_set, "precipitation") ~ "Precipitation",
      TRUE ~ NA_character_
    ),
    scenario = vapply(cluster_set, .parse_scenario_from_cluster_set, character(1)),
    rp_yr = if ("rp" %in% names(raw)) suppressWarnings(as.integer(rp)) else NA_integer_,
    rp_yr = ifelse(is.na(rp_yr), rp_yr_from_folder, rp_yr),
    duration_hr = suppressWarnings(as.numeric(.data[[DURATION_COL]])),
    mean_intensity = suppressWarnings(as.numeric(.data[[INTENSITY_COL]]))
  ) %>%
  filter(!is.na(variable), rp_yr %in% RPS_KEEP, is.finite(duration_hr), is.finite(mean_intensity))

if (nrow(dat) == 0) stop("After filtering, no rows remain. Check column names and file contents.")

# ----------------------------
# Build IFDC curves (binned)
# ----------------------------
curves <- dat %>%
  group_by(variable, rp_yr, scenario) %>%
  group_modify(~{
    out <- .make_ifdc(.x, n_bins = N_BINS)
    if (is.null(out)) return(tibble())
    out
  }) %>%
  ungroup() %>%
  mutate(
    rp_lab = factor(paste0(rp_yr, "-yr"), levels = paste0(RPS_KEEP, "-yr"))
  )

if (nrow(curves) == 0) stop("No curves could be constructed (insufficient events per group?).")

# Optional: define a stable linetype mapping for common scenarios
scenario_levels <- sort(unique(curves$scenario))
linetypes <- rep_len(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), length(scenario_levels))
names(linetypes) <- scenario_levels

# ----------------------------
# Plot
# ----------------------------
p <- ggplot(curves, aes(
  x = duration_hr,
  y = intensity,
  color = rp_lab,
  linetype = scenario,
  group = interaction(variable, rp_lab, scenario)
)) +
  geom_path(linewidth = 0.9, na.rm = TRUE) +
  # Frequency encoded as point size along the curve (events/year within each duration bin)
  geom_point(aes(size = freq_per_year), alpha = 0.75, na.rm = TRUE) +
  facet_wrap(~variable, nrow = 1, scales = "free_y") +
  scale_color_brewer(palette = "Dark2", name = "Return period") +
  scale_linetype_manual(values = linetypes, name = "Scenario (resolution/method)") +
  scale_size_continuous(
    name = "Frequency (events/yr)",
    range = c(1.2, 4.2),
    breaks = pretty_breaks(n = 4)
  ) +
  scale_x_continuous(
    name = "Duration (hours)",
    breaks = pretty_breaks(n = 6)
  ) +
  ylab("Mean intensity (binned)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text = element_text(face = "bold")
  )

# ----------------------------
# Write
# ----------------------------
dir_create(path_dir(OUT_PNG))
ggsave(OUT_PNG, p, width = 16, height = 7, dpi = 200)

message("[OK] Wrote: ", OUT_PNG)