# scripts/fig_04_plot_nca_stm1_stm2_4panel_DEBUGJOIN.R
# =============================================================================
# Fig 04 (NCA map): 4-panel STM1/STM2, heatindex(JJA mean) + precipitation(annual mean)
#
# Panels:
#   (a) heatindex  STM1 (JJA mean)
#   (b) precip     STM1 (annual mean)
#   (c) heatindex  STM2 (JJA mean)
#   (d) precip     STM2 (annual mean)
#
# HARD JOIN DIAGNOSTICS:
#   - Print all region names (shp vs csv, after standardization + mapping)
#   - Print missing region sets in each direction
#   - Validate every CONUS NCA region has finite stm1/stm2 after aggregation
#   - STOP if any region is NA after join
#
# Region mapping (confirmed):
#   Plains == Southern Plains  (ONLY discrepancy)
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(stringr)
  library(readr)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(scales)
  library(grid)
})

sf::sf_use_s2(FALSE)
options(na.print = "NA")

# ----------------------------
# PATHS
# ----------------------------
NCA_SHP <- file.path("V:/users/hquintal/phd2_cluster_conus/data/in/shapefiles/NCA_regions/NCA_regions.shp")

ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"
HEAT_CSV <- file.path(ROOT, "data/out/03_cov_model/nca/heatindex/model/heatindex_daily_optimal_cov_model_nca.csv")
PREC_CSV <- file.path(ROOT, "data/out/03_cov_model/nca/precipitation/model/precipitation_daily_optimal_cov_model_nca.csv")

OUT_DIR <- file.path(ROOT, "figures/fig_04_space_time_metric")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_4P <- file.path(OUT_DIR, "fig_04_nca_stm1_stm2_4panel_sharedlegend_center1.png")

# ----------------------------
# FIGURE SETTINGS
# ----------------------------
CRS_TARGET <- 5070
RES_MAIN   <- 200
WIDTH_4P   <- 5200
HEIGHT_4P  <- 3400

HEAT_MONTHS <- 6:8

DROP_REGIONS <- c("Alaska", "Caribbean", "Hawai'i and Pacific Islands", "Hawaii and Pacific Islands")

# For contrast on 10 regions, use MAX (true range) unless you want p95/p99
LIMIT_METHOD <- "max"   # "max" | "p95" | "p99"

# ----------------------------
# Helpers (match your Fig04 logic style)
# ----------------------------
.parse_day_any <- function(x) {
  suppressWarnings({
    out <- lubridate::mdy(x, quiet = TRUE)
    if (all(is.na(out))) out <- lubridate::ymd(x, quiet = TRUE)
    out
  })
}

.clean_region <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\u2019|\u2018", "'")
  x <- str_replace_all(x, "_", " ")
  x <- str_squish(x)
  x
}

# ONLY mapping you confirmed:
#   Plains == Southern Plains
.std_region_nca <- function(x) {
  x <- .clean_region(x)
  ifelse(x == "Plains", "Southern Plains", x)
}

# Shared limits around 1
.choose_lim_dev <- function(vals, method = "max") {
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(0.01)
  dev <- abs(vals - 1)
  
  method <- tolower(method)
  if (method == "max") {
    m <- max(dev, na.rm = TRUE)
    return(ifelse(is.finite(m) && m > 0, m, 0.01))
  }
  if (method == "p95") {
    q <- as.numeric(stats::quantile(dev, probs = 0.95, na.rm = TRUE, names = FALSE, type = 7))
    return(ifelse(is.finite(q) && q > 0, q, 0.01))
  }
  if (method == "p99") {
    q <- as.numeric(stats::quantile(dev, probs = 0.99, na.rm = TRUE, names = FALSE, type = 7))
    return(ifelse(is.finite(q) && q > 0, q, 0.01))
  }
  q <- as.numeric(stats::quantile(dev, probs = 0.99, na.rm = TRUE, names = FALSE, type = 7))
  ifelse(is.finite(q) && q > 0, q, 0.01)
}

# ----------------------------
# Read shapefile (CONUS only)
# ----------------------------
if (!file.exists(NCA_SHP)) stop("Missing NCA shapefile: ", NCA_SHP)
NCA <- st_read(NCA_SHP, quiet = TRUE) %>% st_make_valid()

if (!("RegionName" %in% names(NCA))) {
  stop("Shapefile missing RegionName. Columns: ", paste(names(NCA), collapse = ", "))
}

NCA <- NCA %>%
  mutate(
    RegionName_raw = as.character(RegionName),
    RegionName_std = .std_region_nca(RegionName_raw)
  ) %>%
  filter(!(RegionName_std %in% DROP_REGIONS)) %>%
  st_transform(CRS_TARGET)

# ----------------------------
# Read STM CSVs (heat + precip)
# ----------------------------
.read_stm_csv <- function(path, expected_variable) {
  if (!file.exists(path)) stop("Missing CSV: ", path)
  df <- readr::read_csv(path, show_col_types = FALSE)
  
  req <- c("variable", "region", "day", "stm1", "stm2")
  miss <- setdiff(req, names(df))
  if (length(miss)) stop("CSV missing columns: ", paste(miss, collapse = ", "), "\nFile: ", path)
  
  df <- df %>%
    mutate(
      variable = as.character(variable),
      region_raw = as.character(region),
      region_std = .std_region_nca(region_raw),
      day = .parse_day_any(as.character(day)),
      month = month(day),
      stm1 = as.numeric(stm1),
      stm2 = as.numeric(stm2)
    ) %>%
    filter(variable == expected_variable) %>%
    filter(!is.na(day))
  
  if (!nrow(df)) stop("No parsed rows for variable=", expected_variable, " in: ", path)
  df
}

heat_df <- .read_stm_csv(HEAT_CSV, "heatindex")
prec_df <- .read_stm_csv(PREC_CSV, "precipitation")

# ----------------------------
# DIAGNOSTICS: region name sets
# ----------------------------
shp_regs <- sort(unique(NCA$RegionName_std))
csv_regs_heat <- sort(unique(heat_df$region_std))
csv_regs_prec <- sort(unique(prec_df$region_std))
csv_regs_all  <- sort(unique(c(csv_regs_heat, csv_regs_prec)))

cat("\n====================\nNCA SHAPEFILE REGIONS (std)\n====================\n")
print.data.frame(data.frame(region_std = shp_regs), row.names = FALSE)

cat("\n====================\nSTM CSV REGIONS (HEAT, std)\n====================\n")
print.data.frame(data.frame(region_std = csv_regs_heat), row.names = FALSE)

cat("\n====================\nSTM CSV REGIONS (PREC, std)\n====================\n")
print.data.frame(data.frame(region_std = csv_regs_prec), row.names = FALSE)

missing_in_csv  <- setdiff(shp_regs, csv_regs_all)
missing_in_shp  <- setdiff(csv_regs_all, shp_regs)

cat("\n====================\nREGIONS IN SHP BUT NOT IN CSVs\n====================\n")
print.data.frame(data.frame(region_std = missing_in_csv), row.names = FALSE)

cat("\n====================\nREGIONS IN CSVs BUT NOT IN SHP\n====================\n")
print.data.frame(data.frame(region_std = missing_in_shp), row.names = FALSE)

# ----------------------------
# DIAGNOSTICS: row counts per region (before/after JJA filter)
# ----------------------------
heat_counts_all <- heat_df %>% count(region_std, name = "n_all_days")
heat_counts_jja <- heat_df %>% filter(month %in% HEAT_MONTHS) %>% count(region_std, name = "n_jja_days")

cat("\n====================\nHEAT rows per region (all vs JJA)\n====================\n")
print.data.frame(
  full_join(heat_counts_all, heat_counts_jja, by = "region_std") %>%
    arrange(region_std),
  row.names = FALSE
)

prec_counts_all <- prec_df %>% count(region_std, name = "n_all_days")

cat("\n====================\nPREC rows per region (annual)\n====================\n")
print.data.frame(prec_counts_all %>% arrange(region_std), row.names = FALSE)

# ----------------------------
# Aggregate STM means (match your intent)
# ----------------------------
heat_sum <- heat_df %>%
  filter(month %in% HEAT_MONTHS) %>%
  group_by(region_std) %>%
  summarise(
    stm1 = mean(stm1, na.rm = TRUE),
    stm2 = mean(stm2, na.rm = TRUE),
    .groups = "drop"
  )

prec_sum <- prec_df %>%
  group_by(region_std) %>%
  summarise(
    stm1 = mean(stm1, na.rm = TRUE),
    stm2 = mean(stm2, na.rm = TRUE),
    .groups = "drop"
  )

# Validate: every region should have finite stm1/stm2
.bad_regions <- function(sum_df, label) {
  bad <- sum_df %>%
    filter(!is.finite(stm1) | !is.finite(stm2) | is.na(stm1) | is.na(stm2)) %>%
    pull(region_std)
  if (length(bad)) {
    message("[WARN] Non-finite aggregated stm in ", label, ": ", paste(bad, collapse = ", "))
  }
  invisible(bad)
}

.bad_regions(heat_sum, "HEAT (JJA)")
.bad_regions(prec_sum, "PREC (annual)")

# ----------------------------
# Join to shapefile
# ----------------------------
nca_heat <- NCA %>% left_join(heat_sum, by = c("RegionName_std" = "region_std"))
nca_prec <- NCA %>% left_join(prec_sum, by = c("RegionName_std" = "region_std"))

# STRICT: stop if any join NA (should not happen)
miss_heat <- nca_heat %>% st_drop_geometry() %>% filter(!is.finite(stm1) | !is.finite(stm2)) %>% pull(RegionName_std)
miss_prec <- nca_prec %>% st_drop_geometry() %>% filter(!is.finite(stm1) | !is.finite(stm2)) %>% pull(RegionName_std)

if (length(miss_heat) || length(miss_prec)) {
  cat("\n====================\nJOIN FAILURE (should be empty)\n====================\n")
  if (length(miss_heat)) cat("Heat missing after join: ", paste(unique(miss_heat), collapse=", "), "\n", sep="")
  if (length(miss_prec)) cat("Precip missing after join: ", paste(unique(miss_prec), collapse=", "), "\n", sep="")
  stop("Stopping because some regions have NA/Non-finite stm after join. See diagnostics above.")
}

# ----------------------------
# Shared limits around 1 (pooled across all 4 panels)
# ----------------------------
pooled_vals <- c(nca_heat$stm1, nca_prec$stm1, nca_heat$stm2, nca_prec$stm2)
LIM_DEV <- .choose_lim_dev(pooled_vals, LIMIT_METHOD)
LOW_LIM  <- 1 - LIM_DEV
HIGH_LIM <- 1 + LIM_DEV
message(sprintf("[LIMITS] method=%s | pooled N=%d | limits=[%.4f, %.4f] centered on 1",
                LIMIT_METHOD, sum(is.finite(pooled_vals)), LOW_LIM, HIGH_LIM))

# ----------------------------
# Plot (single shared legend)
# ----------------------------
bb <- st_bbox(NCA)
XLIM <- c(bb["xmin"], bb["xmax"])
YLIM <- c(bb["ymin"], bb["ymax"])

.cb_guide <- guide_colorbar(
  title.position = "left",
  title.hjust = 0.5,
  barwidth  = unit(7.5, "in"),
  barheight = unit(0.25, "in")
)

theme_map_bottom <- theme_void(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.justification = "center",
    legend.title = element_text(size = 20),
    legend.text  = element_text(size = 20),
    plot.tag = element_text(size=20, face="bold"),
    plot.tag.position = c(0.02, 0.98)
  )

.plot_metric <- function(sf_obj, metric_col, panel_tag, show_legend = FALSE) {
  p <- ggplot() +
    geom_sf(data = sf_obj, aes(fill = .data[[metric_col]]), color = "black", linewidth = 0.55) +
    coord_sf(crs = st_crs(CRS_TARGET), xlim = XLIM, ylim = YLIM, expand = FALSE) +
    scale_fill_gradient2(
      low = "blue",
      mid = "grey80",
      high = "red",
      midpoint = 1,
      limits = c(LOW_LIM, HIGH_LIM),
      oob = scales::squish,
      na.value = "white",
      name = "STM (centered at 1)",
      guide = .cb_guide
    ) +
    labs(tag = panel_tag) +
    theme_map_bottom
  
  if (!isTRUE(show_legend)) p <- p + theme(legend.position = "none")
  p
}

p_a <- .plot_metric(nca_heat, "stm1", "(a)", show_legend = TRUE)
p_b <- .plot_metric(nca_prec, "stm1", "(b)", show_legend = FALSE)
p_c <- .plot_metric(nca_heat, "stm2", "(c)", show_legend = FALSE)
p_d <- .plot_metric(nca_prec, "stm2", "(d)", show_legend = FALSE)

fig_4p <- (p_a + p_b) / (p_c + p_d) +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      legend.position = "bottom",
      legend.box.just = "center",
      legend.justification = "center"
    )
  )

ragg::agg_png(OUT_4P, width = WIDTH_4P, height = HEIGHT_4P, res = RES_MAIN)
print(fig_4p)
dev.off()

message("Wrote: ", OUT_4P)
message("Done.")