# scripts/fig_s03_cov_model_sill_sse.R
# =============================================================================
# Fig S03: Covariance model parameter time series
# 2 x 2 layout:
#   Columns = Sill | SSE total
#   Rows    = Heat Index | Precipitation
#
# Inputs:
#   V:/users/hquintal/phd1_cluster_southeast_v2/data/out/03_cov_model/heatindex/
#     heatindex_daily_optimal_cov_model_southeast.csv
#   V:/users/hquintal/phd1_cluster_southeast_v2/data/out/03_cov_model/precipitation/
#     precipitation_daily_optimal_cov_model_southeast.csv
#
# Output:
#   V:/users/hquintal/phd1_cluster_southeast_v2/figures/fig_s03_cov_model_sill_sse/
#     fig_s03_cov_model_sill_sse.png
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
})

# -----------------------------
# Paths
# -----------------------------
ROOT <- "/proj/sebastianlab/users/hquintal/cisc"

HEAT_FILE <- file.path(
  ROOT,
  "data/out/03_cov_model/heatindex/heatindex_daily_optimal_cov_model_southeast.csv"
)

PREC_FILE <- file.path(
  ROOT,
  "data/out/03_cov_model/precipitation/precipitation_daily_optimal_cov_model_southeast.csv"
)

OUT_DIR <- file.path(ROOT, "figures/fs05_cov_model_sill_sse")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_FIG <- file.path(OUT_DIR, "fs05_cov_model_sill_sse.png")

# -----------------------------
# Helpers
# -----------------------------
parse_day_safe <- function(x) {
  suppressWarnings(parse_date_time(x, orders = c("Y-m-d", "m/d/Y", "m/d/y")) |> as.Date())
}

prepare_cov <- function(file_path, var_name) {
  dt <- read_csv(file_path, show_col_types = FALSE)
  
  req <- c("day", "sill", "sse_total")
  if (!all(req %in% names(dt))) {
    stop(
      "Missing required columns in ", file_path,
      ". Needed: ", paste(req, collapse = ", ")
    )
  }
  
  dt %>%
    mutate(
      date = parse_day_safe(as.character(day)),
      year = year(date),
      variable = var_name,
      sill = as.numeric(sill),
      sse_total = as.numeric(sse_total)
    ) %>%
    filter(!is.na(date))
}

# -----------------------------
# Read + combine
# -----------------------------
heat <- prepare_cov(HEAT_FILE, "HeatIndex")
prec <- prepare_cov(PREC_FILE, "Precipitation")

dat_all <- bind_rows(heat, prec)

# Long format for plotting
dat_long <- bind_rows(
  dat_all %>%
    transmute(
      variable,
      date,
      metric = "Sill",
      value = sill
    ),
  dat_all %>%
    transmute(
      variable,
      date,
      metric = "SSE",
      value = sse_total
    )
) %>%
  filter(is.finite(value))

# Global Sill limits across BOTH variables
sill_limits <- dat_long %>%
  filter(metric == "Sill") %>%
  summarise(
    ymin = min(value, na.rm = TRUE),
    ymax = max(value, na.rm = TRUE)
  )

# Keep global y-limits by metric across both variables
ylims <- dat_long %>%
  group_by(metric) %>%
  summarise(
    ymin = min(value, na.rm = TRUE),
    ymax = max(value, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# Styling
# -----------------------------
COLORS <- c(
  "HeatIndex" = RColorBrewer::brewer.pal(
    n = 9,
    name = "Reds"
  )[6],
  "Precipitation" = RColorBrewer::brewer.pal(
    n = 9,
    name = "Blues"
  )[6]
)

pt_colors <- c(
  "HeatIndex" = scales::alpha(
    COLORS[["HeatIndex"]],
    0.05
  ),
  "Precipitation" = scales::alpha(
    COLORS[["Precipitation"]],
    0.05
  )
)

line_colors <- COLORS

line_colors <- c(
  "HeatIndex"     = "black",
  "Precipitation" = "black"
)

base_theme <- theme_bw() +
  theme(
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.tag = element_text(size = 16, face = "plain"),
    plot.margin = margin(6, 6, 6, 6)
  )

make_panel <- function(data, var_name, metric_name, ylab = NULL,
                       show_x = TRUE, show_y = TRUE, log_y = FALSE) {
  dsub <- data %>%
    filter(variable == var_name, metric == metric_name)
  
  if (metric_name == "SSE") {
    dsub <- dsub %>%
      filter(is.finite(value), value > 0, value <= 100)
  } else {
    dsub <- dsub %>%
      filter(is.finite(value))
  }
  
  p <- ggplot(dsub, aes(x = date, y = value)) +
    geom_point(
      color = pt_colors[[var_name]],
      size = 2.2
    ) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = TRUE,
      color = line_colors[[var_name]],
      fill = scales::alpha(
        line_colors[[var_name]],
        0.20
      ),
      linewidth = 1.1,
      span = 0.35
    ) +
    scale_x_date(
      date_breaks = "20 years",
      date_labels = "%Y",
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      x = if (show_x) "Date" else NULL,
      y = ylab
    ) +
    base_theme
  
  if (metric_name == "SSE") {
    if (log_y) {
      p <- p +
        scale_y_log10(
          limits = c(1e-12, 100),   # match your actual minimum scale
          breaks = c(1e-12, 1e-9, 1e-6, 1e-3, 1e-1, 1, 10, 100),
          labels = scales::label_scientific()
        )
    } else {
      p <- p + coord_cartesian(ylim = c(0, 100))
    }
  }
  
  if (metric_name == "Sill") {
    p <- p + coord_cartesian(
      ylim = c(sill_limits$ymin, sill_limits$ymax)
    )
  }
  
  if (!show_x) {
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  
  if (!show_y) {
    p <- p + theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  
  p
}

# -----------------------------
# Build figure
# Layout:
#   columns = HeatIndex | Precipitation
#   rows    = Sill | SSE
# -----------------------------
p_a <- make_panel(
  data = dat_long,
  var_name = "HeatIndex",
  metric_name = "Sill",
  ylab = "Sill",
  show_x = FALSE,
  show_y = TRUE
)

p_b <- make_panel(
  data = dat_long,
  var_name = "Precipitation",
  metric_name = "Sill",
  ylab = NULL,
  show_x = FALSE,
  show_y = TRUE
)

p_c <- make_panel(
  data = dat_long,
  var_name = "HeatIndex",
  metric_name = "SSE",
  ylab = "SSE",
  show_x = TRUE,
  show_y = TRUE,
  log_y = TRUE
)

p_d <- make_panel(
  data = dat_long,
  var_name = "Precipitation",
  metric_name = "SSE",
  ylab = NULL,
  show_x = TRUE,
  show_y = TRUE,
  log_y = TRUE
)

fig <- (p_a + p_b) / (p_c + p_d) +
  plot_annotation(tag_levels = "a")

# -----------------------------
# Write figure
# -----------------------------
ggsave(
  filename = OUT_FIG,
  plot = fig,
  width = 12,
  height = 8,
  dpi = 300
)

message("[WROTE] ", OUT_FIG)
