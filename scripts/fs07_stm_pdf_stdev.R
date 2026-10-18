#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

# =============================================================================
# Figure FS07: STM probability-density functions with standard deviation
#
# Layout:
#
#                    Heat index        Precipitation
#   gamma[2]         a)                 b)
#   gamma[3]         c)                 d)
#
# Heat index and precipitation use the same plotting format as the STM figures:
#   - vertical STM axis
#   - horizontal probability-density axis
#   - dashed horizontal reference line at STM = 1
#   - ColorBrewer Reds and Blues fills
#
# Summary lines:
#   - solid horizontal line: sample mean
#   - dotted horizontal lines: mean - 1 SD and mean + 1 SD
#
# Heat index:
#   all months listed in HEAT_MONTHS
#
# Precipitation:
#   all months
# =============================================================================

# =============================================================================
# Settings
# =============================================================================

ROOT <- Sys.getenv(
  "CISC_ROOT",
  "/proj/sebastianlab/users/hquintal/cisc"
)

HEAT_FILE <- file.path(
  ROOT,
  "data/out/03_cov_model/heatindex",
  "heatindex_daily_optimal_cov_model_southeast.csv"
)

PREC_FILE <- file.path(
  ROOT,
  "data/out/03_cov_model/precipitation",
  "precipitation_daily_optimal_cov_model_southeast.csv"
)

OUT_DIR <- file.path(
  ROOT,
  "figures",
  "fs07_stm_pdf_stdev"
)

OUT_FILE <- file.path(
  OUT_DIR,
  "fs07_stm_pdf_stdev.png"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

FIG_WIDTH_IN <- 6.5
FIG_HEIGHT_IN <- 5.0
FIG_DPI <- 300

# Retain all heat-index months, matching the current STM plotting script.
# Change this vector if a seasonal restriction is needed later.
HEAT_MONTHS <- 1:12

Y_LIMITS <- c(
  0,
  5
)

# Same ColorBrewer selections used by the STM plotting script.
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

# =============================================================================
# Theme
# =============================================================================

theme_set(
  theme_bw(
    base_size = 8,
    base_family = "Arial"
  ) +
    theme(
      text = element_text(
        family = "Arial",
        size = 8
      ),
      
      axis.text = element_text(
        family = "Arial",
        size = 6,
        color = "black"
      ),
      
      axis.title = element_text(
        family = "Arial",
        size = 8,
        color = "black"
      ),
      
      panel.grid.minor = element_blank(),
      
      plot.margin = margin(
        t = 3,
        r = 3,
        b = 3,
        l = 3,
        unit = "pt"
      )
    )
)

# =============================================================================
# Helpers
# =============================================================================

parse_day_safe <- function(x) {
  
  suppressWarnings(
    parse_date_time(
      x,
      orders = c(
        "Y-m-d",
        "Ymd",
        "m/d/Y",
        "m/d/y"
      )
    ) |>
      as.Date()
  )
}

prepare_stm <- function(
    file_path,
    variable_name,
    metric_col = c(
      "stm1",
      "stm2"
    )
) {
  
  metric_col <- match.arg(
    metric_col
  )
  
  if (!file.exists(file_path)) {
    stop(
      "Input file does not exist:\n",
      file_path
    )
  }
  
  dat <- readr::read_csv(
    file_path,
    show_col_types = FALSE
  )
  
  required_columns <- c(
    "day",
    metric_col
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(dat)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing required column(s) in:\n",
      file_path,
      "\n\nMissing:\n",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  dat %>%
    transmute(
      date = parse_day_safe(
        as.character(day)
      ),
      
      year = lubridate::year(
        date
      ),
      
      month_num = lubridate::month(
        date
      ),
      
      variable = variable_name,
      
      stm_name = metric_col,
      
      stm_val = suppressWarnings(
        as.numeric(
          .data[[metric_col]]
        )
      )
    ) %>%
    filter(
      !is.na(date),
      is.finite(stm_val)
    ) %>%
    arrange(
      date
    )
}

calculate_summary <- function(dat) {
  
  if (nrow(dat) < 2L) {
    stop(
      "At least two finite observations are required to calculate ",
      "a probability density and standard deviation."
    )
  }
  
  dat %>%
    summarise(
      n = n(),
      
      mean = mean(
        stm_val,
        na.rm = TRUE
      ),
      
      sd = sd(
        stm_val,
        na.rm = TRUE
      )
    ) %>%
    mutate(
      lower_sd = mean - sd,
      upper_sd = mean + sd
    )
}

add_reference_line <- function() {
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.35,
    color = "black"
  )
}

add_summary_lines <- function(
    summary_dat
) {
  
  list(
    geom_hline(
      data = summary_dat,
      aes(
        yintercept = mean
      ),
      inherit.aes = FALSE,
      linetype = "solid",
      linewidth = 0.35,
      color = "black"
    ),
    
    geom_hline(
      data = summary_dat,
      aes(
        yintercept = lower_sd
      ),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = 0.35,
      color = "black"
    ),
    
    geom_hline(
      data = summary_dat,
      aes(
        yintercept = upper_sd
      ),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = 0.35,
      color = "black"
    )
  )
}

make_density_panel <- function(
    dat,
    summary_dat,
    fill_color,
    y_label,
    show_y_axis = TRUE
) {
  
  p <- ggplot(
    dat,
    aes(
      y = stm_val
    )
  ) +
    geom_density(
      aes(
        x = after_stat(density)
      ),
      orientation = "y",
      fill = fill_color,
      color = "black",
      linewidth = 0.30,
      alpha = 0.70,
      na.rm = TRUE
    ) +
    add_reference_line() +
    add_summary_lines(
      summary_dat
    ) +
    coord_cartesian(
      ylim = Y_LIMITS
    ) +
    scale_x_continuous(
      expand = expansion(
        mult = c(
          0,
          0.05
        )
      )
    ) +
    labs(
      x = "Probability density",
      y = y_label
    )
  
  if (!show_y_axis) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  p
}

print_summary <- function(
    label,
    summary_dat
) {
  
  message(
    label,
    ": n = ",
    summary_dat$n,
    "; mean = ",
    signif(
      summary_dat$mean,
      6
    ),
    "; SD = ",
    signif(
      summary_dat$sd,
      6
    ),
    "; mean - SD = ",
    signif(
      summary_dat$lower_sd,
      6
    ),
    "; mean + SD = ",
    signif(
      summary_dat$upper_sd,
      6
    )
  )
}

# =============================================================================
# Read and prepare data
# =============================================================================

heat_stm1 <- prepare_stm(
  file_path = HEAT_FILE,
  variable_name = "HeatIndex",
  metric_col = "stm1"
) %>%
  filter(
    month_num %in% HEAT_MONTHS
  )

heat_stm2 <- prepare_stm(
  file_path = HEAT_FILE,
  variable_name = "HeatIndex",
  metric_col = "stm2"
) %>%
  filter(
    month_num %in% HEAT_MONTHS
  )

precip_stm1 <- prepare_stm(
  file_path = PREC_FILE,
  variable_name = "Precipitation",
  metric_col = "stm1"
)

precip_stm2 <- prepare_stm(
  file_path = PREC_FILE,
  variable_name = "Precipitation",
  metric_col = "stm2"
)

# =============================================================================
# Summary statistics
# =============================================================================

heat_stm1_summary <- calculate_summary(
  heat_stm1
)

heat_stm2_summary <- calculate_summary(
  heat_stm2
)

precip_stm1_summary <- calculate_summary(
  precip_stm1
)

precip_stm2_summary <- calculate_summary(
  precip_stm2
)

print_summary(
  "Heat-index gamma[2]",
  heat_stm1_summary
)

print_summary(
  "Precipitation gamma[2]",
  precip_stm1_summary
)

print_summary(
  "Heat-index gamma[3]",
  heat_stm2_summary
)

print_summary(
  "Precipitation gamma[3]",
  precip_stm2_summary
)

# =============================================================================
# Panels
#
# Row 1: gamma[2]
# Row 2: gamma[3]
#
# Column 1: Heat index
# Column 2: Precipitation
# =============================================================================

# 1. Compute the empirical density of your data vector 'r'
# dens <- density(heat_stm1$stm_val)
# max(dens$x)

p_heat_gamma2 <- make_density_panel(
  dat = heat_stm1,
  summary_dat = heat_stm1_summary,
  fill_color = COLORS["HeatIndex"],
  y_label = expression(gamma[2]),
  show_y_axis = TRUE
)

p_heat_gamma3 <- make_density_panel(
  dat = heat_stm2,
  summary_dat = heat_stm2_summary,
  fill_color = COLORS["HeatIndex"],
  y_label = expression(gamma[3]),
  show_y_axis = TRUE
)

p_precip_gamma2 <- make_density_panel(
  dat = precip_stm1,
  summary_dat = precip_stm1_summary,
  fill_color = COLORS["Precipitation"],
  y_label = expression(gamma[2]),
  show_y_axis = TRUE
)

p_precip_gamma3 <- make_density_panel(
  dat = precip_stm2,
  summary_dat = precip_stm2_summary,
  fill_color = COLORS["Precipitation"],
  y_label = expression(gamma[3]),
  show_y_axis = TRUE
)

# =============================================================================
# Combine and save
# =============================================================================

fig <- (
  p_heat_gamma2 |
    p_heat_gamma3
) / (
  p_precip_gamma2 |
    p_precip_gamma3
) +
  plot_layout(
    widths = c(
      1,
      1
    ),
    heights = c(
      1,
      1
    )
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(
        family = "Arial",
        size = 8,
        face = "bold",
        color = "black"
      ),
      
      plot.tag.position = c(
        0.02,
        0.98
      )
    )
  )

ggsave(
  filename = OUT_FILE,
  plot = fig,
  width = FIG_WIDTH_IN,
  height = FIG_HEIGHT_IN,
  units = "in",
  dpi = FIG_DPI,
  bg = "white"
)

message(
  "[WROTE] ",
  OUT_FILE
)