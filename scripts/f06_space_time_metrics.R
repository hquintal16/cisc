#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
})

# =============================================================================
# Space–Time Metric figures: STM1 and STM2
#
# Effective layout: 2 rows x 4 columns
#
#                    Column 1   Column 2   Columns 3-4
#   Heat index       PDF        Month      Annual
#   Precipitation    PDF        Month      Annual
#
# Heat index:
#   May through October only
#
# Precipitation:
#   All months
# =============================================================================

# =============================================================================
# Settings
# =============================================================================

ROOT <- Sys.getenv(
  "CISC_ROOT",
  "/proj/sebastianlab/users/hquintal/cisc"
)

HEAT_MONTHS <- 1:12

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
  "figures/f06_space_time_metric"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

OUT_STM1 <- file.path(
  OUT_DIR,
  "f06_stm1.png"
)

FIG_WIDTH_IN <- 6.5
FIG_HEIGHT_IN <- 3.8
FIG_DPI <- 300

Y_LIMITS <- c(0, 5)

# Use the darkest colors from the maximum-size ColorBrewer palettes.
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
        size = 6
      ),
      
      axis.title = element_text(
        family = "Arial",
        size = 8
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
        "m/d/Y",
        "m/d/y"
      )
    ) |>
      as.Date()
  )
}

prepare_stm <- function(
    file_path,
    var_name,
    metric_col = c("stm1", "stm2")
) {
  
  metric_col <- match.arg(metric_col)
  
  dt <- readr::read_csv(
    file_path,
    show_col_types = FALSE
  )
  
  required_columns <- c(
    "day",
    metric_col
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(dt)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing required column(s) in ",
      file_path,
      ": ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  dt %>%
    mutate(
      date = parse_day_safe(
        as.character(day)
      ),
      
      year = lubridate::year(date),
      
      month_num = lubridate::month(date),
      
      month = factor(
        month.abb[month_num],
        levels = c(
          "Dec",
          "Jan",
          "Feb",
          "Mar",
          "Apr",
          "May",
          "Jun",
          "Jul",
          "Aug",
          "Sep",
          "Oct",
          "Nov"
        )
      ),
      
      variable = var_name,
      
      stm_name = metric_col,
      
      stm_val = as.numeric(
        .data[[metric_col]]
      )
    ) %>%
    filter(
      !is.na(date),
      is.finite(stm_val)
    )
}

theme_hide_y <- theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank()
)

theme_x_rotated <- theme(
  axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5,
    family = "Arial",
    size = 5.5
  )
)

theme_year_axis <- theme(
  axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5,
    family = "Arial",
    size = 4.5
  )
)

add_reference_line <- function() {
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.3
  )
}

# =============================================================================
# Figure function
# =============================================================================

make_stm_figure <- function(
    raw_data,
    metric_col = c("stm1", "stm2")
) {
  
  metric_col <- match.arg(metric_col)
  
  y_lab <- if (metric_col == "stm1") {
    expression(gamma[2])
  } else {
    expression(gamma[3])
  }
  
  dat <- raw_data %>%
    filter(
      stm_name == metric_col
    )
  
  # Heat-index data are restricted to May-October for all three panels.
  dat_heat <- dat %>%
    filter(
      variable == "HeatIndex",
      month_num %in% HEAT_MONTHS
    ) %>%
    mutate(
      month = factor(
        as.character(month),
        levels = month.abb[HEAT_MONTHS]
      )
    )
  
  # Precipitation retains the full year.
  dat_precip <- dat %>%
    filter(
      variable == "Precipitation"
    )
  
  year_levels <- sort(
    unique(dat$year)
  )
  
  dat_heat <- dat_heat %>%
    mutate(
      year_factor = factor(
        year,
        levels = year_levels
      )
    )
  
  dat_precip <- dat_precip %>%
    mutate(
      year_factor = factor(
        year,
        levels = year_levels
      )
    )
  
  # ===========================================================================
  # Heat-index row
  # ===========================================================================
  
  p_hi_daily <- ggplot(
    dat_heat,
    aes(
      y = stm_val
    )
  ) +
    geom_density(
      aes(
        x = after_stat(density)
      ),
      orientation = "y",
      fill = COLORS["HeatIndex"],
      color = "black",
      linewidth = 0.30,
      alpha = 0.70,
      na.rm = TRUE
    ) +
    add_reference_line() +
    coord_cartesian(
      ylim = Y_LIMITS
    ) +
    scale_x_continuous(
      expand = expansion(
        mult = c(0, 0.05)
      )
    ) +
    labs(
      x = "Probability density",
      y = y_lab
    )
  
  p_hi_month <- ggplot(
    dat_heat,
    aes(
      x = month,
      y = stm_val
    )
  ) +
    geom_boxplot(
      fill = COLORS["HeatIndex"],
      outlier.size = 0.20,
      linewidth = 0.30
    ) +
    add_reference_line() +
    coord_cartesian(
      ylim = Y_LIMITS
    ) +
    scale_x_discrete(
      drop = FALSE
    ) +
    labs(
      x = "Month",
      y = NULL
    ) +
    theme_hide_y +
    theme_x_rotated
  
  p_hi_annual <- ggplot(
    dat_heat,
    aes(
      x = year_factor,
      y = stm_val
    )
  ) +
    geom_boxplot(
      fill = COLORS["HeatIndex"],
      outlier.size = 0.12,
      linewidth = 0.22
    ) +
    add_reference_line() +
    coord_cartesian(
      ylim = Y_LIMITS
    ) +
    scale_x_discrete(
      drop = FALSE,
      breaks = as.character(
        seq(
          min(year_levels),
          max(year_levels),
          by = 5
        )
      )
    ) +
    labs(
      x = "Year",
      y = NULL
    ) +
    theme_hide_y +
    theme_year_axis
  
  # ===========================================================================
  # Precipitation row
  # ===========================================================================
  
  p_pr_daily <- ggplot(
    dat_precip,
    aes(
      y = stm_val
    )
  ) +
    geom_density(
      aes(
        x = after_stat(density)
      ),
      orientation = "y",
      fill = COLORS["Precipitation"],
      color = "black",
      linewidth = 0.30,
      alpha = 0.70,
      na.rm = TRUE
    ) +
    add_reference_line() +
    coord_cartesian(
      ylim = Y_LIMITS
    ) +
    scale_x_continuous(
      expand = expansion(
        mult = c(0, 0.05)
      )
    ) +
    labs(
      x = "Probability density",
      y = y_lab
    )
  
  p_pr_month <- ggplot(
    dat_precip,
    aes(
      x = month,
      y = stm_val
    )
  ) +
    geom_boxplot(
      fill = COLORS["Precipitation"],
      outlier.size = 0.20,
      linewidth = 0.30
    ) +
    add_reference_line() +
    coord_cartesian(
      ylim = Y_LIMITS
    ) +
    scale_x_discrete(
      drop = FALSE
    ) +
    labs(
      x = "Month",
      y = NULL
    ) +
    theme_hide_y +
    theme_x_rotated
  
  p_pr_annual <- ggplot(
    dat_precip,
    aes(
      x = year_factor,
      y = stm_val
    )
  ) +
    geom_boxplot(
      fill = COLORS["Precipitation"],
      outlier.size = 0.12,
      linewidth = 0.22
    ) +
    add_reference_line() +
    coord_cartesian(
      ylim = Y_LIMITS
    ) +
    scale_x_discrete(
      drop = FALSE,
      breaks = as.character(
        seq(
          min(year_levels),
          max(year_levels),
          by = 5
        )
      )
    ) +
    labs(
      x = "Year",
      y = NULL
    ) +
    theme_hide_y +
    theme_year_axis
  
  # ===========================================================================
  # Explicit 2-row x 4-column patchwork
  #
  # A B C C
  # D E F F
  #
  # Annual panels span columns 3 and 4.
  # ===========================================================================
  
  design <- "
  ABCC
  DEFF
  "
  
  p_hi_daily +
    p_hi_month +
    p_hi_annual +
    p_pr_daily +
    p_pr_month +
    p_pr_annual +
    plot_layout(
      design = design,
      widths = c(
        1,
        1,
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
      theme = theme(
        plot.tag = element_text(
          family = "Arial",
          size = 8,
          face = "bold"
        ),
        
        plot.tag.position = c(
          0.02,
          0.98
        )
      )
    )
}

# =============================================================================
# Read data
# =============================================================================

stm_heat_1 <- prepare_stm(
  HEAT_FILE,
  "HeatIndex",
  "stm1"
)

stm_heat_2 <- prepare_stm(
  HEAT_FILE,
  "HeatIndex",
  "stm2"
)

stm_prec_1 <- prepare_stm(
  PREC_FILE,
  "Precipitation",
  "stm1"
)

stm_prec_2 <- prepare_stm(
  PREC_FILE,
  "Precipitation",
  "stm2"
)

raw_data <- bind_rows(
  stm_heat_1,
  stm_heat_2,
  stm_prec_1,
  stm_prec_2
)

# =============================================================================
# Summary statistics
# =============================================================================

stm_summary <- raw_data %>%
  filter(
    variable == "Precipitation" |
      (
        variable == "HeatIndex" &
          month_num %in% HEAT_MONTHS
      )
  ) %>%
  group_by(
    variable,
    stm_name
  ) %>%
  summarise(
    n = n(),
    
    mean = mean(
      stm_val,
      na.rm = TRUE
    ),
    
    median = median(
      stm_val,
      na.rm = TRUE
    ),
    
    sd = sd(
      stm_val,
      na.rm = TRUE
    ),
    
    minimum = min(
      stm_val,
      na.rm = TRUE
    ),
    
    maximum = max(
      stm_val,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

print(stm_summary)

# =============================================================================
# Write STM1
# =============================================================================

fig_stm1 <- make_stm_figure(
  raw_data,
  "stm1"
)

ggsave(
  filename = OUT_STM1,
  plot = fig_stm1,
  width = FIG_WIDTH_IN,
  height = FIG_HEIGHT_IN,
  units = "in",
  dpi = FIG_DPI,
  bg = "white"
)

message(
  "[WROTE] ",
  OUT_STM1
)