#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(fs)
  library(scales)
})

# =============================================================================
# Supplemental warning–Ar correlation figures
#
# Heat-index figure:
#   Combined HI + EH warnings
#
# Precipitation figure:
#   Combined FF + TR + SV warnings
#
# Each figure:
#   Left:  daily active-warning count versus Ar
#   Right: daily active-warning area versus Ar
#
# Study period:
#   2002-01-01 through 2023-12-31
# =============================================================================

ROOT <- Sys.getenv(
  "CISC_ROOT",
  "/proj/sebastianlab/users/hquintal/cisc"
)

WARNING_RESULT_DIR <- file.path(
  ROOT,
  "data/out/10_mesonet/stm2_correlations"
)

HEAT_COV_FILE <- file.path(
  ROOT,
  "data/out/03_cov_model/heatindex",
  "heatindex_daily_optimal_cov_model_southeast.csv"
)

PRECIP_COV_FILE <- file.path(
  ROOT,
  "data/out/03_cov_model/precipitation",
  "precipitation_daily_optimal_cov_model_southeast.csv"
)

OUT_DATA_DIR <- file.path(
  ROOT,
  "data/out/10_mesonet/ar_correlations"
)

OUT_FIG_DIR <- file.path(
  ROOT,
  "figures/fig_s10_mesonet_ar"
)

START_DATE <- as.Date("2002-01-01")
END_DATE <- as.Date("2023-12-31")

dir_create(
  OUT_DATA_DIR,
  recurse = TRUE
)

dir_create(
  OUT_FIG_DIR,
  recurse = TRUE
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

identify_column <- function(
    data,
    candidates,
    description
) {
  matched <- candidates[
    candidates %in% names(data)
  ]
  
  if (length(matched) == 0L) {
    stop(
      "Could not identify ",
      description,
      ". Available columns: ",
      paste(
        names(data),
        collapse = ", "
      )
    )
  }
  
  matched[[1]]
}

read_ar <- function(
    file,
    variable_name
) {
  if (!file_exists(file)) {
    stop(
      "Missing ",
      variable_name,
      " covariance-model file: ",
      file
    )
  }
  
  data <- read_csv(
    file,
    show_col_types = FALSE
  )
  
  date_col <- identify_column(
    data,
    c(
      "date",
      "Date",
      "DATE",
      "day"
    ),
    paste0(
      variable_name,
      " date column"
    )
  )
  
  ar_col <- identify_column(
    data,
    c(
      "Ar",
      "AR",
      "ar",
      "A_r"
    ),
    paste0(
      variable_name,
      " Ar column"
    )
  )
  
  data %>%
    transmute(
      date = as.Date(
        .data[[date_col]]
      ),
      Ar = as.numeric(
        .data[[ar_col]]
      )
    ) %>%
    filter(
      date >= START_DATE,
      date <= END_DATE,
      is.finite(Ar)
    ) %>%
    group_by(date) %>%
    summarise(
      Ar = mean(
        Ar,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
}

read_warning_result <- function(
    phenom
) {
  file <- file.path(
    WARNING_RESULT_DIR,
    paste0(
      "warning_stm2_",
      tolower(phenom),
      ".csv"
    )
  )
  
  if (!file_exists(file)) {
    stop(
      "Missing warning-result file: ",
      file
    )
  }
  
  read_csv(
    file,
    show_col_types = FALSE
  ) %>%
    transmute(
      date = as.Date(date),
      event_count = as.numeric(
        event_count
      ),
      area_km2 = as.numeric(
        area_km2
      )
    ) %>%
    filter(
      date >= START_DATE,
      date <= END_DATE
    ) %>%
    group_by(date) %>%
    summarise(
      event_count = sum(
        event_count,
        na.rm = TRUE
      ),
      area_km2 = sum(
        area_km2,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
}

combine_warning_types <- function(
    phenomena
) {
  full_dates <- tibble(
    date = seq.Date(
      START_DATE,
      END_DATE,
      by = "day"
    )
  )
  
  warning_data <- map_dfr(
    phenomena,
    read_warning_result,
    .id = "source"
  ) %>%
    group_by(date) %>%
    summarise(
      event_count = sum(
        event_count,
        na.rm = TRUE
      ),
      area_km2 = sum(
        area_km2,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  full_dates %>%
    left_join(
      warning_data,
      by = "date"
    ) %>%
    mutate(
      event_count = coalesce(
        event_count,
        0
      ),
      area_km2 = coalesce(
        area_km2,
        0
      )
    )
}

correlation_label <- function(
    data,
    response
) {
  complete <- data %>%
    filter(
      is.finite(Ar),
      is.finite(
        .data[[response]]
      )
    )
  
  result <- cor.test(
    complete$Ar,
    complete[[response]],
    method = "pearson"
  )
  
  paste0(
    "r = ",
    number(
      unname(result$estimate),
      accuracy = 0.01
    ),
    "\n",
    "p ",
    pvalue(
      result$p.value,
      accuracy = 0.001,
      add_p = TRUE
    ),
    "\n",
    "n = ",
    comma(
      nrow(complete)
    )
  )
}

make_panel <- function(
    data,
    response,
    x_label,
    y_label,
    panel_label
) {
  ggplot(
    data,
    aes(
      x = Ar,
      y = .data[[response]]
    )
  ) +
    geom_point(
      shape = 16,
      size = 0.6,
      alpha = 0.25,
      na.rm = TRUE
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.7,
      na.rm = TRUE
    ) +
    annotate(
      geom = "label",
      x = -Inf,
      y = Inf,
      label = correlation_label(
        data,
        response
      ),
      hjust = -0.08,
      vjust = 1.08,
      size = 2.5,
      linewidth = 0.25,
      fill = "white",
      color = "black"
    ) +
    labs(
      x = x_label,
      y = y_label,
      tag = panel_label
    ) +
    theme_bw(
      base_size = 8
    ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        linewidth = 0.25
      ),
      axis.title = element_text(
        size = 8
      ),
      axis.text = element_text(
        size = 7
      ),
      plot.tag = element_text(
        size = 9,
        face = "bold"
      ),
      plot.tag.position = c(
        0.01,
        0.99
      )
    )
}

make_figure <- function(
    warning_data,
    ar_data,
    title,
    x_label,
    output_name
) {
  plot_data <- warning_data %>%
    inner_join(
      ar_data,
      by = "date"
    ) %>%
    filter(
      is.finite(Ar),
      event_count > 0,
      area_km2 > 0
    )
  
  write_csv(
    plot_data,
    file.path(
      OUT_DATA_DIR,
      paste0(
        output_name,
        ".csv"
      )
    )
  )
  
  count_plot <- make_panel(
    data = plot_data,
    response = "event_count",
    x_label = x_label,
    y_label = "Active warnings per day",
    panel_label = "a)"
  )
  
  area_plot <- make_panel(
    data = plot_data,
    response = "area_km2",
    x_label = x_label,
    y_label = expression(
      "Active warning area per day (km"^2*")"
    ),
    panel_label = "b)"
  )
  
  figure <- (
    count_plot |
      area_plot
  ) +
    plot_annotation(
      title = title,
      subtitle = "Daily values, 2002–2023",
      theme = theme(
        plot.title = element_text(
          size = 10,
          face = "bold",
          hjust = 0.5
        ),
        plot.subtitle = element_text(
          size = 8,
          hjust = 0.5
        )
      )
    )
  
  output_file <- file.path(
    OUT_FIG_DIR,
    paste0(
      output_name,
      ".png"
    )
  )
  
  ggsave(
    filename = output_file,
    plot = figure,
    width = 6.5,
    height = 3.5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  message(
    "Wrote: ",
    output_file
  )
}

# -----------------------------------------------------------------------------
# Read and combine warning results
# -----------------------------------------------------------------------------

heat_warnings <- combine_warning_types(
  c(
    "HI",
    "EH"
  )
)

precip_warnings <- combine_warning_types(
  c(
    "FF",
    "TR",
    "SV"
  )
)

# -----------------------------------------------------------------------------
# Read Ar
# -----------------------------------------------------------------------------

heat_ar <- read_ar(
  HEAT_COV_FILE,
  "Heat-index"
)

precip_ar <- read_ar(
  PRECIP_COV_FILE,
  "Precipitation"
)

# -----------------------------------------------------------------------------
# Create figures
# -----------------------------------------------------------------------------

make_figure(
  warning_data = heat_warnings,
  ar_data = heat_ar,
  title = "Heat and excessive heat warnings",
  x_label = expression(
    "Heat-index "*A[r]
  ),
  output_name = "fig_s10_mesonet_ar_heatindex"
)

make_figure(
  warning_data = precip_warnings,
  ar_data = precip_ar,
  title = paste(
    "Flash flood, tropical storm,",
    "and severe thunderstorm warnings"
  ),
  x_label = expression(
    "Precipitation "*A[r]
  ),
  output_name = "fig_s10_mesonet_ar_precipitation"
)