#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

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
  "figures"
)

OUT_FIG <- file.path(
  OUT_DIR,
  "fig_s07_stm2_exceedance_days.png"
)

OUT_HEAT_CSV <- file.path(
  OUT_DIR,
  "fig_s07_heatindex_stm2_exceedance_days.csv"
)

OUT_PREC_CSV <- file.path(
  OUT_DIR,
  "fig_s07_precipitation_stm2_exceedance_days.csv"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

parse_day_safe <- function(x) {
  suppressWarnings(
    parse_date_time(
      x,
      orders = c("Y-m-d", "m/d/Y", "m/d/y")
    ) |>
      as.Date()
  )
}

read_stm2 <- function(file, variable, months) {
  
  read_csv(
    file,
    show_col_types = FALSE
  ) %>%
    transmute(
      variable = variable,
      date = parse_day_safe(as.character(day)),
      year = year(date),
      month = month(date),
      stm2 = as.numeric(stm2)
    ) %>%
    filter(
      !is.na(date),
      date >= as.Date("2002-01-01"),
      month %in% months,
      is.finite(stm2)
    )
}

heat <- read_stm2(
  HEAT_FILE,
  variable = "Heat index",
  months = 6:8
)

precip <- read_stm2(
  PREC_FILE,
  variable = "Precipitation",
  months = 1:12
)

find_exceedances <- function(
    dat,
    direction = c("upper", "lower")
) {
  
  direction <- match.arg(direction)
  
  stm2_mean <- mean(
    dat$stm2,
    na.rm = TRUE
  )
  
  stm2_sd <- sd(
    dat$stm2,
    na.rm = TRUE
  )
  
  threshold <- if (direction == "upper") {
    stm2_mean + stm2_sd
  } else {
    stm2_mean - stm2_sd
  }
  
  exceedances <- if (direction == "upper") {
    
    dat %>%
      filter(
        stm2 > threshold
      ) %>%
      mutate(
        stm2_mean = stm2_mean,
        stm2_sd = stm2_sd,
        threshold = threshold,
        threshold_direction = "greater than mean + 1 SD",
        distance_beyond_threshold = stm2 - threshold
      )
    
  } else {
    
    dat %>%
      filter(
        stm2 < threshold
      ) %>%
      mutate(
        stm2_mean = stm2_mean,
        stm2_sd = stm2_sd,
        threshold = threshold,
        threshold_direction = "less than mean - 1 SD",
        distance_beyond_threshold = threshold - stm2
      )
  }
  
  list(
    mean = stm2_mean,
    sd = stm2_sd,
    threshold = threshold,
    direction = direction,
    exceedances = exceedances
  )
}

heat_result <- find_exceedances(
  heat,
  direction = "upper"
)

precip_result <- find_exceedances(
  precip,
  direction = "lower"
)

write_csv(
  heat_result$exceedances,
  OUT_HEAT_CSV
)

write_csv(
  precip_result$exceedances,
  OUT_PREC_CSV
)

message(
  "Heat-index STM2 threshold: ",
  signif(heat_result$threshold, 6),
  " = mean ",
  signif(heat_result$mean, 6),
  " + SD ",
  signif(heat_result$sd, 6)
)

message(
  "Heat-index exceedance days: ",
  nrow(heat_result$exceedances)
)

message(
  "Precipitation STM2 lower threshold: ",
  signif(precip_result$threshold, 6),
  " = mean ",
  signif(precip_result$mean, 6),
  " - SD ",
  signif(precip_result$sd, 6)
)

message(
  "Precipitation days below mean - 1 SD: ",
  nrow(precip_result$exceedances)
)

all_years <- tibble(
  year = seq(
    2002,
    max(c(heat$year, precip$year), na.rm = TRUE)
  )
)

annual_counts <- function(result) {
  
  result$exceedances %>%
    count(
      year,
      name = "exceedance_days"
    ) %>%
    right_join(
      all_years,
      by = "year"
    ) %>%
    mutate(
      exceedance_days = replace_na(
        exceedance_days,
        0L
      )
    ) %>%
    arrange(year)
}

heat_annual <- annual_counts(
  heat_result
)

precip_annual <- annual_counts(
  precip_result
)

y_max <- max(
  heat_annual$exceedance_days,
  precip_annual$exceedance_days,
  na.rm = TRUE
)

y_breaks <- pretty(
  c(0, y_max),
  n = 5
)

theme_set(
  theme_bw(
    base_size = 8,
    base_family = "Calibri"
  ) +
    theme(
      text = element_text(
        family = "Calibri",
        size = 8
      ),
      axis.text = element_text(
        family = "Calibri",
        size = 8,
        color = "black"
      ),
      axis.title = element_text(
        family = "Calibri",
        size = 8,
        color = "black"
      ),
      panel.grid.minor = element_blank(),
      plot.margin = margin(
        4,
        4,
        4,
        4,
        unit = "pt"
      )
    )
)


p_heat <- ggplot(
  heat_annual,
  aes(
    x = year,
    y = exceedance_days
  )
) +
  scale_y_continuous(
    limits = c(0, max(y_breaks)),
    breaks = y_breaks,
    expand = expansion(
      mult = c(0, 0.03)
    )
  ) +
  geom_col(
    fill = "red3",
    width = 0.9
  ) +
  scale_x_continuous(
    breaks = seq(
      2002,
      max(heat_annual$year),
      by = 4
    ),
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  labs(
    x = "Year",
    y = "STM2 exceedance days"
  )

p_precip <- ggplot(
  precip_annual,
  aes(
    x = year,
    y = exceedance_days
  )
) +
  scale_y_continuous(
    limits = c(0, max(y_breaks)),
    breaks = y_breaks,
    expand = expansion(
      mult = c(0, 0.03)
    )
  ) +
  geom_col(
    fill = "blue3",
    width = 0.9
  ) +
  scale_x_continuous(
    breaks = seq(
      2002,
      max(precip_annual$year),
      by = 4
    ),
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  labs(
    x = "Year",
    y = NULL
  )

fig <- (
  p_heat |
    p_precip
) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(
        family = "Calibri",
        size = 8,
        face = "bold"
      ),
      plot.tag.position = c(
        0.02,
        0.98
      )
    )
  )

ggsave(
  filename = OUT_FIG,
  plot = fig,
  width = 6.5,
  height = 3,
  units = "in",
  dpi = 300,
  bg = "white"
)

message("[WROTE] ", OUT_FIG)
message("[WROTE] ", OUT_HEAT_CSV)
message("[WROTE] ", OUT_PREC_CSV)