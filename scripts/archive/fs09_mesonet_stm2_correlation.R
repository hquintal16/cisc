# #!/usr/bin/env Rscript
# 
# suppressPackageStartupMessages({
#   library(readr)
#   library(dplyr)
#   library(lubridate)
#   library(ggplot2)
#   library(patchwork)
#   library(fs)
#   library(scales)
# })
# 
# # =============================================================================
# # Supplemental warning–STM2 correlation figures
# #
# # Heat-index STM2:
# #   HI — Heat
# #   EH — Excessive Heat
# #
# # Precipitation STM2:
# #   FF — Flash Flood
# #   TR — Tropical Storm
# #   SV — Severe Thunderstorm
# #
# # Each figure:
# #   Left:  daily active-warning count versus STM2
# #   Right: daily active-warning area versus STM2
# #
# # Study period:
# #   2002-01-01 through 2023-12-31
# # =============================================================================
# 
# ROOT <- Sys.getenv(
#   "CISC_ROOT",
#   "/proj/sebastianlab/users/hquintal/cisc"
# )
# 
# WARNING_FILE <- file.path(
#   ROOT,
#   "data/out/10_mesonet",
#   "mesonet_warning_daily_southeast_2002_2023.csv"
# )
# 
# HEAT_FILE <- file.path(
#   ROOT,
#   "data/out/03_cov_model/heatindex",
#   "heatindex_daily_optimal_cov_model_southeast.csv"
# )
# 
# PRECIP_FILE <- file.path(
#   ROOT,
#   "data/out/03_cov_model/precipitation",
#   "precipitation_daily_optimal_cov_model_southeast.csv"
# )
# 
# OUT_DATA_DIR <- file.path(
#   ROOT,
#   "data/out/10_mesonet/stm2_correlations"
# )
# 
# OUT_FIG_DIR <- file.path(
#   ROOT,
#   "figures/fig_s09_mesonet_stm2"
# )
# 
# START_DATE <- as.Date("2002-01-01")
# END_DATE <- as.Date("2023-12-31")
# 
# dir_create(
#   OUT_DATA_DIR,
#   recurse = TRUE
# )
# 
# dir_create(
#   OUT_FIG_DIR,
#   recurse = TRUE
# )
# 
# # -----------------------------------------------------------------------------
# # Helpers
# # -----------------------------------------------------------------------------
# 
# identify_column <- function(
#     data,
#     candidates,
#     description
# ) {
#   match <- candidates[
#     candidates %in% names(data)
#   ]
#   
#   if (length(match) == 0L) {
#     stop(
#       "Could not identify ",
#       description,
#       ". Available columns: ",
#       paste(
#         names(data),
#         collapse = ", "
#       )
#     )
#   }
#   
#   match[[1]]
# }
# 
# standardize_stm2 <- function(
#     file,
#     metric_name
# ) {
#   if (!file_exists(file)) {
#     stop(
#       "Missing STM2 file: ",
#       file
#     )
#   }
#   
#   data <- read_csv(
#     file,
#     show_col_types = FALSE
#   )
#   
#   date_col <- identify_column(
#     data,
#     c(
#       "date",
#       "Date",
#       "DATE",
#       "day",
#       "datetime"
#     ),
#     paste0(
#       metric_name,
#       " date column"
#     )
#   )
#   
#   stm2_col <- identify_column(
#     data,
#     c(
#       "stm2",
#       "STM2",
#       "space_time_metric_2",
#       "space_time_metric2",
#       "space_time_metric",
#       "metric"
#     ),
#     paste0(
#       metric_name,
#       " STM2 column"
#     )
#   )
#   
#   data %>%
#     transmute(
#       date = as.Date(
#         .data[[date_col]]
#       ),
#       stm2 = as.numeric(
#         .data[[stm2_col]]
#       )
#     ) %>%
#     filter(
#       date >= START_DATE,
#       date <= END_DATE,
#       is.finite(stm2)
#     ) %>%
#     group_by(date) %>%
#     summarise(
#       stm2 = mean(
#         stm2,
#         na.rm = TRUE
#       ),
#       .groups = "drop"
#     )
# }
# 
# format_correlation <- function(
#     data,
#     response
# ) {
#   complete <- data %>%
#     filter(
#       is.finite(stm2),
#       is.finite(
#         .data[[response]]
#       )
#     )
#   
#   if (nrow(complete) < 3L) {
#     return(
#       "Insufficient data"
#     )
#   }
#   
#   result <- cor.test(
#     complete$stm2,
#     complete[[response]],
#     method = "pearson"
#   )
#   
#   paste0(
#     "r = ",
#     number(
#       unname(result$estimate),
#       accuracy = 0.01
#     ),
#     "\n",
#     "p ",
#     pvalue(
#       result$p.value,
#       accuracy = 0.001,
#       add_p = TRUE
#     ),
#     "\n",
#     "n = ",
#     comma(
#       nrow(complete)
#     )
#   )
# }
# 
# make_panel <- function(
#     data,
#     response,
#     y_label,
#     panel_label
# ) {
#   annotation <- format_correlation(
#     data,
#     response
#   )
#   
#   ggplot(
#     data,
#     aes(
#       x = stm2,
#       y = .data[[response]]
#     )
#   ) +
#     geom_point(
#       shape = 16,
#       size = 0.65,
#       alpha = 0.25,
#       na.rm = TRUE
#     ) +
#     geom_smooth(
#       method = "lm",
#       formula = y ~ x,
#       se = TRUE,
#       linewidth = 0.7,
#       na.rm = TRUE
#     ) +
#     annotate(
#       geom = "label",
#       x = -Inf,
#       y = Inf,
#       label = annotation,
#       hjust = -0.08,
#       vjust = 1.08,
#       size = 2.5,
#       label.size = 0.25,
#       fill = "white",
#       color = "black"
#     ) +
#     labs(
#       x = "Space-time metric 2",
#       y = y_label,
#       tag = panel_label
#     ) +
#     theme_bw(
#       base_size = 8
#     ) +
#     theme(
#       panel.grid.minor = element_blank(),
#       panel.grid.major = element_line(
#         linewidth = 0.25
#       ),
#       axis.title = element_text(
#         size = 8
#       ),
#       axis.text = element_text(
#         size = 7
#       ),
#       plot.tag = element_text(
#         size = 9,
#         face = "bold"
#       ),
#       plot.tag.position = c(
#         0.01,
#         0.99
#       )
#     )
# }
# 
# make_warning_figure <- function(
#     warning_data,
#     stm2_data,
#     phenom,
#     phenom_name,
#     stm2_name,
#     filename
# ) {
#   full_dates <- tibble(
#     date = seq.Date(
#       START_DATE,
#       END_DATE,
#       by = "day"
#     )
#   )
#   
#   daily <- warning_data %>%
#     filter(
#       phenom == !!phenom
#     ) %>%
#     group_by(date) %>%
#     summarise(
#       event_count = sum(
#         event_count,
#         na.rm = TRUE
#       ),
#       area_km2 = sum(
#         area_km2,
#         na.rm = TRUE
#       ),
#       .groups = "drop"
#     )
#   
#   plot_data <- full_dates %>%
#     left_join(
#       daily,
#       by = "date"
#     ) %>%
#     mutate(
#       event_count = coalesce(
#         event_count,
#         0
#       ),
#       area_km2 = coalesce(
#         area_km2,
#         0
#       )
#     ) %>%
#     inner_join(
#       stm2_data,
#       by = "date"
#     )
#   
#   if (nrow(plot_data) == 0L) {
#     warning(
#       "No overlapping warning and STM2 dates for ",
#       phenom
#     )
#     return(
#       invisible(NULL)
#     )
#   }
#   
#   output_csv <- file.path(
#     OUT_DATA_DIR,
#     paste0(
#       "warning_stm2_",
#       tolower(phenom),
#       ".csv"
#     )
#   )
#   
#   write_csv(
#     plot_data,
#     output_csv
#   )
#   
#   count_plot <- make_panel(
#     data = plot_data,
#     response = "event_count",
#     y_label = "Active warnings per day",
#     panel_label = "a)"
#   )
#   
#   area_plot <- make_panel(
#     data = plot_data,
#     response = "area_km2",
#     y_label = expression(
#       "Active warning area per day (km"^2*")"
#     ),
#     panel_label = "b)"
#   )
#   
#   figure <- (
#     count_plot |
#       area_plot
#   ) +
#     plot_annotation(
#       title = paste0(
#         phenom,
#         " — ",
#         phenom_name
#       ),
#       subtitle = paste0(
#         stm2_name,
#         " STM2, 2002–2023"
#       ),
#       theme = theme(
#         plot.title = element_text(
#           size = 10,
#           face = "bold",
#           hjust = 0.5
#         ),
#         plot.subtitle = element_text(
#           size = 8,
#           hjust = 0.5
#         )
#       )
#     )
#   
#   ggsave(
#     filename = filename,
#     plot = figure,
#     width = 6.5,
#     height = 3.5,
#     units = "in",
#     dpi = 300,
#     bg = "white"
#   )
#   
#   message(
#     "Wrote: ",
#     filename
#   )
# }
# 
# # -----------------------------------------------------------------------------
# # Read warning data
# # -----------------------------------------------------------------------------
# 
# if (!file_exists(WARNING_FILE)) {
#   stop(
#     "Missing warning file: ",
#     WARNING_FILE
#   )
# }
# 
# warnings <- read_csv(
#   WARNING_FILE,
#   show_col_types = FALSE
# ) %>%
#   transmute(
#     date = as.Date(date),
#     phenom = toupper(
#       trimws(phenom)
#     ),
#     phenom_name,
#     event_count = as.numeric(
#       event_count
#     ),
#     area_km2 = as.numeric(
#       area_km2
#     )
#   ) %>%
#   filter(
#     date >= START_DATE,
#     date <= END_DATE
#   )
# 
# # -----------------------------------------------------------------------------
# # Read STM2 data
# # -----------------------------------------------------------------------------
# 
# heat_stm2 <- standardize_stm2(
#   HEAT_FILE,
#   "Heat-index"
# )
# 
# precip_stm2 <- standardize_stm2(
#   PRECIP_FILE,
#   "Precipitation"
# )
# 
# # -----------------------------------------------------------------------------
# # Figure definitions
# # -----------------------------------------------------------------------------
# 
# figure_definitions <- tribble(
#   ~phenom, ~phenom_name,          ~stm2_source, ~stm2_name,
#   "HI",    "Heat",                "heat",       "Heat-index",
#   "EH",    "Excessive Heat",      "heat",       "Heat-index",
#   "FF",    "Flash Flood",         "precip",     "Precipitation",
#   "TR",    "Tropical Storm",      "precip",     "Precipitation",
#   "SV",    "Severe Thunderstorm", "precip",     "Precipitation"
# )
# 
# # -----------------------------------------------------------------------------
# # Create figures
# # -----------------------------------------------------------------------------
# 
# for (index in seq_len(nrow(figure_definitions))) {
#   definition <- figure_definitions[
#     index,
#   ]
#   
#   selected_stm2 <- if (
#     definition$stm2_source == "heat"
#   ) {
#     heat_stm2
#   } else {
#     precip_stm2
#   }
#   
#   output_figure <- file.path(
#     OUT_FIG_DIR,
#     paste0(
#       "fig_s09_mesonet_stm2_",
#       tolower(
#         definition$phenom
#       ),
#       ".png"
#     )
#   )
#   
#   make_warning_figure(
#     warning_data = warnings,
#     stm2_data = selected_stm2,
#     phenom = definition$phenom,
#     phenom_name = definition$phenom_name,
#     stm2_name = definition$stm2_name,
#     filename = output_figure
#   )
# }

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lubridate)
  library(ggplot2)
  library(patchwork)
  library(fs)
  library(scales)
})

# =============================================================================
# Supplemental warning–Ar correlation figures
#
# Heat-index warning group:
#   HT — Heat
#   EH — Excessive Heat
#
# Precipitation warning group:
#   FF — Flash Flood
#   SV — Severe Thunderstorm
#   TR — Tropical Storm
#   GL — Gale
#   HI — Inland Hurricane
#   HU — Hurricane
#   SR — Storm
#   TI — Inland Tropical Storm
#
# Warning geometries are:
#   1. Expanded across all active dates from ISSUED through EXPIRED
#   2. Clipped to the Southeast domain
#   3. Unioned within day and warning group
#
# Therefore:
#   count = number of spatially distinct unioned warning components per day
#   area  = total unioned warning area per day
#
# Zero-warning days are excluded.
# =============================================================================

ROOT <- Sys.getenv(
  "CISC_ROOT",
  "/proj/sebastianlab/users/hquintal/cisc"
)

MESONET_DIR <- file.path(
  ROOT,
  "data/in/mesonet"
)

HEAT_FILE <- file.path(
  ROOT,
  "data/out/03_cov_model/heatindex",
  "heatindex_daily_optimal_cov_model_southeast.csv"
)

PRECIP_FILE <- file.path(
  ROOT,
  "data/out/03_cov_model/precipitation",
  "precipitation_daily_optimal_cov_model_southeast.csv"
)

OUT_DATA_DIR <- file.path(
  ROOT,
  "data/out/10_mesonet/ar_correlations_union"
)

CACHE_DIR <- file.path(
  OUT_DATA_DIR,
  "annual"
)

OUT_FIG_DIR <- file.path(
  ROOT,
  "figures/fig_s09_mesonet_ar_union"
)

START_DATE <- as.Date("2002-01-01")
END_DATE <- as.Date("2023-12-31")
YEARS <- 2002:2023

HEAT_TYPES <- c(
  "HT",
  "EH"
)

PRECIP_TYPES <- c(
  "FF",
  "SV",
  "TR",
  "GL",
  "HI",
  "HU",
  "SR",
  "TI"
)

AREA_CRS <- 5070

OVERWRITE_CACHE <- FALSE

dir_create(
  OUT_DATA_DIR,
  recurse = TRUE
)

dir_create(
  CACHE_DIR,
  recurse = TRUE
)

dir_create(
  OUT_FIG_DIR,
  recurse = TRUE
)

# -----------------------------------------------------------------------------
# Southeast domain
# -----------------------------------------------------------------------------

southeast <- st_as_sfc(
  st_bbox(
    c(
      xmin = -95,
      ymin = 25,
      xmax = -75,
      ymax = 40
    ),
    crs = st_crs(4326)
  )
) %>%
  st_as_sf() %>%
  st_transform(AREA_CRS)

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

parse_mesonet_datetime <- function(x) {
  x <- as.character(x)
  x <- sub("\\.0$", "", x)
  
  suppressWarnings(
    ymd_hm(
      x,
      tz = "UTC",
      quiet = TRUE
    )
  )
}

find_year_shapefile <- function(year) {
  year_dir <- file.path(
    MESONET_DIR,
    paste0(
      year,
      "_all"
    )
  )
  
  expected <- file.path(
    year_dir,
    paste0(
      "wwa_",
      year,
      "01010000_",
      year,
      "12312359.shp"
    )
  )
  
  if (file_exists(expected)) {
    return(expected)
  }
  
  candidates <- dir_ls(
    year_dir,
    regexp = "\\.shp$",
    type = "file"
  )
  
  if (length(candidates) != 1L) {
    stop(
      "Expected one shapefile for ",
      year,
      "; found ",
      length(candidates)
    )
  }
  
  candidates[[1]]
}

read_ar <- function(
    file,
    variable_name
) {
  if (!file_exists(file)) {
    stop(
      "Missing ",
      variable_name,
      " file: ",
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
      "day",
      "datetime"
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

safe_make_valid <- function(x) {
  suppressWarnings(
    st_make_valid(x)
  )
}

expand_active_dates <- function(
    issued,
    expired
) {
  start_date <- as.Date(issued)
  end_date <- as.Date(expired)
  
  map2(
    start_date,
    end_date,
    function(start, end) {
      if (
        is.na(start) ||
        is.na(end) ||
        end < start
      ) {
        return(
          as.Date(character())
        )
      }
      
      seq.Date(
        start,
        end,
        by = "day"
      )
    }
  )
}

union_daily_geometry <- function(day_sf) {
  if (nrow(day_sf) == 0L) {
    return(
      tibble(
        event_count = 0L,
        area_km2 = 0
      )
    )
  }
  
  unioned <- suppressWarnings(
    st_union(
      st_geometry(day_sf)
    )
  )
  
  unioned <- safe_make_valid(
    unioned
  )
  
  polygon_parts <- suppressWarnings(
    st_collection_extract(
      unioned,
      "POLYGON"
    )
  )
  
  if (
    length(polygon_parts) == 0L ||
    all(st_is_empty(polygon_parts))
  ) {
    return(
      tibble(
        event_count = 0L,
        area_km2 = 0
      )
    )
  }
  
  polygon_parts <- suppressWarnings(
    st_cast(
      polygon_parts,
      "POLYGON"
    )
  )
  
  polygon_parts <- polygon_parts[
    !st_is_empty(polygon_parts)
  ]
  
  tibble(
    event_count = length(
      polygon_parts
    ),
    area_km2 = sum(
      as.numeric(
        st_area(polygon_parts)
      ),
      na.rm = TRUE
    ) / 1e6
  )
}

# -----------------------------------------------------------------------------
# Annual processing
# -----------------------------------------------------------------------------

process_year <- function(year) {
  cache_file <- file.path(
    CACHE_DIR,
    paste0(
      "mesonet_union_daily_",
      year,
      ".csv"
    )
  )
  
  if (
    file_exists(cache_file) &&
    !OVERWRITE_CACHE
  ) {
    message(
      "Skipping cached year: ",
      year
    )
    
    return(
      read_csv(
        cache_file,
        show_col_types = FALSE
      ) %>%
        mutate(
          date = as.Date(date)
        )
    )
  }
  
  input_file <- find_year_shapefile(
    year
  )
  
  message(
    "Reading ",
    year,
    ": ",
    input_file
  )
  
  warnings <- st_read(
    input_file,
    quiet = TRUE,
    query = paste(
      "SELECT ISSUED, EXPIRED, PHENOM",
      "FROM",
      paste0(
        "'",
        path_ext_remove(
          path_file(input_file)
        ),
        "'"
      )
    )
  )
  
  warnings <- warnings %>%
    mutate(
      PHENOM = toupper(
        trimws(
          as.character(PHENOM)
        )
      ),
      issued_time = parse_mesonet_datetime(
        ISSUED
      ),
      expired_time = parse_mesonet_datetime(
        EXPIRED
      ),
      warning_group = case_when(
        PHENOM %in% HEAT_TYPES ~ "heat",
        PHENOM %in% PRECIP_TYPES ~ "precipitation",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      !is.na(warning_group),
      !is.na(issued_time),
      !is.na(expired_time),
      expired_time >= issued_time
    )
  
  if (nrow(warnings) == 0L) {
    result <- tibble(
      date = as.Date(character()),
      warning_group = character(),
      event_count = integer(),
      area_km2 = double()
    )
    
    write_csv(
      result,
      cache_file
    )
    
    return(result)
  }
  
  if (is.na(st_crs(warnings))) {
    st_crs(warnings) <- 4326
  }
  
  warnings <- warnings %>%
    safe_make_valid() %>%
    st_transform(AREA_CRS)
  
  warnings <- suppressWarnings(
    st_intersection(
      warnings,
      southeast
    )
  )
  
  warnings <- warnings[
    !st_is_empty(warnings),
  ]
  
  warnings <- warnings %>%
    mutate(
      active_date = expand_active_dates(
        issued_time,
        expired_time
      )
    ) %>%
    unnest(active_date) %>%
    filter(
      active_date >= START_DATE,
      active_date <= END_DATE
    )
  
  if (nrow(warnings) == 0L) {
    result <- tibble(
      date = as.Date(character()),
      warning_group = character(),
      event_count = integer(),
      area_km2 = double()
    )
    
    write_csv(
      result,
      cache_file
    )
    
    return(result)
  }
  
  groups <- warnings %>%
    st_drop_geometry() %>%
    distinct(
      date = active_date,
      warning_group
    )
  
  result <- pmap_dfr(
    groups,
    function(
    date,
    warning_group
    ) {
      day_sf <- warnings %>%
        filter(
          active_date == date,
          .data$warning_group == warning_group
        )
      
      union_daily_geometry(
        day_sf
      ) %>%
        mutate(
          date = date,
          warning_group = warning_group,
          .before = 1
        )
    }
  )
  
  write_csv(
    result,
    cache_file
  )
  
  message(
    "Wrote: ",
    cache_file
  )
  
  result
}

# -----------------------------------------------------------------------------
# Process and combine annual warning geometries
# -----------------------------------------------------------------------------

warning_daily <- map_dfr(
  YEARS,
  process_year
) %>%
  mutate(
    date = as.Date(date)
  ) %>%
  group_by(
    date,
    warning_group
  ) %>%
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
  ) %>%
  filter(
    event_count > 0,
    area_km2 > 0
  )

write_csv(
  warning_daily,
  file.path(
    OUT_DATA_DIR,
    "mesonet_union_daily_heat_precipitation_2002_2023.csv"
  )
)

# -----------------------------------------------------------------------------
# Correlation and plotting
# -----------------------------------------------------------------------------

format_correlation <- function(
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
  
  if (nrow(complete) < 3L) {
    return(
      "Insufficient data"
    )
  }
  
  result <- cor.test(
    complete$Ar,
    complete[[response]],
    method = "pearson"
  )
  
  paste0(
    "r = ",
    number(
      unname(
        result$estimate
      ),
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
      size = 0.65,
      alpha = 0.25
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.7
    ) +
    annotate(
      geom = "label",
      x = -Inf,
      y = Inf,
      label = format_correlation(
        data,
        response
      ),
      hjust = -0.08,
      vjust = 1.08,
      size = 2.5,
      linewidth = 0.25,
      fill = "white"
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
    warning_group_name,
    ar_data,
    title,
    x_label,
    output_name
) {
  plot_data <- warning_daily %>%
    filter(
      warning_group == warning_group_name
    ) %>%
    inner_join(
      ar_data,
      by = "date"
    ) %>%
    filter(
      event_count > 0,
      area_km2 > 0,
      is.finite(Ar)
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
    y_label = "Spatially distinct active warnings per day",
    panel_label = "a)"
  )
  
  area_plot <- make_panel(
    data = plot_data,
    response = "area_km2",
    x_label = x_label,
    y_label = expression(
      "Unioned active warning area per day (km"^2*")"
    ),
    panel_label = "b)"
  )
  
  figure <- (
    count_plot |
      area_plot
  ) +
    plot_annotation(
      title = title,
      subtitle = "Daily non-overlapping warning geometries, 2002–2023",
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

heat_ar <- read_ar(
  HEAT_FILE,
  "Heat-index"
)

precip_ar <- read_ar(
  PRECIP_FILE,
  "Precipitation"
)

make_figure(
  warning_group_name = "heat",
  ar_data = heat_ar,
  title = "Heat and excessive heat warnings",
  x_label = expression(
    "Heat-index "*A[r]
  ),
  output_name = "fig_s09_mesonet_ar_heat"
)

make_figure(
  warning_group_name = "precipitation",
  ar_data = precip_ar,
  title = paste(
    "Flash flood, severe thunderstorm, tropical storm,",
    "gale, hurricane, and storm warnings"
  ),
  x_label = expression(
    "Precipitation "*A[r]
  ),
  output_name = "fig_s09_mesonet_ar_precipitation"
)