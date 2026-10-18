#!/usr/bin/env Rscript

# scripts/fs12_fs17_performance_spatial_difference.R
# =============================================================================
# County-level spatial differences between resampled and native model
# performance.
#
# Six 4-column x 4-row figures are written:
#
#   fs12_spatial_difference_heat_f1.png
#   fs13_spatial_difference_precip_f1.png
#   fs14_spatial_difference_heat_recall.png
#   fs15_spatial_difference_precip_recall.png
#   fs16_spatial_difference_heat_precision.png
#   fs17_spatial_difference_precip_precision.png
#
# Each county value is:
#
#   resampled metric - native metric
#
# Rows:
#   5-year, 10-year, 25-year, and 50-year return periods
#
# Columns:
#   STM1 bilinear
#   STM1 max
#   STM2 bilinear
#   STM2 max
#
# Fill:
#   negative difference = resampled performs worse than native
#   zero                = no difference
#   positive difference = resampled performs better than native
#
# All panels within each figure use the same symmetric color scale centered
# on zero.
#
# Default warning combinations:
#   heat          : EH_W
#   precipitation : FF_W
#
# Output:
#   figures/fs12_fs17_performance_spatial/
# =============================================================================


suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(grid)
})


# =============================================================================
# CONFIGURATION
# =============================================================================

ROOT <- Sys.getenv(
  "CISC_ROOT",
  "/proj/sebastianlab/users/hquintal/cisc"
)

PERF_BASE <- Sys.getenv(
  "CISC_PERF_DIR",
  file.path(
    ROOT,
    "data",
    "out",
    "20260802_model_performance"
  )
)

HEAT_COMBINATION <- Sys.getenv(
  "CISC_HEAT_COMBINATION",
  "EH_W"
)

PRECIPITATION_COMBINATION <- Sys.getenv(
  "CISC_PRECIPITATION_COMBINATION",
  "FF_W"
)

COUNTY_SHP_ENV <- trimws(
  Sys.getenv(
    "CISC_COUNTY_SHP",
    ""
  )
)

COUNTY_CACHE <- file.path(
  PERF_BASE,
  "_shared",
  "control",
  "counties_2022_cb.gpkg"
)

COUNTY_ZIP_CACHE <- file.path(
  PERF_BASE,
  "_shared",
  "control",
  "cb_2022_us_county_500k.zip"
)

COUNTY_EXTRACT_DIR <- file.path(
  PERF_BASE,
  "_shared",
  "control",
  "cb_2022_us_county_500k"
)

COUNTY_URL <- Sys.getenv(
  "CISC_COUNTY_URL",
  paste0(
    "https://www2.census.gov/geo/tiger/GENZ2022/shp/",
    "cb_2022_us_county_500k.zip"
  )
)

OUT_DIR <- file.path(
  ROOT,
  "figures",
  "fs12_fs17_performance_spatial_difference"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

WINDOW_TAG <- "20020101_20231231"

RETURN_PERIODS <- c(
  5L,
  10L,
  25L,
  50L
)

BBOX <- c(
  xmin = -95,
  ymin = 25,
  xmax = -75,
  ymax = 40
)

NA_COLOR <- "grey75"

FIGURE_WIDTH_IN <- 13.5
FIGURE_HEIGHT_IN <- 11.0
FIGURE_DPI <- 600


# =============================================================================
# DATASET CONFIGURATION
# =============================================================================

DATASET_CONFIG <- data.table(
  var_group = c(
    rep("heatindex", 4L),
    rep("precipitation", 4L)
  ),
  
  native_tag = c(
    rep("heatindex_0.25", 4L),
    rep("precipitation_0.25", 4L)
  ),
  
  resampled_tag = c(
    "heatindex_0.30_bilinear",
    "heatindex_0.30_max",
    "heatindex_0.36_bilinear",
    "heatindex_0.36_max",
    
    "precipitation_0.11_bilinear",
    "precipitation_0.11_max",
    "precipitation_0.14_bilinear",
    "precipitation_0.14_max"
  ),
  
  strategy = rep(
    c(
      "STM1 bilinear",
      "STM1 max",
      "STM2 bilinear",
      "STM2 max"
    ),
    2L
  ),
  
  strategy_order = rep(
    seq_len(4L),
    2L
  )
)

STRATEGY_LEVELS <- c(
  "STM1 bilinear",
  "STM1 max",
  "STM2 bilinear",
  "STM2 max"
)

RP_LEVELS <- sprintf(
  "%d-year",
  RETURN_PERIODS
)


# =============================================================================
# LOGGING
# =============================================================================

.log <- function(...) {
  cat(
    sprintf(
      "[%s] ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    sprintf(...),
    "\n",
    sep = ""
  )
}


# =============================================================================
# GENERIC HELPERS
# =============================================================================

.normalize_geoid <- function(x) {
  x <- as.character(x)
  
  x <- gsub(
    "[^0-9]",
    "",
    x
  )
  
  output <- rep(
    NA_character_,
    length(x)
  )
  
  valid <- nzchar(x) &
    nchar(x) <= 5L
  
  output[valid] <- sprintf(
    "%05d",
    as.integer(x[valid])
  )
  
  output
}


.metric_label <- function(metric_name) {
  switch(
    metric_name,
    f1 = "F1 difference",
    recall = "Recall difference",
    precision = "Precision difference",
    stop(
      "Unsupported metric: ",
      metric_name
    )
  )
}


# =============================================================================
# READ PERFORMANCE RESULTS
# =============================================================================

.read_metrics <- function(
    variable_dir,
    combination_tag
) {
  input_dir <- file.path(
    PERF_BASE,
    variable_dir,
    combination_tag
  )
  
  csv_file <- file.path(
    input_dir,
    "metrics_county.csv"
  )
  
  parquet_file <- file.path(
    input_dir,
    "metrics_county.parquet"
  )
  
  if (file.exists(csv_file)) {
    .log(
      "Reading %s",
      csv_file
    )
    
    return(
      fread(
        csv_file,
        na.strings = c(
          "",
          "NA",
          "NaN"
        )
      )
    )
  }
  
  if (file.exists(parquet_file)) {
    if (!requireNamespace(
      "arrow",
      quietly = TRUE
    )) {
      stop(
        "Found only the Parquet metrics file, but package ",
        "'arrow' is unavailable:\n",
        parquet_file
      )
    }
    
    .log(
      "Reading %s",
      parquet_file
    )
    
    return(
      as.data.table(
        arrow::read_parquet(
          parquet_file
        )
      )
    )
  }
  
  stop(
    "Missing metrics_county.csv and metrics_county.parquet under:\n",
    input_dir
  )
}


.validate_metrics <- function(
    dt,
    expected_var_group,
    expected_combination
) {
  required_columns <- c(
    "var_group",
    "combination_tag",
    "dataset_tag",
    "return_period",
    "window_tag",
    "GEOID",
    "f1",
    "recall",
    "precision"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(dt)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Metrics table is missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  output <- copy(dt)
  
  character_columns <- intersect(
    c(
      "variable",
      "var_group",
      "combination_tag",
      "warnings",
      "dataset_tag",
      "window_tag",
      "GEOID"
    ),
    names(output)
  )
  
  for (column_name in character_columns) {
    set(
      output,
      j = column_name,
      value = as.character(
        output[[column_name]]
      )
    )
  }
  
  output[
    ,
    var_group := tolower(
      var_group
    )
  ]
  
  output[
    ,
    return_period := as.integer(
      return_period
    )
  ]
  
  output[
    ,
    GEOID := .normalize_geoid(
      GEOID
    )
  ]
  
  for (metric_name in c(
    "f1",
    "recall",
    "precision"
  )) {
    set(
      output,
      j = metric_name,
      value = as.numeric(
        output[[metric_name]]
      )
    )
  }
  
  valid_dataset_tags <- unique(
    c(
      DATASET_CONFIG[
        var_group == expected_var_group,
        native_tag
      ],
      DATASET_CONFIG[
        var_group == expected_var_group,
        resampled_tag
      ]
    )
  )
  
  keep <-
    output[["var_group"]] ==
    expected_var_group &
    output[["combination_tag"]] ==
    expected_combination &
    output[["window_tag"]] ==
    WINDOW_TAG &
    output[["return_period"]] %in%
    RETURN_PERIODS &
    output[["dataset_tag"]] %in%
    valid_dataset_tags
  
  keep[
    is.na(keep)
  ] <- FALSE
  
  output <- output[
    which(keep)
  ]
  
  if (nrow(output) == 0L) {
    stop(
      "No rows remain after filtering metrics for:\n",
      "  var_group: ",
      expected_var_group,
      "\n",
      "  combination_tag: ",
      expected_combination,
      "\n",
      "  window_tag: ",
      WINDOW_TAG
    )
  }
  
  duplicate_rows <- output[
    ,
    .N,
    by = .(
      return_period,
      dataset_tag,
      GEOID
    )
  ][
    N > 1L
  ]
  
  if (nrow(duplicate_rows) > 0L) {
    stop(
      "Metrics contain duplicate county/dataset/return-period rows:\n",
      paste(
        capture.output(
          print(
            head(
              duplicate_rows,
              20L
            )
          )
        ),
        collapse = "\n"
      )
    )
  }
  
  output[]
}


# =============================================================================
# COUNTY GEOMETRY
# =============================================================================

.resolve_county_source <- function() {
  if (nzchar(COUNTY_SHP_ENV)) {
    if (!file.exists(COUNTY_SHP_ENV)) {
      stop(
        "CISC_COUNTY_SHP points to a missing file:\n",
        COUNTY_SHP_ENV
      )
    }
    
    return(
      COUNTY_SHP_ENV
    )
  }
  
  if (file.exists(COUNTY_CACHE)) {
    return(
      COUNTY_CACHE
    )
  }
  
  dir.create(
    dirname(COUNTY_ZIP_CACHE),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if (!file.exists(COUNTY_ZIP_CACHE)) {
    .log(
      "Downloading county boundaries from %s",
      COUNTY_URL
    )
    
    download.file(
      url = COUNTY_URL,
      destfile = COUNTY_ZIP_CACHE,
      mode = "wb",
      quiet = FALSE
    )
  }
  
  dir.create(
    COUNTY_EXTRACT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  unzip(
    COUNTY_ZIP_CACHE,
    exdir = COUNTY_EXTRACT_DIR
  )
  
  shapefiles <- list.files(
    COUNTY_EXTRACT_DIR,
    pattern = "\\.shp$",
    full.names = TRUE
  )
  
  if (length(shapefiles) == 0L) {
    stop(
      "Downloaded county ZIP contained no shapefile:\n",
      COUNTY_ZIP_CACHE
    )
  }
  
  county_sf <- st_read(
    shapefiles[[1L]],
    quiet = TRUE
  )
  
  st_write(
    county_sf,
    COUNTY_CACHE,
    driver = "GPKG",
    append = FALSE,
    quiet = TRUE
  )
  
  COUNTY_CACHE
}


.read_counties <- function() {
  county_source <- .resolve_county_source()
  
  .log(
    "Reading county geometry from %s",
    county_source
  )
  
  counties <- st_read(
    county_source,
    quiet = TRUE
  )
  
  required_columns <- c(
    "GEOID",
    "STATEFP"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(counties)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "County geometry is missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  counties$GEOID <- .normalize_geoid(
    counties$GEOID
  )
  
  counties$STATEFP <- as.character(
    counties$STATEFP
  )
  
  excluded_state_codes <- c(
    "02",
    "15",
    "60",
    "66",
    "69",
    "72",
    "78"
  )
  
  counties <- counties[
    !counties$STATEFP %in%
      excluded_state_codes &
      !is.na(counties$GEOID),
  ]
  
  counties <- st_transform(
    counties,
    4326
  )
  
  counties <- st_make_valid(
    counties
  )
  
  bbox_geometry <- st_as_sfc(
    st_bbox(
      BBOX,
      crs = st_crs(4326)
    )
  )
  
  intersects_bbox <- lengths(
    st_intersects(
      counties,
      bbox_geometry
    )
  ) > 0L
  
  counties <- counties[
    intersects_bbox,
    c(
      "GEOID",
      "STATEFP"
    )
  ]
  
  counties <- counties[
    order(
      counties$GEOID
    ),
  ]
  
  if (anyDuplicated(
    counties$GEOID
  )) {
    stop(
      "County geometry contains duplicate GEOIDs after filtering"
    )
  }
  
  states <- aggregate(
    counties[
      ,
      "STATEFP",
      drop = FALSE
    ],
    by = list(
      STATEFP = counties$STATEFP
    ),
    FUN = function(x) {
      x[[1L]]
    },
    do_union = TRUE
  )
  
  list(
    counties = counties,
    states = states
  )
}


# =============================================================================
# BUILD RESAMPLED-MINUS-NATIVE DIFFERENCES
# =============================================================================

.prepare_difference_data <- function(
    metrics,
    expected_var_group,
    metric_name
) {
  configuration <- DATASET_CONFIG[
    var_group == expected_var_group
  ]
  
  native_tags <- unique(
    configuration$native_tag
  )
  
  if (length(native_tags) != 1L) {
    stop(
      "Expected exactly one native dataset for ",
      expected_var_group
    )
  }
  
  native_tag <- native_tags[[1L]]
  
  native_dt <- metrics[
    dataset_tag == native_tag,
    .(
      GEOID,
      return_period,
      native_value = get(
        metric_name
      )
    )
  ]
  
  if (nrow(native_dt) == 0L) {
    stop(
      "No native rows found for dataset ",
      native_tag
    )
  }
  
  resampled_dt <- metrics[
    dataset_tag %in%
      configuration$resampled_tag,
    .(
      GEOID,
      return_period,
      resampled_tag = dataset_tag,
      resampled_value = get(
        metric_name
      )
    )
  ]
  
  resampled_dt <- merge(
    resampled_dt,
    configuration[
      ,
      .(
        resampled_tag,
        strategy,
        strategy_order
      )
    ],
    by = "resampled_tag",
    all.x = TRUE,
    sort = FALSE
  )
  
  output <- merge(
    resampled_dt,
    native_dt,
    by = c(
      "GEOID",
      "return_period"
    ),
    all = FALSE,
    sort = FALSE
  )
  
  output[
    ,
    difference := resampled_value -
      native_value
  ]
  
  # A difference requires finite values for both native and resampled
  # performance. Undefined comparisons remain NA and plot in grey.
  output[
    !is.finite(native_value) |
      !is.finite(resampled_value),
    difference := NA_real_
  ]
  
  output[
    ,
    rp_label := factor(
      sprintf(
        "%d-year",
        return_period
      ),
      levels = RP_LEVELS
    )
  ]
  
  output[
    ,
    strategy := factor(
      strategy,
      levels = STRATEGY_LEVELS
    )
  ]
  
  expected_panels <- CJ(
    return_period = RETURN_PERIODS,
    resampled_tag = configuration$resampled_tag,
    unique = TRUE
  )
  
  observed_panels <- unique(
    output[
      ,
      .(
        return_period,
        resampled_tag
      )
    ]
  )
  
  missing_panels <- fsetdiff(
    expected_panels,
    observed_panels
  )
  
  if (nrow(missing_panels) > 0L) {
    stop(
      "Missing native/resampled comparisons for these panels:\n",
      paste(
        capture.output(
          print(
            missing_panels
          )
        ),
        collapse = "\n"
      )
    )
  }
  
  output[]
}


# =============================================================================
# COLOR-SCALE HELPERS
# =============================================================================

.round_up_limit <- function(
    maximum_absolute_value
) {
  if (
    !is.finite(maximum_absolute_value) ||
    maximum_absolute_value <= 0
  ) {
    return(
      0.1
    )
  }
  
  candidate_steps <- c(
    0.01,
    0.02,
    0.025,
    0.05,
    0.10,
    0.20,
    0.25,
    0.50,
    1.00
  )
  
  target_step <- maximum_absolute_value / 4
  
  step <- candidate_steps[
    which(
      candidate_steps >= target_step
    )[1L]
  ]
  
  if (is.na(step)) {
    step <- 1
  }
  
  ceiling(
    maximum_absolute_value /
      step
  ) * step
}


.make_difference_breaks <- function(
    scale_limit
) {
  pretty(
    c(
      -scale_limit,
      scale_limit
    ),
    n = 5L
  )
}


# =============================================================================
# THEME
# =============================================================================

.theme_spatial <- function(
    base_size = 20
) {
  theme_void(
    base_size = base_size,
    base_family = "Arial"
  ) +
    theme(
      text = element_text(
        family = "Arial",
        size = base_size,
        color = "black"
      ),
      
      strip.background = element_rect(
        fill = "white",
        color = "grey30",
        linewidth = 0.65
      ),
      
      strip.text.x = element_text(
        family = "Arial",
        face = "bold",
        size = 18,
        margin = margin(
          4,
          4,
          4,
          4
        )
      ),
      
      strip.text.y = element_text(
        family = "Arial",
        face = "bold",
        size = 18,
        angle = -90,
        margin = margin(
          4,
          4,
          4,
          4
        )
      ),
      
      panel.spacing.x = unit(
        0,
        "mm"
      ),
      
      panel.spacing.y = unit(
        0,
        "mm"
      ),
      
      legend.position = "bottom",
      
      legend.direction = "horizontal",
      
      legend.title = element_text(
        family = "Arial",
        size = 20,
        hjust = 0.5
      ),
      
      legend.text = element_text(
        family = "Arial",
        size = 15
      ),
      
      legend.margin = margin(
        0,
        0,
        0,
        0
      ),
      
      legend.box.margin = margin(
        0,
        0,
        0,
        0
      ),
      
      plot.margin = margin(
        6,
        6,
        6,
        6
      )
    )
}


# =============================================================================
# FIGURE BUILDER
# =============================================================================

.build_difference_figure <- function(
    metrics,
    county_geometry,
    state_geometry,
    expected_var_group,
    metric_name,
    output_stem
) {
  difference_dt <- .prepare_difference_data(
    metrics = metrics,
    expected_var_group = expected_var_group,
    metric_name = metric_name
  )
  
  spatial <- merge(
    county_geometry,
    difference_dt,
    by = "GEOID",
    all.x = FALSE,
    sort = FALSE
  )
  
  if (nrow(spatial) == 0L) {
    stop(
      "County join produced no spatial rows for ",
      expected_var_group,
      " / ",
      metric_name
    )
  }
  
  finite_differences <- difference_dt[
    is.finite(difference),
    difference
  ]
  
  if (length(finite_differences) == 0L) {
    stop(
      "No finite spatial differences were found for ",
      expected_var_group,
      " / ",
      metric_name
    )
  }
  
  raw_limit <- max(
    abs(
      finite_differences
    )
  )
  
  scale_limit <- .round_up_limit(
    raw_limit
  )
  
  legend_breaks <- .make_difference_breaks(
    scale_limit
  )
  
  legend_breaks <- legend_breaks[
    legend_breaks >= -scale_limit &
      legend_breaks <= scale_limit
  ]
  
  if (!any(
    abs(
      legend_breaks
    ) < .Machine$double.eps^0.5
  )) {
    legend_breaks <- sort(
      unique(
        c(
          legend_breaks,
          0
        )
      )
    )
  }
  
  fill_colors <- colorRampPalette(
    rev(
      brewer.pal(
        11L,
        "RdBu"
      )
    ),
    space = "Lab"
  )(
    256L
  )
  
  p <- ggplot() +
    geom_sf(
      data = spatial,
      aes(
        fill = difference
      ),
      color = "grey72",
      linewidth = 0.04
    ) +
    geom_sf(
      data = state_geometry,
      fill = NA,
      color = "black",
      linewidth = 0.35
    ) +
    facet_grid(
      rows = vars(
        rp_label
      ),
      cols = vars(
        strategy
      ),
      drop = FALSE
    ) +
    scale_fill_gradientn(
      colours = fill_colors,
      values = c(
        0,
        0.5,
        1
      ),
      limits = c(
        -scale_limit,
        scale_limit
      ),
      breaks = legend_breaks,
      labels = number(
        legend_breaks,
        accuracy = 0.01,
        trim = TRUE
      ),
      na.value = NA_COLOR,
      oob = squish,
      name = .metric_label(
        metric_name
      ),
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        direction = "horizontal",
        barwidth = unit(
          8.5,
          "in"
        ),
        barheight = unit(
          0.18,
          "in"
        ),
        ticks = TRUE,
        frame.colour = "black",
        frame.linewidth = 0.4
      )
    ) +
    coord_sf(
      xlim = BBOX[
        c(
          "xmin",
          "xmax"
        )
      ],
      ylim = BBOX[
        c(
          "ymin",
          "ymax"
        )
      ],
      expand = FALSE,
      datum = NA,
      clip = "on"
    ) +
    .theme_spatial()
  
  output_file <- file.path(
    OUT_DIR,
    paste0(
      output_stem,
      ".png"
    )
  )
  
  ggsave(
    filename = output_file,
    plot = p,
    width = FIGURE_WIDTH_IN,
    height = FIGURE_HEIGHT_IN,
    units = "in",
    dpi = FIGURE_DPI,
    bg = "white",
    limitsize = FALSE
  )
  
  panel_summary <- difference_dt[
    ,
    .(
      n_counties = uniqueN(
        GEOID
      ),
      n_finite = sum(
        is.finite(
          difference
        )
      ),
      n_undefined = sum(
        !is.finite(
          difference
        )
      ),
      mean_difference = mean(
        difference,
        na.rm = TRUE
      ),
      median_difference = median(
        difference,
        na.rm = TRUE
      ),
      minimum_difference = min(
        difference,
        na.rm = TRUE
      ),
      maximum_difference = max(
        difference,
        na.rm = TRUE
      )
    ),
    by = .(
      return_period,
      rp_label,
      strategy,
      resampled_tag
    )
  ]
  
  panel_summary[
    !is.finite(mean_difference),
    mean_difference := NA_real_
  ]
  
  panel_summary[
    !is.finite(median_difference),
    median_difference := NA_real_
  ]
  
  panel_summary[
    !is.finite(minimum_difference),
    minimum_difference := NA_real_
  ]
  
  panel_summary[
    !is.finite(maximum_difference),
    maximum_difference := NA_real_
  ]
  
  panel_summary[
    ,
    `:=`(
      raw_absolute_limit = raw_limit,
      plotted_absolute_limit = scale_limit
    )
  ]
  
  summary_file <- file.path(
    OUT_DIR,
    paste0(
      output_stem,
      "_panel_summary.csv"
    )
  )
  
  fwrite(
    panel_summary,
    summary_file
  )
  
  .log(
    "%s / %s: raw maximum absolute difference = %.6f; plotted limit = ±%.3f",
    expected_var_group,
    metric_name,
    raw_limit,
    scale_limit
  )
  
  .log(
    "Wrote %s",
    output_file
  )
  
  .log(
    "Wrote %s",
    summary_file
  )
  
  output_file
}


# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  .log(
    "PERF_BASE=%s",
    PERF_BASE
  )
  
  .log(
    "OUT_DIR=%s",
    OUT_DIR
  )
  
  geometry <- .read_counties()
  
  heat_metrics <- .read_metrics(
    variable_dir = "heat",
    combination_tag = HEAT_COMBINATION
  )
  
  heat_metrics <- .validate_metrics(
    dt = heat_metrics,
    expected_var_group = "heatindex",
    expected_combination = HEAT_COMBINATION
  )
  
  precipitation_metrics <- .read_metrics(
    variable_dir = "precipitation",
    combination_tag = PRECIPITATION_COMBINATION
  )
  
  precipitation_metrics <- .validate_metrics(
    dt = precipitation_metrics,
    expected_var_group = "precipitation",
    expected_combination = PRECIPITATION_COMBINATION
  )
  
  outputs <- c(
    .build_difference_figure(
      metrics = heat_metrics,
      county_geometry = geometry$counties,
      state_geometry = geometry$states,
      expected_var_group = "heatindex",
      metric_name = "f1",
      output_stem = "fs12_spatial_difference_heat_f1"
    ),
    
    .build_difference_figure(
      metrics = precipitation_metrics,
      county_geometry = geometry$counties,
      state_geometry = geometry$states,
      expected_var_group = "precipitation",
      metric_name = "f1",
      output_stem = "fs13_spatial_difference_precip_f1"
    ),
    
    .build_difference_figure(
      metrics = heat_metrics,
      county_geometry = geometry$counties,
      state_geometry = geometry$states,
      expected_var_group = "heatindex",
      metric_name = "recall",
      output_stem = "fs14_spatial_difference_heat_recall"
    ),
    
    .build_difference_figure(
      metrics = precipitation_metrics,
      county_geometry = geometry$counties,
      state_geometry = geometry$states,
      expected_var_group = "precipitation",
      metric_name = "recall",
      output_stem = "fs15_spatial_difference_precip_recall"
    ),
    
    .build_difference_figure(
      metrics = heat_metrics,
      county_geometry = geometry$counties,
      state_geometry = geometry$states,
      expected_var_group = "heatindex",
      metric_name = "precision",
      output_stem = "fs16_spatial_difference_heat_precision"
    ),
    
    .build_difference_figure(
      metrics = precipitation_metrics,
      county_geometry = geometry$counties,
      state_geometry = geometry$states,
      expected_var_group = "precipitation",
      metric_name = "precision",
      output_stem = "fs17_spatial_difference_precip_precision"
    )
  )
  
  manifest <- data.table(
    figure = basename(
      outputs
    ),
    path = outputs
  )
  
  manifest_file <- file.path(
    OUT_DIR,
    "figure_manifest.csv"
  )
  
  fwrite(
    manifest,
    manifest_file
  )
  
  cat(
    "\nSaved figures:\n"
  )
  
  cat(
    outputs,
    sep = "\n"
  )
  
  cat(
    "\n\nManifest:\n",
    manifest_file,
    "\n",
    sep = ""
  )
}


main()