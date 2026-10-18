#!/usr/bin/env Rscript

# scripts/dissertation_02_performance_spatial.R
# =============================================================================
# Dissertation county-level spatial model-performance figures
#
# Native datasets only:
#   Heat index:
#     heatindex_0.25
#
#   Precipitation:
#     precipitation_0.25
#
# Return period:
#   5-year only
#
# Each variable produces one 1 x 3 figure:
#   F1 | Recall | Precision
#
# Metric values:
#   > 0        = plotted using the variable-specific color scale
#   0 or NA    = grey
#
# Warning combinations:
#   Heat          = EH_W
#   Precipitation = FF_W
#
# Geographic display:
#   - EPSG:4326
#   - Southeast domain: 95W–75W, 25N–40N
#   - longitude/latitude axes and graticules
#   - county boundaries
#   - state boundaries
#
# Outputs:
#   /proj/sebastianlab/users/hquintal/cisc/figures/dissertation/
#     dissertation_02_performance_spatial_heat.png
#     dissertation_02_performance_spatial_precipitation.png
#
# =============================================================================


suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(grid)
  library(tigris)
  library(ragg)
})


options(
  tigris_use_cache = TRUE,
  tigris_class = "sf"
)

sf::sf_use_s2(FALSE)


# =============================================================================
# USER SETTINGS
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

OUT_DIR <- file.path(
  ROOT,
  "figures",
  "dissertation"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# ANALYSIS SETTINGS
# =============================================================================

WINDOW_TAG <- "20020101_20231231"

RETURN_PERIOD <- 5L

DATASETS <- data.table(
  variable = c(
    "heatindex",
    "precipitation"
  ),
  variable_dir = c(
    "heat",
    "precipitation"
  ),
  dataset_tag = c(
    "heatindex_0.25",
    "precipitation_0.25"
  ),
  combination_tag = c(
    HEAT_COMBINATION,
    PRECIPITATION_COMBINATION
  ),
  palette = c(
    "Reds",
    "Blues"
  ),
  output_name = c(
    "dissertation_02_performance_spatial_heat.png",
    "dissertation_02_performance_spatial_precipitation.png"
  )
)


# =============================================================================
# SPATIAL SETTINGS
# =============================================================================

CRS_PLOT <- 4326

BBOX <- c(
  xmin = -95,
  xmax = -75,
  ymin = 25,
  ymax = 40
)

METRIC_LIMITS <- c(
  0,
  1
)

METRIC_BREAKS <- seq(
  0,
  1,
  by = 0.2
)

ZERO_AND_NA_COLOR <- "grey75"


# =============================================================================
# FIGURE SETTINGS
# =============================================================================

FIG_WIDTH_IN <- 15.0
FIG_HEIGHT_IN <- 5.7
DPI <- 300

BASE_SIZE <- 18


# =============================================================================
# LOGGING
# =============================================================================

.log <- function(fmt, ...) {
  message(
    sprintf(
      "[%s] %s",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      ),
      sprintf(
        fmt,
        ...
      )
    )
  )
}


# =============================================================================
# HELPERS
# =============================================================================

.normalize_geoid <- function(x) {
  
  x <- as.character(
    x
  )
  
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
    as.integer(
      x[valid]
    )
  )
  
  output
}


.metric_palette <- function(
    palette_name,
    n_colors = 256L
) {
  
  colorRampPalette(
    RColorBrewer::brewer.pal(
      9L,
      palette_name
    ),
    space = "Lab"
  )(
    n_colors
  )
}


# =============================================================================
# READ METRICS
# =============================================================================

.read_metrics <- function(
    variable_dir,
    combination_tag,
    dataset_tag,
    variable
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
    
    x <- data.table::fread(
      csv_file,
      na.strings = c(
        "",
        "NA",
        "NaN"
      )
    )
    
  } else if (file.exists(parquet_file)) {
    
    if (!requireNamespace(
      "arrow",
      quietly = TRUE
    )) {
      stop(
        "Found only the Parquet metrics file, but package ",
        "'arrow' is unavailable:\n",
        parquet_file,
        call. = FALSE
      )
    }
    
    .log(
      "Reading %s",
      parquet_file
    )
    
    x <- as.data.table(
      arrow::read_parquet(
        parquet_file
      )
    )
    
  } else {
    
    stop(
      "Missing metrics_county.csv and metrics_county.parquet under:\n",
      input_dir,
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Validate required columns
  # ---------------------------------------------------------------------------
  
  required <- c(
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
  
  missing_cols <- setdiff(
    required,
    names(x)
  )
  
  if (length(missing_cols) > 0L) {
    stop(
      "Metrics table is missing required columns: ",
      paste(
        missing_cols,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Normalize fields
  # ---------------------------------------------------------------------------
  
  x[
    ,
    GEOID := .normalize_geoid(
      GEOID
    )
  ]
  
  x[
    ,
    var_group := tolower(
      as.character(
        var_group
      )
    )
  ]
  
  x[
    ,
    combination_tag := as.character(
      combination_tag
    )
  ]
  
  x[
    ,
    dataset_tag := as.character(
      dataset_tag
    )
  ]
  
  x[
    ,
    window_tag := as.character(
      window_tag
    )
  ]
  
  x[
    ,
    return_period := as.integer(
      return_period
    )
  ]
  
  for (metric_name in c(
    "f1",
    "recall",
    "precision"
  )) {
    
    set(
      x,
      j = metric_name,
      value = as.numeric(
        x[[metric_name]]
      )
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # Native RP5 only
  # ---------------------------------------------------------------------------
  
  expected_variable <- variable
  expected_combination <- combination_tag
  expected_dataset <- dataset_tag
  expected_rp <- RETURN_PERIOD
  expected_window <- WINDOW_TAG
  
  keep <-
    x[["var_group"]] == expected_variable &
    x[["combination_tag"]] == expected_combination &
    x[["dataset_tag"]] == expected_dataset &
    x[["return_period"]] == expected_rp &
    x[["window_tag"]] == expected_window
  
  keep[is.na(keep)] <- FALSE
  
  x <- x[
    which(keep)
  ]
  
  # ---------------------------------------------------------------------------
  # Require one row per county
  # ---------------------------------------------------------------------------
  
  duplicates <- x[
    ,
    .N,
    by = GEOID
  ][
    N > 1L
  ]
  
  if (nrow(duplicates) > 0L) {
    stop(
      "Duplicate county rows remain after filtering.\n",
      paste(
        capture.output(
          print(
            head(
              duplicates,
              20L
            )
          )
        ),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }
  
  
  x[]
}


# =============================================================================
# LOAD STATES AND COUNTIES
# =============================================================================

.load_states <- function() {
  
  .log(
    "Loading states."
  )
  
  states <- tigris::states(
    cb = TRUE,
    year = 2023,
    class = "sf"
  )
  
  states <- states[
    !states$STUSPS %in% c(
      "AK",
      "HI",
      "PR",
      "VI",
      "GU",
      "MP",
      "AS"
    ),
  ]
  
  states <- sf::st_make_valid(
    states
  )
  
  states <- sf::st_transform(
    states,
    CRS_PLOT
  )
  
  states
}


.load_counties <- function() {
  
  .log(
    "Loading counties."
  )
  
  counties <- tigris::counties(
    cb = TRUE,
    year = 2023,
    class = "sf"
  )
  
  counties <- counties[
    !counties$STATEFP %in% c(
      "02",
      "15",
      "60",
      "66",
      "69",
      "72",
      "78"
    ),
  ]
  
  counties <- sf::st_make_valid(
    counties
  )
  
  counties <- sf::st_transform(
    counties,
    CRS_PLOT
  )
  
  bbox_polygon <- sf::st_as_sfc(
    sf::st_bbox(
      c(
        xmin = BBOX[["xmin"]],
        ymin = BBOX[["ymin"]],
        xmax = BBOX[["xmax"]],
        ymax = BBOX[["ymax"]]
      ),
      crs = sf::st_crs(
        CRS_PLOT
      )
    )
  )
  
  keep <- lengths(
    sf::st_intersects(
      counties,
      bbox_polygon
    )
  ) > 0L
  
  counties <- counties[
    keep,
  ]
  
  counties$GEOID <- .normalize_geoid(
    counties$GEOID
  )
  
  counties <- counties[
    !is.na(
      counties$GEOID
    ),
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
      "County geometry contains duplicate GEOIDs.",
      call. = FALSE
    )
  }
  
  counties
}


# =============================================================================
# PREPARE METRICS FOR PLOTTING
# =============================================================================

.prepare_spatial_metrics <- function(
    metrics,
    counties
) {
  
  # ---------------------------------------------------------------------------
  # Convert F1, recall, and precision from columns to rows so they can be
  # displayed as three map panels with a common 0-1 scale.
  # ---------------------------------------------------------------------------
  
  long <- data.table::melt(
    metrics[
      ,
      .(
        GEOID,
        f1,
        recall,
        precision
      )
    ],
    id.vars = "GEOID",
    measure.vars = c(
      "f1",
      "recall",
      "precision"
    ),
    variable.name = "metric",
    value.name = "metric_value"
  )
  
  long[
    ,
    metric := factor(
      metric,
      levels = c(
        "recall",
        "precision",
        "f1"
      ),
      labels = c(
        "Recall",
        "Precision",
        "F1"
      )
    )
  ]
  
  
  # Match the existing performance-spatial workflow:
  # exact zero and undefined values plot in grey.
  
  long[
    ,
    plot_value := data.table::fifelse(
      is.finite(
        metric_value
      ) &
        metric_value > 0,
      metric_value,
      NA_real_
    )
  ]
  
  
  # ---------------------------------------------------------------------------
  # Join metrics to county polygons
  # ---------------------------------------------------------------------------
  
  spatial <- merge(
    counties,
    long,
    by = "GEOID",
    all.x = FALSE,
    sort = FALSE
  )
  
  if (nrow(spatial) == 0L) {
    stop(
      "County/metric join produced zero rows.",
      call. = FALSE
    )
  }
  
  spatial
}


# =============================================================================
# BUILD FIGURE
# =============================================================================

.build_figure <- function(
    variable,
    palette_name,
    metrics,
    counties,
    states
) {
  
  spatial <- .prepare_spatial_metrics(
    metrics = metrics,
    counties = counties
  )
  
  fill_colors <- .metric_palette(
    palette_name = palette_name
  )
  
  
  p <- ggplot2::ggplot() +
    
    # -------------------------------------------------------------------------
  # White state background
  # -------------------------------------------------------------------------
  
  ggplot2::geom_sf(
    data = states,
    fill = "white",
    color = NA,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
    
    # -------------------------------------------------------------------------
  # County metric values
  # -------------------------------------------------------------------------
  
  ggplot2::geom_sf(
    data = spatial,
    mapping = ggplot2::aes(
      fill = plot_value
    ),
    color = ZERO_AND_NA_COLOR,
    linewidth = 0.05,
    show.legend = TRUE
  ) +
    
    # -------------------------------------------------------------------------
  # State boundaries
  # -------------------------------------------------------------------------
  
  ggplot2::geom_sf(
    data = states,
    fill = NA,
    color = "black",
    linewidth = 0.45,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
    
    # -------------------------------------------------------------------------
  # Three native RP5 metric panels
  # -------------------------------------------------------------------------
  
  ggplot2::facet_wrap(
    ~ metric,
    nrow = 1,
    drop = FALSE
  ) +
    
    # -------------------------------------------------------------------------
  # Common metric scale
  # -------------------------------------------------------------------------
  
  ggplot2::scale_fill_gradientn(
    colours = fill_colors,
    limits = METRIC_LIMITS,
    breaks = METRIC_BREAKS,
    labels = scales::number(
      METRIC_BREAKS,
      accuracy = 0.1
    ),
    na.value = ZERO_AND_NA_COLOR,
    oob = scales::squish,
    name = NULL,
    guide = ggplot2::guide_colorbar(
      direction = "horizontal",
      barwidth = grid::unit(
        8.0,
        "in"
      ),
      barheight = grid::unit(
        0.18,
        "in"
      ),
      ticks = TRUE,
      frame.colour = "black",
      frame.linewidth = 0.4
    )
  ) +
    
    # -------------------------------------------------------------------------
  # EPSG:4326 with longitude/latitude borders and labels
  # -------------------------------------------------------------------------
  
  ggplot2::coord_sf(
    crs = sf::st_crs(
      CRS_PLOT
    ),
    datum = sf::st_crs(
      CRS_PLOT
    ),
    xlim = c(
      BBOX[["xmin"]],
      BBOX[["xmax"]]
    ),
    ylim = c(
      BBOX[["ymin"]],
      BBOX[["ymax"]]
    ),
    expand = FALSE,
    label_axes = list(
      bottom = "E",
      left = "N"
    )
  ) +
    
    ggplot2::labs(
      x = NULL,
      y = NULL
    ) +
    
    # -------------------------------------------------------------------------
  # Dissertation map styling
  # -------------------------------------------------------------------------
  
  ggplot2::theme_bw(
    base_size = BASE_SIZE,
    base_family = "Arial"
  ) +
    
    ggplot2::theme(
      
      text = ggplot2::element_text(
        family = "Arial",
        size = BASE_SIZE,
        color = "black"
      ),
      
      # -----------------------------------------------------------------------
      # Latitude / longitude graticules
      # -----------------------------------------------------------------------
      
      panel.grid.major = ggplot2::element_line(
        color = "grey82",
        linewidth = 0.35
      ),
      
      panel.grid.minor = ggplot2::element_blank(),
      
      panel.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      
      # -----------------------------------------------------------------------
      # Longitude / latitude labels
      # -----------------------------------------------------------------------
      
      axis.title = ggplot2::element_blank(),
      
      axis.text = ggplot2::element_text(
        family = "Arial",
        size = 14,
        color = "grey25"
      ),
      
      axis.ticks = ggplot2::element_line(
        color = "grey25",
        linewidth = 0.3
      ),
      
      # -----------------------------------------------------------------------
      # Black panel border
      # -----------------------------------------------------------------------
      
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.55
      ),
      
      # -----------------------------------------------------------------------
      # Metric labels
      # -----------------------------------------------------------------------
      
      strip.background = ggplot2::element_blank(),
      
      strip.text = ggplot2::element_text(
        family = "Arial",
        face = "bold",
        size = 18,
        color = "black",
        margin = ggplot2::margin(
          b = 6
        )
      ),
      
      panel.spacing.x = grid::unit(
        5,
        "mm"
      ),
      
      # -----------------------------------------------------------------------
      # Legend
      # -----------------------------------------------------------------------
      
      legend.position = "bottom",
      
      legend.direction = "horizontal",
      
      legend.text = ggplot2::element_text(
        family = "Arial",
        size = 14
      ),
      
      legend.margin = ggplot2::margin(
        t = 2,
        r = 0,
        b = 0,
        l = 0
      ),
      
      legend.box.margin = ggplot2::margin(
        0,
        0,
        0,
        0
      ),
      
      # -----------------------------------------------------------------------
      # Figure margins
      # -----------------------------------------------------------------------
      
      plot.margin = ggplot2::margin(
        8,
        8,
        8,
        8
      )
    )
  
  p
}


# =============================================================================
# SAVE FIGURE
# =============================================================================

.save_png <- function(
    plot,
    path
) {
  
  ragg::agg_png(
    filename = path,
    width = FIG_WIDTH_IN,
    height = FIG_HEIGHT_IN,
    units = "in",
    res = DPI,
    background = "white"
  )
  
  print(
    plot
  )
  
  grDevices::dev.off()
  
  .log(
    "Wrote: %s",
    path
  )
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
  
  
  # ---------------------------------------------------------------------------
  # Common geometry
  # ---------------------------------------------------------------------------
  
  states <- .load_states()
  
  counties <- .load_counties()
  
  .log(
    "Evaluation domain contains %d counties.",
    nrow(counties)
  )
  
  
  # ---------------------------------------------------------------------------
  # Heat and precipitation
  # ---------------------------------------------------------------------------
  
  for (i in seq_len(
    nrow(DATASETS)
  )) {
    
    config <- DATASETS[
      i
    ]
    
    .log(
      "============================================================"
    )
    
    .log(
      "Processing %s / %s / RP%d",
      config$variable,
      config$dataset_tag,
      RETURN_PERIOD
    )
    
    
    metrics <- .read_metrics(
      variable_dir = config$variable_dir,
      combination_tag = config$combination_tag,
      dataset_tag = config$dataset_tag,
      variable = config$variable
    )
    
    
    .log(
      "Selected %d county metric rows.",
      nrow(metrics)
    )
    
    
    p <- .build_figure(
      variable = config$variable,
      palette_name = config$palette,
      metrics = metrics,
      counties = counties,
      states = states
    )
    
    
    output_file <- file.path(
      OUT_DIR,
      config$output_name
    )
    
    
    .save_png(
      plot = p,
      path = output_file
    )
  }
  
  
  .log(
    "============================================================"
  )
  
  .log(
    "Complete."
  )
  
  .log(
    "Outputs:"
  )
  
  .log(
    "%s",
    file.path(
      OUT_DIR,
      "dissertation_02_performance_spatial_heat.png"
    )
  )
  
  .log(
    "%s",
    file.path(
      OUT_DIR,
      "dissertation_02_performance_spatial_precipitation.png"
    )
  )
}


main()