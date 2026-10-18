#!/usr/bin/env Rscript

# scripts/fig_01_wwa_heat_precip_two_panel.R
# =============================================================================
# Build a two-panel county map of warning-day frequency:
#   (a) Excessive Heat Warning county-days (EH.W)
#   (b) Flash Flood Warning county-days (FF.W)
#
# The county-day counts are calculated directly from the raw Mesonet/IEM CSV
# files immediately before plotting.
#
# Each panel has:
#   - its own fill scale;
#   - its own maximum warning-day value;
#   - its own sequential color palette;
#   - its own centered legend beneath the panel.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(lubridate)
  library(ggplot2)
  library(ragg)
  library(tigris)
  library(fs)
  library(patchwork)
  library(scales)
  library(grid)
  library(RColorBrewer)
})

options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)

# =============================================================================
# User settings
# =============================================================================

ROOT <- Sys.getenv(
  "CISC_ROOT",
  "/proj/sebastianlab/users/hquintal/cisc"
)

MESONET_BASE <- Sys.getenv(
  "CISC_MESONET_DIR",
  file.path(
    ROOT,
    "data/in/mesonet"
  )
)

ZONE_SHP <- Sys.getenv(
  "CISC_ZONE_SHP",
  file.path(
    ROOT,
    "data/in/shapefile/nws_public_forecast_zone/z_18mr25.shp"
  )
)

COUNTY_SHP <- Sys.getenv(
  "CISC_COUNTY_SHP",
  ""
)

OUT_DIR <- file.path(
  ROOT,
  "figures/f04_mesonet_warning"
)

OUT_PNG <- file.path(
  OUT_DIR,
  "f04_mesonet_warning.png"
)

OUT_PDF <- file.path(
  OUT_DIR,
  "wwa_freq_panels_heat_vs_precip_county.pdf"
)

WRITE_PDF <- FALSE

YEAR_START <- 2002L
YEAR_END <- 2023L

START_DATE <- as.Date(
  sprintf(
    "%d-01-01",
    YEAR_START
  )
)

END_DATE <- as.Date(
  sprintf(
    "%d-12-31",
    YEAR_END
  )
)

# Southeast plotting/evaluation extent.
BBOX_LONLAT <- c(
  xmin = -95,
  xmax = -75,
  ymin = 25,
  ymax = 40
)

CRS_PLOT <- 4326

HEAT_PHENOM <- "EH"
HEAT_SIG <- "W"

PRECIP_PHENOM <- "FF"
PRECIP_SIG <- "W"

RES_MAIN <- 300L
WIDTH_IN <- 14
HEIGHT_IN <- 7.0

# Number of colors used in each stepped scale.
N_FILL_COLOURS <- 9L

# Horizontal legend dimensions for each individual panel.
LEGEND_BAR_WIDTH_IN <- 4.9
LEGEND_BAR_HEIGHT_IN <- 0.24

# =============================================================================
# Helpers
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

.parse_datetime <- function(x) {
  x <- trimws(
    as.character(x)
  )
  
  x[x == ""] <- NA_character_
  
  parsed <- suppressWarnings(
    lubridate::ymd_hms(
      x,
      quiet = TRUE,
      tz = "UTC",
      truncated = 3
    )
  )
  
  missing <- is.na(parsed) & !is.na(x)
  
  if (any(missing)) {
    parsed[missing] <- suppressWarnings(
      lubridate::ymd_hm(
        x[missing],
        quiet = TRUE,
        tz = "UTC"
      )
    )
  }
  
  parsed
}

.expand_days <- function(
    start_date,
    end_date
) {
  if (
    is.na(start_date) ||
    is.na(end_date) ||
    end_date < start_date
  ) {
    return(
      as.Date(character())
    )
  }
  
  seq.Date(
    start_date,
    end_date,
    by = "day"
  )
}

# Convert UGC identifiers to the STATE_ZONE form used by the forecast-zone
# shapefile:
#   ALZ019 -> AL019
#   OHC027 -> OH027
.drop_third_character <- function(x) {
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  x[x == ""] <- NA_character_
  
  output <- x
  
  valid <- !is.na(x) &
    nchar(x) >= 4L
  
  output[valid] <- paste0(
    substr(
      x[valid],
      1L,
      2L
    ),
    substr(
      x[valid],
      4L,
      nchar(x[valid])
    )
  )
  
  output
}

.normalize_geoid <- function(x) {
  x <- gsub(
    "[^0-9]",
    "",
    as.character(x)
  )
  
  x[nchar(x) == 0L] <- NA_character_
  
  output <- rep(
    NA_character_,
    length(x)
  )
  
  valid <- !is.na(x)
  
  output[valid] <- sprintf(
    "%05d",
    as.integer(x[valid])
  )
  
  output
}

.read_counties <- function() {
  if (nzchar(COUNTY_SHP)) {
    if (!file.exists(COUNTY_SHP)) {
      stop(
        "Missing CISC_COUNTY_SHP: ",
        COUNTY_SHP
      )
    }
    
    counties <- sf::st_read(
      COUNTY_SHP,
      quiet = TRUE
    )
  } else {
    counties <- tigris::counties(
      cb = TRUE,
      year = 2022,
      class = "sf"
    )
  }
  
  required <- c(
    "GEOID",
    "STATEFP"
  )
  
  missing <- setdiff(
    required,
    names(counties)
  )
  
  if (length(missing)) {
    stop(
      "County geometry is missing fields: ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }
  
  drop_statefp <- c(
    "02",
    "15",
    "60",
    "66",
    "69",
    "72",
    "78"
  )
  
  counties <- counties[
    !(as.character(counties$STATEFP) %in% drop_statefp),
  ]
  
  counties$GEOID <- .normalize_geoid(
    counties$GEOID
  )
  
  counties <- counties[
    !is.na(counties$GEOID),
  ]
  
  counties <- sf::st_make_valid(
    counties
  )
  
  counties <- sf::st_transform(
    counties,
    4326
  )
  
  bbox_geometry <- sf::st_as_sfc(
    sf::st_bbox(
      BBOX_LONLAT,
      crs = sf::st_crs(4326)
    )
  )
  
  counties <- counties[
    lengths(
      sf::st_intersects(
        counties,
        bbox_geometry
      )
    ) > 0L,
  ]
  
  counties <- counties[
    ,
    c(
      "GEOID",
      "STATEFP"
    )
  ]
  
  sf::st_transform(
    counties,
    CRS_PLOT
  )
}

.read_zones <- function() {
  if (!file.exists(ZONE_SHP)) {
    stop(
      "Missing forecast-zone shapefile: ",
      ZONE_SHP
    )
  }
  
  zones <- sf::st_read(
    ZONE_SHP,
    quiet = TRUE
  )
  
  if (!("STATE_ZONE" %in% names(zones))) {
    stop(
      "Forecast-zone shapefile is missing STATE_ZONE."
    )
  }
  
  zones$STATE_ZONE <- toupper(
    trimws(
      as.character(
        zones$STATE_ZONE
      )
    )
  )
  
  zones <- zones[
    !is.na(zones$STATE_ZONE) &
      zones$STATE_ZONE != "",
  ]
  
  zones <- sf::st_make_valid(
    zones
  )
  
  zones <- sf::st_transform(
    zones,
    4326
  )
  
  bbox_geometry <- sf::st_as_sfc(
    sf::st_bbox(
      BBOX_LONLAT,
      crs = sf::st_crs(4326)
    )
  )
  
  zones <- zones[
    lengths(
      sf::st_intersects(
        zones,
        bbox_geometry
      )
    ) > 0L,
  ]
  
  zones <- zones[
    ,
    "STATE_ZONE"
  ]
  
  sf::st_transform(
    zones,
    CRS_PLOT
  )
}

.build_zone_county_lookup <- function(
    zones,
    counties
) {
  .log(
    "Building forecast-zone-to-county intersection lookup"
  )
  
  intersections <- suppressWarnings(
    sf::st_join(
      zones,
      counties[, "GEOID"],
      join = sf::st_intersects,
      left = FALSE
    )
  )
  
  lookup <- data.table::as.data.table(
    sf::st_drop_geometry(
      intersections
    )
  )[
    ,
    .(
      STATE_ZONE = toupper(
        trimws(
          as.character(STATE_ZONE)
        )
      ),
      GEOID = .normalize_geoid(GEOID)
    )
  ]
  
  lookup <- unique(
    lookup[
      !is.na(STATE_ZONE) &
        STATE_ZONE != "" &
        !is.na(GEOID)
    ]
  )
  
  if (!nrow(lookup)) {
    stop(
      "The zone-to-county lookup contains no rows."
    )
  }
  
  data.table::setkey(
    lookup,
    STATE_ZONE,
    GEOID
  )
  
  .log(
    "Zone-to-county lookup: %s unique zone-county pairs",
    format(
      nrow(lookup),
      big.mark = ","
    )
  )
  
  lookup
}

.discover_mesonet_csvs <- function() {
  if (!fs::dir_exists(MESONET_BASE)) {
    stop(
      "Missing MESONET_BASE: ",
      MESONET_BASE
    )
  }
  
  files <- fs::dir_ls(
    MESONET_BASE,
    recurse = TRUE,
    type = "file",
    glob = "*.csv"
  )
  
  if (!length(files)) {
    stop(
      "No Mesonet CSV files found under: ",
      MESONET_BASE
    )
  }
  
  # Retain files from the requested period when a four-digit year is present
  # in the path. Files without an identifiable year remain eligible because
  # their records are filtered by ISSUED/EXPIRED dates after reading.
  path_year <- suppressWarnings(
    as.integer(
      sub(
        ".*(?:^|[^0-9])((?:19|20)[0-9]{2})(?:[^0-9]|$).*",
        "\\1",
        as.character(files),
        perl = TRUE
      )
    )
  )
  
  has_year <- !is.na(path_year)
  
  keep <- !has_year |
    (
      path_year >= YEAR_START &
        path_year <= YEAR_END
    )
  
  files <- files[keep]
  
  .log(
    "Mesonet CSV files selected: %s",
    format(
      length(files),
      big.mark = ","
    )
  )
  
  files
}

.read_warning_days_one_file <- function(
    file,
    wanted_pairs,
    valid_zone_keys
) {
  frame <- tryCatch(
    data.table::fread(
      file,
      select = c(
        "PHENOM",
        "SIG",
        "ISSUED",
        "EXPIRED",
        "NWS_UGC"
      ),
      showProgress = FALSE
    ),
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(frame)) {
    return(NULL)
  }
  
  required <- c(
    "PHENOM",
    "SIG",
    "ISSUED",
    "EXPIRED",
    "NWS_UGC"
  )
  
  if (!all(required %in% names(frame))) {
    return(NULL)
  }
  
  frame[
    ,
    `:=`(
      PHENOM = toupper(
        trimws(
          as.character(PHENOM)
        )
      ),
      SIG = toupper(
        trimws(
          as.character(SIG)
        )
      ),
      NWS_UGC = toupper(
        trimws(
          as.character(NWS_UGC)
        )
      )
    )
  ]
  
  frame <- frame[
    wanted_pairs,
    on = .(
      PHENOM,
      SIG
    ),
    nomatch = 0L
  ]
  
  frame <- frame[
    !is.na(NWS_UGC) &
      NWS_UGC != ""
  ]
  
  if (!nrow(frame)) {
    return(NULL)
  }
  
  frame[
    ,
    `:=`(
      issued_datetime = .parse_datetime(
        ISSUED
      ),
      expired_datetime = .parse_datetime(
        EXPIRED
      )
    )
  ]
  
  frame <- frame[
    !is.na(issued_datetime) &
      !is.na(expired_datetime) &
      expired_datetime >= issued_datetime
  ]
  
  if (!nrow(frame)) {
    return(NULL)
  }
  
  frame[
    ,
    `:=`(
      issued_date = as.Date(
        issued_datetime
      ),
      expired_date = as.Date(
        expired_datetime
      )
    )
  ]
  
  frame[
    ,
    `:=`(
      issued_date = pmax(
        issued_date,
        START_DATE
      ),
      expired_date = pmin(
        expired_date,
        END_DATE
      )
    )
  ]
  
  frame <- frame[
    expired_date >= issued_date
  ]
  
  if (!nrow(frame)) {
    return(NULL)
  }
  
  day_lists <- Map(
    .expand_days,
    frame$issued_date,
    frame$expired_date
  )
  
  output <- data.table::data.table(
    PHENOM = rep(
      frame$PHENOM,
      lengths(day_lists)
    ),
    SIG = rep(
      frame$SIG,
      lengths(day_lists)
    ),
    NWS_UGC = rep(
      frame$NWS_UGC,
      lengths(day_lists)
    ),
    date = as.Date(
      unlist(
        day_lists,
        use.names = FALSE
      ),
      origin = "1970-01-01"
    )
  )
  
  if (!nrow(output)) {
    return(NULL)
  }
  
  output[
    ,
    STATE_ZONE := .drop_third_character(
      NWS_UGC
    )
  ]
  
  output <- output[
    !is.na(STATE_ZONE) &
      STATE_ZONE %in% valid_zone_keys
  ]
  
  if (!nrow(output)) {
    return(NULL)
  }
  
  unique(
    output[
      ,
      .(
        PHENOM,
        SIG,
        STATE_ZONE,
        date
      )
    ]
  )
}

.build_warning_zone_days <- function(
    zone_county_lookup
) {
  files <- .discover_mesonet_csvs()
  
  wanted_pairs <- data.table::data.table(
    PHENOM = c(
      HEAT_PHENOM,
      PRECIP_PHENOM
    ),
    SIG = c(
      HEAT_SIG,
      PRECIP_SIG
    )
  )
  
  valid_zone_keys <- unique(
    zone_county_lookup$STATE_ZONE
  )
  
  parts <- vector(
    "list",
    length(files)
  )
  
  for (index in seq_along(files)) {
    if (
      index == 1L ||
      index %% 100L == 0L ||
      index == length(files)
    ) {
      .log(
        "Reading Mesonet file %s of %s",
        format(
          index,
          big.mark = ","
        ),
        format(
          length(files),
          big.mark = ","
        )
      )
    }
    
    parts[[index]] <- .read_warning_days_one_file(
      files[[index]],
      wanted_pairs,
      valid_zone_keys
    )
  }
  
  parts <- Filter(
    Negate(is.null),
    parts
  )
  
  if (!length(parts)) {
    stop(
      paste0(
        "No EH.W or FF.W zone-day records were produced ",
        "from the Mesonet files."
      )
    )
  }
  
  warning_zone_days <- data.table::rbindlist(
    parts,
    use.names = TRUE,
    fill = TRUE
  )
  
  warning_zone_days <- unique(
    warning_zone_days[
      date >= START_DATE &
        date <= END_DATE
    ]
  )
  
  data.table::setkey(
    warning_zone_days,
    PHENOM,
    SIG,
    STATE_ZONE,
    date
  )
  
  .log(
    "Unique warning zone-days: %s",
    format(
      nrow(warning_zone_days),
      big.mark = ","
    )
  )
  
  warning_zone_days
}

.count_warning_county_days <- function(
    warning_zone_days,
    zone_county_lookup
) {
  # A warning zone may intersect several counties, and a county may intersect
  # several forecast zones. Merge first and then deduplicate at
  # PHENOM/SIG/GEOID/date so each county contributes at most one count per day.
  county_days <- merge(
    warning_zone_days,
    zone_county_lookup,
    by = "STATE_ZONE",
    allow.cartesian = TRUE,
    sort = FALSE
  )
  
  county_days <- unique(
    county_days[
      ,
      .(
        PHENOM,
        SIG,
        GEOID,
        date
      )
    ]
  )
  
  counts <- county_days[
    ,
    .(
      warning_days = uniqueN(date)
    ),
    by = .(
      PHENOM,
      SIG,
      GEOID
    )
  ]
  
  expected <- data.table::data.table(
    PHENOM = c(
      HEAT_PHENOM,
      PRECIP_PHENOM
    ),
    SIG = c(
      HEAT_SIG,
      PRECIP_SIG
    )
  )
  
  missing_pairs <- expected[
    !unique(
      counts[
        ,
        .(
          PHENOM,
          SIG
        )
      ]
    ),
    on = .(
      PHENOM,
      SIG
    )
  ]
  
  if (nrow(missing_pairs)) {
    stop(
      "No county-day counts were found for: ",
      paste(
        paste0(
          missing_pairs$PHENOM,
          ".",
          missing_pairs$SIG
        ),
        collapse = ", "
      )
    )
  }
  
  .log(
    "Unique warning county-days: %s",
    format(
      nrow(county_days),
      big.mark = ","
    )
  )
  
  list(
    county_days = county_days,
    counts = counts
  )
}

.pretty_breaks_integer <- function(
    limits,
    n = 7L
) {
  breaks <- scales::breaks_pretty(
    n = n
  )(limits)
  
  breaks <- unique(
    as.integer(
      round(
        breaks[
          breaks >= limits[[1L]] &
            breaks <= limits[[2L]]
        ]
      )
    )
  )
  
  if (
    !length(breaks) ||
    breaks[[1L]] != limits[[1L]]
  ) {
    breaks <- c(
      as.integer(limits[[1L]]),
      breaks
    )
  }
  
  if (
    tail(
      breaks,
      1L
    ) != limits[[2L]]
  ) {
    breaks <- c(
      breaks,
      as.integer(limits[[2L]])
    )
  }
  
  unique(breaks)
}

.theme_map <- function() {
  ggplot2::theme_minimal(
    base_size = 18,
    base_family = "Arial"
  ) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = "Arial",
        size = 18
      ),
      
      panel.grid.major = ggplot2::element_line(
        color = "grey75",
        linewidth = 0.25
      ),
      
      panel.grid.minor = ggplot2::element_blank(),
      
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.5
      ),
      
      axis.title = ggplot2::element_blank(),
      
      axis.text = ggplot2::element_text(
        color = "grey25",
        size = 14
      ),
      
      axis.ticks = ggplot2::element_line(
        color = "grey25",
        linewidth = 0.3
      ),
      
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.box.just = "center",
      legend.box = "vertical",
      
      legend.title = ggplot2::element_text(
        size = 17,
        hjust = 0.5
      ),
      
      legend.text = ggplot2::element_text(
        size = 14
      ),
      
      legend.margin = ggplot2::margin(
        t = 2,
        r = 0,
        b = 0,
        l = 0
      ),
      
      legend.box.margin = ggplot2::margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      ),
      
      plot.tag = ggplot2::element_text(
        size = 18,
        face = "bold"
      ),
      
      # Lower than the original 0.985 setting.
      plot.tag.position = c(
        0.015,
        0.965
      ),
      
      plot.margin = ggplot2::margin(
        t = 3,
        r = 5,
        b = 3,
        l = 5
      )
    )
}

.make_panel <- function(
    county_sf,
    state_sf,
    fill_limits,
    fill_breaks,
    fill_colours,
    legend_title,
    panel_tag
) {
  county_sf$plot_warning_days <- ifelse(
    county_sf$warning_days > 0,
    county_sf$warning_days,
    NA_real_
  )
  
  ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = county_sf,
      ggplot2::aes(
        fill = plot_warning_days
      ),
      color = "grey70",
      linewidth = 0.08
    ) +
    
    ggplot2::geom_sf(
      data = state_sf,
      fill = NA,
      color = "black",
      linewidth = 0.45
    ) +
    
    ggplot2::scale_fill_stepsn(
      colours = fill_colours,
      limits = fill_limits,
      breaks = fill_breaks,
      labels = scales::label_number(
        accuracy = 1,
        big.mark = ","
      ),
      oob = scales::squish,
      na.value = "grey85",
      name = legend_title,
      guide = ggplot2::guide_colorsteps(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        barwidth = grid::unit(
          LEGEND_BAR_WIDTH_IN,
          "in"
        ),
        barheight = grid::unit(
          LEGEND_BAR_HEIGHT_IN,
          "in"
        ),
        even.steps = TRUE,
        show.limits = TRUE
      )
    ) +
    
    ggplot2::coord_sf(
      crs = sf::st_crs(4326),
      datum = sf::st_crs(4326),
      xlim = c(
        BBOX_LONLAT[["xmin"]],
        BBOX_LONLAT[["xmax"]]
      ),
      ylim = c(
        BBOX_LONLAT[["ymin"]],
        BBOX_LONLAT[["ymax"]]
      ),
      expand = FALSE,
      label_axes = list(
        bottom = "E",
        left = "N"
      )
    ) +
    
    ggplot2::labs(
      tag = panel_tag
    ) +
    
    .theme_map()
}

# =============================================================================
# Build data immediately before plotting
# =============================================================================

fs::dir_create(
  OUT_DIR,
  recurse = TRUE
)

.log(
  "Building county warning-day maps for %d-%d",
  YEAR_START,
  YEAR_END
)

counties <- .read_counties()

zones <- .read_zones()

zone_county_lookup <- .build_zone_county_lookup(
  zones,
  counties
)

warning_zone_days <- .build_warning_zone_days(
  zone_county_lookup
)

warning_results <- .count_warning_county_days(
  warning_zone_days,
  zone_county_lookup
)

warning_counts <- warning_results$counts

heat_counts <- warning_counts[
  PHENOM == HEAT_PHENOM &
    SIG == HEAT_SIG,
  .(
    GEOID,
    warning_days
  )
]

precip_counts <- warning_counts[
  PHENOM == PRECIP_PHENOM &
    SIG == PRECIP_SIG,
  .(
    GEOID,
    warning_days
  )
]

county_heat <- merge(
  counties,
  heat_counts,
  by = "GEOID",
  all.x = TRUE,
  sort = FALSE
)

county_precip <- merge(
  counties,
  precip_counts,
  by = "GEOID",
  all.x = TRUE,
  sort = FALSE
)

county_heat$warning_days[
  is.na(
    county_heat$warning_days
  )
] <- 0L

county_precip$warning_days[
  is.na(
    county_precip$warning_days
  )
] <- 0L

states <- counties |>
  dplyr::group_by(
    STATEFP
  ) |>
  dplyr::summarise(
    geometry = sf::st_union(geometry),
    .groups = "drop"
  )

# =============================================================================
# Independent panel scales
# =============================================================================

heat_positive <- county_heat$warning_days[
  county_heat$warning_days > 0
]

precip_positive <- county_precip$warning_days[
  county_precip$warning_days > 0
]

if (!length(heat_positive)) {
  stop(
    "No counties have Excessive Heat Warning days greater than zero."
  )
}

if (!length(precip_positive)) {
  stop(
    "No counties have Flash Flood Warning days greater than zero."
  )
}

max_heat <- max(
  heat_positive,
  na.rm = TRUE
)

max_precip <- max(
  precip_positive,
  na.rm = TRUE
)

# Both scales begin at one warning day, because zero-warning counties are shown
# as white rather than assigned to the lowest color class.
heat_fill_limits <- c(
  1L,
  max_heat
)

precip_fill_limits <- c(
  1L,
  max_precip
)

# Heat panel: 25-day bins
heat_fill_breaks <- unique(c(
  1L,
  seq(
    25L,
    ceiling(max_heat / 25) * 25L,
    by = 25L
  ),
  max_heat
))

# Remove any break above the maximum
heat_fill_breaks <- heat_fill_breaks[
  heat_fill_breaks <= max_heat
]

# Ensure the maximum is always included
if (tail(heat_fill_breaks, 1) != max_heat) {
  heat_fill_breaks <- c(
    heat_fill_breaks,
    max_heat
  )
}

precip_fill_breaks <- .pretty_breaks_integer(
  precip_fill_limits,
  n = 7L
)

# Reverse these palettes so lower warning-day counts are pale and higher counts
# are dark.
heat_colours <- RColorBrewer::brewer.pal(
  n = 9,
  name = "Reds"
)

precip_colours <- RColorBrewer::brewer.pal(
  n = 9,
  name = "Blues"
)

.log(
  "Heat warning days: maximum county count = %s",
  format(
    max_heat,
    big.mark = ","
  )
)

.log(
  "Heat legend limits: %s to %s",
  format(
    heat_fill_limits[[1L]],
    big.mark = ","
  ),
  format(
    heat_fill_limits[[2L]],
    big.mark = ","
  )
)

.log(
  "Precipitation warning days: maximum county count = %s",
  format(
    max_precip,
    big.mark = ","
  )
)

.log(
  "Precipitation legend limits: %s to %s",
  format(
    precip_fill_limits[[1L]],
    big.mark = ","
  ),
  format(
    precip_fill_limits[[2L]],
    big.mark = ","
  )
)

county_heat <- sf::st_transform(
  county_heat,
  4326
)

county_precip <- sf::st_transform(
  county_precip,
  4326
)

states <- sf::st_transform(
  states,
  4326
)

# =============================================================================
# Plot
# =============================================================================

p_heat <- .make_panel(
  county_sf = county_heat,
  state_sf = states,
  fill_limits = heat_fill_limits,
  fill_breaks = heat_fill_breaks,
  fill_colours = heat_colours,
  legend_title = "Excessive Heat Warning days",
  panel_tag = "(a)"
)

p_precip <- .make_panel(
  county_sf = county_precip,
  state_sf = states,
  fill_limits = precip_fill_limits,
  fill_breaks = precip_fill_breaks,
  fill_colours = precip_colours,
  legend_title = "Flash Flood Warning days",
  panel_tag = "(b)"
)

# Do not collect guides. Each panel retains its own independently centered
# legend beneath its own map.
figure <- p_heat + p_precip +
  patchwork::plot_layout(
    ncol = 2L,
    widths = c(
      1,
      1
    ),
    guides = "keep"
  ) +
  patchwork::plot_annotation(
    theme = ggplot2::theme(
      text = ggplot2::element_text(
        family = "Arial",
        size = 18
      )
    )
  )

ragg::agg_png(
  filename = OUT_PNG,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  units = "in",
  res = RES_MAIN,
  background = "white"
)

print(figure)

grDevices::dev.off()

.log(
  "Wrote: %s",
  OUT_PNG
)

if (WRITE_PDF) {
  grDevices::cairo_pdf(
    filename = OUT_PDF,
    width = WIDTH_IN,
    height = HEIGHT_IN,
    family = "Arial",
    onefile = TRUE
  )
  
  print(figure)
  
  grDevices::dev.off()
  
  .log(
    "Wrote: %s",
    OUT_PDF
  )
}

.log(
  "Done."
)