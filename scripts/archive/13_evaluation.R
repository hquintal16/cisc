suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(fs)
  library(sf)
  library(terra)
  library(tigris)
  library(future)
  library(future.apply)
  library(progressr)
  library(lubridate)
  library(stringr)
})

# =============================================================================
# County-level WWA evaluation workflow for clustered ERA5 heat index + precipitation
# =============================================================================
# WHAT THIS SCRIPT DOES
#   1) Preflight: discovers existing modeled, reference-grid, and observed WWA files
#   2) Builds / updates county lookup tables for every dataset grid used in evaluation
#   3) Builds definitive WWA county-day observed products from forecast-zone day products
#   4) Builds modeled county-day products from clustered outputs
#   5) Computes county confusion counts + metrics across multiple windows
#
# IMPORTANT
#   - This script assumes your current WWA workflow has already written sparse
#     zone-day .qs files like: wwa_freq_<PHENOM>_<SIG>_<...>.qs
#   - Evaluation is against WWA at COUNTY resolution, not NOAA Storm Events.
#   - Modeled source is the clustered outputs in data/out/06_clusters.
#   - Lookup tables are built from the reference-grid TIFF for each dataset_tag.
# =============================================================================

# ----------------------------
# USER SETTINGS
# ----------------------------
ROOT <- "V:/users/hquintal/phd1_cluster_southeast_v2"

# Current observed WWA products from your existing script
WWA_DIR <- file.path(ROOT, "figures/fig_01_wwa")
ZONE_SHP <- file.path(ROOT, "data/in/shapefile/nws_public_forecast_zone/z_18mr25.shp")

# Optional local county shapefile. If NA, uses tigris::counties().
COUNTY_SHP <- NA_character_
COUNTY_YEAR <- 2022L
COUNTY_CB <- TRUE

# Bounding box used throughout the Southeast workflow
BBOX_LONLAT <- c(xmin = -95, xmax = -75, ymin = 22, ymax = 40)

# Clustered modeled events
CLUST_BASE <- file.path(ROOT, "data/out/06_clusters")
SCENARIO_PATTERN <- "^hi\\d+_pr\\d+$"

# Scenario folders to evaluate. Set NULL to auto-discover all.
SCENARIOS_KEEP <- NULL

# Current reference grids already created upstream
REF_HEAT_DIR <- file.path(ROOT, "data/out/03_cov_model/heatindex")
REF_PREC_DIR <- file.path(ROOT, "data/out/03_cov_model/precipitation")

# Performance output base
PERF_BASE <- file.path(ROOT, "data/out/09_performance")
LOOKUP_XY_TO_COUNTY_DIR <- file.path(PERF_BASE, "lookup_xy_to_county")
LOOKUP_ZONE_TO_COUNTY_DIR <- file.path(PERF_BASE, "lookup_zone_to_county")
OBS_DIR <- file.path(PERF_BASE, "wwa_county_days_definitive")
MODELED_DIR <- file.path(PERF_BASE, "modeled_county_daily")
MET_DIR <- file.path(PERF_BASE, "metrics_county")
PART_DIR <- file.path(MET_DIR, "parts")
QA_DIR <- file.path(PERF_BASE, "qa_wwa_eval")

dir_create(c(LOOKUP_XY_TO_COUNTY_DIR, LOOKUP_ZONE_TO_COUNTY_DIR, OBS_DIR, MODELED_DIR, MET_DIR, PART_DIR, QA_DIR), recurse = TRUE)

# Study date domain for evaluation
YEAR_START <- 2002L
YEAR_END <- 2023L
GLOBAL_START <- as.Date(sprintf("%d-01-01", YEAR_START))
GLOBAL_END <- as.Date(sprintf("%d-12-31", YEAR_END))

FIXED_WINDOWS <- list(
  list(tag = "20020101_20231231", start = as.Date("2002-01-01"), end = as.Date("2023-12-31")),
  list(tag = "20140101_20231231", start = as.Date("2014-01-01"), end = as.Date("2023-12-31")),
  list(tag = "20190101_20231231", start = as.Date("2019-01-01"), end = as.Date("2023-12-31"))
)

# WWA grouping rules used for evaluation keys
# Edit these if you want stricter or broader definitions.
# HCQ EDIT ----
WWA_GROUPS <- list(
  all_heat = c("EH_W"),
  all_precip = c("SV_W",
                 "FF_W", "HI_W", "HU_W", "TI_W", "TR_W")
)

# Optional: also evaluate all individual WWA tags found in WWA_DIR
EVALUATE_INDIVIDUAL_WWA_TAGS <- FALSE

# Workers
WORKERS_LOOKUP <- 8L
WORKERS_MODEL <- 16L
WORKERS_EVAL <- 16L

# Behavior
OVERWRITE_LOOKUPS <- FALSE
OVERWRITE_OBSERVED <- FALSE
OVERWRITE_MODELED <- FALSE
OVERWRITE_METRICS <- TRUE
RUN_PREFLIGHT_ONLY <- FALSE

# ----------------------------
# HELPERS
# ----------------------------
.log <- function(fmt, ...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))

.write_atomic_qs <- function(obj, out_path, preset = "high") {
  dir_create(path_dir(out_path), recurse = TRUE)
  tmp <- paste0(out_path, ".tmp")
  if (file_exists(tmp)) file_delete(tmp)
  qsave(obj, tmp, preset = preset)
  file_move(tmp, out_path)
  invisible(out_path)
}

.write_atomic_csv <- function(dt, out_path) {
  dir_create(path_dir(out_path), recurse = TRUE)
  tmp <- paste0(out_path, ".tmp")
  if (file_exists(tmp)) file_delete(tmp)
  fwrite(dt, tmp)
  file_move(tmp, out_path)
  invisible(out_path)
}

.norm_geoid5 <- function(x) {
  x <- as.character(x)
  x <- str_extract(x, "\\d+")
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x) & nchar(x) <= 5
  out[ok] <- str_pad(x[ok], width = 5, pad = "0")
  out
}

.slug <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

.event_tag_norm <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

.safe_match_first <- function(x, patterns) {
  hit <- patterns[patterns %in% names(x)]
  if (!length(hit)) return(NA_character_)
  hit[1]
}

.parse_scenario_floors <- function(scenario) {
  m <- str_match(scenario, "^hi(\\d+)_pr(\\d+)$")
  if (any(is.na(m))) return(list(hi_floor = NA_real_, pr_floor = NA_real_))
  list(hi_floor = as.numeric(m[1, 2]), pr_floor = as.numeric(m[1, 3]))
}

.build_windows <- function() {
  yearly <- rbindlist(lapply(YEAR_START:YEAR_END, function(y) {
    data.table(
      window_tag = sprintf("%d", y),
      window_type = "year",
      year = y,
      start = as.Date(sprintf("%d-01-01", y)),
      end = as.Date(sprintf("%d-12-31", y))
    )
  }))
  
  fixed <- rbindlist(lapply(FIXED_WINDOWS, function(w) {
    data.table(
      window_tag = w$tag,
      window_type = "fixed",
      year = NA_integer_,
      start = w$start,
      end = w$end
    )
  }))
  
  out <- rbindlist(list(yearly, fixed), use.names = TRUE, fill = TRUE)
  out[, total_days := as.integer(end - start) + 1L]
  out[]
}

WINDOWS <- .build_windows()

# ----------------------------
# DISCOVERY HELPERS
# ----------------------------
.discover_scenarios <- function() {
  dirs <- dir_ls(CLUST_BASE, type = "directory", fail = FALSE)
  dirs <- dirs[grepl(SCENARIO_PATTERN, path_file(dirs))]
  dirs <- sort(as.character(dirs))
  if (!is.null(SCENARIOS_KEEP)) {
    dirs <- dirs[path_file(dirs) %in% SCENARIOS_KEEP]
  }
  dirs
}

.discover_dataset_dirs <- function(scenario_dir) {
  ds <- dir_ls(scenario_dir, type = "directory", fail = FALSE)
  sort(as.character(ds))
}

.discover_cluster_files <- function(dataset_dir) {
  files <- dir_ls(dataset_dir, type = "file", glob = "*.qs", fail = FALSE)
  sort(as.character(files))
}

# ----------------------------
# REFERENCE GRID DISCOVERY
# ----------------------------
.find_ref_tif_for_dataset <- function(dataset_tag) {
  var_group <- if (grepl("^heatindex", dataset_tag, ignore.case = TRUE)) "heatindex" else "precipitation"
  ref_base <- if (var_group == "heatindex") REF_HEAT_DIR else REF_PREC_DIR
  
  # preferred folder names based on user's current structure
  cand1 <- file.path(ref_base, dataset_tag, paste0("ref_grid_", dataset_tag, ".tif"))
  cand2 <- file.path(ref_base, dataset_tag, paste0(dataset_tag, ".tif"))
  cand3 <- file.path(ref_base, dataset_tag, paste0("ref_grid_", .slug(dataset_tag), ".tif"))
  
  cands <- c(cand1, cand2, cand3)
  hit <- cands[file_exists(cands)]
  if (length(hit)) return(hit[1])
  
  # broader recursive fallback
  all_tifs <- dir_ls(ref_base, recurse = TRUE, glob = "*.tif", type = "file", fail = FALSE)
  bn <- path_file(all_tifs)
  keep <- grepl(gsub("([.])", "\\\\1", dataset_tag), bn, fixed = FALSE)
  hit <- all_tifs[keep]
  if (length(hit)) return(as.character(hit[1]))
  
  NA_character_
}

.discover_dataset_registry <- function() {
  sc_dirs <- .discover_scenarios()
  if (!length(sc_dirs)) stop("No scenario directories found in: ", CLUST_BASE)
  
  reg <- rbindlist(lapply(sc_dirs, function(sc_dir) {
    scenario <- path_file(sc_dir)
    ds_dirs <- .discover_dataset_dirs(sc_dir)
    if (!length(ds_dirs)) return(NULL)
    
    rbindlist(lapply(ds_dirs, function(ds_dir) {
      dataset_tag <- path_file(ds_dir)
      floors <- .parse_scenario_floors(scenario)
      data.table(
        scenario = scenario,
        dataset_tag = dataset_tag,
        var_group = if (grepl("^heatindex", dataset_tag, ignore.case = TRUE)) "heatindex" else "precipitation",
        scenario_dir = sc_dir,
        dataset_dir = ds_dir,
        ref_tif = .find_ref_tif_for_dataset(dataset_tag),
        hi_floor = floors$hi_floor,
        pr_floor = floors$pr_floor
      )
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE)
  
  reg[, ref_exists := file_exists(ref_tif)]
  reg[]
}

# ----------------------------
# COUNTY / ZONE GEOMETRY
# ----------------------------
.read_counties_sf <- function() {
  if (!is.na(COUNTY_SHP) && file_exists(COUNTY_SHP)) {
    cty <- st_read(COUNTY_SHP, quiet = TRUE)
  } else {
    cty <- tigris::counties(cb = COUNTY_CB, year = COUNTY_YEAR, class = "sf")
  }
  
  stopifnot(all(c("GEOID", "STATEFP") %in% names(cty)))
  
  # Drop AK/HI + territories first
  drop_statefp <- c("02", "15", "60", "66", "69", "72", "78")
  cty <- cty[!(cty$STATEFP %in% drop_statefp), ]
  
  # Clean geometry
  cty <- st_make_valid(cty)
  cty <- st_transform(cty, 4326)
  
  # Build lon/lat bbox polygon
  bb <- st_as_sfc(
    st_bbox(
      c(
        xmin = unname(BBOX_LONLAT["xmin"]),
        xmax = unname(BBOX_LONLAT["xmax"]),
        ymin = unname(BBOX_LONLAT["ymin"]),
        ymax = unname(BBOX_LONLAT["ymax"])
      ),
      crs = st_crs(4326)
    )
  )
  
  # Use crop first, then intersection fallback if needed
  cty_crop <- tryCatch(
    st_crop(cty, st_bbox(bb)),
    error = function(e) NULL
  )
  
  if (is.null(cty_crop) || nrow(cty_crop) == 0) {
    hit <- tryCatch(
      st_intersects(st_geometry(cty), st_geometry(bb), sparse = FALSE)[, 1],
      error = function(e) NULL
    )
    
    if (is.null(hit)) {
      stop("Failed to subset counties to bbox. Check county geometries and bbox.")
    }
    
    cty_crop <- cty[hit, ]
  }
  
  # Final cleanup
  cty_crop <- cty_crop[!is.na(cty_crop$GEOID) & nzchar(cty_crop$GEOID), ]
  cty_crop
}

.read_zones_sf <- function() {
  if (!file_exists(ZONE_SHP)) stop("Missing ZONE_SHP: ", ZONE_SHP)
  
  z <- st_read(ZONE_SHP, quiet = TRUE)
  if (!("STATE_ZONE" %in% names(z))) stop("ZONE_SHP missing STATE_ZONE field.")
  
  z$STATE_ZONE <- toupper(trimws(as.character(z$STATE_ZONE)))
  z <- st_make_valid(z)
  z <- st_transform(z, 4326)
  
  bb <- st_as_sfc(
    st_bbox(
      c(
        xmin = unname(BBOX_LONLAT["xmin"]),
        xmax = unname(BBOX_LONLAT["xmax"]),
        ymin = unname(BBOX_LONLAT["ymin"]),
        ymax = unname(BBOX_LONLAT["ymax"])
      ),
      crs = st_crs(4326)
    )
  )
  
  z_crop <- tryCatch(
    st_crop(z, st_bbox(bb)),
    error = function(e) NULL
  )
  
  if (is.null(z_crop) || nrow(z_crop) == 0) {
    hit <- tryCatch(
      st_intersects(st_geometry(z), st_geometry(bb), sparse = FALSE)[, 1],
      error = function(e) NULL
    )
    
    if (is.null(hit)) {
      stop("Failed to subset zones to bbox. Check zone geometries and bbox.")
    }
    
    z_crop <- z[hit, ]
  }
  
  z_crop <- z_crop[!is.na(z_crop$STATE_ZONE) & nzchar(z_crop$STATE_ZONE), ]
  z_crop
}

COUNTIES_SF <- .read_counties_sf()
COUNTY_UNIVERSE <- sort(unique(as.character(COUNTIES_SF$GEOID)))
ZONES_SF <- .read_zones_sf()

# ----------------------------
# PREFLIGHT
# ----------------------------
run_preflight <- function(dataset_registry) {
  .log("Running preflight checks...")
  
  # WWA observed products
  wwa_qs <- dir_ls(WWA_DIR, type = "file", glob = "wwa_freq_*.qs", fail = FALSE)
  wwa_qs <- sort(as.character(wwa_qs))
  
  # Existing modeled county outputs
  modeled_qs <- dir_ls(MODELED_DIR, recurse = TRUE, type = "file", glob = "*.qs", fail = FALSE)
  modeled_qs <- sort(as.character(modeled_qs))
  
  # Existing XY->county lookups
  xy_lookup_qs <- dir_ls(LOOKUP_XY_TO_COUNTY_DIR, type = "file", glob = "lookup_xy_to_county_*.qs", fail = FALSE)
  xy_lookup_qs <- sort(as.character(xy_lookup_qs))
  
  # Existing zone->county crosswalk
  zone_lookup_qs <- dir_ls(LOOKUP_ZONE_TO_COUNTY_DIR, type = "file", glob = "*.qs", fail = FALSE)
  zone_lookup_qs <- sort(as.character(zone_lookup_qs))
  
  # Existing observed county-day files
  obs_qs <- dir_ls(OBS_DIR, type = "file", glob = "*.qs", fail = FALSE)
  obs_qs <- sort(as.character(obs_qs))
  
  # Existing metrics
  met_qs <- dir_ls(MET_DIR, type = "file", glob = "metrics_*.qs", fail = FALSE)
  met_qs <- sort(as.character(met_qs))
  
  # Dataset-level needs
  need_lookup <- dataset_registry[, .(dataset_tag, ref_tif, ref_exists)]
  need_lookup[, lookup_qs := file.path(LOOKUP_XY_TO_COUNTY_DIR, paste0("lookup_xy_to_county_", dataset_tag, ".qs"))]
  need_lookup[, lookup_exists := file_exists(lookup_qs)]
  
  # QA on lookup completeness if file exists
  qa <- rbindlist(lapply(seq_len(nrow(need_lookup)), function(i) {
    row <- need_lookup[i]
    if (!row$lookup_exists || !row$ref_exists) {
      return(data.table(
        dataset_tag = row$dataset_tag,
        ref_exists = row$ref_exists,
        lookup_exists = row$lookup_exists,
        n_ref_cells = NA_integer_,
        n_lookup_rows = NA_integer_,
        n_lookup_unique_xy = NA_integer_,
        n_missing_assignments = NA_integer_
      ))
    }
    
    r <- terra::rast(row$ref_tif)
    df <- as.data.table(as.data.frame(r, xy = TRUE, cells = TRUE, na.rm = FALSE))
    setnames(df, old = c("x", "y", "cell"), new = c("x", "y", "cell_id"), skip_absent = TRUE)
    df <- df[, .(cell_id = as.integer(cell_id), x = as.numeric(x), y = as.numeric(y))]
    
    lk <- as.data.table(qread(row$lookup_qs))
    if (!all(c("cell_id", "GEOID") %in% names(lk))) {
      return(data.table(
        dataset_tag = row$dataset_tag,
        ref_exists = TRUE,
        lookup_exists = TRUE,
        n_ref_cells = nrow(df),
        n_lookup_rows = nrow(lk),
        n_lookup_unique_xy = NA_integer_,
        n_missing_assignments = NA_integer_
      ))
    }
    
    lk[, cell_id := as.integer(cell_id)]
    merged <- merge(df, lk[, .(cell_id, GEOID)], by = "cell_id", all.x = TRUE)
    
    data.table(
      dataset_tag = row$dataset_tag,
      ref_exists = TRUE,
      lookup_exists = TRUE,
      n_ref_cells = nrow(df),
      n_lookup_rows = nrow(lk),
      n_lookup_unique_xy = uniqueN(lk$cell_id),
      n_missing_assignments = sum(is.na(merged$GEOID) | !nzchar(as.character(merged$GEOID)))
    )
  }), use.names = TRUE, fill = TRUE)
  
  .write_atomic_csv(need_lookup, file.path(QA_DIR, "preflight_dataset_lookup_status.csv"))
  .write_atomic_csv(qa, file.path(QA_DIR, "preflight_lookup_completeness.csv"))
  
  summary_txt <- c(
    sprintf("WWA sparse zone-day files found: %d", length(wwa_qs)),
    sprintf("Scenario x dataset rows discovered: %d", nrow(dataset_registry)),
    sprintf("Distinct dataset tags discovered: %d", uniqueN(dataset_registry$dataset_tag)),
    sprintf("Reference TIFFs found: %d / %d", sum(dataset_registry$ref_exists), nrow(dataset_registry)),
    sprintf("Existing XY->county lookups found: %d", length(xy_lookup_qs)),
    sprintf("Existing zone->county lookup files found: %d", length(zone_lookup_qs)),
    sprintf("Existing observed county-day qs files found: %d", length(obs_qs)),
    sprintf("Existing modeled county-day qs files found: %d", length(modeled_qs)),
    sprintf("Existing metrics qs files found: %d", length(met_qs))
  )
  writeLines(summary_txt, file.path(QA_DIR, "preflight_summary.txt"))
  
  invisible(list(
    need_lookup = need_lookup,
    lookup_qa = qa,
    wwa_qs = wwa_qs,
    modeled_qs = modeled_qs,
    obs_qs = obs_qs,
    met_qs = met_qs
  ))
}

# ----------------------------
# LOOKUPS
# ----------------------------
build_zone_to_county_lookup <- function(out_qs, overwrite = FALSE) {
  if (file_exists(out_qs) && !overwrite) return(out_qs)
  
  z <- st_make_valid(ZONES_SF)
  cty <- st_make_valid(COUNTIES_SF)
  z <- st_transform(z, 5070)
  cty <- st_transform(cty, 5070)
  
  inter <- st_intersection(z[, c("STATE_ZONE")], cty[, c("GEOID")])
  out <- as.data.table(st_drop_geometry(inter))[, .(STATE_ZONE, GEOID)]
  out[, STATE_ZONE := toupper(trimws(as.character(STATE_ZONE)))]
  out[, GEOID := .norm_geoid5(GEOID)]
  out <- out[!is.na(STATE_ZONE) & !is.na(GEOID)]
  out <- unique(out)
  setorder(out, STATE_ZONE, GEOID)
  
  .write_atomic_qs(out, out_qs)
  .write_atomic_csv(out, sub("\\.qs$", ".csv", out_qs))
  out_qs
}

build_lookup_xy_to_county_one <- function(dataset_tag, ref_tif, out_qs, overwrite = FALSE) {
  if (file_exists(out_qs) && !overwrite) return(out_qs)
  if (!file_exists(ref_tif)) stop("Missing ref_tif for ", dataset_tag, ": ", ref_tif)
  
  r <- terra::rast(ref_tif)
  df <- as.data.table(as.data.frame(r, xy = TRUE, cells = TRUE, na.rm = FALSE))
  setnames(df, old = c("x", "y", "cell"), new = c("x", "y", "cell_id"), skip_absent = TRUE)
  df <- df[, .(cell_id = as.integer(cell_id), x = as.numeric(x), y = as.numeric(y))]
  
  pts <- terra::vect(df[, .(x, y)], geom = c("x", "y"), crs = "EPSG:4326")
  cty_v <- terra::vect(st_transform(COUNTIES_SF, 4326))
  ex <- terra::extract(cty_v, pts)
  if (!("GEOID" %in% names(ex))) stop("terra::extract did not return GEOID for dataset_tag=", dataset_tag)
  
  out <- data.table(
    cell_id = df$cell_id,
    x = df$x,
    y = df$y,
    GEOID = .norm_geoid5(ex$GEOID)
  )
  out <- out[!is.na(GEOID)]
  
  .write_atomic_qs(out, out_qs)
  .write_atomic_csv(out, sub("\\.qs$", ".csv", out_qs))
  out_qs
}

# ----------------------------
# OBSERVED WWA COUNTY-DAYS
# ----------------------------
.discover_wwa_tag_files <- function() {
  if (!dir_exists(WWA_DIR)) stop("WWA_DIR does not exist: ", WWA_DIR)
  
  files <- dir_ls(
    WWA_DIR,
    type = "file",
    glob = "*.qs",
    fail = FALSE
  )
  
  files <- as.character(files)
  files <- files[grepl("^wwa_freq_[A-Z]{2}_[A-Z]_.+\\.qs$", basename(files))]
  files <- sort(files)
  
  .log("[WWA] files found in %s: %d", WWA_DIR, length(files))
  if (length(files)) {
    .log("[WWA] example files: %s", paste(head(basename(files), 5), collapse = ", "))
  }
  
  files
}

.parse_wwa_tag_from_filename <- function(path) {
  bn <- path_file(path)
  m <- str_match(bn, "^wwa_freq_([A-Z]{2})_([A-Z])_.+\\.qs$")
  if (any(is.na(m))) return(NA_character_)
  tolower(paste0(m[1, 2], "_", m[1, 3]))
}

.build_wwa_tag_map <- function() {
  files <- .discover_wwa_tag_files()
  if (!length(files)) stop("No WWA sparse .qs files found in: ", WWA_DIR)
  
  out <- data.table(
    path = files,
    event_tag = vapply(files, .parse_wwa_tag_from_filename, character(1))
  )
  
  out <- out[!is.na(event_tag)]
  out <- unique(out, by = c("event_tag", "path"))
  
  .log("[WWA] parsed tag files: %d", nrow(out))
  if (nrow(out)) print(out[1:min(.N, 10)])
  
  out
}

.build_wwa_individual_county_days <- function(zone_to_county_qs, overwrite = FALSE) {
  tag_map <- .build_wwa_tag_map()
  z2c <- as.data.table(qread(zone_to_county_qs))
  stopifnot(all(c("STATE_ZONE", "GEOID") %in% names(z2c)))
  
  z2c[, STATE_ZONE := toupper(trimws(as.character(STATE_ZONE)))]
  z2c[, GEOID := .norm_geoid5(GEOID)]
  z2c <- z2c[!is.na(STATE_ZONE) & !is.na(GEOID)]
  setkey(z2c, STATE_ZONE)
  
  out_index <- vector("list", nrow(tag_map))
  
  for (i in seq_len(nrow(tag_map))) {
    tag_i <- tolower(tag_map$event_tag[i])
    in_qs <- tag_map$path[i]
    
    out_qs  <- file.path(OBS_DIR, sprintf("wwa_county_days_%s.qs", tag_i))
    out_csv <- file.path(OBS_DIR, sprintf("wwa_county_days_%s.csv", tag_i))
    
    if (file.exists(out_qs) && !overwrite) {
      out_index[[i]] <- data.table(obs_key = tag_i, path = out_qs, type = "tag")
      next
    }
    
    dt <- as.data.table(qread(in_qs))
    
    if (!all(c("STATE_ZONE", "date") %in% names(dt))) {
      warning("Skipping malformed WWA file: ", in_qs)
      out_index[[i]] <- data.table(obs_key = tag_i, path = out_qs, type = "tag")
      next
    }
    
    dt[, STATE_ZONE := toupper(trimws(as.character(STATE_ZONE)))]
    dt[, date := as.Date(date)]
    dt <- dt[!is.na(STATE_ZONE) & !is.na(date)]
    dt <- dt[date >= GLOBAL_START & date <= GLOBAL_END]
    
    if (!nrow(dt)) {
      empty_dt <- data.table(
        event_tag = character(),
        GEOID = character(),
        date = as.Date(character()),
        observed = integer()
      )
      .write_atomic_qs(empty_dt, out_qs)
      .write_atomic_csv(empty_dt, out_csv)
      out_index[[i]] <- data.table(obs_key = tag_i, path = out_qs, type = "tag")
      next
    }
    
    x <- merge(
      dt[, .(STATE_ZONE, date)],
      z2c,
      by = "STATE_ZONE",
      all.x = FALSE,
      allow.cartesian = TRUE
    )
    
    if (!nrow(x)) {
      empty_dt <- data.table(
        event_tag = character(),
        GEOID = character(),
        date = as.Date(character()),
        observed = integer()
      )
      .write_atomic_qs(empty_dt, out_qs)
      .write_atomic_csv(empty_dt, out_csv)
      out_index[[i]] <- data.table(obs_key = tag_i, path = out_qs, type = "tag")
      next
    }
    
    x <- unique(
      x[, .(event_tag = tag_i, GEOID, date)],
      by = c("GEOID", "date")
    )
    x[, observed := 1L]
    setorder(x, date, GEOID)
    
    .write_atomic_qs(x, out_qs)
    .write_atomic_csv(x, out_csv)
    
    out_index[[i]] <- data.table(obs_key = tag_i, path = out_qs, type = "tag")
  }
  
  rbindlist(out_index, use.names = TRUE, fill = TRUE)
}

.build_wwa_group_county_days <- function(tag_index, overwrite = FALSE) {
  out_idx <- list()
  
  for (grp in names(WWA_GROUPS)) {
    obs_key <- tolower(grp)
    tags <- tolower(.event_tag_norm(WWA_GROUPS[[grp]]))
    paths <- tag_index[obs_key %chin% tags, path]
    paths <- paths[file_exists(paths)]
    
    out_qs <- file.path(OBS_DIR, sprintf("wwa_county_days__%s.qs", obs_key))
    out_csv <- file.path(OBS_DIR, sprintf("wwa_county_days__%s.csv", obs_key))
    
    if (file_exists(out_qs) && !overwrite) {
      out_idx[[grp]] <- data.table(obs_key = obs_key, path = out_qs, type = "group")
      next
    }
    
    if (!length(paths)) {
      .write_atomic_qs(
        data.table(event_tag = character(), GEOID = character(), date = as.Date(character()), observed = integer()),
        out_qs
      )
      .write_atomic_csv(
        data.table(event_tag = character(), GEOID = character(), date = as.Date(character()), observed = integer()),
        out_csv
      )
      out_idx[[grp]] <- data.table(obs_key = obs_key, path = out_qs, type = "group")
      next
    }
    
    x <- rbindlist(lapply(paths, qread), use.names = TRUE, fill = TRUE)
    u <- unique(as.data.table(x)[, .(GEOID, date)], by = c("GEOID", "date"))
    u[, `:=`(event_tag = NA_character_, event_group = obs_key, observed = 1L)]
    setorder(u, date, GEOID)
    
    .write_atomic_qs(u, out_qs)
    .write_atomic_csv(u, out_csv)
    out_idx[[grp]] <- data.table(obs_key = obs_key, path = out_qs, type = "group")
  }
  
  rbindlist(out_idx, use.names = TRUE, fill = TRUE)
}

.build_wwa_lookup <- function(tag_index, grp_index, overwrite = FALSE) {
  out_qs  <- file.path(OBS_DIR, "wwa_county_day_lookup.qs")
  out_csv <- file.path(OBS_DIR, "wwa_county_day_lookup.csv")
  
  if (file.exists(out_qs) && !overwrite) return(out_qs)
  
  parts <- list()
  
  # Individual tag rows
  if (nrow(tag_index)) {
    tag_parts <- rbindlist(lapply(seq_len(nrow(tag_index)), function(i) {
      obs_key_i <- tag_index$obs_key[i]
      p <- tag_index$path[i]
      
      if (!file.exists(p)) return(NULL)
      
      x <- as.data.table(qread(p))
      if (!nrow(x)) return(NULL)
      
      data.table(
        obs_key = obs_key_i,
        event_tag = obs_key_i,
        event_group = NA_character_,
        GEOID = .norm_geoid5(x$GEOID),
        date = as.Date(x$date),
        observed = 1L
      )
    }), use.names = TRUE, fill = TRUE)
    
    parts[[length(parts) + 1L]] <- tag_parts
  }
  
  # Group rows
  if (nrow(grp_index)) {
    grp_parts <- rbindlist(lapply(seq_len(nrow(grp_index)), function(i) {
      obs_key_i <- grp_index$obs_key[i]
      p <- grp_index$path[i]
      
      if (!file.exists(p)) return(NULL)
      
      x <- as.data.table(qread(p))
      if (!nrow(x)) return(NULL)
      
      data.table(
        obs_key = obs_key_i,
        event_tag = NA_character_,
        event_group = obs_key_i,
        GEOID = .norm_geoid5(x$GEOID),
        date = as.Date(x$date),
        observed = 1L
      )
    }), use.names = TRUE, fill = TRUE)
    
    parts[[length(parts) + 1L]] <- grp_parts
  }
  
  out <- rbindlist(parts, use.names = TRUE, fill = TRUE)
  if (!nrow(out)) {
    out <- data.table(
      key = character(),
      event_tag = character(),
      event_group = character(),
      GEOID = character(),
      date = as.Date(character()),
      observed = integer()
    )
    .write_atomic_qs(out, out_qs)
    .write_atomic_csv(out, out_csv)
    return(out_qs)
  }
  
  out <- out[!is.na(GEOID) & !is.na(date)]
  out <- unique(out, by = c("obs_key", "GEOID", "date", "event_tag", "event_group"))
  
  # Rename only at the end
  setnames(out, "obs_key", "key")
  setorder(out, key, GEOID, date)
  
  .write_atomic_qs(out, out_qs)
  .write_atomic_csv(out, out_csv)
  out_qs
}

# ----------------------------
# MODELED COUNTY-DAYS
# ----------------------------
.read_clustered_points <- function(qs_path) {
  dt <- as.data.table(qread(qs_path))
  if (!nrow(dt)) return(dt[0])
  
  need <- c("x", "y", "time_utc", "cluster", "cluster_id")
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("Missing columns in ", qs_path, ": ", paste(miss, collapse = ", "))
  
  dt <- dt[is.finite(x) & is.finite(y)]
  dt <- dt[is.finite(cluster) & cluster > 0]
  dt <- dt[!is.na(cluster_id) & nzchar(as.character(cluster_id))]
  if (!nrow(dt)) return(dt[0])
  
  dt[, date := as.Date(time_utc, tz = "UTC")]
  dt <- dt[!is.na(date)]
  dt <- dt[date >= GLOBAL_START & date <= GLOBAL_END]
  dt[, .(x, y, date)]
}

.points_to_modeled_daily_one_file <- function(qs_path, r_ref, lookup_dt) {
  pts <- .read_clustered_points(qs_path)
  if (!nrow(pts)) return(NULL)
  
  pts[, cell_id := terra::cellFromXY(r_ref, cbind(x, y))]
  pts <- pts[!is.na(cell_id)]
  if (!nrow(pts)) return(NULL)
  
  pts <- lookup_dt[pts, on = "cell_id", nomatch = 0L]
  pts <- pts[!is.na(GEOID) & nzchar(GEOID)]
  if (!nrow(pts)) return(NULL)
  
  unique(pts[, .(GEOID, date)][, modeled := 1L], by = c("GEOID", "date"))
}

.build_modeled_one <- function(dataset_row) {
  out_qs <- file.path(MODELED_DIR, dataset_row$scenario, sprintf("modeled_county_daily_%s.qs", dataset_row$dataset_tag))
  if (file_exists(out_qs) && !OVERWRITE_MODELED) return(out_qs)
  
  files <- .discover_cluster_files(dataset_row$dataset_dir)
  if (!length(files)) {
    .write_atomic_qs(data.table(), out_qs)
    return(out_qs)
  }
  
  lookup_qs <- file.path(LOOKUP_XY_TO_COUNTY_DIR, paste0("lookup_xy_to_county_", dataset_row$dataset_tag, ".qs"))
  if (!file_exists(lookup_qs)) stop("Missing lookup for dataset_tag=", dataset_row$dataset_tag, ": ", lookup_qs)
  
  r_ref <- terra::rast(dataset_row$ref_tif)
  lookup_dt <- as.data.table(qread(lookup_qs))
  lookup_dt[, cell_id := as.integer(cell_id)]
  lookup_dt[, GEOID := .norm_geoid5(GEOID)]
  lookup_dt <- lookup_dt[!is.na(GEOID)]
  setkey(lookup_dt, cell_id)
  
  parts <- lapply(files, function(f) {
    tryCatch(.points_to_modeled_daily_one_file(f, r_ref, lookup_dt), error = function(e) {
      .log("[MODEL][ERR] %s | %s", path_file(f), e$message)
      NULL
    })
  })
  
  out <- rbindlist(Filter(Negate(is.null), parts), use.names = TRUE, fill = TRUE)
  if (!nrow(out)) {
    .write_atomic_qs(data.table(), out_qs)
    return(out_qs)
  }
  
  out <- out[, .(modeled = 1L), by = .(GEOID, date)]
  setorder(out, date, GEOID)
  .write_atomic_qs(out, out_qs)
  out_qs
}

# ----------------------------
# EVALUATION ENGINE
# ----------------------------
.discover_modeled_jobs <- function(dataset_registry) {
  out <- copy(dataset_registry)
  out[, mod_path := file.path(MODELED_DIR, scenario, sprintf("modeled_county_daily_%s.qs", dataset_tag))]
  out
}

.discover_wwa_keys_from_lookup <- function(lookup_qs) {
  lk <- as.data.table(qread(lookup_qs))
  sort(unique(as.character(lk$key)))
}

.read_obs_from_lookup <- function(key_value, lookup_qs) {
  lk <- as.data.table(qread(lookup_qs))
  lk <- lk[key == key_value, .(GEOID, date, observed)]
  lk[, GEOID := .norm_geoid5(GEOID)]
  lk[, date := as.Date(date)]
  lk <- lk[!is.na(GEOID) & !is.na(date)]
  lk <- lk[date >= GLOBAL_START & date <= GLOBAL_END]
  lk <- lk[observed > 0L]
  unique(lk[, .(GEOID, date)])
}

.key_to_var_group <- function(k) {
  k <- tolower(as.character(k))
  if (k == "all_heat") return("heatindex")
  if (k == "all_precip") return("precipitation")
  if (grepl("^(eh|ht)_", k)) return("heatindex")
  "precipitation"
}

.read_mod_sparse <- function(path) {
  if (!file_exists(path)) return(data.table(GEOID = character(), date = as.Date(character())))
  dt <- as.data.table(qread(path))
  if (!nrow(dt)) return(data.table(GEOID = character(), date = as.Date(character())))
  stopifnot(all(c("GEOID", "date", "modeled") %in% names(dt)))
  dt[, GEOID := .norm_geoid5(GEOID)]
  dt[, date := as.Date(date)]
  dt[, modeled := as.integer(modeled)]
  dt <- dt[!is.na(GEOID) & !is.na(date)]
  dt <- dt[date >= GLOBAL_START & date <= GLOBAL_END]
  dt <- dt[modeled > 0L]
  unique(dt[, .(GEOID, date)])
}

.metrics_from_counts <- function(TP, FP, FN, TN) {
  N <- TP + FP + FN + TN
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  accuracy <- if (N > 0) (TP + TN) / N else NA_real_
  fpr <- if ((FP + TN) > 0) FP / (FP + TN) else NA_real_
  fnr <- if ((FN + TP) > 0) FN / (FN + TP) else NA_real_
  f1 <- if (is.finite(precision) && is.finite(recall) && (precision + recall) > 0) 2 * precision * recall / (precision + recall) else NA_real_
  beta <- 2
  f2 <- if (is.finite(precision) && is.finite(recall) && (beta^2 * precision + recall) > 0) (1 + beta^2) * precision * recall / (beta^2 * precision + recall) else NA_real_
  auc_balanced <- if (is.finite(recall) && is.finite(specificity)) (recall + specificity) / 2 else NA_real_
  
  list(TP = TP, FP = FP, FN = FN, TN = TN, N = N,
       precision = precision, recall = recall,
       specificity = specificity, accuracy = accuracy,
       f1 = f1, f2 = f2, fpr = fpr, fnr = fnr, auc = auc_balanced)
}

.compute_counts_all_windows <- function(obs_set, mod_set) {
  keys <- CJ(GEOID = COUNTY_UNIVERSE, window_tag = WINDOWS$window_tag, unique = TRUE)
  setkey(keys, GEOID, window_tag)
  
  all_dates <- data.table(date = seq.Date(GLOBAL_START, GLOBAL_END, by = "day"))
  all_dates[, year := year(date)]
  date_windows <- all_dates[year >= YEAR_START & year <= YEAR_END, .(date, window_tag_year = as.character(year))]
  for (w in FIXED_WINDOWS) {
    tag <- w$tag
    date_windows[, (tag) := as.integer(date >= w$start & date <= w$end)]
  }
  
  .attach_windows <- function(set_dt) {
    if (!nrow(set_dt)) return(data.table(GEOID = character(), date = as.Date(character()), window_tag = character()))
    x <- merge(set_dt, date_windows, by = "date", all.x = TRUE)
    out <- list()
    out[[1]] <- x[!is.na(window_tag_year), .(GEOID, date, window_tag = window_tag_year)]
    for (w in FIXED_WINDOWS) {
      tag <- w$tag
      out[[length(out) + 1L]] <- x[get(tag) == 1L, .(GEOID, date, window_tag = tag)]
    }
    rbindlist(out, use.names = TRUE, fill = TRUE)
  }
  
  obs_w <- .attach_windows(obs_set)
  mod_w <- .attach_windows(mod_set)
  
  obs_n <- if (nrow(obs_w)) obs_w[, .(obs_n = .N), by = .(GEOID, window_tag)] else data.table(GEOID = character(), window_tag = character(), obs_n = integer())
  mod_n <- if (nrow(mod_w)) mod_w[, .(mod_n = .N), by = .(GEOID, window_tag)] else data.table(GEOID = character(), window_tag = character(), mod_n = integer())
  
  tp_set <- if (nrow(obs_set) && nrow(mod_set)) {
    setkey(obs_set, GEOID, date)
    setkey(mod_set, GEOID, date)
    inter <- obs_set[mod_set, nomatch = 0L]
    .attach_windows(inter)
  } else data.table(GEOID = character(), date = as.Date(character()), window_tag = character())
  tp_n <- if (nrow(tp_set)) tp_set[, .(tp_n = .N), by = .(GEOID, window_tag)] else data.table(GEOID = character(), window_tag = character(), tp_n = integer())
  
  setkey(obs_n, GEOID, window_tag)
  setkey(mod_n, GEOID, window_tag)
  setkey(tp_n, GEOID, window_tag)
  
  out <- copy(keys)
  out <- out[obs_n]
  out <- out[mod_n]
  out <- out[tp_n]
  
  out[is.na(obs_n), obs_n := 0L]
  out[is.na(mod_n), mod_n := 0L]
  out[is.na(tp_n), tp_n := 0L]
  
  win_days <- WINDOWS[, .(window_tag, total_days, window_type, year, start, end)]
  out <- merge(out, win_days, by = "window_tag", all.x = TRUE)
  
  out[, FP := pmax.int(0L, mod_n - tp_n)]
  out[, FN := pmax.int(0L, obs_n - tp_n)]
  out[, TP := tp_n]
  out[, TN := pmax.int(0L, total_days - TP - FP - FN)]
  out[]
}

.part_path <- function(var_group_value, event_name, scenario, dataset_tag) {
  safe_scn <- gsub("[^A-Za-z0-9_\\-]+", "_", scenario)
  safe_dsn <- gsub("[^A-Za-z0-9_\\-]+", "_", dataset_tag)
  file.path(PART_DIR, sprintf("part_%s_%s__%s__%s.qs", var_group_value, event_name, safe_scn, safe_dsn))
}

.eval_one_job_part <- function(obs_set, var_group_value, event_name, job_row) {
  out_part <- .part_path(var_group_value, event_name, job_row$scenario, job_row$dataset_tag)
  if (file_exists(out_part) && !OVERWRITE_METRICS) return(out_part)
  
  mod_set <- .read_mod_sparse(job_row$mod_path)
  cnt <- .compute_counts_all_windows(obs_set, mod_set)
  if (!nrow(cnt)) return(NA_character_)
  
  met <- cnt[, as.list(.metrics_from_counts(TP, FP, FN, TN)),
             by = .(GEOID, window_tag, window_type, year, start, end)]
  met[, `:=`(
    scenario = job_row$scenario,
    dataset_tag = job_row$dataset_tag,
    var_group = var_group_value,
    event = event_name
  )]
  
  .write_atomic_qs(met, out_part)
  out_part
}

.eval_one_event_from_set <- function(var_group_value, event_name, obs_set, jobs_dt) {
  out_qs <- file.path(MET_DIR, sprintf("metrics_%s_%s_multiwindow.qs", var_group_value, event_name))
  out_csv <- file.path(MET_DIR, sprintf("metrics_%s_%s_multiwindow.csv", var_group_value, event_name))
  
  jobs <- jobs_dt[var_group == var_group_value]
  if (!nrow(jobs)) return(invisible(FALSE))
  
  .log("[EVAL] %s | %s | modeled_jobs=%d | obs_rows=%d", var_group_value, event_name, nrow(jobs), nrow(obs_set))
  
  plan(multisession, workers = min(WORKERS_EVAL, nrow(jobs)))
  handlers(global = TRUE)
  handlers("txtprogressbar")
  
  part_files <- with_progress({
    p <- progressor(steps = nrow(jobs))
    future_lapply(seq_len(nrow(jobs)), function(i) {
      row <- jobs[i]
      p(sprintf("%s | %s", row$scenario, row$dataset_tag))
      .eval_one_job_part(obs_set, var_group_value, event_name, row)
    }, future.seed = TRUE, future.packages = c("data.table", "qs", "fs", "lubridate", "stringr"))
  })
  part_files <- unlist(part_files, use.names = FALSE)
  part_files <- part_files[file_exists(part_files)]
  if (!length(part_files)) return(invisible(FALSE))
  
  parts <- lapply(part_files, function(pth) as.data.table(qread(pth)))
  met_all <- rbindlist(parts, use.names = TRUE, fill = TRUE)
  
  setcolorder(met_all, c("scenario", "dataset_tag", "var_group", "event", "window_type", "window_tag", "year", "start", "end", "GEOID",
                         "TP", "FP", "FN", "TN", "N", "precision", "recall", "specificity", "accuracy", "f1", "f2", "fpr", "fnr", "auc"))
  
  .write_atomic_qs(met_all, out_qs)
  .write_atomic_csv(met_all, out_csv)
  .log("[EVAL] wrote: %s (rows=%d)", path_file(out_csv), nrow(met_all))
  
  invisible(TRUE)
}

# ----------------------------
# MAIN
# ----------------------------
dataset_registry <- .discover_dataset_registry()
preflight <- run_preflight(dataset_registry)

if (RUN_PREFLIGHT_ONLY) {
  .log("RUN_PREFLIGHT_ONLY=TRUE. Stopping after preflight.")
  quit(save = "no")
}

# 1) zone -> county lookup
zone_to_county_qs <- file.path(LOOKUP_ZONE_TO_COUNTY_DIR, "lookup_zone_to_county_bbox.qs")
build_zone_to_county_lookup(zone_to_county_qs, overwrite = OVERWRITE_LOOKUPS)

# 2) dataset XY -> county lookups
need_lookup <- unique(dataset_registry[ref_exists == TRUE, .(dataset_tag, ref_tif)])
if (nrow(need_lookup)) {
  plan(multisession, workers = min(WORKERS_LOOKUP, nrow(need_lookup)))
  handlers(global = TRUE)
  handlers("txtprogressbar")
  
  with_progress({
    p <- progressor(steps = nrow(need_lookup))
    future_lapply(seq_len(nrow(need_lookup)), function(i) {
      row <- need_lookup[i]
      p(row$dataset_tag)
      out_qs <- file.path(LOOKUP_XY_TO_COUNTY_DIR, paste0("lookup_xy_to_county_", row$dataset_tag, ".qs"))
      build_lookup_xy_to_county_one(row$dataset_tag, row$ref_tif, out_qs, overwrite = OVERWRITE_LOOKUPS)
    }, future.seed = TRUE, future.packages = c("data.table", "qs", "fs", "terra", "sf", "tigris"))
  })
}

# 3) build definitive observed WWA county-day lookup
obs_tag_index <- .build_wwa_individual_county_days(zone_to_county_qs, overwrite = OVERWRITE_OBSERVED)
obs_grp_index <- .build_wwa_group_county_days(obs_tag_index, overwrite = OVERWRITE_OBSERVED)
wwa_lookup_qs <- .build_wwa_lookup(obs_tag_index, obs_grp_index, overwrite = OVERWRITE_OBSERVED)

# 4) build modeled county-day caches
plan(multisession, workers = min(WORKERS_MODEL, nrow(dataset_registry)))
handlers(global = TRUE)
handlers("txtprogressbar")
with_progress({
  p <- progressor(steps = nrow(dataset_registry))
  future_lapply(seq_len(nrow(dataset_registry)), function(i) {
    row <- dataset_registry[i]
    p(sprintf("%s | %s", row$scenario, row$dataset_tag))
    .build_modeled_one(row)
  }, future.seed = TRUE, future.packages = c("data.table", "qs", "fs", "terra", "lubridate"))
})
plan(sequential)

# 5) evaluate
# HCQ EDIT ----
HEAT_WWA_PHENOM <- c("EH")

PRECIP_WWA_PHENOM <- c(
  "SV", 
  "FF", "HI", "HU", "TI", "TR"
)

VALID_INDIVIDUAL_WWA_KEYS <- tolower(c(
  # heat
  "EH_A", "EH_W", "HT_A", "HT_W",
  
  # precipitation
  "BZ_A", "BZ_W",
  "CF_A", "CF_W",
  "DF_A", "DF_W",
  "FA_A", "FA_W",
  "FF_A", "FF_W",
  "FL_A", "FL_W",
  "HI_A", "HI_W",
  "HS_A", "HS_W",
  "HU_A", "HU_W",
  "SR_A", "SR_W",
  "SV_A", "SV_W",
  "TR_A", "TR_W",
  "TY_A", "TY_W",
  "WS_A", "WS_W"
))

.get_wwa_phenom_from_key <- function(key) {
  key <- toupper(as.character(key))
  m <- stringr::str_match(key, "^([A-Z]{2})_[A-Z]$")
  out <- m[, 2]
  out[is.na(out) | out == ""] <- NA_character_
  out
}

.key_to_var_group <- function(key) {
  k <- tolower(as.character(key))
  
  # grouped products
  if (k == "all_heat")   return("heatindex")
  if (k == "all_precip") return("precipitation")
  
  # individual WWA keys
  phenom <- .get_wwa_phenom_from_key(k)
  if (is.na(phenom)) return(NA_character_)
  
  if (phenom %in% HEAT_WWA_PHENOM)   return("heatindex")
  if (phenom %in% PRECIP_WWA_PHENOM) return("precipitation")
  
  NA_character_
}

jobs_dt <- .discover_modeled_jobs(dataset_registry)
keys <- .discover_wwa_keys_from_lookup(wwa_lookup_qs)
.log("WWA keys found in lookup (raw): %d", length(keys))

keys_keep <- c("all_heat", "all_precip", VALID_INDIVIDUAL_WWA_KEYS)

if (isTRUE(EVALUATE_INDIVIDUAL_WWA_TAGS)) {
  keys <- sort(unique(keys[keys %in% keys_keep]))
} else {
  keys <- sort(unique(keys[keys %in% c("all_heat", "all_precip")]))
}

.log("WWA keys retained for evaluation: %d", length(keys))
.log("Retained keys: %s", paste(keys, collapse = ", "))

for (k in keys) {
  vg <- .key_to_var_group(k)
  evt <- k
  
  if (is.na(vg)) {
    .log("[EVAL] skipping key with no physical hazard mapping: %s", k)
    next
  }
  
  obs_set <- .read_obs_from_lookup(k, wwa_lookup_qs)
  if (!nrow(obs_set)) {
    .log("[EVAL] skipping empty observed set: %s", k)
    next
  }
  
  .log("[EVAL] key=%s routed_to=%s | obs_rows=%d", k, vg, nrow(obs_set))
  .eval_one_event_from_set(vg, evt, obs_set, jobs_dt)
}

# 6) quick overall summary
met_files <- dir_ls(MET_DIR, glob = "metrics_*.csv", type = "file", fail = FALSE)
if (length(met_files)) {
  all_met <- rbindlist(lapply(met_files, fread), use.names = TRUE, fill = TRUE)
  sum_dt <- all_met[, .(
    n = .N,
    mean_f1 = mean(f1, na.rm = TRUE),
    mean_auc = mean(auc, na.rm = TRUE),
    mean_precision = mean(precision, na.rm = TRUE),
    mean_recall = mean(recall, na.rm = TRUE)
  ), by = .(var_group, event, scenario, dataset_tag, window_tag)]
  .write_atomic_csv(sum_dt, file.path(MET_DIR, "metrics_summary_all.csv"))
}

.log("✅ WWA county evaluation workflow complete")
