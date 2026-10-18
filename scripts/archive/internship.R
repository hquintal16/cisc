# =============================================================================
# DFO ZIP -> plot ALL raster locations in solid blue
# (Any non-NA value is blue)
# Fixed CONUS extent
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(tigris)
  library(fs)
})

options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)

# ----------------------------
# USER SETTINGS
# ----------------------------
DL_DIR <- "C:/Users/hquintal.AD/Downloads"
ZIP_REGEX <- "(?i)^DFO_[0-9]+_From_[0-9]{8}_to_[0-9]{8}\\.zip$"

# ----------------------------
# Helpers
# ----------------------------
.list_zip_tifs <- function(zip_path) {
  info <- utils::unzip(zip_path, list = TRUE)
  info$Name[grepl("\\.tif(f)?$", info$Name, ignore.case = TRUE)]
}

.read_tif_from_zip <- function(zip_path, internal_name) {
  zp <- gsub("\\\\", "/", zip_path)
  iname <- gsub("\\\\", "/", internal_name)
  terra::rast(paste0("/vsizip/", zp, "/", iname))
}

.ext_overlap <- function(e1, e2) {
  !(e1$xmax <= e2$xmin || e1$xmin >= e2$xmax ||
      e1$ymax <= e2$ymin || e1$ymin >= e2$ymax)
}

# ----------------------------
# Fixed CONUS extent (EPSG:4326)
# ----------------------------
states_sf <- tigris::states(cb = TRUE, year = 2024, class = "sf")
drop <- c("AK","HI","PR","GU","VI","MP","AS")
states_sf <- states_sf[!(states_sf$STUSPS %in% drop), ]

states_sf_ll <- sf::st_transform(states_sf, 4326)
states_v_ll  <- terra::vect(states_sf_ll)
conus_ext_ll <- terra::ext(states_v_ll)

# ----------------------------
# Find ZIPs
# ----------------------------
all_zips  <- fs::dir_ls(DL_DIR, regexp = "(?i)\\.zip$", type = "file")
zip_files <- all_zips[grepl(ZIP_REGEX, basename(all_zips), perl = TRUE)]
if (length(zip_files) == 0) stop("No matching DFO ZIPs found.")

zip_tifs <- lapply(zip_files, function(z) {
  tifs <- tryCatch(.list_zip_tifs(z), error = function(e) character(0))
  list(zip = z, tifs = tifs)
})
zip_tifs <- Filter(function(x) length(x$tifs) > 0, zip_tifs)

# ----------------------------
# Plot loop
# ----------------------------
for (zx in zip_tifs) {
  zip_path <- zx$zip
  
  for (internal_tif in zx$tifs) {
    
    r <- tryCatch(.read_tif_from_zip(zip_path, internal_tif),
                  error = function(e) NULL)
    if (is.null(r)) next
    
    # Use first layer
    r1 <- r[[1]]
    
    # Project to lon/lat
    r_ll <- tryCatch(
      terra::project(r1, "EPSG:4326", method = "near"),
      error = function(e) NULL
    )
    if (is.null(r_ll)) next
    
    if (!.ext_overlap(terra::ext(r_ll), conus_ext_ll)) next
    
    r_crop <- terra::crop(r_ll, conus_ext_ll)
    
    # Convert all non-NA cells to 1 (presence everywhere data exists)
    r_blue <- !is.na(r_crop)
    r_blue[r_blue == 0] <- NA
    
    ttl <- basename(internal_tif)
    
    # 1) plot states first (fixed extent)
    terra::plot(
      states_v_ll,
      main = ttl,
      col = "white",
      border = "black",
      xlim = c(conus_ext_ll$xmin, conus_ext_ll$xmax),
      ylim = c(conus_ext_ll$ymin, conus_ext_ll$ymax)
    )
    
    # 2) plot raster as solid blue
    terra::plot(
      r_blue,
      add = TRUE,
      col = "blue",
      legend = FALSE
    )
    
    # 3) redraw state outlines
    terra::lines(states_v_ll, col = "black", lwd = 0.8)
    
    message("[PLOT] ", ttl)
  }
}