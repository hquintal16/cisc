# suppressPackageStartupMessages({
#   library(data.table)
#   library(fs)
#   library(qs)
#   library(sf)
#   library(terra)
# })
# 
# source("V:/users/hquintal/phd1_cluster_southeast_v2/scripts/14_eval_config.R")
# 
# .read_cwa_polygons_terra <- function() {
#   stopifnot(file_exists(CWA_SHP))
#   wfo_sf <- sf::st_read(CWA_SHP, quiet = TRUE)
#   if (!("CWA" %in% names(wfo_sf))) {
#     stop("CWA field not found. Have: ", paste(names(wfo_sf), collapse = ", "))
#   }
#   wfo_sf$CWA <- toupper(trimws(as.character(wfo_sf$CWA)))
#   
#   # Convert to terra vector (terra handles many geometry issues more gracefully)
#   wfo_v <- terra::vect(wfo_sf)
#   wfo_v
# }
# 
# build_lookup_xy_to_cwa <- function(ref_cells_qs, out_qs, wfo_v, overwrite = FALSE) {
#   if (file_exists(out_qs) && !overwrite) return(out_qs)
#   
#   cells <- as.data.table(qread(ref_cells_qs)) # cell_id, x, y
#   stopifnot(all(c("cell_id","x","y") %in% names(cells)))
#   
#   # Points in lon/lat WGS84 (your ref grids are lon/lat)
#   pts_v <- terra::vect(cells[, .(x, y)], geom = c("x","y"), crs = "EPSG:4326")
#   
#   # If WFO is NAD83, project points to match WFO CRS
#   # (terra::project is safe; use WFO CRS as target)
#   pts_v <- terra::project(pts_v, terra::crs(wfo_v))
#   
#   ex <- terra::extract(wfo_v, pts_v)
#   # ex has ID column + polygon attributes
#   if (!("CWA" %in% names(ex))) {
#     stop("terra::extract did not return CWA. Returned cols: ", paste(names(ex), collapse=", "))
#   }
#   
#   out <- data.table(
#     cell_id = cells$cell_id,
#     x = cells$x,
#     y = cells$y,
#     CWA = toupper(trimws(as.character(ex$CWA)))
#   )
#   
#   # Drop outside polygons (ocean, etc.)
#   out <- out[!is.na(CWA) & nzchar(CWA)]
#   
#   dir_create(path_dir(out_qs), recurse = TRUE)
#   tmp <- paste0(out_qs, ".tmp")
#   if (file_exists(tmp)) file_delete(tmp)
#   qsave(out, tmp); file_move(tmp, out_qs)
#   
#   out_qs
# }
# 
# # ---- run ----
# ref_map <- as.data.table(qread(file.path(REF_DIR, "ref_map.qs")))
# # ref_map <- ref_map[c(2,5),]
# wfo_v <- .read_cwa_polygons_terra()
# 
# for (i in seq_len(nrow(ref_map))) {
#   tag <- ref_map$dataset_tag[i]
#   out_qs <- file.path(LOOKUP_DIR, paste0("lookup_xy_to_cwa_", tag, ".qs"))
#   build_lookup_xy_to_cwa(ref_map$ref_cells_qs[i], out_qs, wfo_v, overwrite = OVERWRITE_LOOKUP)
#   .log("[LOOKUP] %s", out_qs)
# }
# 
# # # check ----
# # library(qs); library(data.table)
# # lk <- as.data.table(qread("V:/users/hquintal/phd1_cluster_southeast_v2/data/out/09_performance/lookup_xy_to_cwa/lookup_xy_to_cwa_heatindex_0.25.qs"))
# # lk[, .(n=.N, n_unique_cwa=uniqueN(CWA), frac_na=mean(is.na(CWA) | !nzchar(CWA)))]

# 20260224 ----
suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(qs)
  library(sf)
  library(terra)
  library(tigris)
})

source("V:/users/hquintal/phd1_cluster_southeast_v2/scripts/14_eval_config.R")

options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)

.read_counties_terra <- function() {
  cty_sf <- tigris::counties(cb = TRUE, year = 2022, class = "sf")
  stopifnot(all(c("GEOID","STATEFP") %in% names(cty_sf)))
  cty_sf <- cty_sf[!(cty_sf$STATEFP %in% DROP_STATEFP), ]
  cty_sf <- st_make_valid(cty_sf)
  cty_sf <- st_transform(cty_sf, 4326)
  terra::vect(cty_sf)
}

build_lookup_xy_to_county <- function(ref_cells_qs, out_qs, cty_v, overwrite = FALSE) {
  if (file_exists(out_qs) && !overwrite) return(out_qs)
  
  cells <- as.data.table(qread(ref_cells_qs)) # cell_id, x, y
  stopifnot(all(c("cell_id","x","y") %in% names(cells)))
  
  pts_v <- terra::vect(cells[, .(x, y)], geom = c("x","y"), crs = "EPSG:4326")
  pts_v <- terra::project(pts_v, terra::crs(cty_v))
  
  ex <- terra::extract(cty_v, pts_v)
  if (!("GEOID" %in% names(ex))) {
    stop("terra::extract did not return GEOID. Returned cols: ", paste(names(ex), collapse=", "))
  }
  
  out <- data.table(
    cell_id = as.integer(cells$cell_id),
    x = cells$x,
    y = cells$y,
    GEOID = as.character(ex$GEOID)
  )
  
  out <- out[!is.na(GEOID) & nzchar(GEOID)]
  
  dir_create(path_dir(out_qs), recurse = TRUE)
  tmp <- paste0(out_qs, ".tmp")
  if (file_exists(tmp)) file_delete(tmp)
  qsave(out, tmp); file_move(tmp, out_qs)
  
  out_qs
}

# ---- run (one per resolution) ----
ref_map <- as.data.table(qread(file.path(REF_DIR, "ref_map_by_resolution.qs")))
stopifnot(all(c("ref_grid_tag","ref_cells_qs") %in% names(ref_map)))

cty_v <- .read_counties_terra()

for (i in seq_len(nrow(ref_map))) {
  ref_grid_tag <- ref_map$ref_grid_tag[i]
  out_qs <- file.path(LOOKUP_DIR_COUNTY, paste0("lookup_xy_to_county_", ref_grid_tag, ".qs"))
  
  build_lookup_xy_to_county(
    ref_cells_qs = ref_map$ref_cells_qs[i],
    out_qs = out_qs,
    cty_v = cty_v,
    overwrite = OVERWRITE_LOOKUP
  )
  
  .log("[LOOKUP COUNTY] %s", out_qs)
}

# test <- qs::qread('V:/users/hquintal/phd1_cluster_southeast_v2/data/out/09_performance/lookup_xy_to_county/lookup_xy_to_county_res_0.12.qs')
# test
