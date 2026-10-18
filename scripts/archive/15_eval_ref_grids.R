# suppressPackageStartupMessages({
#   library(data.table); library(fs); library(qs); library(terra); library(lubridate)
# })
# source("V:/users/hquintal/phd1_cluster_southeast_v2/scripts/14_eval_config.R")
# 
# get_stm_resolution <- function(cov_csv, months_keep, base_res=0.25) {
#   dt <- fread(cov_csv)
#   dt[, day := as.Date(day)]
#   if (all(is.na(dt$day))) dt[, day := as.Date(day, format="%m/%d/%Y")]
#   if (all(is.na(dt$day))) stop("Could not parse day in: ", cov_csv)
#   dt <- dt[month(day) %in% months_keep]
#   stm_col <- if ("stm1" %in% names(dt)) "stm1" else if ("stm" %in% names(dt)) "stm" else NA_character_
#   if (is.na(stm_col)) stop("No stm/stm2 in: ", cov_csv)
#   base_res * mean(dt[[stm_col]], na.rm=TRUE)
# }
# 
# build_ref_raster_from_mask <- function(mask_tif, res_deg, out_tif, overwrite=FALSE) {
#   if (file_exists(out_tif) && !overwrite) return(out_tif)
#   r_mask <- rast(mask_tif)
#   r_ref <- rast(extent = ext(r_mask), resolution = res_deg, crs = crs(r_mask))
#   values(r_ref) <- 1
#   dir_create(path_dir(out_tif), recurse = TRUE)
#   writeRaster(r_ref, out_tif, overwrite = TRUE)
#   out_tif
# }
# 
# write_ref_cells_qs <- function(ref_tif, out_qs, overwrite=FALSE) {
#   if (file_exists(out_qs) && !overwrite) return(out_qs)
#   r <- rast(ref_tif)
#   df <- as.data.frame(r, xy=TRUE, cells=TRUE, na.rm=FALSE)
#   dt <- as.data.table(df)[, .(cell_id = cell, x, y)]
#   dir_create(path_dir(out_qs), recurse = TRUE)
#   qsave(dt, out_qs)
#   out_qs
# }
# 
# stopifnot(file_exists(MASK_TIF), file_exists(HI_COV_CSV), file_exists(PR_COV_CSV))
# 
# hi_res <- get_stm_resolution(HI_COV_CSV, HI_STM_MONTHS, BASE_RES)
# pr_res <- get_stm_resolution(PR_COV_CSV, PR_STM_MONTHS, BASE_RES)
# 
# hi_tag <- sprintf("%.2f", hi_res)  # e.g. 0.44
# pr_tag <- sprintf("%.2f", pr_res)  # e.g. 0.17
# 
# .log("[REF] hi_res=%.6f tag=%s | pr_res=%.6f tag=%s", hi_res, hi_tag, pr_res, pr_tag)
# 
# grid_specs <- list(
#   list(tag="heatindex_0.25", res=BASE_RES),
#   list(tag=paste0("heatindex_", hi_tag, "_bilinear"), res=hi_res),
#   list(tag=paste0("heatindex_", hi_tag, "_max"),      res=hi_res),
#   list(tag="precipitation_0.25", res=BASE_RES),
#   list(tag=paste0("precipitation_", pr_tag, "_bilinear"), res=pr_res),
#   list(tag=paste0("precipitation_", pr_tag, "_max"),      res=pr_res)
# )
# 
# ref_map <- rbindlist(lapply(grid_specs, function(s) {
#   tag <- s$tag; res <- s$res
#   ref_tif <- file.path(REF_DIR, tag, paste0("ref_grid_", tag, ".tif"))
#   ref_qs  <- file.path(REF_DIR, tag, paste0("ref_grid_cells_", tag, ".qs"))
#   build_ref_raster_from_mask(MASK_TIF, res, ref_tif, overwrite=OVERWRITE_REF)
#   write_ref_cells_qs(ref_tif, ref_qs, overwrite=OVERWRITE_REF)
#   data.table(dataset_tag=tag, res_deg=res, ref_tif=ref_tif, ref_cells_qs=ref_qs)
# }))
# 
# qsave(ref_map, file.path(REF_DIR, "ref_map.qs"))
# .log("Wrote ref_map.qs (rows=%d)", nrow(ref_map))

# 20260224 ----
suppressPackageStartupMessages({
  library(data.table); library(fs); library(qs); library(terra); library(lubridate)
})

# County pipeline config (swap to 14_eval_config.R if staying WFO)
source("V:/users/hquintal/phd1_cluster_southeast_v2/scripts/14_eval_config.R")

# ----------------------------
# STM diagnostics (stm1 + stm2 + stm if present)
# ----------------------------
.get_stm_resolution <- function(cov_csv, months_keep, base_res=0.25, stm_col=c("stm1","stm2","stm")) {
  stm_col <- match.arg(stm_col)
  dt <- fread(cov_csv)
  dt[, day := as.Date(day)]
  if (all(is.na(dt$day))) dt[, day := as.Date(day, format="%m/%d/%Y")]
  if (all(is.na(dt$day))) stop("Could not parse day in: ", cov_csv)
  
  dt <- dt[month(day) %in% months_keep]
  if (!nrow(dt)) return(NA_real_)
  if (!(stm_col %in% names(dt))) return(NA_real_)
  
  base_res * mean(dt[[stm_col]], na.rm=TRUE)
}

# ----------------------------
# Parse resolution from dataset_tag
# - supports: heatindex_0.44_bilinear, precipitation_0.12_max, heatindex_0.25
# ----------------------------
.parse_res_from_tag <- function(tag) {
  # first numeric token after underscore
  m <- regexpr("_[0-9]+\\.[0-9]+", tag)
  if (m[1] == -1) return(NA_real_)
  as.numeric(sub("^_", "", regmatches(tag, m)))
}

# ----------------------------
# Ref grid builders
# ----------------------------
build_ref_raster_from_mask <- function(mask_tif, res_deg, out_tif, overwrite=FALSE) {
  if (file_exists(out_tif) && !overwrite) return(out_tif)
  r_mask <- rast(mask_tif)
  r_ref  <- rast(extent = ext(r_mask), resolution = res_deg, crs = crs(r_mask))
  values(r_ref) <- 1
  dir_create(path_dir(out_tif), recurse = TRUE)
  writeRaster(r_ref, out_tif, overwrite = TRUE)
  out_tif
}

write_ref_cells_qs <- function(ref_tif, out_qs, overwrite=FALSE) {
  if (file_exists(out_qs) && !overwrite) return(out_qs)
  r  <- rast(ref_tif)
  df <- as.data.frame(r, xy=TRUE, cells=TRUE, na.rm=FALSE)
  dt <- as.data.table(df)[, .(cell_id = as.integer(cell), x, y)]
  dir_create(path_dir(out_qs), recurse = TRUE)
  qsave(dt, out_qs)
  out_qs
}

stopifnot(file_exists(MASK_TIF), file_exists(HI_COV_CSV), file_exists(PR_COV_CSV))

# ----------------------------
# STM diagnostics log (so you "see" stm1 vs stm2)
# ----------------------------
hi_stm1 <- .get_stm_resolution(HI_COV_CSV, HI_STM_MONTHS, BASE_RES, "stm1")
hi_stm2 <- .get_stm_resolution(HI_COV_CSV, HI_STM_MONTHS, BASE_RES, "stm2")
pr_stm1 <- .get_stm_resolution(PR_COV_CSV, PR_STM_MONTHS, BASE_RES, "stm1")
pr_stm2 <- .get_stm_resolution(PR_COV_CSV, PR_STM_MONTHS, BASE_RES, "stm2")

.log("[STM DIAG] heatindex: stm1=%s | stm2=%s",
     ifelse(is.na(hi_stm1), "NA", sprintf("%.4f", hi_stm1)),
     ifelse(is.na(hi_stm2), "NA", sprintf("%.4f", hi_stm2)))
.log("[STM DIAG] precip:   stm1=%s | stm2=%s",
     ifelse(is.na(pr_stm1), "NA", sprintf("%.4f", pr_stm1)),
     ifelse(is.na(pr_stm2), "NA", sprintf("%.4f", pr_stm2)))

# ----------------------------
# Determine required unique resolutions
# - authoritative: DATASET_TAGS
# - include BASE_RES always
# ----------------------------
dt_tags <- data.table(dataset_tag = DATASET_TAGS)
dt_tags[, res_deg := vapply(dataset_tag, .parse_res_from_tag, numeric(1))]
dt_tags[is.na(res_deg), res_deg := BASE_RES]  # fallback safety

# round to match your naming convention (0.12, 0.17, 0.34, 0.44, 0.25, ...)
dt_tags[, res_tag := sprintf("%.2f", res_deg)]
dt_tags[, ref_grid_tag := paste0("res_", res_tag)]

# unique resolutions only (across variables)
needed_res <- unique(dt_tags[, .(res_deg = unique(res_deg))])
needed_res[, res_tag := sprintf("%.2f", res_deg)]
needed_res[, ref_grid_tag := paste0("res_", res_tag)]
setorder(needed_res, res_deg)

.log("[REF] unique resolutions needed (n=%d): %s",
     nrow(needed_res),
     paste(needed_res$res_tag, collapse=", "))

# ----------------------------
# Build ref grid per unique resolution
# ----------------------------
ref_map <- rbindlist(lapply(seq_len(nrow(needed_res)), function(i) {
  res_deg <- needed_res$res_deg[i]
  ref_grid_tag <- needed_res$ref_grid_tag[i]
  
  ref_tif <- file.path(REF_DIR, ref_grid_tag, paste0("ref_grid_", ref_grid_tag, ".tif"))
  ref_qs  <- file.path(REF_DIR, ref_grid_tag, paste0("ref_grid_cells_", ref_grid_tag, ".qs"))
  
  build_ref_raster_from_mask(MASK_TIF, res_deg, ref_tif, overwrite=OVERWRITE_REF)
  write_ref_cells_qs(ref_tif, ref_qs, overwrite=OVERWRITE_REF)
  
  data.table(
    ref_grid_tag = ref_grid_tag,
    res_deg = res_deg,
    ref_tif = ref_tif,
    ref_cells_qs = ref_qs
  )
}), use.names=TRUE, fill=TRUE)

setkey(ref_map, ref_grid_tag)

# ----------------------------
# Build dataset_tag -> ref_grid mapping (shared lookups by resolution)
# ----------------------------
dataset_map <- merge(
  dt_tags[, .(dataset_tag, res_deg, res_tag, ref_grid_tag,
              var_group = fifelse(grepl("^heatindex", dataset_tag, ignore.case=TRUE), "heatindex", "precipitation"))],
  ref_map,
  by = "ref_grid_tag",
  all.x = TRUE
)

if (any(is.na(dataset_map$ref_tif))) {
  stop("dataset_map has missing ref_tif for some tags. Check DATASET_TAGS parsing.")
}
setkey(dataset_map, dataset_tag)

# write outputs
qsave(ref_map,     file.path(REF_DIR, "ref_map_by_resolution.qs"))
qsave(dataset_map, file.path(REF_DIR, "dataset_to_refgrid.qs"))

.log("Wrote ref_map_by_resolution.qs (rows=%d)", nrow(ref_map))
.log("Wrote dataset_to_refgrid.qs (rows=%d)", nrow(dataset_map))