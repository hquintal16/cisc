# Heat index
test <- read.csv("/proj/sebastianlab/users/hquintal/cisc/data/out/03_cov_model/heatindex/heatindex_daily_optimal_cov_model_southeast.csv")
head(test)

sd(test$stm2, na.rm = TRUE)
density(test$stm2, na.rm = TRUE)

# Precipitation
test <- read.csv("/proj/sebastianlab/users/hquintal/cisc/data/out/03_cov_model/precipitation/precipitation_daily_optimal_cov_model_southeast.csv")
head(test)

sd(test$stm2, na.rm = TRUE)
density(test$stm2, na.rm = TRUE)

#####
test <- read.csv(
  "/proj/sebastianlab/users/hquintal/cisc/data/out/03_cov_model/heatindex/heatindex_daily_optimal_cov_model_southeast.csv"
)

test$day <- as.Date(test$day)

# Keep May through October
test_warm <- test[
  as.integer(format(test$day, "%m")) %in% 5:10,
]

# Standard deviation of stm2
sd(test_warm$stm2, na.rm = TRUE)

# Density of stm2
stm2_values <- test_warm$stm2[
  is.finite(test_warm$stm2)
]

test2 <- density(stm2_values)

plot(
  test2,
  main = "Heat-index stm2 density, May–October",
  xlab = "stm2"
)

c(
  n = length(stm2_values),
  mean = mean(stm2_values),
  sd = sd(stm2_values),
  min = min(stm2_values),
  max = max(stm2_values)
)

# correlate mesonet counts per day and stm daily values ----
library(arrow)
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/10_mesonet/annual_combined/mesonet_warning_combined_daily_2002.parquet')
head(test)

table(test$PHENOM)

# check clusters ----
library(arrow)

# heat index ----
# 0.25
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/heatindex_0.25/cluster_rp5_2023.parquet')
length(unique(test$cluster_id))
# 0.11 bil
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/heatindex_0.30_bilinear/cluster_rp5_2023.parquet')
length(unique(test$cluster_id))
# 0.11 max
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/heatindex_0.30_max/cluster_rp5_2023.parquet')
length(unique(test$cluster_id))
# 0.14 bil
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/heatindex_0.36_bilinear/cluster_rp5_2023.parquet')
length(unique(test$cluster_id))
# 0.14 max
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/heatindex_0.36_max/cluster_rp5_2023.parquet')
length(unique(test$cluster_id))
# precipitation ----
# 0.25
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/precipitation_0.25/cluster_rp5_2023.parquet')
length(unique(test$cluster_id))
# 0.11 bil
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/precipitation_0.11_bilinear/cluster_rp5_2023.parquet')
length(unique(test$cluster_id))
# 0.11 max
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/precipitation_0.11_max/cluster_rp5_2023.parquet')
length(unique(test$cluster_id))
# 0.14 bil
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/precipitation_0.14_bilinear/cluster_rp5_2023.parquet')
length(unique(test$cluster_id))
# 0.14 max
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/precipitation_0.14_max/cluster_rp5_2023.parquet')
length(unique(test$cluster_id))

# check tifs ----
library(rast)
test <- terra::rast("/proj/sebastianlab/users/hquintal/cisc/data/out/04_rp/heatindex_idf_hi1_return_periods_bbox.nc")
plot(test[[1]])

# check wwa files ----
library(data.table)
library(qs)

cache_file <- paste0(
  "/proj/sebastianlab/users/hquintal/cisc/",
  "figures/fig_02_wwa/",
  "wwa_cache_ugcday_bbox_2002_2023.qs"
)

x <- as.data.table(qread(cache_file))

x[
  PHENOM %in% c("HT", "DF", "SR", "TY"),
  .N,
  by = .(PHENOM, SIG)
][order(PHENOM, SIG)]

# check mesonet info ----
library(sf)
test <- sf::st_read("/proj/sebastianlab/users/hquintal/cisc/data/in/mesonet/2023_all/wwa_202301010000_202312312359.shp")
head(test)
test <- test %>%
  filter(test$PHENOM == 'EH' & test$SIG == 'W')

library(dplyr)
library(lubridate)

month_range <- test %>%
  st_drop_geometry() %>%
  mutate(
    issued_datetime  = ymd_hm(as.character(ISSUED), tz = "UTC"),
    expired_datetime = ymd_hm(as.character(EXPIRED), tz = "UTC")
  ) %>%
  summarise(
    earliest_issued_month = min(month(issued_datetime), na.rm = TRUE),
    latest_expired_month  = max(month(expired_datetime), na.rm = TRUE)
  )

month_range

# check cluster info ----
# heat index ----
# 0.25
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/heatindex_0.25/cluster_rp5_2023.parquet')
head(test)
table(test$cluster_id)

# heat index ----
# 0.25
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/heatindex_0.25/cluster_rp5_2023.parquet')
head(test)
test <- test %>%
  filter(cluster_id == "heatindex-heatindex_0.25-2023-06-15-001")
range(test$time_utc)

# 0.30 bilinear
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/heatindex_0.30_bilinear/cluster_rp5_2023.parquet')
test <- test %>%
  filter(cluster_id == "heatindex-heatindex_0.30_bilinear-2023-06-15-002")
range(test$time_utc)

# 0.30 max
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/heatindex_0.30_max/cluster_rp5_2023.parquet')
test <- test %>%
  filter(cluster_id == "heatindex-heatindex_0.30_max-2023-06-15-001")
range(test$time_utc)

# 0.36 bilinear
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/heatindex_0.36_bilinear/cluster_rp5_2023.parquet')
test <- test %>%
  filter(cluster_id == "heatindex-heatindex_0.36_bilinear-2023-06-15-001")
range(test$time_utc)

# 0.36 max
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/heatindex_0.36_max/cluster_rp5_2023.parquet')
test <- test %>%
  filter(cluster_id == "heatindex-heatindex_0.36_max-2023-06-15-001")
range(test$time_utc)

# precipitation ----
# 0.25
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/precipitation_0.25/cluster_rp5_2023.parquet')
test <- test %>%
  filter(cluster_id == "precipitation-precipitation_0.25-2023-08-30-001")
range(test$time_utc)

# 0.11 bilinear
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/precipitation_0.11_bilinear/cluster_rp5_2023.parquet')
test <- test %>%
  filter(cluster_id == "precipitation-precipitation_0.11_bilinear-2023-08-30-002")
range(test$time_utc)

# 0.11 max
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/precipitation_0.11_max/cluster_rp5_2023.parquet')
test <- test %>%
  filter(cluster_id == "precipitation-precipitation_0.11_max-2023-08-30-002")
range(test$time_utc)

# 0.14 bilinear
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/precipitation_0.14_bilinear/cluster_rp5_2023.parquet')
test <- test %>%
  filter(cluster_id == "precipitation-precipitation_0.14_bilinear-2023-08-30-002")
range(test$time_utc)

# 0.14 max
test <- arrow::read_parquet('/proj/sebastianlab/users/hquintal/cisc/data/out/06_clusters/filtered_rp5/precipitation_0.14_max/cluster_rp5_2023.parquet')
test <- test %>%
  filter(cluster_id == "precipitation-precipitation_0.14_max-2023-08-30-001")
range(test$time_utc)

# performance ----
test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/20260802_model_performance/heat/EH_W/confusion_matrices.csv')
test <- test %>%
  select(dataset_tag,return_period,f1)
test

test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/20260802_model_performance/precipitation/FF_W/confusion_matrices.csv')
test <- test %>%
  select(dataset_tag,return_period,f1)
test

# improvement ----
test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/20260802_model_performance/heat/EH_W/f1_comparison_vs_native.csv')
test

test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/20260802_model_performance/precipitation/FF_W/f1_comparison_vs_native.csv')
test <- test %>%
  select(dataset_tag,return_period,t_p_value) %>%
  filter(t_p_value < 0.05)
test

# cov models ----
test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/03_cov_model/heatindex/heatindex_daily_optimal_cov_model_southeast.csv')
head(test)
class(test$day)
test <- test %>%
  filter()

# Confusion matrix ----
test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/20260802_model_performance/warning_combination_summary/all_warning_combination_confusion_matrices.csv')
test <- test %>%
  filter(var_group == "precipitation" & return_period == 50) %>%
  select(dataset_tag,return_period,TP,FP,FN,TN)
test

test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/20260802_model_performance/warning_combination_summary/all_warning_combination_confusion_matrices.csv')
test <- test %>%
  filter(var_group == "precipitation" & return_period == 50) %>%
  select(dataset_tag,return_period,precision)
test

# characteristics ----
test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/07_characteristics/filtered_rp5/heatindex_0.25/cluster_characteristics.csv')
test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/07_characteristics/filtered_rp5/heatindex_0.30_bilinear/cluster_characteristics.csv')
test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/07_characteristics/filtered_rp5/heatindex_0.30_max/cluster_characteristics.csv')
test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/07_characteristics/filtered_rp5/heatindex_0.36_bilinear/cluster_characteristics.csv')
test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/07_characteristics/filtered_rp5/heatindex_0.36_max/cluster_characteristics.csv')

test %>%
  mutate(
    start_month = as.integer(format(as.Date(start_time), "%m")),
    end_month   = as.integer(format(as.Date(end_time), "%m"))
  ) %>%
  summarise(
    earliest_month = min(start_month),
    latest_month   = max(end_month)
  )

# performance ----
test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/20260802_model_performance/heat/EH_W/metrics_county.csv')
test

test <- read.csv('/proj/sebastianlab/users/hquintal/cisc/data/out/20260802_model_performance/heat/EH_W/confusion_matrices.csv')
test <- test %>%
  select()

# =============================================================================
# Minimum and maximum county recall west of 85°W
# 5-year return period, native 0.25° heat-index dataset only
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(sf)
  library(tigris)
  library(stringr)
})

options(tigris_use_cache = TRUE)

METRICS_FILE <- paste0(
  "/proj/sebastianlab/users/hquintal/cisc/data/out/",
  "20260802_model_performance/heat/EH_W/metrics_county.csv"
)

# -----------------------------------------------------------------------------
# Read and filter performance data
# -----------------------------------------------------------------------------

metrics_native_rp5 <- read_csv(
  METRICS_FILE,
  col_types = cols(
    GEOID = col_character(),
    .default = col_guess()
  )
) %>%
  mutate(
    GEOID = str_pad(GEOID, width = 5, side = "left", pad = "0")
  ) %>%
  filter(
    return_period == 5,
    dataset_tag == "heatindex_0.25"
  ) %>%
  select(
    GEOID,
    recall
  )

# Confirm there is only one recall value per county
duplicate_geoids <- metrics_native_rp5 %>%
  count(GEOID) %>%
  filter(n > 1)

if (nrow(duplicate_geoids) > 0L) {
  stop(
    "Some GEOIDs still have more than one row after filtering:\n",
    paste(duplicate_geoids$GEOID, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Identify counties west of 85°W
# -----------------------------------------------------------------------------

counties_west_85 <- tigris::counties(
  cb = TRUE,
  resolution = "20m",
  year = 2023,
  class = "sf"
) %>%
  st_transform(4326) %>%
  select(
    GEOID,
    county = NAME,
    state = STUSPS,
    geometry
  ) %>%
  st_point_on_surface()

coords <- st_coordinates(counties_west_85)

counties_west_85 <- counties_west_85 %>%
  st_drop_geometry() %>%
  mutate(
    longitude = coords[, "X"],
    latitude  = coords[, "Y"]
  ) %>%
  filter(
    longitude < -85#,
    # latitude > 35
  )

# -----------------------------------------------------------------------------
# Join and calculate overall minimum and maximum recall
# -----------------------------------------------------------------------------

county_recall_west_85 <- metrics_native_rp5 %>%
  inner_join(
    counties_west_85,
    by = "GEOID"
  ) %>%
  filter(
    is.finite(recall)
  )

recall_min_max_west_85 <- county_recall_west_85 %>%
  summarise(
    n_counties = n_distinct(GEOID),
    minimum_recall = min(recall),
    maximum_recall = max(recall),
    mean_recall = mean(recall)
  )

print(recall_min_max_west_85)


# -----------------------------------------------------------------------------
# Florida counties only
# -----------------------------------------------------------------------------

florida_counties <- tigris::counties(
  cb = TRUE,
  resolution = "20m",
  year = 2023,
  class = "sf"
) %>%
  st_drop_geometry() %>%
  filter(STUSPS == "FL") %>%
  select(
    GEOID,
    county = NAME
  )

florida_recall <- metrics_native_rp5 %>%
  inner_join(
    florida_counties,
    by = "GEOID"
  )

# -----------------------------------------------------------------------------
# Mean recall across Florida counties
# -----------------------------------------------------------------------------

florida_summary <- florida_recall %>%
  summarise(
    n_counties = n(),
    mean_recall = mean(recall, na.rm = TRUE),
    minimum_recall = min(recall, na.rm = TRUE),
    maximum_recall = max(recall, na.rm = TRUE)
  )

print(florida_summary)

# =============================================================================
# Minimum and maximum county recall west of 85°W
# 5-year return period, native 0.25° heat-index dataset only
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(sf)
  library(tigris)
  library(stringr)
})

options(tigris_use_cache = TRUE)

METRICS_FILE <- paste0(
  "/proj/sebastianlab/users/hquintal/cisc/data/out/",
  "20260802_model_performance/heat/EH_W/metrics_county.csv"
)

# -----------------------------------------------------------------------------
# Read and filter performance data
# -----------------------------------------------------------------------------

metrics_native_rp5 <- read_csv(
  METRICS_FILE,
  col_types = cols(
    GEOID = col_character(),
    .default = col_guess()
  )
) %>%
  mutate(
    GEOID = str_pad(GEOID, width = 5, side = "left", pad = "0")
  ) %>%
  filter(
    return_period == 5,
    dataset_tag == "heatindex_0.25"
  ) %>%
  select(
    GEOID,
    precision
  )

# Confirm there is only one recall value per county
duplicate_geoids <- metrics_native_rp5 %>%
  count(GEOID) %>%
  filter(n > 1)

if (nrow(duplicate_geoids) > 0L) {
  stop(
    "Some GEOIDs still have more than one row after filtering:\n",
    paste(duplicate_geoids$GEOID, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Identify counties north of 35°N
# -----------------------------------------------------------------------------

counties_west_85 <- tigris::counties(
  cb = TRUE,
  resolution = "20m",
  year = 2023,
  class = "sf"
) %>%
  st_transform(4326) %>%
  select(
    GEOID,
    county = NAME,
    state = STUSPS,
    geometry
  ) %>%
  st_point_on_surface()

coords <- st_coordinates(counties_west_85)

counties_west_85 <- counties_west_85 %>%
  st_drop_geometry() %>%
  mutate(
    longitude = coords[, "X"],
    latitude  = coords[, "Y"]
  ) %>%
  filter(
    # longitude < -85#,
    latitude < 35
  )

# -----------------------------------------------------------------------------
# Join and calculate overall minimum and maximum recall
# -----------------------------------------------------------------------------

county_recall_west_85 <- metrics_native_rp5 %>%
  inner_join(
    counties_west_85,
    by = "GEOID"
  ) %>%
  filter(
    is.finite(precision)
  )

recall_min_max_west_85 <- county_recall_west_85 %>%
  summarise(
    n_counties = n_distinct(GEOID),
    mean_precision = mean(precision)
  )

print(recall_min_max_west_85)


# =============================================================================
# Minimum and maximum county recall west of 85°W
# 5-year return period, native 0.25° heat-index dataset only
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(sf)
  library(tigris)
  library(stringr)
})

options(tigris_use_cache = TRUE)

METRICS_FILE <- paste0(
  "/proj/sebastianlab/users/hquintal/cisc/data/out/",
  "20260802_model_performance/precipitation/FF_W/metrics_county.csv"
)

# -----------------------------------------------------------------------------
# Read and filter performance data
# -----------------------------------------------------------------------------

metrics_native_rp5 <- read_csv(
  METRICS_FILE,
  col_types = cols(
    GEOID = col_character(),
    .default = col_guess()
  )
) %>%
  mutate(
    GEOID = str_pad(GEOID, width = 5, side = "left", pad = "0")
  ) %>%
  filter(
    return_period == 5,
    dataset_tag == "precipitation_0.25"
  ) %>%
  select(
    GEOID,
    precision
  )

# Confirm there is only one recall value per county
duplicate_geoids <- metrics_native_rp5 %>%
  count(GEOID) %>%
  filter(n > 1)

if (nrow(duplicate_geoids) > 0L) {
  stop(
    "Some GEOIDs still have more than one row after filtering:\n",
    paste(duplicate_geoids$GEOID, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Identify counties north of 35°N
# -----------------------------------------------------------------------------

counties_west_85 <- tigris::counties(
  cb = TRUE,
  resolution = "20m",
  year = 2023,
  class = "sf"
) %>%
  st_transform(4326) %>%
  select(
    GEOID,
    county = NAME,
    state = STUSPS,
    geometry
  ) %>%
  st_point_on_surface()

coords <- st_coordinates(counties_west_85)

counties_west_85 <- counties_west_85 %>%
  st_drop_geometry() %>%
  mutate(
    longitude = coords[, "X"],
    latitude  = coords[, "Y"]
  ) %>%
  filter(
    longitude < -85#,
    # latitude < 35
  )

# -----------------------------------------------------------------------------
# Join and calculate overall minimum and maximum recall
# -----------------------------------------------------------------------------

county_recall_west_85 <- metrics_native_rp5 %>%
  inner_join(
    counties_west_85,
    by = "GEOID"
  ) %>%
  filter(
    is.finite(precision)
  )

recall_min_max_west_85 <- county_recall_west_85 %>%
  summarise(
    n_counties = n_distinct(GEOID),
    mean_precision = mean(precision)
  )

print(recall_min_max_west_85)
