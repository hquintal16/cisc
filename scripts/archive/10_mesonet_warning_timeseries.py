# #!/usr/bin/env python3
# """
# Process Iowa Environmental Mesonet warning shapefiles for the southeastern
# United States, 2002-2023.
# 
# Stages
# ------
# annual:
#     Process one annual shapefile, clip warning polygons to the Southeast,
#     calculate clipped area, expand warnings across active UTC dates, and write
#     a daily summary by PHENOM code.
# 
# combine:
#     Combine annual summaries and create one two-row figure per PHENOM code.
# 
# Spatial domain
# --------------
# Longitude: -95 to -75
# Latitude:   25 to  40
# """
# 
# from __future__ import annotations
# 
# import argparse
# import gc
# import logging
# import math
# import os
# import sys
# from dataclasses import dataclass
# from pathlib import Path
# 
# import geopandas as gpd
# import matplotlib
# 
# matplotlib.use("Agg")
# 
# import matplotlib.dates as mdates
# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
# from pyproj import CRS
# import shapely
# from shapely.geometry import box
# 
# 
# # =============================================================================
# # Configuration
# # =============================================================================
# 
# DEFAULT_ROOT = Path(
#     os.environ.get(
#         "CISC_ROOT",
#         "/proj/sebastianlab/users/hquintal/cisc",
#     )
# )
# 
# START_YEAR = 2002
# END_YEAR = 2023
# YEARS = tuple(range(START_YEAR, END_YEAR + 1))
# 
# DOMAIN_WEST = -95.0
# DOMAIN_EAST = -75.0
# DOMAIN_SOUTH = 25.0
# DOMAIN_NORTH = 40.0
# 
# SOURCE_CRS = CRS.from_epsg(4326)
# AREA_CRS = CRS.from_epsg(5070)
# 
# SOUTHEAST_WGS84 = box(
#     DOMAIN_WEST,
#     DOMAIN_SOUTH,
#     DOMAIN_EAST,
#     DOMAIN_NORTH,
# )
# 
# SOUTHEAST_GDF = gpd.GeoDataFrame(
#     {"region": ["Southeast"]},
#     geometry=[SOUTHEAST_WGS84],
#     crs=SOURCE_CRS,
# )
# 
# REQUIRED_COLUMNS = [
#     "WFO",
#     "ISSUED",
#     "EXPIRED",
#     "PHENOM",
#     "SIG",
#     "ETN",
#     "STATUS",
#     "PROD_ID",
# ]
# 
# OUTPUT_COLUMNS = [
#     "date",
#     "year",
#     "phenom",
#     "phenom_name",
#     "event_count",
#     "area_km2",
# ]
# 
# VTEC_PHENOMENA = {
#     "AF": "Ashfall",
#     "AS": "Air Stagnation",
#     "BH": "Beach Hazard",
#     "BS": "Blowing Snow",
#     "BW": "Brisk Wind",
#     "BZ": "Blizzard",
#     "CF": "Coastal Flood",
#     "CW": "Cold Weather",
#     "DF": "Debris Flow",
#     "DS": "Dust Storm",
#     "DU": "Blowing Dust",
#     "EC": "Extreme Cold",
#     "EH": "Excessive Heat",
#     "EW": "Extreme Wind",
#     "FA": "Flood",
#     "FF": "Flash Flood",
#     "FG": "Dense Fog",
#     "FL": "Flood",
#     "FR": "Frost",
#     "FW": "Red Flag",
#     "FZ": "Freeze",
#     "GL": "Gale",
#     "HF": "Hurricane Force Wind",
#     "HI": "Inland Hurricane",
#     "HS": "Heavy Snow",
#     "HT": "Heat",
#     "HU": "Hurricane",
#     "HW": "High Wind",
#     "HY": "Hydrologic",
#     "HZ": "Hard Freeze",
#     "IP": "Sleet",
#     "IS": "Ice Storm",
#     "LB": "Lake Effect Snow and Blowing Snow",
#     "LE": "Lake Effect Snow",
#     "LO": "Low Water",
#     "LS": "Lakeshore Flood",
#     "LW": "Lake Wind",
#     "MA": "Marine",
#     "MF": "Marine Dense Fog",
#     "MH": "Marine Ashfall",
#     "MS": "Marine Dense Smoke",
#     "RB": "Small Craft for Rough Bar",
#     "RP": "Rip Current",
#     "SB": "Snow and Blowing Snow",
#     "SC": "Small Craft",
#     "SE": "Hazardous Seas",
#     "SI": "Small Craft for Winds",
#     "SM": "Dense Smoke",
#     "SN": "Snow",
#     "SQ": "Snow Squall",
#     "SR": "Storm",
#     "SS": "Storm Surge",
#     "SU": "High Surf",
#     "SV": "Severe Thunderstorm",
#     "SW": "Small Craft for Hazardous Seas",
#     "TI": "Inland Tropical Storm",
#     "TO": "Tornado",
#     "TR": "Tropical Storm",
#     "TS": "Tsunami",
#     "TY": "Typhoon",
#     "UP": "Freezing Spray",
#     "WC": "Wind Chill",
#     "WI": "Wind",
#     "WS": "Winter Storm",
#     "WW": "Winter Weather",
#     "XH": "Extreme Heat",
#     "ZF": "Freezing Fog",
#     "ZR": "Freezing Rain",
# }
# 
# 
# # =============================================================================
# # Paths
# # =============================================================================
# 
# 
# @dataclass(frozen=True)
# class WorkflowPaths:
#     root: Path
# 
#     @property
#     def input_dir(self) -> Path:
#         return self.root / "data" / "in" / "mesonet"
# 
#     @property
#     def output_dir(self) -> Path:
#         return self.root / "data" / "out" / "10_mesonet"
# 
#     @property
#     def annual_dir(self) -> Path:
#         return self.output_dir / "annual"
# 
#     @property
#     def figure_dir(self) -> Path:
#         return self.root / "figures" / "fig_s08_mesonet"
# 
#     @property
#     def combined_parquet(self) -> Path:
#         return (
#             self.output_dir
#             / "mesonet_warning_daily_southeast_2002_2023.parquet"
#         )
# 
#     @property
#     def combined_csv(self) -> Path:
#         return (
#             self.output_dir
#             / "mesonet_warning_daily_southeast_2002_2023.csv"
#         )
# 
#     @property
#     def lookup_csv(self) -> Path:
#         return self.output_dir / "iem_vtec_phenomenon_lookup.csv"
# 
#     @property
#     def figure_manifest_csv(self) -> Path:
#         return self.output_dir / "fig_s08_mesonet_manifest.csv"
# 
#     def source_year_dir(self, year: int) -> Path:
#         return self.input_dir / f"{year}_all"
# 
#     def annual_output(self, year: int) -> Path:
#         return (
#             self.annual_dir
#             / f"mesonet_warning_daily_{year}.parquet"
#         )
# 
# 
# # =============================================================================
# # Logging
# # =============================================================================
# 
# 
# def configure_logging() -> None:
#     logging.basicConfig(
#         level=logging.INFO,
#         format=(
#             "%(asctime)s | %(processName)s | "
#             "%(levelname)s | %(message)s"
#         ),
#         datefmt="%Y-%m-%d %H:%M:%S",
#         force=True,
#     )
# 
# 
# # =============================================================================
# # Helpers
# # =============================================================================
# 
# 
# def phenomenon_name(code: str) -> str:
#     return VTEC_PHENOMENA.get(
#         code,
#         f"Unknown phenomenon {code}",
#     )
# 
# 
# def empty_annual_result() -> pd.DataFrame:
#     return pd.DataFrame(
#         {
#             "date": pd.Series(dtype="datetime64[ns]"),
#             "year": pd.Series(dtype="int16"),
#             "phenom": pd.Series(dtype="string"),
#             "phenom_name": pd.Series(dtype="string"),
#             "event_count": pd.Series(dtype="int64"),
#             "area_km2": pd.Series(dtype="float64"),
#         }
#     )
# 
# 
# def parse_iem_datetime(values: pd.Series) -> pd.Series:
#     cleaned = (
#         values.astype("string")
#         .str.strip()
#         .str.replace(r"\.0$", "", regex=True)
#         .replace(
#             {
#                 "": pd.NA,
#                 "NA": pd.NA,
#                 "None": pd.NA,
#                 "<NA>": pd.NA,
#             }
#         )
#     )
# 
#     return pd.to_datetime(
#         cleaned,
#         format="%Y%m%d%H%M",
#         errors="coerce",
#     )
# 
# 
# def write_parquet_atomic(
#     frame: pd.DataFrame,
#     output_file: Path,
# ) -> None:
#     output_file.parent.mkdir(
#         parents=True,
#         exist_ok=True,
#     )
# 
#     temporary_file = output_file.with_name(
#         output_file.name + ".tmp"
#     )
# 
#     frame.to_parquet(
#         temporary_file,
#         index=False,
#         compression="zstd",
#     )
# 
#     temporary_file.replace(output_file)
# 
# 
# def valid_annual_cache(output_file: Path) -> bool:
#     if not output_file.is_file():
#         return False
# 
#     try:
#         cached = pd.read_parquet(output_file)
#     except Exception:
#         return False
# 
#     return set(OUTPUT_COLUMNS).issubset(cached.columns)
# 
# 
# def find_annual_shapefile(
#     paths: WorkflowPaths,
#     year: int,
# ) -> Path:
#     year_dir = paths.source_year_dir(year)
# 
#     if not year_dir.is_dir():
#         raise FileNotFoundError(
#             f"Missing annual input directory: {year_dir}"
#         )
# 
#     expected = (
#         year_dir
#         / f"wwa_{year}01010000_{year}12312359.shp"
#     )
# 
#     if expected.is_file():
#         return expected
# 
#     candidates = sorted(year_dir.glob("wwa_*.shp"))
# 
#     if len(candidates) == 1:
#         return candidates[0]
# 
#     if not candidates:
#         raise FileNotFoundError(
#             f"No wwa_*.shp file found in {year_dir}"
#         )
# 
#     raise RuntimeError(
#         f"Multiple shapefiles found for {year}: "
#         + ", ".join(str(path) for path in candidates)
#     )
# 
# 
# # =============================================================================
# # Geometry processing
# # =============================================================================
# 
# 
# def repair_geometries(
#     frame: gpd.GeoDataFrame,
# ) -> gpd.GeoDataFrame:
#     frame = frame.loc[
#         frame.geometry.notna()
#         & ~frame.geometry.is_empty
#     ].copy()
# 
#     if frame.empty:
#         return frame
# 
#     invalid = ~frame.geometry.is_valid
# 
#     if invalid.any():
#         logging.info(
#             "Repairing %s invalid geometries",
#             f"{int(invalid.sum()):,}",
#         )
# 
#         repaired = shapely.make_valid(
#             frame.loc[invalid, "geometry"].array
#         )
# 
#         frame.loc[invalid, "geometry"] = gpd.GeoSeries(
#             repaired,
#             index=frame.index[invalid],
#             crs=frame.crs,
#         )
# 
#     return frame.loc[
#         frame.geometry.notna()
#         & ~frame.geometry.is_empty
#     ].copy()
# 
# 
# def clip_and_calculate_area(
#     frame: gpd.GeoDataFrame,
# ) -> gpd.GeoDataFrame:
#     if frame.empty:
#         frame["clipped_area_km2"] = pd.Series(
#             dtype="float64"
#         )
#         return frame
# 
#     if frame.crs is None:
#         logging.warning(
#             "Input CRS is missing; assuming EPSG:4326"
#         )
#         frame = frame.set_crs(SOURCE_CRS)
#     elif CRS.from_user_input(frame.crs) != SOURCE_CRS:
#         frame = frame.to_crs(SOURCE_CRS)
# 
#     frame = repair_geometries(frame)
# 
#     if frame.empty:
#         frame["clipped_area_km2"] = pd.Series(
#             dtype="float64"
#         )
#         return frame
# 
#     projected = frame.to_crs(AREA_CRS)
# 
#     projected_domain = (
#         SOUTHEAST_GDF
#         .to_crs(AREA_CRS)
#         .geometry
#         .iloc[0]
#     )
# 
#     dropped = 0
# 
#     def safe_intersection(geometry):
#         nonlocal dropped
# 
#         if geometry is None or geometry.is_empty:
#             dropped += 1
#             return None
# 
#         try:
#             repaired = shapely.make_valid(geometry)
# 
#             repaired = shapely.set_precision(
#                 repaired,
#                 grid_size=1.0,
#             )
# 
#             clipped = shapely.intersection(
#                 repaired,
#                 projected_domain,
#                 grid_size=1.0,
#             )
# 
#             if clipped is None or clipped.is_empty:
#                 return None
# 
#             return clipped
# 
#         except (
#             shapely.GEOSException,
#             ValueError,
#             TypeError,
#         ):
#             dropped += 1
#             return None
# 
#     projected["geometry"] = projected.geometry.apply(
#         safe_intersection
#     )
# 
#     projected = projected.loc[
#         projected.geometry.notna()
#         & ~projected.geometry.is_empty
#     ].copy()
# 
#     if dropped > 0:
#         logging.warning(
#             "Dropped %s unrecoverable geometries",
#             f"{dropped:,}",
#         )
# 
#     if projected.empty:
#         projected["clipped_area_km2"] = pd.Series(
#             dtype="float64"
#         )
#         return projected
# 
#     polygonal = projected.geometry.geom_type.isin(
#         [
#             "Polygon",
#             "MultiPolygon",
#         ]
#     )
# 
#     projected = projected.loc[
#         polygonal
#     ].copy()
# 
#     if projected.empty:
#         projected["clipped_area_km2"] = pd.Series(
#             dtype="float64"
#         )
#         return projected
# 
#     projected["clipped_area_km2"] = (
#         projected.geometry.area / 1_000_000.0
#     )
# 
#     projected = projected.loc[
#         np.isfinite(projected["clipped_area_km2"])
#         & (projected["clipped_area_km2"] > 0)
#     ].copy()
# 
#     return projected
# 
# 
# # =============================================================================
# # Daily expansion
# # =============================================================================
# 
# 
# def expand_records_to_daily(
#     frame: pd.DataFrame,
# ) -> pd.DataFrame:
#     if frame.empty:
#         return empty_annual_result()
# 
#     frame = frame.copy()
# 
#     frame["start_date"] = (
#         frame["issued_time"].dt.normalize()
#     )
# 
#     frame["end_date"] = (
#         frame["expired_time"].dt.normalize()
#     )
# 
#     valid = (
#         frame["start_date"].notna()
#         & frame["end_date"].notna()
#         & (frame["end_date"] >= frame["start_date"])
#     )
# 
#     frame = frame.loc[valid].copy()
# 
#     if frame.empty:
#         return empty_annual_result()
# 
#     analysis_start = pd.Timestamp(
#         f"{START_YEAR}-01-01"
#     )
# 
#     analysis_end = pd.Timestamp(
#         f"{END_YEAR}-12-31"
#     )
# 
#     frame["start_date"] = frame["start_date"].clip(
#         lower=analysis_start,
#         upper=analysis_end,
#     )
# 
#     frame["end_date"] = frame["end_date"].clip(
#         lower=analysis_start,
#         upper=analysis_end,
#     )
# 
#     duration_days = (
#         frame["end_date"] - frame["start_date"]
#     ).dt.days + 1
# 
#     valid_duration = duration_days.between(
#         1,
#         730,
#     )
# 
#     if (~valid_duration).any():
#         logging.warning(
#             "Dropping %s records with active intervals "
#             "outside 1-730 days",
#             f"{int((~valid_duration).sum()):,}",
#         )
# 
#     frame = frame.loc[valid_duration].copy()
#     duration_days = duration_days.loc[valid_duration]
# 
#     if frame.empty:
#         return empty_annual_result()
# 
#     repeated_rows = np.repeat(
#         np.arange(len(frame), dtype=np.int64),
#         duration_days.to_numpy(dtype=np.int64),
#     )
# 
#     day_offsets = np.concatenate(
#         [
#             np.arange(
#                 number_days,
#                 dtype="timedelta64[D]",
#             )
#             for number_days in duration_days.to_numpy(
#                 dtype=np.int64
#             )
#         ]
#     )
# 
#     expanded = pd.DataFrame(
#         {
#             "date": (
#                 frame["start_date"]
#                 .to_numpy(dtype="datetime64[D]")[
#                     repeated_rows
#                 ]
#                 + day_offsets
#             ),
#             "phenom": (
#                 frame["PHENOM"]
#                 .to_numpy(dtype=object)[repeated_rows]
#             ),
#             "area_km2": (
#                 frame["clipped_area_km2"]
#                 .to_numpy(dtype=np.float64)[
#                     repeated_rows
#                 ]
#             ),
#         }
#     )
# 
#     expanded["date"] = pd.to_datetime(
#         expanded["date"]
#     )
# 
#     expanded["event_count"] = 1
# 
#     daily = (
#         expanded.groupby(
#             [
#                 "date",
#                 "phenom",
#             ],
#             as_index=False,
#             observed=True,
#         )
#         .agg(
#             event_count=("event_count", "sum"),
#             area_km2=("area_km2", "sum"),
#         )
#         .sort_values(
#             [
#                 "phenom",
#                 "date",
#             ]
#         )
#         .reset_index(drop=True)
#     )
# 
#     daily["year"] = (
#         daily["date"].dt.year.astype("int16")
#     )
# 
#     daily["phenom_name"] = (
#         daily["phenom"]
#         .astype(str)
#         .map(phenomenon_name)
#     )
# 
#     return daily[OUTPUT_COLUMNS]
# 
# 
# # =============================================================================
# # Annual stage
# # =============================================================================
# 
# 
# def process_annual_stage(
#     paths: WorkflowPaths,
#     year: int,
#     overwrite: bool,
# ) -> None:
#     if year not in YEARS:
#         raise ValueError(
#             f"Year must be between {START_YEAR} "
#             f"and {END_YEAR}: {year}"
#         )
# 
#     output_file = paths.annual_output(year)
# 
#     if (
#         not overwrite
#         and valid_annual_cache(output_file)
#     ):
#         logging.info(
#             "Year %d is already complete: %s",
#             year,
#             output_file,
#         )
#         return
# 
#     shapefile = find_annual_shapefile(
#         paths,
#         year,
#     )
# 
#     logging.info(
#         "Reading year %d: %s",
#         year,
#         shapefile,
#     )
# 
#     frame = gpd.read_file(
#         shapefile,
#         engine="pyogrio",
#         columns=REQUIRED_COLUMNS,
#         bbox=(
#             DOMAIN_WEST,
#             DOMAIN_SOUTH,
#             DOMAIN_EAST,
#             DOMAIN_NORTH,
#         ),
#         use_arrow=True,
#     )
# 
#     logging.info(
#         "Year %d: %s source records intersect bounding box",
#         year,
#         f"{len(frame):,}",
#     )
# 
#     if frame.empty:
#         write_parquet_atomic(
#             empty_annual_result(),
#             output_file,
#         )
#         return
# 
#     frame.columns = [
#         (
#             column.upper()
#             if column != "geometry"
#             else column
#         )
#         for column in frame.columns
#     ]
# 
#     missing = set(REQUIRED_COLUMNS).difference(
#         frame.columns
#     )
# 
#     if missing:
#         raise KeyError(
#             f"Year {year} is missing columns: "
#             f"{sorted(missing)}"
#         )
# 
#     frame["PHENOM"] = (
#         frame["PHENOM"]
#         .astype("string")
#         .str.strip()
#         .str.upper()
#     )
# 
#     frame["SIG"] = (
#         frame["SIG"]
#         .astype("string")
#         .str.strip()
#         .str.upper()
#     )
# 
#     frame = frame.loc[
#         frame["PHENOM"].notna()
#         & frame["PHENOM"].ne("")
#     ].copy()
# 
#     frame["issued_time"] = parse_iem_datetime(
#         frame["ISSUED"]
#     )
# 
#     frame["expired_time"] = parse_iem_datetime(
#         frame["EXPIRED"]
#     )
# 
#     valid_time = (
#         frame["issued_time"].notna()
#         & frame["expired_time"].notna()
#         & (
#             frame["expired_time"]
#             >= frame["issued_time"]
#         )
#     )
# 
#     if (~valid_time).any():
#         logging.warning(
#             "Year %d: dropping %s records with invalid "
#             "ISSUED or EXPIRED values",
#             year,
#             f"{int((~valid_time).sum()):,}",
#         )
# 
#     frame = frame.loc[
#         valid_time
#     ].copy()
# 
#     frame = clip_and_calculate_area(
#         frame
#     )
# 
#     logging.info(
#         "Year %d: %s clipped polygons retained",
#         year,
#         f"{len(frame):,}",
#     )
# 
#     daily = expand_records_to_daily(
#         frame
#     )
# 
#     write_parquet_atomic(
#         daily,
#         output_file,
#     )
# 
#     logging.info(
#         "Year %d: wrote %s daily phenomenon rows: %s",
#         year,
#         f"{len(daily):,}",
#         output_file,
#     )
# 
#     del frame
#     del daily
#     gc.collect()
# 
# 
# # =============================================================================
# # Combine stage
# # =============================================================================
# 
# 
# def combine_annual_results(
#     paths: WorkflowPaths,
# ) -> pd.DataFrame:
#     missing_files = [
#         paths.annual_output(year)
#         for year in YEARS
#         if not valid_annual_cache(
#             paths.annual_output(year)
#         )
#     ]
# 
#     if missing_files:
#         raise FileNotFoundError(
#             "Annual outputs are missing or invalid:\n"
#             + "\n".join(
#                 str(path)
#                 for path in missing_files
#             )
#         )
# 
#     annual_frames: list[pd.DataFrame] = []
# 
#     for year in YEARS:
#         annual_file = paths.annual_output(
#             year
#         )
# 
#         logging.info(
#             "Reading annual output: %s",
#             annual_file,
#         )
# 
#         annual = pd.read_parquet(
#             annual_file
#         )
# 
#         if not annual.empty:
#             annual_frames.append(
#                 annual
#             )
# 
#     if not annual_frames:
#         raise RuntimeError(
#             "All annual outputs are empty."
#         )
# 
#     combined = pd.concat(
#         annual_frames,
#         ignore_index=True,
#     )
# 
#     combined["date"] = pd.to_datetime(
#         combined["date"]
#     )
# 
#     combined["phenom"] = (
#         combined["phenom"]
#         .astype("string")
#         .str.upper()
#     )
# 
#     combined["phenom_name"] = (
#         combined["phenom"]
#         .astype(str)
#         .map(phenomenon_name)
#     )
# 
#     combined = (
#         combined.groupby(
#             [
#                 "date",
#                 "phenom",
#                 "phenom_name",
#             ],
#             as_index=False,
#             observed=True,
#         )
#         .agg(
#             event_count=("event_count", "sum"),
#             area_km2=("area_km2", "sum"),
#         )
#         .sort_values(
#             [
#                 "phenom",
#                 "date",
#             ]
#         )
#         .reset_index(drop=True)
#     )
# 
#     combined = combined.loc[
#         combined["date"].between(
#             pd.Timestamp(
#                 f"{START_YEAR}-01-01"
#             ),
#             pd.Timestamp(
#                 f"{END_YEAR}-12-31"
#             ),
#         )
#     ].copy()
# 
#     combined["year"] = (
#         combined["date"].dt.year.astype("int16")
#     )
# 
#     combined = combined[
#         OUTPUT_COLUMNS
#     ]
# 
#     write_parquet_atomic(
#         combined,
#         paths.combined_parquet,
#     )
# 
#     combined.to_csv(
#         paths.combined_csv,
#         index=False,
#         date_format="%Y-%m-%d",
#     )
# 
#     logging.info(
#         "Wrote combined Parquet: %s",
#         paths.combined_parquet,
#     )
# 
#     logging.info(
#         "Wrote combined CSV: %s",
#         paths.combined_csv,
#     )
# 
#     return combined
# 
# 
# # =============================================================================
# # Lookup table
# # =============================================================================
# 
# 
# def write_lookup_table(
#     paths: WorkflowPaths,
#     observed_codes: list[str],
# ) -> None:
#     observed_set = set(
#         observed_codes
#     )
# 
#     all_codes = sorted(
#         set(VTEC_PHENOMENA).union(
#             observed_set
#         )
#     )
# 
#     lookup = pd.DataFrame(
#         {
#             "phenom": all_codes,
#             "phenom_name": [
#                 phenomenon_name(code)
#                 for code in all_codes
#             ],
#             "observed_2002_2023": [
#                 code in observed_set
#                 for code in all_codes
#             ],
#         }
#     )
# 
#     lookup.to_csv(
#         paths.lookup_csv,
#         index=False,
#     )
# 
#     logging.info(
#         "Wrote phenomenon lookup: %s",
#         paths.lookup_csv,
#     )
# 
# 
# # =============================================================================
# # Plotting
# # =============================================================================
# 
# 
# def complete_daily_series(
#     phenomenon_data: pd.DataFrame,
# ) -> pd.DataFrame:
#     full_dates = pd.date_range(
#         start=f"{START_YEAR}-01-01",
#         end=f"{END_YEAR}-12-31",
#         freq="D",
#     )
# 
#     daily = (
#         phenomenon_data.groupby(
#             "date",
#             as_index=True,
#         )
#         .agg(
#             event_count=("event_count", "sum"),
#             area_km2=("area_km2", "sum"),
#         )
#         .reindex(
#             full_dates,
#             fill_value=0,
#         )
#         .rename_axis("date")
#         .reset_index()
#     )
# 
#     daily["event_count"] = (
#         daily["event_count"]
#         .fillna(0)
#         .astype("int64")
#     )
# 
#     daily["area_km2"] = (
#         daily["area_km2"]
#         .fillna(0.0)
#         .astype("float64")
#     )
# 
#     return daily
# 
# 
# def rounded_axis_upper(
#     maximum: float,
# ) -> float:
#     if (
#         not np.isfinite(maximum)
#         or maximum <= 0
#     ):
#         return 1.0
# 
#     exponent = math.floor(
#         math.log10(maximum)
#     )
# 
#     magnitude = 10**exponent
#     normalized = maximum / magnitude
# 
#     for candidate in (
#         1,
#         2,
#         2.5,
#         5,
#         10,
#     ):
#         if normalized <= candidate:
#             return candidate * magnitude
# 
#     return 10 * magnitude
# 
# 
# def plot_phenomenon(
#     paths: WorkflowPaths,
#     combined: pd.DataFrame,
#     phenom: str,
# ) -> dict[str, object]:
#     subset = combined.loc[
#         combined["phenom"].eq(phenom)
#     ].copy()
# 
#     if subset.empty:
#         raise ValueError(
#             f"No records found for {phenom}"
#         )
# 
#     name = phenomenon_name(
#         phenom
#     )
# 
#     daily = complete_daily_series(
#         subset
#     )
# 
#     output_file = (
#         paths.figure_dir
#         / f"fig_s08_mesonet_{phenom.lower()}.png"
#     )
# 
#     figure, axes = plt.subplots(
#         nrows=2,
#         ncols=1,
#         figsize=(6.5, 5.6),
#         sharex=True,
#         constrained_layout=True,
#     )
# 
#     count_axis = axes[0]
#     area_axis = axes[1]
# 
#     count_axis.plot(
#         daily["date"],
#         daily["event_count"],
#         linewidth=0.45,
#     )
# 
#     area_axis.plot(
#         daily["date"],
#         daily["area_km2"],
#         linewidth=0.45,
#     )
# 
#     count_axis.set_ylabel(
#         "Active warnings per day",
#         fontsize=8,
#     )
# 
#     area_axis.set_ylabel(
#         "Active warning area\nper day (km²)",
#         fontsize=8,
#     )
# 
#     area_axis.set_xlabel(
#         "Year",
#         fontsize=8,
#     )
# 
#     count_axis.set_ylim(
#         0,
#         rounded_axis_upper(
#             float(
#                 daily["event_count"].max()
#             )
#         ),
#     )
# 
#     area_axis.set_ylim(
#         0,
#         rounded_axis_upper(
#             float(
#                 daily["area_km2"].max()
#             )
#         ),
#     )
# 
#     for axis in axes:
#         axis.grid(
#             axis="y",
#             linewidth=0.35,
#             alpha=0.35,
#         )
# 
#         axis.tick_params(
#             axis="both",
#             labelsize=7,
#             width=0.6,
#             length=3,
#         )
# 
#         axis.spines["top"].set_visible(
#             False
#         )
# 
#         axis.spines["right"].set_visible(
#             False
#         )
# 
#     area_axis.xaxis.set_major_locator(
#         mdates.YearLocator(base=2)
#     )
# 
#     area_axis.xaxis.set_major_formatter(
#         mdates.DateFormatter("%Y")
#     )
# 
#     area_axis.set_xlim(
#         pd.Timestamp(
#             f"{START_YEAR}-01-01"
#         ),
#         pd.Timestamp(
#             f"{END_YEAR}-12-31"
#         ),
#     )
# 
#     figure.suptitle(
#         (
#             f"{phenom} — {name}\n"
#             "Daily warning count and warning area "
#             "in the Southeast, 2002–2023"
#         ),
#         fontsize=10,
#         fontweight="bold",
#     )
# 
#     figure.savefig(
#         output_file,
#         dpi=300,
#         bbox_inches="tight",
#         facecolor="white",
#     )
# 
#     plt.close(
#         figure
#     )
# 
#     active_dates = daily.loc[
#         daily["event_count"] > 0,
#         "date",
#     ]
# 
#     logging.info(
#         "Wrote figure: %s",
#         output_file,
#     )
# 
#     return {
#         "phenom": phenom,
#         "phenom_name": name,
#         "figure_file": str(output_file),
#         "first_active_date": (
#             active_dates.min()
#             .date()
#             .isoformat()
#             if not active_dates.empty
#             else None
#         ),
#         "last_active_date": (
#             active_dates.max()
#             .date()
#             .isoformat()
#             if not active_dates.empty
#             else None
#         ),
#         "number_active_dates": int(
#             (
#                 daily["event_count"] > 0
#             ).sum()
#         ),
#         "maximum_daily_count": int(
#             daily["event_count"].max()
#         ),
#         "maximum_daily_area_km2": float(
#             daily["area_km2"].max()
#         ),
#     }
# 
# 
# def plot_all_phenomena(
#     paths: WorkflowPaths,
#     combined: pd.DataFrame,
# ) -> None:
#     paths.figure_dir.mkdir(
#         parents=True,
#         exist_ok=True,
#     )
# 
#     phenomena = sorted(
#         combined["phenom"]
#         .dropna()
#         .astype(str)
#         .unique()
#         .tolist()
#     )
# 
#     manifest_rows: list[
#         dict[str, object]
#     ] = []
# 
#     for index, phenom in enumerate(
#         phenomena,
#         start=1,
#     ):
#         logging.info(
#             "Plotting %02d/%02d: %s — %s",
#             index,
#             len(phenomena),
#             phenom,
#             phenomenon_name(phenom),
#         )
# 
#         manifest_rows.append(
#             plot_phenomenon(
#                 paths,
#                 combined,
#                 phenom,
#             )
#         )
# 
#     manifest = pd.DataFrame(
#         manifest_rows
#     )
# 
#     manifest.to_csv(
#         paths.figure_manifest_csv,
#         index=False,
#     )
# 
#     logging.info(
#         "Wrote figure manifest: %s",
#         paths.figure_manifest_csv,
#     )
# 
# 
# def process_combine_stage(
#     paths: WorkflowPaths,
# ) -> None:
#     combined = combine_annual_results(
#         paths
#     )
# 
#     observed_codes = sorted(
#         combined["phenom"]
#         .dropna()
#         .astype(str)
#         .unique()
#         .tolist()
#     )
# 
#     write_lookup_table(
#         paths,
#         observed_codes,
#     )
# 
#     plot_all_phenomena(
#         paths,
#         combined,
#     )
# 
#     logging.info(
#         "Combine and plotting stage complete."
#     )
# 
# 
# # =============================================================================
# # Command line
# # =============================================================================
# 
# 
# def parse_arguments() -> argparse.Namespace:
#     parser = argparse.ArgumentParser(
#         description=(
#             "Process and plot IEM warning data "
#             "for the Southeast."
#         )
#     )
# 
#     parser.add_argument(
#         "--root",
#         type=Path,
#         default=DEFAULT_ROOT,
#         help="CISC project root.",
#     )
# 
#     parser.add_argument(
#         "--stage",
#         required=True,
#         choices=[
#             "annual",
#             "combine",
#         ],
#         help="Workflow stage.",
#     )
# 
#     parser.add_argument(
#         "--year",
#         type=int,
#         default=None,
#         help=(
#             "Year to process during the annual stage."
#         ),
#     )
# 
#     parser.add_argument(
#         "--overwrite",
#         action="store_true",
#         help=(
#             "Overwrite an existing valid annual output."
#         ),
#     )
# 
#     return parser.parse_args()
# 
# 
# def main() -> int:
#     configure_logging()
#     arguments = parse_arguments()
# 
#     paths = WorkflowPaths(
#         arguments.root
#         .expanduser()
#         .resolve()
#     )
# 
#     paths.output_dir.mkdir(
#         parents=True,
#         exist_ok=True,
#     )
# 
#     paths.annual_dir.mkdir(
#         parents=True,
#         exist_ok=True,
#     )
# 
#     paths.figure_dir.mkdir(
#         parents=True,
#         exist_ok=True,
#     )
# 
#     logging.info(
#         "CISC root: %s",
#         paths.root,
#     )
# 
#     logging.info(
#         "Stage: %s",
#         arguments.stage,
#     )
# 
#     if arguments.stage == "annual":
#         if arguments.year is None:
#             raise ValueError(
#                 "--year is required for "
#                 "--stage annual"
#             )
# 
#         process_annual_stage(
#             paths=paths,
#             year=arguments.year,
#             overwrite=arguments.overwrite,
#         )
# 
#     elif arguments.stage == "combine":
#         process_combine_stage(
#             paths
#         )
# 
#     return 0
# 
# 
# if __name__ == "__main__":
#     try:
#         sys.exit(
#             main()
#         )
#     except Exception:
#         logging.exception(
#             "Mesonet workflow failed"
#         )
#         sys.exit(1)

#!/usr/bin/env python3
"""
Process Iowa Environmental Mesonet warning polygons for the southeastern
United States, 2002-2023.

Warning groups
--------------
Heat:
    HT — Heat
    EH — Excessive Heat

Precipitation:
    FF — Flash Flood
    SV — Severe Thunderstorm
    TR — Tropical Storm
    GL — Gale
    HI — Inland Hurricane
    HU — Hurricane
    SR — Storm
    TI — Inland Tropical Storm

Annual stage
------------
For one year:

1. Read only the selected warning types.
2. Restrict polygons to 25-40 N and 95-75 W.
3. Expand each warning from ISSUED through EXPIRED to every active UTC date.
4. Within each date and warning group, union all active polygons.
5. Calculate:
       event_count = number of spatially distinct unioned polygon components
       area_km2    = total unioned polygon area
6. Write one annual Parquet file.

Combining geometries before calculating count and area prevents spatially and
temporally overlapping warning records and warning types from being counted
more than once on a date.

Combine stage
-------------
1. Combine all annual files.
2. Read daily heat-index and precipitation Ar.
3. Exclude dates with no active warnings.
4. Create two figures:
       heat warnings versus heat-index Ar
       precipitation warnings versus precipitation Ar

Each figure contains:
    left:  Ar versus spatially distinct active-warning count
    right: Ar versus unioned active-warning area
"""

from __future__ import annotations

import argparse
import gc
import logging
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import CRS
from scipy.stats import pearsonr
import shapely
from shapely.geometry import box


# =============================================================================
# Configuration
# =============================================================================

DEFAULT_ROOT = Path(
    os.environ.get(
        "CISC_ROOT",
        "/proj/sebastianlab/users/hquintal/cisc",
    )
)

START_YEAR = 2002
END_YEAR = 2023
YEARS = tuple(range(START_YEAR, END_YEAR + 1))

START_DATE = pd.Timestamp(f"{START_YEAR}-01-01")
END_DATE = pd.Timestamp(f"{END_YEAR}-12-31")

DOMAIN_WEST = -95.0
DOMAIN_EAST = -75.0
DOMAIN_SOUTH = 25.0
DOMAIN_NORTH = 40.0

SOURCE_CRS = CRS.from_epsg(4326)
AREA_CRS = CRS.from_epsg(5070)

SOUTHEAST_WGS84 = box(
    DOMAIN_WEST,
    DOMAIN_SOUTH,
    DOMAIN_EAST,
    DOMAIN_NORTH,
)

SOUTHEAST_GDF = gpd.GeoDataFrame(
    {"region": ["Southeast"]},
    geometry=[SOUTHEAST_WGS84],
    crs=SOURCE_CRS,
)

HEAT_TYPES = {
    "HT",
    "EH",
}

PRECIPITATION_TYPES = {
    "FF",
    "SV",
    "TR",
    "HI",
    "HU",
    "SR",
    "TI",
}

SELECTED_TYPES = HEAT_TYPES | PRECIPITATION_TYPES

GROUP_CODES = {
    "heat": sorted(HEAT_TYPES),
    "precipitation": sorted(PRECIPITATION_TYPES),
}

GROUP_TITLES = {
    "heat": "Heat and excessive heat warnings",
    "precipitation": (
        "Flash flood, severe thunderstorm, tropical storm, "
        "hurricane, and storm warnings"
    ),
}

REQUIRED_COLUMNS = [
    "ISSUED",
    "EXPIRED",
    "PHENOM",
]

ANNUAL_OUTPUT_COLUMNS = [
    "date",
    "year",
    "warning_group",
    "event_count",
    "area_km2",
]

CORRELATION_OUTPUT_COLUMNS = [
    "date",
    "year",
    "warning_group",
    "Ar",
    "event_count",
    "area_km2",
]


# =============================================================================
# Paths
# =============================================================================


@dataclass(frozen=True)
class WorkflowPaths:
    root: Path

    @property
    def input_dir(self) -> Path:
        return self.root / "data" / "in" / "mesonet"

    @property
    def output_dir(self) -> Path:
        return self.root / "data" / "out" / "10_mesonet"

    @property
    def annual_dir(self) -> Path:
        return self.output_dir / "annual_combined"

    @property
    def correlation_dir(self) -> Path:
        return self.output_dir / "ar_correlations"

    @property
    def figure_dir(self) -> Path:
        return self.root / "figures" / "fig_s09_mesonet_ar"

    @property
    def combined_parquet(self) -> Path:
        return (
            self.output_dir
            / "mesonet_warning_combined_daily_southeast_2002_2023.parquet"
        )

    @property
    def combined_csv(self) -> Path:
        return (
            self.output_dir
            / "mesonet_warning_combined_daily_southeast_2002_2023.csv"
        )

    @property
    def group_lookup_csv(self) -> Path:
        return self.output_dir / "mesonet_warning_group_lookup.csv"

    @property
    def correlation_summary_csv(self) -> Path:
        return self.correlation_dir / "mesonet_ar_correlation_summary.csv"

    @property
    def heat_ar_file(self) -> Path:
        return (
            self.root
            / "data"
            / "out"
            / "03_cov_model"
            / "heatindex"
            / "heatindex_daily_optimal_cov_model_southeast.csv"
        )

    @property
    def precipitation_ar_file(self) -> Path:
        return (
            self.root
            / "data"
            / "out"
            / "03_cov_model"
            / "precipitation"
            / "precipitation_daily_optimal_cov_model_southeast.csv"
        )

    def source_year_dir(self, year: int) -> Path:
        return self.input_dir / f"{year}_all"

    def annual_output(self, year: int) -> Path:
        return (
            self.annual_dir
            / f"mesonet_warning_combined_daily_{year}.parquet"
        )

    def correlation_output(self, warning_group: str) -> Path:
        return (
            self.correlation_dir
            / f"mesonet_ar_{warning_group}_2002_2023.csv"
        )

    def figure_output(self, warning_group: str) -> Path:
        return (
            self.figure_dir
            / f"fig_s09_mesonet_ar_{warning_group}.png"
        )


# =============================================================================
# Logging
# =============================================================================


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format=(
            "%(asctime)s | %(processName)s | "
            "%(levelname)s | %(message)s"
        ),
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )


# =============================================================================
# General helpers
# =============================================================================


def empty_annual_result() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "date": pd.Series(dtype="datetime64[ns]"),
            "year": pd.Series(dtype="int16"),
            "warning_group": pd.Series(dtype="string"),
            "event_count": pd.Series(dtype="int64"),
            "area_km2": pd.Series(dtype="float64"),
        }
    )


def parse_iem_datetime(values: pd.Series) -> pd.Series:
    cleaned = (
        values.astype("string")
        .str.strip()
        .str.replace(r"\.0$", "", regex=True)
        .replace(
            {
                "": pd.NA,
                "NA": pd.NA,
                "None": pd.NA,
                "<NA>": pd.NA,
            }
        )
    )

    return pd.to_datetime(
        cleaned,
        format="%Y%m%d%H%M",
        errors="coerce",
    )


def warning_group_from_code(code: str) -> str | None:
    if code in HEAT_TYPES:
        return "heat"

    if code in PRECIPITATION_TYPES:
        return "precipitation"

    return None


def write_parquet_atomic(
    frame: pd.DataFrame,
    output_file: Path,
) -> None:
    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary_file = output_file.with_name(
        output_file.name + ".tmp"
    )

    frame.to_parquet(
        temporary_file,
        index=False,
        compression="zstd",
    )

    temporary_file.replace(output_file)


def valid_annual_cache(output_file: Path) -> bool:
    if not output_file.is_file():
        return False

    try:
        cached = pd.read_parquet(output_file)
    except Exception:
        return False

    return set(ANNUAL_OUTPUT_COLUMNS).issubset(
        cached.columns
    )


def find_annual_shapefile(
    paths: WorkflowPaths,
    year: int,
) -> Path:
    year_dir = paths.source_year_dir(year)

    if not year_dir.is_dir():
        raise FileNotFoundError(
            f"Missing annual input directory: {year_dir}"
        )

    expected = (
        year_dir
        / f"wwa_{year}01010000_{year}12312359.shp"
    )

    if expected.is_file():
        return expected

    candidates = sorted(year_dir.glob("wwa_*.shp"))

    if len(candidates) == 1:
        return candidates[0]

    if not candidates:
        raise FileNotFoundError(
            f"No wwa_*.shp file found in {year_dir}"
        )

    raise RuntimeError(
        f"Multiple shapefiles found for {year}: "
        + ", ".join(str(path) for path in candidates)
    )


# =============================================================================
# Geometry processing
# =============================================================================


def safe_clip_geometry(
    geometry,
    domain,
):
    if geometry is None or geometry.is_empty:
        return None

    try:
        repaired = shapely.make_valid(geometry)

        repaired = shapely.set_precision(
            repaired,
            grid_size=1.0,
        )

        clipped = shapely.intersection(
            repaired,
            domain,
            grid_size=1.0,
        )

        if clipped is None or clipped.is_empty:
            return None

        return clipped

    except (
        shapely.GEOSException,
        ValueError,
        TypeError,
    ):
        return None


def clip_warning_geometries(
    frame: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    if frame.empty:
        return frame

    frame = frame.loc[
        frame.geometry.notna()
        & ~frame.geometry.is_empty
    ].copy()

    if frame.empty:
        return frame

    if frame.crs is None:
        logging.warning(
            "Input CRS is missing; assuming EPSG:4326"
        )
        frame = frame.set_crs(SOURCE_CRS)
    elif CRS.from_user_input(frame.crs) != SOURCE_CRS:
        frame = frame.to_crs(SOURCE_CRS)

    projected = frame.to_crs(AREA_CRS)

    projected_domain = (
        SOUTHEAST_GDF
        .to_crs(AREA_CRS)
        .geometry
        .iloc[0]
    )

    before = len(projected)

    projected["geometry"] = projected.geometry.apply(
        safe_clip_geometry,
        domain=projected_domain,
    )

    projected = projected.loc[
        projected.geometry.notna()
        & ~projected.geometry.is_empty
    ].copy()

    dropped = before - len(projected)

    if dropped:
        logging.warning(
            "Dropped %s empty or unrecoverable geometries",
            f"{dropped:,}",
        )

    return projected


def iter_polygon_parts(geometry) -> Iterable:
    if geometry is None or geometry.is_empty:
        return

    geometry_type = geometry.geom_type

    if geometry_type == "Polygon":
        yield geometry
        return

    if geometry_type in {
        "MultiPolygon",
        "GeometryCollection",
    }:
        for part in geometry.geoms:
            yield from iter_polygon_parts(part)


def safe_union_metrics(
    geometries: Iterable,
) -> tuple[int, float]:
    geometry_list = [
        geometry
        for geometry in geometries
        if geometry is not None
        and not geometry.is_empty
    ]

    if not geometry_list:
        return 0, 0.0

    try:
        repaired = [
            shapely.set_precision(
                shapely.make_valid(geometry),
                grid_size=1.0,
            )
            for geometry in geometry_list
        ]

        unioned = shapely.union_all(
            repaired,
            grid_size=1.0,
        )

    except (
        shapely.GEOSException,
        ValueError,
        TypeError,
    ):
        unioned = None

        for geometry in geometry_list:
            try:
                repaired_geometry = shapely.set_precision(
                    shapely.make_valid(geometry),
                    grid_size=1.0,
                )

                if unioned is None:
                    unioned = repaired_geometry
                else:
                    unioned = shapely.union(
                        unioned,
                        repaired_geometry,
                        grid_size=1.0,
                    )

            except (
                shapely.GEOSException,
                ValueError,
                TypeError,
            ):
                continue

    polygon_parts = list(
        iter_polygon_parts(unioned)
    )

    if not polygon_parts:
        return 0, 0.0

    event_count = len(polygon_parts)

    area_km2 = (
        sum(
            float(shapely.area(part))
            for part in polygon_parts
        )
        / 1_000_000.0
    )

    return event_count, area_km2


# =============================================================================
# Active-date expansion and daily union
# =============================================================================


def expand_active_dates(
    frame: gpd.GeoDataFrame,
    year: int,
) -> gpd.GeoDataFrame:
    if frame.empty:
        return frame

    frame = frame.copy()

    frame["start_date"] = (
        frame["issued_time"].dt.normalize()
    )

    frame["end_date"] = (
        frame["expired_time"].dt.normalize()
    )

    duration_days = (
        frame["end_date"] - frame["start_date"]
    ).dt.days + 1

    valid_duration = duration_days.between(
        1,
        730,
    )

    if (~valid_duration).any():
        logging.warning(
            "Dropping %s records with active intervals "
            "outside 1-730 days",
            f"{int((~valid_duration).sum()):,}",
        )

    frame = frame.loc[valid_duration].copy()
    duration_days = duration_days.loc[valid_duration]

    if frame.empty:
        return frame

    repeated_rows = np.repeat(
        np.arange(len(frame), dtype=np.int64),
        duration_days.to_numpy(dtype=np.int64),
    )

    offsets = np.concatenate(
        [
            np.arange(
                number_days,
                dtype="timedelta64[D]",
            )
            for number_days in duration_days.to_numpy(
                dtype=np.int64
            )
        ]
    )

    active_dates = (
        frame["start_date"]
        .to_numpy(dtype="datetime64[D]")[repeated_rows]
        + offsets
    )

    expanded = frame.iloc[
        repeated_rows
    ].copy()

    expanded["date"] = pd.to_datetime(
        active_dates
    )

    year_start = pd.Timestamp(f"{year}-01-01")
    year_end = pd.Timestamp(f"{year}-12-31")

    expanded = expanded.loc[
        expanded["date"].between(
            year_start,
            year_end,
        )
    ].copy()

    return expanded


def aggregate_daily_unions(
    expanded: gpd.GeoDataFrame,
) -> pd.DataFrame:
    if expanded.empty:
        return empty_annual_result()

    output_rows: list[dict[str, object]] = []

    grouped = expanded.groupby(
        [
            "date",
            "warning_group",
        ],
        sort=True,
        observed=True,
    )

    number_groups = grouped.ngroups

    for index, (
        (active_date, warning_group),
        daily_group,
    ) in enumerate(
        grouped,
        start=1,
    ):
        event_count, area_km2 = safe_union_metrics(
            daily_group.geometry
        )

        if event_count <= 0 or area_km2 <= 0:
            continue

        output_rows.append(
            {
                "date": active_date,
                "year": int(active_date.year),
                "warning_group": warning_group,
                "event_count": event_count,
                "area_km2": area_km2,
            }
        )

        if (
            index == 1
            or index % 100 == 0
            or index == number_groups
        ):
            logging.info(
                "Daily union progress: %s/%s",
                f"{index:,}",
                f"{number_groups:,}",
            )

    if not output_rows:
        return empty_annual_result()

    result = pd.DataFrame(output_rows)

    result["date"] = pd.to_datetime(
        result["date"]
    )

    result["year"] = result["year"].astype(
        "int16"
    )

    result["warning_group"] = result[
        "warning_group"
    ].astype("string")

    result["event_count"] = result[
        "event_count"
    ].astype("int64")

    result["area_km2"] = result[
        "area_km2"
    ].astype("float64")

    return result[ANNUAL_OUTPUT_COLUMNS]


# =============================================================================
# Annual stage
# =============================================================================


def process_annual_stage(
    paths: WorkflowPaths,
    year: int,
    overwrite: bool,
) -> None:
    if year not in YEARS:
        raise ValueError(
            f"Year must be between {START_YEAR} "
            f"and {END_YEAR}: {year}"
        )

    output_file = paths.annual_output(year)

    if (
        not overwrite
        and valid_annual_cache(output_file)
    ):
        logging.info(
            "Year %d is already complete: %s",
            year,
            output_file,
        )
        return

    shapefile = find_annual_shapefile(
        paths,
        year,
    )

    logging.info(
        "Reading year %d: %s",
        year,
        shapefile,
    )

    frame = gpd.read_file(
        shapefile,
        engine="pyogrio",
        columns=REQUIRED_COLUMNS,
        bbox=(
            DOMAIN_WEST,
            DOMAIN_SOUTH,
            DOMAIN_EAST,
            DOMAIN_NORTH,
        ),
        use_arrow=True,
    )

    frame.columns = [
        (
            column.upper()
            if column != "geometry"
            else column
        )
        for column in frame.columns
    ]

    missing_columns = set(
        REQUIRED_COLUMNS
    ).difference(frame.columns)

    if missing_columns:
        raise KeyError(
            f"Year {year} is missing columns: "
            f"{sorted(missing_columns)}"
        )

    frame["PHENOM"] = (
        frame["PHENOM"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    frame = frame.loc[
        frame["PHENOM"].isin(
            SELECTED_TYPES
        )
    ].copy()

    logging.info(
        "Year %d: %s selected source records",
        year,
        f"{len(frame):,}",
    )

    if frame.empty:
        write_parquet_atomic(
            empty_annual_result(),
            output_file,
        )
        return

    frame["warning_group"] = (
        frame["PHENOM"]
        .astype(str)
        .map(warning_group_from_code)
    )

    frame["issued_time"] = parse_iem_datetime(
        frame["ISSUED"]
    )

    frame["expired_time"] = parse_iem_datetime(
        frame["EXPIRED"]
    )

    valid_time = (
        frame["issued_time"].notna()
        & frame["expired_time"].notna()
        & (
            frame["expired_time"]
            >= frame["issued_time"]
        )
    )

    if (~valid_time).any():
        logging.warning(
            "Year %d: dropping %s records with invalid "
            "ISSUED or EXPIRED values",
            year,
            f"{int((~valid_time).sum()):,}",
        )

    frame = frame.loc[
        valid_time
    ].copy()

    frame = clip_warning_geometries(
        frame
    )

    logging.info(
        "Year %d: %s clipped polygons retained",
        year,
        f"{len(frame):,}",
    )

    expanded = expand_active_dates(
        frame,
        year,
    )

    logging.info(
        "Year %d: %s active-date polygon records",
        year,
        f"{len(expanded):,}",
    )

    daily = aggregate_daily_unions(
        expanded
    )

    write_parquet_atomic(
        daily,
        output_file,
    )

    logging.info(
        "Year %d: wrote %s daily group rows: %s",
        year,
        f"{len(daily):,}",
        output_file,
    )

    del frame
    del expanded
    del daily
    gc.collect()


# =============================================================================
# Combine annual results
# =============================================================================


def combine_annual_results(
    paths: WorkflowPaths,
) -> pd.DataFrame:
    missing_files = [
        paths.annual_output(year)
        for year in YEARS
        if not valid_annual_cache(
            paths.annual_output(year)
        )
    ]

    if missing_files:
        raise FileNotFoundError(
            "Annual outputs are missing or invalid:\n"
            + "\n".join(
                str(path)
                for path in missing_files
            )
        )

    annual_frames: list[pd.DataFrame] = []

    for year in YEARS:
        annual_file = paths.annual_output(
            year
        )

        logging.info(
            "Reading annual output: %s",
            annual_file,
        )

        annual = pd.read_parquet(
            annual_file
        )

        if not annual.empty:
            annual_frames.append(
                annual
            )

    if not annual_frames:
        raise RuntimeError(
            "All annual outputs are empty."
        )

    combined = pd.concat(
        annual_frames,
        ignore_index=True,
    )

    combined["date"] = pd.to_datetime(
        combined["date"]
    )

    combined = combined.loc[
        combined["date"].between(
            START_DATE,
            END_DATE,
        )
        & (combined["event_count"] > 0)
        & (combined["area_km2"] > 0)
    ].copy()

    combined = combined.sort_values(
        [
            "warning_group",
            "date",
        ]
    ).reset_index(drop=True)

    combined = combined[
        ANNUAL_OUTPUT_COLUMNS
    ]

    write_parquet_atomic(
        combined,
        paths.combined_parquet,
    )

    combined.to_csv(
        paths.combined_csv,
        index=False,
        date_format="%Y-%m-%d",
    )

    lookup = pd.DataFrame(
        [
            {
                "warning_group": warning_group,
                "phenom": code,
            }
            for warning_group, codes in GROUP_CODES.items()
            for code in codes
        ]
    )

    lookup.to_csv(
        paths.group_lookup_csv,
        index=False,
    )

    logging.info(
        "Wrote combined Parquet: %s",
        paths.combined_parquet,
    )

    logging.info(
        "Wrote combined CSV: %s",
        paths.combined_csv,
    )

    return combined


# =============================================================================
# Ar processing
# =============================================================================


def identify_column(
    frame: pd.DataFrame,
    candidates: list[str],
    description: str,
) -> str:
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate

    raise KeyError(
        f"Could not identify {description}. "
        f"Available columns: {list(frame.columns)}"
    )


def read_daily_ar(
    input_file: Path,
    variable_name: str,
) -> pd.DataFrame:
    if not input_file.is_file():
        raise FileNotFoundError(
            f"Missing {variable_name} file: {input_file}"
        )

    frame = pd.read_csv(
        input_file
    )

    date_column = identify_column(
        frame,
        [
            "date",
            "Date",
            "DATE",
            "day",
            "datetime",
        ],
        f"{variable_name} date column",
    )

    ar_column = identify_column(
        frame,
        [
            "Ar",
            "AR",
            "ar",
            "A_r",
        ],
        f"{variable_name} Ar column",
    )

    output = pd.DataFrame(
        {
            "date": pd.to_datetime(
                frame[date_column],
                errors="coerce",
            ),
            "Ar": pd.to_numeric(
                frame[ar_column],
                errors="coerce",
            ),
        }
    )

    output = output.loc[
        output["date"].notna()
        & np.isfinite(output["Ar"])
        & output["date"].between(
            START_DATE,
            END_DATE,
        )
    ].copy()

    output = (
        output.groupby(
            "date",
            as_index=False,
        )
        .agg(
            Ar=("Ar", "mean"),
        )
        .sort_values("date")
        .reset_index(drop=True)
    )

    return output


# =============================================================================
# Correlations and figures
# =============================================================================


def correlation_statistics(
    frame: pd.DataFrame,
    response: str,
) -> dict[str, float | int]:
    complete = frame.loc[
        np.isfinite(frame["Ar"])
        & np.isfinite(frame[response])
    ].copy()

    if len(complete) < 3:
        return {
            "n": len(complete),
            "r": np.nan,
            "p_value": np.nan,
        }

    result = pearsonr(
        complete["Ar"],
        complete[response],
    )

    return {
        "n": len(complete),
        "r": float(result.statistic),
        "p_value": float(result.pvalue),
    }


def format_p_value(value: float) -> str:
    if not np.isfinite(value):
        return "NA"

    if value < 0.001:
        return "< 0.001"

    return f"= {value:.3f}"


def annotation_text(
    statistics: dict[str, float | int],
) -> str:
    correlation = statistics["r"]

    correlation_text = (
        "NA"
        if not np.isfinite(correlation)
        else f"{correlation:.2f}"
    )

    return (
        f"r = {correlation_text}\n"
        f"p {format_p_value(statistics['p_value'])}\n"
        f"n = {statistics['n']:,}"
    )


def rounded_axis_upper(maximum: float) -> float:
    if not np.isfinite(maximum) or maximum <= 0:
        return 1.0

    exponent = math.floor(
        math.log10(maximum)
    )

    magnitude = 10**exponent
    normalized = maximum / magnitude

    for candidate in (
        1,
        2,
        2.5,
        5,
        10,
    ):
        if normalized <= candidate:
            return candidate * magnitude

    return 10 * magnitude


def add_regression(
    axis,
    x_values: np.ndarray,
    y_values: np.ndarray,
) -> None:
    valid = (
        np.isfinite(x_values)
        & np.isfinite(y_values)
    )

    x_values = x_values[valid]
    y_values = y_values[valid]

    if len(x_values) < 2:
        return

    coefficients = np.polyfit(
        x_values,
        y_values,
        deg=1,
    )

    x_line = np.linspace(
        x_values.min(),
        x_values.max(),
        200,
    )

    y_line = np.polyval(
        coefficients,
        x_line,
    )

    axis.plot(
        x_line,
        y_line,
        linewidth=0.8,
    )


def create_correlation_figure(
    paths: WorkflowPaths,
    warnings: pd.DataFrame,
    ar_data: pd.DataFrame,
    warning_group: str,
    x_label: str,
) -> list[dict[str, object]]:
    plot_data = (
        warnings.loc[
            warnings["warning_group"].eq(
                warning_group
            )
        ]
        .merge(
            ar_data,
            on="date",
            how="inner",
            validate="one_to_one",
        )
    )

    plot_data = plot_data.loc[
        (plot_data["event_count"] > 0)
        & (plot_data["area_km2"] > 0)
        & np.isfinite(plot_data["Ar"])
    ].copy()

    if plot_data.empty:
        raise RuntimeError(
            f"No overlapping nonzero warning and Ar dates "
            f"for {warning_group}"
        )

    plot_data = plot_data[
        CORRELATION_OUTPUT_COLUMNS
    ]

    plot_data.to_csv(
        paths.correlation_output(
            warning_group
        ),
        index=False,
        date_format="%Y-%m-%d",
    )

    count_statistics = correlation_statistics(
        plot_data,
        "event_count",
    )

    area_statistics = correlation_statistics(
        plot_data,
        "area_km2",
    )

    figure, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(6.5, 3.4),
        constrained_layout=True,
    )

    panel_definitions = [
        (
            axes[0],
            "event_count",
            "Spatially distinct active warnings per day",
            count_statistics,
        ),
        (
            axes[1],
            "area_km2",
            "Unioned active warning area per day (km²)",
            area_statistics,
        ),
    ]

    for panel_index, (
        axis,
        response,
        y_label,
        statistics,
    ) in enumerate(
        panel_definitions
    ):
        axis.scatter(
            plot_data["Ar"],
            plot_data[response],
            s=8,
            alpha=0.25,
            linewidths=0,
        )

        add_regression(
            axis,
            plot_data["Ar"].to_numpy(
                dtype=float
            ),
            plot_data[response].to_numpy(
                dtype=float
            ),
        )

        axis.text(
            0.03,
            0.97,
            annotation_text(statistics),
            transform=axis.transAxes,
            horizontalalignment="left",
            verticalalignment="top",
            fontsize=7,
            bbox={
                "facecolor": "white",
                "edgecolor": "black",
                "linewidth": 0.4,
                "pad": 3,
            },
        )

        axis.text(
            0.01,
            0.99,
            f"{chr(97 + panel_index)})",
            transform=axis.transAxes,
            horizontalalignment="left",
            verticalalignment="top",
            fontsize=9,
            fontweight="bold",
        )

        axis.set_xlabel(
            x_label,
            fontsize=8,
        )

        axis.set_ylabel(
            y_label,
            fontsize=8,
        )

        axis.set_ylim(
            0,
            rounded_axis_upper(
                float(
                    plot_data[response].max()
                )
            ),
        )

        axis.grid(
            axis="both",
            linewidth=0.3,
            alpha=0.3,
        )

        axis.tick_params(
            axis="both",
            labelsize=7,
            width=0.6,
            length=3,
        )

        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    figure.suptitle(
        (
            f"{GROUP_TITLES[warning_group]}\n"
            "Daily non-overlapping warning geometries, 2002–2023"
        ),
        fontsize=10,
        fontweight="bold",
    )

    output_figure = paths.figure_output(
        warning_group
    )

    figure.savefig(
        output_figure,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(figure)

    logging.info(
        "Wrote figure: %s",
        output_figure,
    )

    return [
        {
            "warning_group": warning_group,
            "response": "event_count",
            **count_statistics,
        },
        {
            "warning_group": warning_group,
            "response": "area_km2",
            **area_statistics,
        },
    ]


def process_combine_stage(
    paths: WorkflowPaths,
) -> None:
    warnings = combine_annual_results(
        paths
    )

    heat_ar = read_daily_ar(
        paths.heat_ar_file,
        "heat-index",
    )

    precipitation_ar = read_daily_ar(
        paths.precipitation_ar_file,
        "precipitation",
    )

    summary_rows = []

    summary_rows.extend(
        create_correlation_figure(
            paths=paths,
            warnings=warnings,
            ar_data=heat_ar,
            warning_group="heat",
            x_label="Heat-index Ar",
        )
    )

    summary_rows.extend(
        create_correlation_figure(
            paths=paths,
            warnings=warnings,
            ar_data=precipitation_ar,
            warning_group="precipitation",
            x_label="Precipitation Ar",
        )
    )

    pd.DataFrame(
        summary_rows
    ).to_csv(
        paths.correlation_summary_csv,
        index=False,
    )

    logging.info(
        "Combine and plotting stage complete."
    )


# =============================================================================
# Command line
# =============================================================================


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Process selected IEM warning types, union daily geometries, "
            "and compare warning activity with Ar."
        )
    )

    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="CISC project root.",
    )

    parser.add_argument(
        "--stage",
        required=True,
        choices=[
            "annual",
            "combine",
        ],
        help="Workflow stage.",
    )

    parser.add_argument(
        "--year",
        type=int,
        default=None,
        help="Year for the annual stage.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing valid annual output.",
    )

    return parser.parse_args()


def main() -> int:
    configure_logging()
    arguments = parse_arguments()

    paths = WorkflowPaths(
        arguments.root
        .expanduser()
        .resolve()
    )

    paths.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    paths.annual_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    paths.correlation_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    paths.figure_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    logging.info(
        "CISC root: %s",
        paths.root,
    )

    logging.info(
        "Stage: %s",
        arguments.stage,
    )

    if arguments.stage == "annual":
        if arguments.year is None:
            raise ValueError(
                "--year is required for --stage annual"
            )

        process_annual_stage(
            paths=paths,
            year=arguments.year,
            overwrite=arguments.overwrite,
        )

    elif arguments.stage == "combine":
        process_combine_stage(
            paths
        )

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        logging.exception(
            "Mesonet workflow failed"
        )
        sys.exit(1)
