#!/usr/bin/env python3
"""Write one single-layer maximum-return-period GeoTIFF and PNG per 2023 cluster.

The workflow reads yearly cluster Parquets only from:

    data/out/06_clusters/filtered_rp25/<dataset_tag>/

For every available heat-index and precipitation dataset, it:

1. infers the exact grid-center lattice from the Parquet x/y coordinates;
2. builds a full Southeast grid covering the reference raster domain;
3. assigns each source coordinate to its nearest grid-cell center;
4. collapses repeated space-time rows using maximum rp_max_years per cell;
5. writes one GeoTIFF and one title-free PNG for each cluster_id.

Outputs
-------
figures/cluster_tifs_pngs/
    geotiffs/<dataset_tag>/<cluster_id>.tif
    images/<dataset_tag>/<cluster_id>.png
    reference_grids/<dataset_tag>.tif
    cluster_geotiff_png_manifest_2023.csv
"""

from __future__ import annotations

import argparse
import logging
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import rasterio
from rasterio.transform import Affine
from rasterio.warp import Resampling, reproject


LOG = logging.getLogger("cisc08")

DATASET_RE = re.compile(
    r"^(?P<variable>heatindex|precipitation)_"
    r"(?P<resolution>\d+(?:\.\d+)?)"
    r"(?:_(?P<method>bilinear|max))?$"
)

RP_BREAKS = np.array(
    [5, 10, 25, 50, 100, 200, 500, 1000],
    dtype=float,
)
RP_BOUNDARIES = np.array(
    [5, 10, 25, 50, 100, 200, 500, 1000, 1000.0001],
    dtype=float,
)
RP_COLORS = [
    "#d9d9d9",
    "#fff5f0",
    "#fee0d2",
    "#fcbba1",
    "#fc9272",
    "#fb6a4a",
    "#ef3b2c",
    "#cb181d",
    "#99000d",
]


@dataclass(frozen=True)
class DatasetSpec:
    tag: str
    variable: str
    nominal_resolution: float
    method: str
    parquet_path: Path


@dataclass(frozen=True)
class GridSpec:
    width: int
    height: int
    resolution_x: float
    resolution_y: float
    x_center0: float
    y_center0: float
    transform: Affine
    crs_wkt: str
    nodata: float


@dataclass(frozen=True)
class ClusterJob:
    dataset: DatasetSpec
    cluster_id: str
    frame: pd.DataFrame
    grid: GridSpec
    geotiff_dir: Path
    image_dir: Path
    overwrite: bool
    state_shapefile: str | None


def parse_args() -> argparse.Namespace:
    root_default = os.environ.get(
        "CISC_ROOT",
        "/proj/sebastianlab/users/hquintal/cisc",
    )

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(root_default))
    parser.add_argument("--year", type=int, default=2016)
    parser.add_argument(
        "--reference-raster",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", "10")),
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--datasets", nargs="*", default=None)
    parser.add_argument(
        "--state-shapefile",
        default=os.environ.get("CISC_STATE_SHP"),
    )
    return parser.parse_args()


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def sanitize_cluster_id(value: object) -> str:
    text = str(value).strip()
    if not text or text.lower() in {"nan", "<na>"}:
        raise ValueError("Blank cluster_id encountered")
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text)


def discover_year_parquet(folder: Path, year: int) -> Path | None:
    preferred = folder / f"cluster_rp5_{year}.parquet"
    if preferred.exists():
        return preferred

    candidates = sorted(folder.glob(f"*{year}*.parquet"))
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]

    exact = [
        path for path in candidates
        if re.search(rf"(?:^|_){year}(?:_|$)", path.stem)
    ]
    if len(exact) == 1:
        return exact[0]

    raise RuntimeError(
        f"Multiple candidate {year} Parquets in {folder}: "
        + ", ".join(str(path) for path in candidates)
    )


def discover_datasets(
    cluster_root: Path,
    year: int,
    requested: list[str] | None,
) -> list[DatasetSpec]:
    if not cluster_root.exists():
        raise FileNotFoundError(
            f"filtered_rp25 directory does not exist: {cluster_root}"
        )

    requested_set = set(requested or [])
    specs: list[DatasetSpec] = []

    for folder in sorted(path for path in cluster_root.iterdir() if path.is_dir()):
        match = DATASET_RE.fullmatch(folder.name)
        if match is None:
            continue
        if requested_set and folder.name not in requested_set:
            continue

        parquet_path = discover_year_parquet(folder, year)
        if parquet_path is None:
            LOG.warning("No %d Parquet found in %s", year, folder)
            continue

        specs.append(
            DatasetSpec(
                tag=folder.name,
                variable=match.group("variable"),
                nominal_resolution=float(match.group("resolution")),
                method=match.group("method") or "native",
                parquet_path=parquet_path,
            )
        )

    if requested_set:
        found = {spec.tag for spec in specs}
        missing = sorted(requested_set - found)
        if missing:
            raise FileNotFoundError(
                f"Requested datasets were not found for {year}: {missing}"
            )

    if not specs:
        raise FileNotFoundError(
            f"No valid {year} heatindex/precipitation datasets found under "
            f"{cluster_root}"
        )

    return specs


def read_reference(
    path: Path,
) -> tuple[rasterio.coords.BoundingBox, str, Affine]:
    if not path.exists():
        raise FileNotFoundError(f"Reference raster does not exist: {path}")

    with rasterio.open(path) as src:
        if src.crs is None:
            raise ValueError(f"Reference raster has no CRS: {path}")
        return src.bounds, src.crs.to_wkt(), src.transform


def robust_spacing(values: np.ndarray, fallback: float) -> float:
    """Infer lattice spacing from unique coordinate centers."""
    values = values[np.isfinite(values)]
    if values.size < 2:
        return fallback

    unique = np.unique(np.round(values, 10))
    if unique.size < 2:
        return fallback

    diffs = np.diff(unique)
    diffs = diffs[diffs > 1e-8]
    if diffs.size == 0:
        return fallback

    # The lower quantiles are robust when some intermediate lattice values
    # are absent from the filtered event table.
    estimate = float(np.median(diffs[diffs <= np.quantile(diffs, 0.35)]))
    if not np.isfinite(estimate) or estimate <= 0:
        return fallback
    return estimate


def lattice_phase(values: np.ndarray, spacing: float) -> float:
    """Estimate the coordinate-center phase modulo spacing."""
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("Cannot infer grid phase without finite coordinates")

    phases = np.mod(values, spacing)
    angles = phases / spacing * 2.0 * np.pi
    mean_angle = math.atan2(
        float(np.mean(np.sin(angles))),
        float(np.mean(np.cos(angles))),
    )
    if mean_angle < 0:
        mean_angle += 2.0 * np.pi
    return mean_angle / (2.0 * np.pi) * spacing


def first_center_at_or_below(bound: float, phase: float, spacing: float) -> float:
    return phase + math.floor((bound - phase) / spacing) * spacing


def first_center_at_or_above(bound: float, phase: float, spacing: float) -> float:
    return phase + math.ceil((bound - phase) / spacing) * spacing


def build_dataset_grid(
    bounds: rasterio.coords.BoundingBox,
    crs_wkt: str,
    frame: pd.DataFrame,
    nominal_resolution: float,
) -> GridSpec:
    x = pd.to_numeric(frame["x"], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(frame["y"], errors="coerce").to_numpy(dtype=float)

    dx = robust_spacing(x, nominal_resolution)
    dy = robust_spacing(y, nominal_resolution)

    # Dataset tags are rounded, so keep the measured lattice spacing when it is
    # plausible; otherwise fall back to the declared resolution.
    if not (0.5 * nominal_resolution <= dx <= 1.5 * nominal_resolution):
        LOG.warning(
            "Measured x spacing %.10f is inconsistent with nominal %.10f; "
            "using nominal resolution",
            dx,
            nominal_resolution,
        )
        dx = nominal_resolution
    if not (0.5 * nominal_resolution <= dy <= 1.5 * nominal_resolution):
        LOG.warning(
            "Measured y spacing %.10f is inconsistent with nominal %.10f; "
            "using nominal resolution",
            dy,
            nominal_resolution,
        )
        dy = nominal_resolution

    x_phase = lattice_phase(x, dx)
    y_phase = lattice_phase(y, dy)

    left_center = first_center_at_or_below(bounds.left, x_phase, dx)
    right_center = first_center_at_or_above(bounds.right, x_phase, dx)
    bottom_center = first_center_at_or_below(bounds.bottom, y_phase, dy)
    top_center = first_center_at_or_above(bounds.top, y_phase, dy)

    width = int(round((right_center - left_center) / dx)) + 1
    height = int(round((top_center - bottom_center) / dy)) + 1

    transform = Affine(
        dx,
        0.0,
        left_center - dx / 2.0,
        0.0,
        -dy,
        top_center + dy / 2.0,
    )

    return GridSpec(
        width=width,
        height=height,
        resolution_x=dx,
        resolution_y=dy,
        x_center0=left_center,
        y_center0=top_center,
        transform=transform,
        crs_wkt=crs_wkt,
        nodata=-9999.0,
    )


def map_coordinates_to_nearest_cells(
    x: np.ndarray,
    y: np.ndarray,
    grid: GridSpec,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Map every finite point to its nearest target-grid cell center."""
    finite = np.isfinite(x) & np.isfinite(y)

    rows = np.full(x.shape, -1, dtype=np.int64)
    cols = np.full(x.shape, -1, dtype=np.int64)
    x_error = np.full(x.shape, np.nan, dtype=float)
    y_error = np.full(x.shape, np.nan, dtype=float)

    cols[finite] = np.rint(
        (x[finite] - grid.x_center0) / grid.resolution_x
    ).astype(np.int64)
    rows[finite] = np.rint(
        (grid.y_center0 - y[finite]) / grid.resolution_y
    ).astype(np.int64)

    mapped_x = grid.x_center0 + cols[finite] * grid.resolution_x
    mapped_y = grid.y_center0 - rows[finite] * grid.resolution_y
    x_error[finite] = np.abs(x[finite] - mapped_x)
    y_error[finite] = np.abs(y[finite] - mapped_y)

    valid = (
        finite
        & (rows >= 0)
        & (rows < grid.height)
        & (cols >= 0)
        & (cols < grid.width)
    )

    return rows, cols, valid, x_error, y_error


def aggregate_cluster(
    frame: pd.DataFrame,
    grid: GridSpec,
) -> tuple[np.ndarray, int, int, float, float]:
    work = frame.loc[:, ["x", "y", "rp_max_years"]].copy()
    work["x"] = pd.to_numeric(work["x"], errors="coerce")
    work["y"] = pd.to_numeric(work["y"], errors="coerce")
    work["rp_max_years"] = pd.to_numeric(
        work["rp_max_years"],
        errors="coerce",
    )

    x = work["x"].to_numpy(dtype=float)
    y = work["y"].to_numpy(dtype=float)
    rp = work["rp_max_years"].to_numpy(dtype=float)

    rows, cols, valid_xy, x_error, y_error = map_coordinates_to_nearest_cells(
        x,
        y,
        grid,
    )
    valid = valid_xy & np.isfinite(rp)

    n_rejected = int((~valid).sum())
    if not valid.any():
        raise ValueError("No finite cluster rows map inside the dataset grid")

    mapped = pd.DataFrame(
        {
            "row": rows[valid],
            "col": cols[valid],
            "rp_max_years": rp[valid],
        }
    )

    grouped = (
        mapped.groupby(["row", "col"], sort=False, as_index=False)
        ["rp_max_years"]
        .max()
    )

    arr = np.full((grid.height, grid.width), np.nan, dtype=np.float32)
    arr[
        grouped["row"].to_numpy(dtype=np.int64),
        grouped["col"].to_numpy(dtype=np.int64),
    ] = grouped["rp_max_years"].to_numpy(dtype=np.float32)

    max_x_error = float(np.nanmax(x_error[valid])) if valid.any() else np.nan
    max_y_error = float(np.nanmax(y_error[valid])) if valid.any() else np.nan

    return (
        arr,
        int(len(grouped)),
        n_rejected,
        max_x_error,
        max_y_error,
    )


def raster_profile(grid: GridSpec) -> dict[str, object]:
    return {
        "driver": "GTiff",
        "height": grid.height,
        "width": grid.width,
        "count": 1,
        "dtype": "float32",
        "crs": grid.crs_wkt,
        "transform": grid.transform,
        "nodata": grid.nodata,
        "compress": "deflate",
        "predictor": 2,
        "zlevel": 4,
        "tiled": True,
        "BIGTIFF": "IF_SAFER",
    }


def write_reference_grid(
    path: Path,
    grid: GridSpec,
    reference_path: Path,
    method: str,
    overwrite: bool,
) -> None:
    if path.exists() and not overwrite:
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    resampling = {
        "native": Resampling.nearest,
        "bilinear": Resampling.bilinear,
        "max": Resampling.max,
    }[method]

    with rasterio.open(reference_path) as src:
        source = src.read(1)
        destination = np.full(
            (grid.height, grid.width),
            grid.nodata,
            dtype=np.float32,
        )

        reproject(
            source=source,
            destination=destination,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=src.nodata,
            dst_transform=grid.transform,
            dst_crs=grid.crs_wkt,
            dst_nodata=grid.nodata,
            resampling=resampling,
        )

    with rasterio.open(path, "w", **raster_profile(grid)) as dst:
        dst.write(destination, 1)
        dst.set_band_description(1, "resampled_reference_band_1")
        dst.update_tags(
            source_reference=str(reference_path),
            resample_method=method,
            resolution_x_deg=f"{grid.resolution_x:.12g}",
            resolution_y_deg=f"{grid.resolution_y:.12g}",
        )


def write_cluster_tif(
    path: Path,
    arr: np.ndarray,
    grid: GridSpec,
    tags: dict[str, str],
    overwrite: bool,
) -> None:
    if path.exists() and not overwrite:
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    output = np.where(np.isfinite(arr), arr, grid.nodata).astype(np.float32)

    with rasterio.open(path, "w", **raster_profile(grid)) as dst:
        dst.write(output, 1)
        dst.set_band_description(1, "max_rp_years")
        dst.update_tags(**tags)


def load_state_boundaries(path: str | None, target_crs: str):
    if not path:
        return None

    try:
        import geopandas as gpd

        states = gpd.read_file(path)
        if states.crs is None:
            raise ValueError("State-boundary file has no CRS")
        return states.to_crs(target_crs)
    except Exception as exc:
        LOG.warning("Could not read state boundaries %s: %s", path, exc)
        return None


def write_cluster_png(
    path: Path,
    arr: np.ndarray,
    grid: GridSpec,
    state_shapefile: str | None,
    overwrite: bool,
) -> None:
    if path.exists() and not overwrite:
        return

    path.parent.mkdir(parents=True, exist_ok=True)

    cmap = ListedColormap(RP_COLORS)
    cmap.set_bad("white")
    norm = BoundaryNorm(RP_BOUNDARIES, cmap.N, clip=True)
    masked = np.ma.masked_invalid(arr)

    left = grid.transform.c
    right = left + grid.width * grid.resolution_x
    top = grid.transform.f
    bottom = top - grid.height * grid.resolution_y

    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    image = ax.imshow(
        masked,
        extent=(left, right, bottom, top),
        origin="upper",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
    )

    states = load_state_boundaries(state_shapefile, grid.crs_wkt)
    if states is not None:
        states.boundary.plot(
            ax=ax,
            color="black",
            linewidth=0.45,
            zorder=3,
        )

    # Intentionally no title, subtitle, annotation, or suptitle.
    ax.set_title("")
    fig.suptitle("")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("equal", adjustable="box")

    colorbar = fig.colorbar(
        image,
        ax=ax,
        orientation="horizontal",
        pad=0.10,
        fraction=0.06,
        ticks=RP_BREAKS,
        boundaries=RP_BOUNDARIES,
        spacing="proportional",
    )
    colorbar.set_label("Maximum return period (yr)")
    colorbar.ax.set_xticklabels([str(int(value)) for value in RP_BREAKS])
    colorbar.outline.set_linewidth(0.6)

    fig.savefig(
        path,
        dpi=220,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def read_dataset_frame(spec: DatasetSpec, year: int) -> pd.DataFrame:
    required = {
        "x",
        "y",
        "time_utc",
        "rp_max_years",
        "cluster_id",
    }

    schema_names = set(pq.read_schema(spec.parquet_path).names)
    missing = sorted(required - schema_names)
    if missing:
        raise ValueError(
            f"{spec.parquet_path} is missing required columns: {missing}"
        )

    optional = {
        "resolution_deg",
        "resample_method",
        "variable",
        "dataset_tag",
        "cluster_max_rp_years",
        "event_filter_rp_years",
    }
    columns = sorted(required | (optional & schema_names))

    frame = pq.read_table(
        spec.parquet_path,
        columns=columns,
    ).to_pandas()

    frame["time_utc"] = pd.to_datetime(
        frame["time_utc"],
        utc=True,
        errors="coerce",
    )
    frame = frame.loc[frame["time_utc"].dt.year.eq(year)].copy()

    frame["cluster_id"] = (
        frame["cluster_id"]
        .astype("string")
        .str.strip()
    )
    frame = frame.loc[
        frame["cluster_id"].notna()
        & frame["cluster_id"].ne("")
    ].copy()

    return frame


def process_cluster(job: ClusterJob) -> dict[str, object]:
    dataset = job.dataset
    cluster_id = sanitize_cluster_id(job.cluster_id)

    frame = job.frame.copy()
    frame["cluster_id"] = (
        frame["cluster_id"]
        .astype("string")
        .str.strip()
    )
    frame = frame.loc[frame["cluster_id"].eq(job.cluster_id)].copy()

    unique_ids = frame["cluster_id"].dropna().unique()
    if len(unique_ids) != 1 or str(unique_ids[0]) != job.cluster_id:
        raise ValueError(
            f"Expected only cluster_id={job.cluster_id}; "
            f"found {unique_ids.tolist()}"
        )

    tif_path = job.geotiff_dir / dataset.tag / f"{cluster_id}.tif"
    png_path = job.image_dir / dataset.tag / f"{cluster_id}.png"

    times = pd.to_datetime(frame["time_utc"], utc=True, errors="coerce")
    start = times.min()
    end = times.max()

    (
        arr,
        n_cells,
        n_rejected,
        max_x_error,
        max_y_error,
    ) = aggregate_cluster(frame, job.grid)

    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        raise ValueError("Cluster raster contains no populated cells")

    cluster_max = float(np.max(finite))

    tags = {
        "cluster_id": cluster_id,
        "dataset_tag": dataset.tag,
        "variable": dataset.variable,
        "nominal_resolution_deg": f"{dataset.nominal_resolution:g}",
        "actual_resolution_x_deg": f"{job.grid.resolution_x:.12g}",
        "actual_resolution_y_deg": f"{job.grid.resolution_y:.12g}",
        "resample_method": dataset.method,
        "start_utc": start.isoformat() if pd.notna(start) else "",
        "end_utc": end.isoformat() if pd.notna(end) else "",
        "aggregation": (
            "maximum rp_max_years per nearest target-grid cell "
            "across cluster time steps"
        ),
    }

    write_cluster_tif(
        tif_path,
        arr,
        job.grid,
        tags,
        job.overwrite,
    )
    write_cluster_png(
        png_path,
        arr,
        job.grid,
        job.state_shapefile,
        job.overwrite,
    )

    return {
        "dataset_tag": dataset.tag,
        "variable": dataset.variable,
        "nominal_resolution_deg": dataset.nominal_resolution,
        "actual_resolution_x_deg": job.grid.resolution_x,
        "actual_resolution_y_deg": job.grid.resolution_y,
        "resample_method": dataset.method,
        "source_parquet": str(dataset.parquet_path),
        "cluster_id": cluster_id,
        "start_utc": start.isoformat() if pd.notna(start) else None,
        "end_utc": end.isoformat() if pd.notna(end) else None,
        "n_source_rows": int(len(frame)),
        "n_raster_cells": n_cells,
        "n_rejected_rows": n_rejected,
        "max_x_mapping_error_deg": max_x_error,
        "max_y_mapping_error_deg": max_y_error,
        "cluster_max_rp_years": cluster_max,
        "geotiff": str(tif_path),
        "png": str(png_path),
        "status": "ok",
    }


def main() -> int:
    args = parse_args()
    setup_logging()

    root = args.root.resolve()
    reference = (
        args.reference_raster
        or root / "data/out/04_rp/heatindex_idf_hi1_return_periods_bbox.nc"
    )

    # This workflow is intentionally fixed to filtered_rp25.
    cluster_root = root / "data/out/06_clusters/filtered_rp25"

    output_root = root / "figures/cluster_tifs_pngs"
    geotiff_dir = output_root / "geotiffs"
    image_dir = output_root / "images"
    reference_grid_dir = output_root / "reference_grids"
    manifest_path = (
        output_root
        / f"cluster_geotiff_png_manifest_{args.year}.csv"
    )
    output_root.mkdir(parents=True, exist_ok=True)

    bounds, crs_wkt, reference_transform = read_reference(reference)
    LOG.info("Reference raster: %s", reference)
    LOG.info("Reference bounds: %s", bounds)
    LOG.info(
        "Reference resolution: x=%.10f, y=%.10f",
        abs(reference_transform.a),
        abs(reference_transform.e),
    )
    LOG.info("Cluster root: %s", cluster_root)

    datasets = discover_datasets(
        cluster_root,
        args.year,
        args.datasets,
    )
    LOG.info(
        "Discovered %d datasets: %s",
        len(datasets),
        ", ".join(spec.tag for spec in datasets),
    )

    jobs: list[ClusterJob] = []
    failures: list[dict[str, object]] = []

    for spec in datasets:
        LOG.info("Reading %s", spec.parquet_path)

        try:
            frame = read_dataset_frame(spec, args.year)
            if frame.empty:
                LOG.warning(
                    "No %d cluster rows in %s",
                    args.year,
                    spec.parquet_path,
                )
                continue

            grid = build_dataset_grid(
                bounds=bounds,
                crs_wkt=crs_wkt,
                frame=frame,
                nominal_resolution=spec.nominal_resolution,
            )

            LOG.info(
                "%s grid: %d columns x %d rows; dx=%.10f; dy=%.10f",
                spec.tag,
                grid.width,
                grid.height,
                grid.resolution_x,
                grid.resolution_y,
            )

            write_reference_grid(
                reference_grid_dir / f"{spec.tag}.tif",
                grid,
                reference,
                spec.method,
                args.overwrite,
            )

            cluster_ids = sorted(
                str(value)
                for value in frame["cluster_id"].dropna().unique()
            )
            LOG.info("%s: %d clusters", spec.tag, len(cluster_ids))

            for cluster_id in cluster_ids:
                group = frame.loc[
                    frame["cluster_id"].eq(cluster_id)
                ].copy()

                jobs.append(
                    ClusterJob(
                        dataset=spec,
                        cluster_id=cluster_id,
                        frame=group,
                        grid=grid,
                        geotiff_dir=geotiff_dir,
                        image_dir=image_dir,
                        overwrite=args.overwrite,
                        state_shapefile=args.state_shapefile,
                    )
                )

        except Exception as exc:
            LOG.exception("Dataset failed: %s", spec.tag)
            failures.append(
                {
                    "dataset_tag": spec.tag,
                    "source_parquet": str(spec.parquet_path),
                    "cluster_id": None,
                    "status": f"dataset error: {exc}",
                }
            )

    worker_count = max(1, min(args.workers, len(jobs) or 1))
    LOG.info(
        "Processing %d clusters with %d workers",
        len(jobs),
        worker_count,
    )

    results: list[dict[str, object]] = []

    if worker_count == 1:
        for index, job in enumerate(jobs, start=1):
            try:
                results.append(process_cluster(job))
            except Exception as exc:
                LOG.exception(
                    "Cluster failed: %s / %s",
                    job.dataset.tag,
                    job.cluster_id,
                )
                results.append(
                    {
                        "dataset_tag": job.dataset.tag,
                        "variable": job.dataset.variable,
                        "nominal_resolution_deg": (
                            job.dataset.nominal_resolution
                        ),
                        "resample_method": job.dataset.method,
                        "source_parquet": str(
                            job.dataset.parquet_path
                        ),
                        "cluster_id": job.cluster_id,
                        "status": f"error: {exc}",
                    }
                )

            if index % 25 == 0 or index == len(jobs):
                LOG.info(
                    "Completed %d/%d clusters",
                    index,
                    len(jobs),
                )

    else:
        with ProcessPoolExecutor(max_workers=worker_count) as pool:
            future_map = {
                pool.submit(process_cluster, job): job
                for job in jobs
            }

            for index, future in enumerate(
                as_completed(future_map),
                start=1,
            ):
                job = future_map[future]

                try:
                    results.append(future.result())
                except Exception as exc:
                    LOG.exception(
                        "Cluster failed: %s / %s",
                        job.dataset.tag,
                        job.cluster_id,
                    )
                    results.append(
                        {
                            "dataset_tag": job.dataset.tag,
                            "variable": job.dataset.variable,
                            "nominal_resolution_deg": (
                                job.dataset.nominal_resolution
                            ),
                            "resample_method": job.dataset.method,
                            "source_parquet": str(
                                job.dataset.parquet_path
                            ),
                            "cluster_id": job.cluster_id,
                            "status": f"error: {exc}",
                        }
                    )

                if (
                    index % 25 == 0
                    or index == len(future_map)
                ):
                    LOG.info(
                        "Completed %d/%d clusters",
                        index,
                        len(future_map),
                    )

    manifest = pd.DataFrame(results + failures)

    if not manifest.empty:
        sort_columns = [
            column
            for column in [
                "variable",
                "nominal_resolution_deg",
                "resample_method",
                "cluster_id",
            ]
            if column in manifest.columns
        ]
        manifest = manifest.sort_values(
            sort_columns,
            na_position="last",
        )

    manifest.to_csv(manifest_path, index=False)

    status = manifest.get(
        "status",
        pd.Series(dtype="string"),
    )
    n_ok = int(status.eq("ok").sum())
    n_error = int(len(manifest) - n_ok)

    LOG.info("Wrote manifest: %s", manifest_path)
    LOG.info(
        "Complete: %d successful cluster outputs; %d errors",
        n_ok,
        n_error,
    )

    return 1 if n_error else 0


if __name__ == "__main__":
    sys.exit(main())
