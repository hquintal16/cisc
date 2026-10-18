#!/usr/bin/env python3
"""Cluster RP >= 5-year exceedances once, then filter whole clusters by peak RP.

This CISC workflow processes heat-index HI1 and precipitation PR24 exceedance
files by dataset and year. DBSCAN is run only on points with rp_max_years >= 5.
Each positive cluster is then assigned its maximum return period. Four nested
event products retain the complete RP >= 5 footprint of clusters whose maximum
return period reaches 5, 10, 25, or 50 years.

Modes
-----
build_manifest
    Discover annual exceedance files and write one manifest row per
    variable/dataset/year.
cluster_worker
    Process the strided subset assigned to the current SLURM array worker.
summarize
    Combine task JSON records into global and dataset summaries.

Input layouts
-------------
<data_out>/05_exceedances/heatindex/<dataset_suffix>/
    hi1_exceedances_maxrp_<year>.qs|parquet
<data_out>/05_exceedances/precipitation/pr24/<dataset_suffix>/
    pr24_exceedances_maxrp_<year>.qs|parquet

Output layout
-------------
<data_out>/11_clusters/
    clustered_rp5/<dataset_tag>/cluster_rp5_<year>.parquet
    filtered_rp5/<dataset_tag>/cluster_rp5_<year>.parquet
    filtered_rp10/<dataset_tag>/cluster_rp10_<year>.parquet
    filtered_rp25/<dataset_tag>/cluster_rp25_<year>.parquet
    _control/...

QS inputs are converted to temporary Parquet files with Rscript using the R
packages qs and arrow. Parquet inputs are read directly.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN


# ── CONFIG ────────────────────────────────────────────────────────────────

ROOT = Path(
    os.environ.get(
        "CISC_ROOT",
        "/proj/sebastianlab/users/hquintal/cisc",
    )
)

DATA_OUT = Path(
    os.environ.get(
        "CISC_DATA_OUT",
        str(ROOT / "data/out"),
    )
)

IN_ROOT = Path(
    os.environ.get(
        "CISC_EXCEEDANCE_DIR",
        str(DATA_OUT / "05_exceedances"),
    )
)

OUT_ROOT = Path(
    os.environ.get(
        "CISC_CLUSTER_DIR",
        str(DATA_OUT / "11_clusters"),
    )
)

CONTROL_DIR = OUT_ROOT / "_control"
RESULT_DIR = CONTROL_DIR / "task_results"
TEMP_DIR = Path(
    os.environ.get(
        "CISC_CLUSTER_TEMP_DIR",
        str(DATA_OUT / "tmp_python/11_cluster_rp5_filter"),
    )
)

MANIFEST_PATH = Path(
    os.environ.get(
        "CISC_CLUSTER_MANIFEST",
        str(CONTROL_DIR / "cluster_manifest.csv"),
    )
)

MODE = os.environ.get(
    "CISC_CLUSTER_MODE",
    "build_manifest",
).strip().lower()

WORKER_INDEX = int(
    os.environ.get(
        "SLURM_ARRAY_TASK_ID",
        os.environ.get("CISC_CLUSTER_WORKER_INDEX", "0"),
    )
)

WORKER_COUNT = int(
    os.environ.get(
        "CISC_CLUSTER_WORKER_COUNT",
        "4",
    )
)

EPS = float(os.environ.get("CISC_DBSCAN_EPS", "4.0"))
MIN_SAMPLES = int(os.environ.get("CISC_DBSCAN_MIN_SAMPLES", "4"))
TIME_STEP_HOURS = float(os.environ.get("CISC_CLUSTER_TIME_STEP_HOURS", "24"))
CLUSTER_FLOOR = float(os.environ.get("CISC_CLUSTER_FLOOR_YEARS", "5"))
FILTER_THRESHOLDS = tuple(
    sorted(
        {
            float(value.strip())
            for value in os.environ.get(
                "CISC_CLUSTER_FILTER_RPS",
                "5,10,25,50",
            ).split(",")
            if value.strip()
        }
    )
)

OVERWRITE = os.environ.get(
    "CISC_OVERWRITE",
    "false",
).strip().lower() in {"true", "1", "yes"}

RSCRIPT = os.environ.get("CISC_RSCRIPT", "Rscript")

YEAR_RE = re.compile(r"((?:19|20)\d{2})")

INPUT_SPECS = (
    {
        "variable": "heatindex",
        "base_dir": IN_ROOT / "heatindex",
        "prefix": "hi1_exceedances_maxrp_",
        "dataset_prefix": "heatindex_",
    },
    {
        "variable": "precipitation",
        "base_dir": IN_ROOT / "precipitation/pr24",
        "prefix": "pr24_exceedances_maxrp_",
        "dataset_prefix": "precipitation_",
    },
)


# ── GENERAL HELPERS ───────────────────────────────────────────────────────


def log(message: str) -> None:
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {message}", flush=True)


def threshold_label(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else str(value).replace(".", "p")


def atomic_write_csv(frame: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_suffix(".tmp.csv")
    frame.to_csv(tmp_path, index=False)
    os.replace(tmp_path, output_path)


def atomic_write_json(payload: dict[str, Any], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_suffix(".tmp.json")
    tmp_path.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    os.replace(tmp_path, output_path)


def atomic_write_parquet(frame: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    TEMP_DIR.mkdir(parents=True, exist_ok=True)

    fd, tmp_name = tempfile.mkstemp(
        prefix=f"{output_path.stem}.",
        suffix=".tmp.parquet",
        dir=TEMP_DIR,
    )
    os.close(fd)
    tmp_path = Path(tmp_name)

    try:
        frame.to_parquet(
            tmp_path,
            index=False,
            engine="pyarrow",
            compression="zstd",
        )
        os.replace(tmp_path, output_path)
    finally:
        tmp_path.unlink(missing_ok=True)


def parse_year(path: Path) -> int:
    match = YEAR_RE.search(path.name)
    if match is None:
        raise ValueError(f"Could not parse year from {path.name}")
    return int(match.group(1))


def output_paths(dataset_tag: str, year: int) -> dict[str, Path]:
    paths = {
        "clustered": (
            OUT_ROOT
            / "clustered_rp5"
            / dataset_tag
            / f"cluster_rp5_{year}.parquet"
        )
    }
    for threshold in FILTER_THRESHOLDS:
        label = threshold_label(threshold)
        paths[f"rp{label}"] = (
            OUT_ROOT
            / f"filtered_rp{label}"
            / dataset_tag
            / f"cluster_rp{label}_{year}.parquet"
        )
    return paths


# ── INPUT DISCOVERY AND QS BRIDGE ─────────────────────────────────────────


def discover_inputs() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    for spec in INPUT_SPECS:
        base_dir = Path(spec["base_dir"])
        if not base_dir.exists():
            log(f"Input base not found; skipping: {base_dir}")
            continue

        dataset_dirs = sorted(path for path in base_dir.iterdir() if path.is_dir())
        for dataset_dir in dataset_dirs:
            dataset_name = dataset_dir.name

            if dataset_name.startswith(
                spec["dataset_prefix"]
            ):
                dataset_tag = dataset_name
            else:
                dataset_tag = (
                    f"{spec['dataset_prefix']}"
                    f"{dataset_name}"
                )
            candidates = sorted(
                list(dataset_dir.glob(f"{spec['prefix']}*.qs"))
                + list(dataset_dir.glob(f"{spec['prefix']}*.parquet"))
            )

            by_year: dict[int, Path] = {}
            for path in candidates:
                year = parse_year(path)
                current = by_year.get(year)
                # Prefer Parquet when both formats exist.
                if current is None or path.suffix.lower() == ".parquet":
                    by_year[year] = path

            for year, input_path in sorted(by_year.items()):
                paths = output_paths(dataset_tag, year)
                task_index = len(rows)
                rows.append(
                    {
                        "task_index": task_index,
                        "variable": spec["variable"],
                        "dataset_tag": dataset_tag,
                        "year": year,
                        "input_path": str(input_path),
                        "clustered_output_path": str(paths["clustered"]),
                        "result_path": str(RESULT_DIR / f"task_{task_index:06d}.json"),
                    }
                )

    return rows


def convert_qs_to_parquet(input_path: Path, output_path: Path) -> None:
    r_code = r'''
args <- commandArgs(trailingOnly = TRUE)
infile <- args[[1]]
outfile <- args[[2]]
suppressPackageStartupMessages({
  library(qs)
  library(arrow)
  library(data.table)
})
x <- as.data.table(qread(infile))
arrow::write_parquet(x, outfile, compression = "zstd")
'''
    command = [RSCRIPT, "-e", r_code, str(input_path), str(output_path)]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            "QS-to-Parquet conversion failed for "
            f"{input_path}:\nSTDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )


def read_input(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path, engine="pyarrow")
    if suffix != ".qs":
        raise ValueError(f"Unsupported input format: {path}")

    TEMP_DIR.mkdir(parents=True, exist_ok=True)
    temp_path = TEMP_DIR / f"qs_{os.getpid()}_{path.stem}.parquet"
    try:
        convert_qs_to_parquet(path, temp_path)
        return pd.read_parquet(temp_path, engine="pyarrow")
    finally:
        temp_path.unlink(missing_ok=True)


# ── STANDARDIZATION ───────────────────────────────────────────────────────


def first_present(columns: pd.Index, candidates: tuple[str, ...]) -> str | None:
    return next((column for column in candidates if column in columns), None)


def standardize_input(frame: pd.DataFrame, variable: str) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()

    columns = frame.columns

    x_col = first_present(columns, ("x", "lon", "longitude"))
    y_col = first_present(columns, ("y", "lat", "latitude"))
    time_col = first_present(
        columns,
        ("time_utc", "time", "datetime_utc", "datetime", "day", "date"),
    )
    rp_col = first_present(columns, ("rp_max_years", "rp"))

    if variable == "heatindex":
        intensity_candidates = (
            "hi1_intensity",
            "hi3_intensity",
            "hi_intensity",
            "heatindex",
            "heat_index",
            "HeatIndex",
            "intensity_C",
            "value",
            "intensity",
        )
    else:
        intensity_candidates = (
            "pr24_intensity_mm",
            "intensity_mm",
            "mm",
            "precip_mm",
            "prcp_mm",
            "precipitation",
            "precip",
            "pr",
            "ppt",
            "value",
            "intensity",
        )

    intensity_col = first_present(columns, intensity_candidates)

    required = {
        "x/lon": x_col,
        "y/lat": y_col,
        "time": time_col,
        "return period": rp_col,
        "intensity": intensity_col,
    }
    missing = [name for name, value in required.items() if value is None]
    if missing:
        raise ValueError("Missing required fields: " + ", ".join(missing))

    rename = {
        str(x_col): "x",
        str(y_col): "y",
        str(time_col): "time_utc",
        str(rp_col): "rp_max_years",
        str(intensity_col): "intensity",
    }
    frame = frame.rename(columns=rename).copy()

    for column in ("x", "y", "rp_max_years", "intensity"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")

    frame["time_utc"] = pd.to_datetime(frame["time_utc"], errors="coerce", utc=True)

    frame = frame[
        frame["x"].notna()
        & frame["y"].notna()
        & frame["time_utc"].notna()
        & frame["rp_max_years"].notna()
        & frame["intensity"].notna()
        & np.isfinite(frame["x"])
        & np.isfinite(frame["y"])
        & np.isfinite(frame["rp_max_years"])
        & np.isfinite(frame["intensity"])
        & (frame["rp_max_years"] >= CLUSTER_FLOOR)
    ].copy()

    if frame.empty:
        return frame

    # One record per grid cell and native timestep. Keep the greatest RP, then
    # greatest intensity, if duplicate rows occur.
    frame.sort_values(
        ["x", "y", "time_utc", "rp_max_years", "intensity"],
        ascending=[True, True, True, False, False],
        kind="mergesort",
        inplace=True,
    )
    frame.drop_duplicates(
        subset=["x", "y", "time_utc"],
        keep="first",
        inplace=True,
    )
    frame.sort_values(["time_utc", "y", "x"], kind="mergesort", inplace=True)
    frame.reset_index(drop=True, inplace=True)
    return frame


# ── DBSCAN CUBE AND IDS ───────────────────────────────────────────────────


def minimum_positive_spacing(
    values: np.ndarray,
) -> float:
    unique_values = np.unique(
        values.astype(np.float64)
    )

    if unique_values.size < 5:
        return 1.0

    differences = np.diff(unique_values)
    positive = differences[differences > 0]

    if positive.size == 0:
        return 1.0

    return float(np.min(positive))


def prepare_cube(
    frame: pd.DataFrame,
) -> np.ndarray:
    x = frame["x"].to_numpy(
        dtype=np.float64
    )

    y = frame["y"].to_numpy(
        dtype=np.float64
    )

    dx = minimum_positive_spacing(x)
    dy = minimum_positive_spacing(y)

    x0 = float(np.min(x))
    y0 = float(np.min(y))

    col_ix = np.rint(
        (x - x0) / dx
    ).astype(np.int64)

    row_ix = np.rint(
        (y - y0) / dy
    ).astype(np.int64)

    t0 = frame["time_utc"].min()

    delta_hours = (
        (
            frame["time_utc"] - t0
        )
        .dt.total_seconds()
        .to_numpy(dtype=np.float64)
        / 3600.0
    )

    t_ix = np.rint(
        delta_hours / TIME_STEP_HOURS
    ).astype(np.int64)

    return np.column_stack(
        (
            row_ix,
            col_ix,
            t_ix,
        )
    ).astype(
        np.float64,
        copy=False,
    )


def add_stable_cluster_ids(
    frame: pd.DataFrame,
    variable: str,
    dataset_tag: str,
) -> pd.DataFrame:
    frame["cluster_id"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    positive = frame[frame["cluster"] > 0]
    if positive.empty:
        return frame

    metadata_rows: list[dict[str, Any]] = []
    for cluster_number, group in positive.groupby("cluster", sort=True):
        first_time = group["time_utc"].min()
        first = group[group["time_utc"] == first_time]
        metadata_rows.append(
            {
                "cluster": int(cluster_number),
                "first_time": first_time,
                "first_y": float(first["y"].min()),
                "first_x": float(first["x"].min()),
            }
        )

    metadata = pd.DataFrame(metadata_rows).sort_values(
        ["first_time", "first_y", "first_x", "cluster"],
        kind="mergesort",
    )
    metadata["start_date"] = metadata["first_time"].dt.strftime("%Y-%m-%d")
    metadata["daily_sequence"] = metadata.groupby("start_date").cumcount() + 1
    metadata["cluster_id"] = metadata.apply(
        lambda row: (
            f"{variable}-{dataset_tag}-{row['start_date']}-"
            f"{int(row['daily_sequence']):03d}"
        ),
        axis=1,
    )

    id_map = metadata.set_index("cluster")["cluster_id"]
    frame.loc[frame["cluster"] > 0, "cluster_id"] = (
        frame.loc[frame["cluster"] > 0, "cluster"].map(id_map).astype("string")
    )
    return frame


# ── MANIFEST ──────────────────────────────────────────────────────────────


def build_manifest() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    rows = discover_inputs()
    if not rows:
        raise RuntimeError(f"No annual exceedance files found beneath {IN_ROOT}")

    manifest = pd.DataFrame(rows)
    atomic_write_csv(manifest, MANIFEST_PATH)
    (CONTROL_DIR / "cluster_manifest_count.txt").write_text(f"{len(manifest)}\n")

    log(f"Wrote manifest: {MANIFEST_PATH}")
    log(f"Manifest tasks: {len(manifest)}")
    log(f"Datasets: {manifest['dataset_tag'].nunique()}")
    log(f"Years: {manifest['year'].min()}-{manifest['year'].max()}")


# ── ONE TASK ──────────────────────────────────────────────────────────────


def run_manifest_row(row: pd.Series) -> str:
    task_index = int(row["task_index"])
    variable = str(row["variable"])
    dataset_tag = str(row["dataset_tag"])
    year = int(row["year"])
    input_path = Path(row["input_path"])
    clustered_path = Path(row["clustered_output_path"])
    result_path = Path(row["result_path"])
    paths = output_paths(dataset_tag, year)

    expected_outputs = list(paths.values())
    if all(path.exists() for path in expected_outputs) and not OVERWRITE:
        result = {
            "task_index": task_index,
            "worker_index": WORKER_INDEX,
            "variable": variable,
            "dataset_tag": dataset_tag,
            "year": year,
            "input_path": str(input_path),
            "status": "existing",
        }
        atomic_write_json(result, result_path)
        log(f"Task {task_index}: all outputs exist; skipped")
        return "existing"

    if not input_path.exists():
        raise FileNotFoundError(f"Missing input: {input_path}")

    log(f"Task {task_index}: {variable}/{dataset_tag}/{year}")
    raw = read_input(input_path)
    n_input = len(raw)
    frame = standardize_input(raw, variable)
    n_cluster_input = len(frame)

    frame["cluster"] = pd.Series(dtype="int32")
    frame["cluster_id"] = pd.Series(dtype="string")
    frame["cluster_max_rp_years"] = pd.Series(dtype="float64")

    if n_cluster_input > 0:
        cube = prepare_cube(frame)
        labels = DBSCAN(
            eps=EPS,
            min_samples=MIN_SAMPLES,
            metric="euclidean",
            algorithm="auto",
            n_jobs=1,
        ).fit_predict(cube)

        frame["cluster"] = np.where(labels < 0, 0, labels + 1).astype(np.int32)
        frame = add_stable_cluster_ids(frame, variable, dataset_tag)

        cluster_max = (
            frame.loc[frame["cluster"] > 0]
            .groupby("cluster", sort=False)["rp_max_years"]
            .max()
        )
        frame["cluster_max_rp_years"] = frame["cluster"].map(cluster_max)

    frame["variable"] = variable
    frame["dataset_tag"] = dataset_tag
    frame["cluster_floor_years"] = CLUSTER_FLOOR
    frame["eps"] = EPS
    frame["min_samples"] = MIN_SAMPLES
    frame["time_step_hours"] = TIME_STEP_HOURS

    atomic_write_parquet(frame, clustered_path)

    positive = frame[frame["cluster"] > 0].copy()
    threshold_counts: dict[str, int] = {}
    threshold_cluster_counts: dict[str, int] = {}

    for threshold in FILTER_THRESHOLDS:
        label = threshold_label(threshold)
        subset = positive[
            positive["cluster_max_rp_years"] >= threshold
        ].copy()
        subset["event_filter_rp_years"] = threshold
        atomic_write_parquet(subset, paths[f"rp{label}"])
        threshold_counts[f"n_rows_rp{label}"] = int(len(subset))
        threshold_cluster_counts[f"n_clusters_rp{label}"] = int(
            subset["cluster"].nunique() if not subset.empty else 0
        )

    n_noise = int((frame["cluster"] == 0).sum()) if not frame.empty else 0
    n_clusters = int(positive["cluster"].nunique()) if not positive.empty else 0

    result = {
        "task_index": task_index,
        "worker_index": WORKER_INDEX,
        "variable": variable,
        "dataset_tag": dataset_tag,
        "year": year,
        "input_path": str(input_path),
        "clustered_output_path": str(clustered_path),
        "n_input": n_input,
        "n_after_rp5_filter": n_cluster_input,
        "n_noise": n_noise,
        "n_clusters_rp5_base": n_clusters,
        **threshold_counts,
        **threshold_cluster_counts,
        "status": "written",
    }
    atomic_write_json(result, result_path)
    log(
        f"Task {task_index} complete: rows={n_cluster_input}; "
        f"clusters={n_clusters}; noise={n_noise}"
    )
    return "written"


# ── ARRAY WORKER ──────────────────────────────────────────────────────────


def run_cluster_worker() -> None:
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(f"Missing manifest: {MANIFEST_PATH}")
    if not 0 <= WORKER_INDEX < WORKER_COUNT:
        raise ValueError(
            f"Worker index {WORKER_INDEX} outside 0-{WORKER_COUNT - 1}"
        )

    manifest = pd.read_csv(MANIFEST_PATH)
    assigned = list(range(WORKER_INDEX, len(manifest), WORKER_COUNT))
    log(
        f"Worker {WORKER_INDEX} assigned {len(assigned)} of "
        f"{len(manifest)} tasks"
    )

    status_counts: dict[str, int] = {}
    for sequence, position in enumerate(assigned, start=1):
        log(f"Worker {WORKER_INDEX} progress {sequence}/{len(assigned)}")
        status = run_manifest_row(manifest.iloc[position])
        status_counts[status] = status_counts.get(status, 0) + 1

    log(f"Worker {WORKER_INDEX} finished: {status_counts}")


# ── SUMMARIZE ─────────────────────────────────────────────────────────────


def summarize() -> None:
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(f"Missing manifest: {MANIFEST_PATH}")

    manifest = pd.read_csv(MANIFEST_PATH)
    records: list[dict[str, Any]] = []
    missing: list[int] = []

    for row in manifest.itertuples(index=False):
        path = Path(row.result_path)
        if not path.exists():
            missing.append(int(row.task_index))
        else:
            records.append(json.loads(path.read_text()))

    if missing:
        raise RuntimeError(
            f"Missing {len(missing)} task result files; first: {missing[:20]}"
        )

    summary = pd.DataFrame(records).sort_values(
        ["variable", "dataset_tag", "year"],
        kind="mergesort",
    )
    summary.reset_index(drop=True, inplace=True)
    atomic_write_csv(summary, OUT_ROOT / "cluster_summary_all_datasets.csv")

    for dataset_tag, subset in summary.groupby("dataset_tag", sort=True):
        dataset_dir = OUT_ROOT / "clustered_rp5" / str(dataset_tag)
        atomic_write_csv(
            subset.copy(),
            dataset_dir / f"cluster_summary_{dataset_tag}.csv",
        )

    log("Wrote global and per-dataset summaries")


# ── VALIDATION AND MAIN ───────────────────────────────────────────────────


def validate() -> None:
    try:
        import pyarrow  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "Python package pyarrow is required"
        ) from exc

    qs_inputs_exist = any(
        spec["base_dir"].exists()
        and any(
            spec["base_dir"].rglob("*.qs")
        )
        for spec in INPUT_SPECS
    )

    if qs_inputs_exist:
        try:
            completed = subprocess.run(
                [
                    RSCRIPT,
                    "-e",
                    (
                        'pkgs <- c("qs", "arrow", "data.table"); '
                        'missing <- pkgs['
                        '!vapply(pkgs, requireNamespace, '
                        'logical(1), quietly=TRUE)]; '
                        'if (length(missing)) '
                        'stop(paste("Missing R packages:", '
                        'paste(missing, collapse=", ")))'
                    ),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
        except FileNotFoundError as exc:
            raise RuntimeError(
                f"Rscript executable not found: {RSCRIPT}"
            ) from exc

        if completed.returncode != 0:
            raise RuntimeError(
                "R dependencies needed for .qs inputs "
                "are unavailable.\n"
                f"STDOUT:\n{completed.stdout}\n"
                f"STDERR:\n{completed.stderr}"
            )

    if EPS <= 0:
        raise ValueError(
            "CISC_DBSCAN_EPS must be > 0"
        )

    if MIN_SAMPLES < 1:
        raise ValueError(
            "CISC_DBSCAN_MIN_SAMPLES must be >= 1"
        )

    if TIME_STEP_HOURS <= 0:
        raise ValueError(
            "CISC_CLUSTER_TIME_STEP_HOURS must be > 0"
        )

    if CLUSTER_FLOOR < 1:
        raise ValueError(
            "CISC_CLUSTER_FLOOR_YEARS must be >= 1"
        )

    if WORKER_COUNT < 1:
        raise ValueError(
            "CISC_CLUSTER_WORKER_COUNT must be >= 1"
        )

    OUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    CONTROL_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    RESULT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    TEMP_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )


def main() -> None:
    validate()
    log("==== 11_cluster_rp5_filter.py ====")
    log(f"MODE={MODE}")
    log(f"DATA_OUT={DATA_OUT}")
    log(f"IN_ROOT={IN_ROOT}")
    log(f"OUT_ROOT={OUT_ROOT}")
    log(f"EPS={EPS}; MIN_SAMPLES={MIN_SAMPLES}")
    log(f"CLUSTER_FLOOR={CLUSTER_FLOOR:g}")
    log(f"FILTER_THRESHOLDS={FILTER_THRESHOLDS}")
    log(f"WORKER_INDEX={WORKER_INDEX}; WORKER_COUNT={WORKER_COUNT}")
    log(f"OVERWRITE={OVERWRITE}")

    if MODE == "build_manifest":
        build_manifest()
    elif MODE == "cluster_worker":
        run_cluster_worker()
    elif MODE == "summarize":
        summarize()
    else:
        raise ValueError(
            "CISC_CLUSTER_MODE must be build_manifest, cluster_worker, or summarize"
        )

    log("Done.")


if __name__ == "__main__":
    main()
