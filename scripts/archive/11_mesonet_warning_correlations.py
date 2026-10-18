#!/usr/bin/env python3
"""
Correlate completed daily Mesonet warning metrics with daily covariance-model
Ar and sill values for heat index and precipitation.

Mesonet warning groups are read from:
    data/out/10_mesonet/
    mesonet_warning_combined_daily_southeast_2002_2023.parquet

Figures
-------
One two-column figure for each variable and covariance metric:

    Heat index Ar
    Heat index sill
    Precipitation Ar
    Precipitation sill

Left panel:
    Daily spatially distinct active-warning count

Right panel:
    Daily unioned active-warning area

Only dates with positive warning count and positive warning area are retained.
"""

from __future__ import annotations

import argparse
import logging
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr


# =============================================================================
# Configuration
# =============================================================================

DEFAULT_ROOT = Path(
    os.environ.get(
        "CISC_ROOT",
        "/proj/sebastianlab/users/hquintal/cisc",
    )
)

START_DATE = pd.Timestamp("2002-01-01")
END_DATE = pd.Timestamp("2023-12-31")

WARNING_GROUP_TITLES = {
    "heat": "Heat and excessive heat warnings",
    "precipitation": (
        "Flash flood, severe thunderstorm, tropical storm, "
        "hurricane, and storm warnings"
    ),
}

VARIABLE_LABELS = {
    "heat": "Heat index",
    "precipitation": "Precipitation",
}

METRIC_LABELS = {
    "Ar": r"$A_r$",
    "sill": "Sill",
    "stm1": "STM1",
    "stm2": "STM2",
}


# =============================================================================
# Paths
# =============================================================================


@dataclass(frozen=True)
class WorkflowPaths:
    root: Path

    @property
    def warning_file(self) -> Path:
        return (
            self.root
            / "data"
            / "out"
            / "10_mesonet"
            / "mesonet_warning_combined_daily_southeast_2002_2023.parquet"
        )

    @property
    def heat_file(self) -> Path:
        return (
            self.root
            / "data"
            / "out"
            / "03_cov_model"
            / "heatindex"
            / "heatindex_daily_optimal_cov_model_southeast.csv"
        )

    @property
    def precipitation_file(self) -> Path:
        return (
            self.root
            / "data"
            / "out"
            / "03_cov_model"
            / "precipitation"
            / "precipitation_daily_optimal_cov_model_southeast.csv"
        )

    @property
    def output_dir(self) -> Path:
        return (
            self.root
            / "data"
            / "out"
            / "10_mesonet"
            / "covariance_correlations"
        )

    @property
    def figure_dir(self) -> Path:
        return (
            self.root
            / "figures"
            / "fig_s09_mesonet_covariance"
        )

    @property
    def summary_file(self) -> Path:
        return (
            self.output_dir
            / "mesonet_covariance_correlation_summary.csv"
        )

    def covariance_file(self, warning_group: str) -> Path:
        if warning_group == "heat":
            return self.heat_file

        if warning_group == "precipitation":
            return self.precipitation_file

        raise ValueError(
            f"Unknown warning group: {warning_group}"
        )

    def merged_output(
        self,
        warning_group: str,
        metric: str,
    ) -> Path:
        return (
            self.output_dir
            / (
                f"mesonet_{warning_group}_"
                f"{metric.lower()}_2002_2023.csv"
            )
        )

    def figure_output(
        self,
        warning_group: str,
        metric: str,
    ) -> Path:
        return (
            self.figure_dir
            / (
                f"fig_s09_mesonet_"
                f"{warning_group}_{metric.lower()}.png"
            )
        )


# =============================================================================
# Logging
# =============================================================================


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format=(
            "%(asctime)s | %(levelname)s | %(message)s"
        ),
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )


# =============================================================================
# Input processing
# =============================================================================


def identify_column(
    frame: pd.DataFrame,
    candidates: list[str],
    description: str,
) -> str:
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate

    lower_lookup = {
        str(column).lower(): column
        for column in frame.columns
    }

    for candidate in candidates:
        match = lower_lookup.get(
            candidate.lower()
        )

        if match is not None:
            return match

    raise KeyError(
        f"Could not identify {description}. "
        f"Available columns: {list(frame.columns)}"
    )


def read_warnings(
    paths: WorkflowPaths,
) -> pd.DataFrame:
    if not paths.warning_file.is_file():
        raise FileNotFoundError(
            f"Missing warning file: {paths.warning_file}"
        )

    warnings = pd.read_parquet(
        paths.warning_file
    )

    required = {
        "date",
        "year",
        "warning_group",
        "event_count",
        "area_km2",
    }

    missing = required.difference(
        warnings.columns
    )

    if missing:
        raise KeyError(
            "Warning file is missing columns: "
            f"{sorted(missing)}"
        )

    warnings["date"] = pd.to_datetime(
        warnings["date"],
        errors="coerce",
    )

    warnings["event_count"] = pd.to_numeric(
        warnings["event_count"],
        errors="coerce",
    )

    warnings["area_km2"] = pd.to_numeric(
        warnings["area_km2"],
        errors="coerce",
    )

    warnings = warnings.loc[
        warnings["date"].notna()
        & warnings["date"].between(
            START_DATE,
            END_DATE,
        )
        & warnings["warning_group"].isin(
            [
                "heat",
                "precipitation",
            ]
        )
        & np.isfinite(
            warnings["event_count"]
        )
        & np.isfinite(
            warnings["area_km2"]
        )
        & (warnings["event_count"] > 0)
        & (warnings["area_km2"] > 0)
    ].copy()

    # Protect against any remaining duplicated date/group rows.
    warnings = (
        warnings.groupby(
            [
                "date",
                "warning_group",
            ],
            as_index=False,
            observed=True,
        )
        .agg(
            event_count=(
                "event_count",
                "sum",
            ),
            area_km2=(
                "area_km2",
                "sum",
            ),
        )
        .sort_values(
            [
                "warning_group",
                "date",
            ]
        )
        .reset_index(drop=True)
    )

    warnings["year"] = (
        warnings["date"].dt.year.astype(
            "int16"
        )
    )

    return warnings


def read_covariance_metrics(
    input_file: Path,
    variable_name: str,
) -> pd.DataFrame:
    if not input_file.is_file():
        raise FileNotFoundError(
            f"Missing {variable_name} file: "
            f"{input_file}"
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

    sill_column = identify_column(
        frame,
        [
            "sill",
            "Sill",
            "SILL",
            "partial_sill",
            "partial.sill",
            "psill",
        ],
        f"{variable_name} sill column",
    )
    
    stm1_column = identify_column(
        frame,
        [
            "stm1",
            "STM1",
            "space_time_metric_1",
            "space_time_metric1",
        ],
        f"{variable_name} STM1 column",
    )

    stm2_column = identify_column(
        frame,
        [
            "stm2",
            "STM2",
            "space_time_metric_2",
            "space_time_metric2",
        ],
        f"{variable_name} STM2 column",
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
            "sill": pd.to_numeric(
                frame[sill_column],
                errors="coerce",
            ),
            "stm1": pd.to_numeric(
                frame[stm1_column],
                errors="coerce",
            ),
            "stm2": pd.to_numeric(
                frame[stm2_column],
                errors="coerce",
            ),
        }
    )

    output = output.loc[
        output["date"].notna()
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
            sill=("sill", "mean"),
            stm1=("stm1", "mean"),
            stm2=("stm2", "mean"),
        )
        .sort_values("date")
        .reset_index(drop=True)
    )

    return output


# =============================================================================
# Correlation helpers
# =============================================================================


def correlation_statistics(
    frame: pd.DataFrame,
    metric: str,
    response: str,
) -> dict[str, float | int]:
    complete = frame.loc[
        np.isfinite(frame[metric])
        & np.isfinite(frame[response])
    ].copy()

    if len(complete) < 3:
        return {
            "n": len(complete),
            "r": np.nan,
            "p_value": np.nan,
        }

    result = pearsonr(
        complete[metric],
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


def rounded_axis_upper(
    maximum: float,
) -> float:
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

    axis.plot(
        x_line,
        np.polyval(
            coefficients,
            x_line,
        ),
        linewidth=0.8,
    )


# =============================================================================
# Plotting
# =============================================================================


def create_correlation_figure(
    paths: WorkflowPaths,
    warnings: pd.DataFrame,
    covariance: pd.DataFrame,
    warning_group: str,
    metric: str,
) -> list[dict[str, object]]:
    plot_data = (
        warnings.loc[
            warnings["warning_group"].eq(
                warning_group
            )
        ]
        .merge(
            covariance[
                [
                    "date",
                    metric,
                ]
            ],
            on="date",
            how="inner",
            validate="one_to_one",
        )
    )

    plot_data = plot_data.loc[
        np.isfinite(
            plot_data[metric]
        )
        & (plot_data["event_count"] > 0)
        & (plot_data["area_km2"] > 0)
    ].copy()

    if plot_data.empty:
        raise RuntimeError(
            f"No overlapping {warning_group} and "
            f"{metric} observations"
        )

    plot_data.to_csv(
        paths.merged_output(
            warning_group,
            metric,
        ),
        index=False,
        date_format="%Y-%m-%d",
    )

    count_statistics = correlation_statistics(
        plot_data,
        metric,
        "event_count",
    )

    area_statistics = correlation_statistics(
        plot_data,
        metric,
        "area_km2",
    )

    figure, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(6.5, 3.4),
        constrained_layout=True,
    )

    panels = [
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

    x_label = (
        f"{VARIABLE_LABELS[warning_group]} "
        f"{METRIC_LABELS[metric]}"
    )

    for panel_index, (
        axis,
        response,
        y_label,
        statistics,
    ) in enumerate(panels):
        axis.scatter(
            plot_data[metric],
            plot_data[response],
            s=8,
            alpha=0.25,
            linewidths=0,
        )

        add_regression(
            axis,
            plot_data[metric].to_numpy(
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
            linewidth=0.3,
            alpha=0.3,
        )

        axis.tick_params(
            axis="both",
            labelsize=7,
            width=0.6,
            length=3,
        )

        axis.spines["top"].set_visible(
            False
        )

        axis.spines["right"].set_visible(
            False
        )

    figure.suptitle(
        (
            f"{WARNING_GROUP_TITLES[warning_group]}\n"
            f"{VARIABLE_LABELS[warning_group]} "
            f"{METRIC_LABELS[metric]}, 2002–2023"
        ),
        fontsize=10,
        fontweight="bold",
    )

    output_file = paths.figure_output(
        warning_group,
        metric,
    )

    figure.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(figure)

    logging.info(
        "Wrote figure: %s",
        output_file,
    )

    return [
        {
            "warning_group": warning_group,
            "covariance_metric": metric,
            "response": "event_count",
            **count_statistics,
        },
        {
            "warning_group": warning_group,
            "covariance_metric": metric,
            "response": "area_km2",
            **area_statistics,
        },
    ]


# =============================================================================
# Main workflow
# =============================================================================


def run_correlations(
    paths: WorkflowPaths,
) -> None:
    warnings = read_warnings(
        paths
    )

    covariance_by_group = {
        "heat": read_covariance_metrics(
            paths.heat_file,
            "heat index",
        ),
        "precipitation": read_covariance_metrics(
            paths.precipitation_file,
            "precipitation",
        ),
    }

    summary_rows: list[
        dict[str, object]
    ] = []

    for warning_group in (
        "heat",
        "precipitation",
    ):
        for metric in (
            "Ar",
            "sill",
            "stm1",
            "stm2",
        ):
            summary_rows.extend(
                create_correlation_figure(
                    paths=paths,
                    warnings=warnings,
                    covariance=(
                        covariance_by_group[
                            warning_group
                        ]
                    ),
                    warning_group=warning_group,
                    metric=metric,
                )
            )

    pd.DataFrame(
        summary_rows
    ).to_csv(
        paths.summary_file,
        index=False,
    )

    logging.info(
        "Wrote correlation summary: %s",
        paths.summary_file,
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Correlate completed daily Mesonet warning "
            "metrics with daily Ar and sill."
        )
    )

    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="CISC project root.",
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

    paths.figure_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    run_correlations(
        paths
    )

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        logging.exception(
            "Mesonet correlation workflow failed"
        )
        sys.exit(1)
