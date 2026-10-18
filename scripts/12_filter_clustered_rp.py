#!/usr/bin/env python3

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pandas as pd


ROOT = Path(
    os.environ.get(
        "CISC_ROOT",
        "/proj/sebastianlab/users/hquintal/cisc",
    )
)

CLUSTER_ROOT = ROOT / "data/out/12_clusters"
INPUT_ROOT = CLUSTER_ROOT / "clustered_rp5"

THRESHOLDS = (10.0, 25.0, 50.0)

WORKER_INDEX = int(
    os.environ.get(
        "SLURM_ARRAY_TASK_ID",
        "0",
    )
)

WORKER_COUNT = int(
    os.environ.get(
        "CISC_FILTER_WORKER_COUNT",
        "4",
    )
)


def threshold_label(value: float) -> str:
    if float(value).is_integer():
        return str(int(value))

    return str(value).replace(".", "p")


def atomic_write_parquet(
    frame: pd.DataFrame,
    output_path: Path,
) -> None:
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fd, temporary_name = tempfile.mkstemp(
        prefix=f"{output_path.stem}.",
        suffix=".tmp.parquet",
        dir=output_path.parent,
    )

    os.close(fd)
    temporary_path = Path(temporary_name)

    try:
        frame.to_parquet(
            temporary_path,
            index=False,
            engine="pyarrow",
            compression="zstd",
        )

        os.replace(
            temporary_path,
            output_path,
        )

    finally:
        temporary_path.unlink(
            missing_ok=True
        )


def main() -> None:
    files = sorted(
        INPUT_ROOT.glob(
            "*/cluster_rp5_*.parquet"
        )
    )

    if not files:
        raise RuntimeError(
            f"No clustered RP5 files found under {INPUT_ROOT}"
        )

    assigned = files[
        WORKER_INDEX::WORKER_COUNT
    ]

    print(
        f"Worker {WORKER_INDEX}: "
        f"{len(assigned)} of {len(files)} files",
        flush=True,
    )

    for sequence, input_path in enumerate(
        assigned,
        start=1,
    ):
        dataset_tag = input_path.parent.name
        year = input_path.stem.rsplit("_", 1)[-1]

        frame = pd.read_parquet(
            input_path
        )

        required = {
            "cluster",
            "cluster_id",
            "cluster_max_rp_years",
        }

        missing = required.difference(
            frame.columns
        )

        if missing:
            raise ValueError(
                f"{input_path} is missing: "
                + ", ".join(sorted(missing))
            )

        positive = frame[
            frame["cluster"] > 0
        ].copy()

        for threshold in THRESHOLDS:
            label = threshold_label(
                threshold
            )

            subset = positive[
                positive[
                    "cluster_max_rp_years"
                ]
                >= threshold
            ].copy()

            subset[
                "event_filter_rp_years"
            ] = threshold

            output_path = (
                CLUSTER_ROOT
                / f"filtered_rp{label}"
                / dataset_tag
                / f"cluster_rp{label}_{year}.parquet"
            )

            atomic_write_parquet(
                subset,
                output_path,
            )

        if (
            sequence == 1
            or sequence % 20 == 0
            or sequence == len(assigned)
        ):
            print(
                f"Worker {WORKER_INDEX}: "
                f"{sequence}/{len(assigned)}",
                flush=True,
            )

    print(
        f"Worker {WORKER_INDEX} complete",
        flush=True,
    )


if __name__ == "__main__":
    main()
