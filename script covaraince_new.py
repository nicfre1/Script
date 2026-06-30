#script updated:

#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def build_arg_parser():
    p = argparse.ArgumentParser(
        description="Compute group structural covariance from thickness TSV files."
    )
    p.add_argument(
        "in_rh_stats",
        type=Path,
        help="Right hemisphere",
    )
    p.add_argument(
        "in_lh_stats",
        type=Path,
        help="Left hemisphere",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="djdj",
    )

    return p


def main():
    args = build_arg_parser().parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    subject_data = {}

    for tsv_file in [args.in_rh_stats, args.in_lh_stats]:
        df = pd.read_csv(tsv_file, sep="\t")

        if "sample" not in df.columns:
            raise ValueError(f"Missing sample column in {tsv_file}")

        thickness_cols = [
            col
            for col in df.columns
            if col.endswith("_thickness") and "MeanThickness" not in col
        ]

        for _, row in df.iterrows():
            subject_id = str(row["sample"]).strip()
            subject_id = subject_id.replace("_fs", "").replace("_ses-baseline", "")

            values = pd.to_numeric(row[thickness_cols], errors="coerce").dropna()

            subject_data.setdefault(subject_id, {})
            subject_data[subject_id].update(values.to_dict())

    if not subject_data:
        raise ValueError("No thickness data found in input TSV files.")

    thickness = pd.DataFrame.from_dict(subject_data, orient="index")
    thickness.index.name = "subject"
    thickness = thickness.sort_index().dropna(axis=1)

    if thickness.shape[0] < 2:
        raise ValueError("Group covariance requires at least 2 subjects.")

    regions = list(thickness.columns)
    x = thickness.to_numpy(dtype=float)

    x_z = (x - np.mean(x, axis=0)) / np.std(x, axis=0)
    group_covariance = np.corrcoef(x_z, rowvar=False)

    thickness.to_csv(output_dir / "subject_region_thickness_matrix.csv")

    pd.DataFrame(
        group_covariance,
        index=regions,
        columns=regions,
    ).to_csv(output_dir / "group_covariance.csv")


if __name__ == "__main__":
    main()