#try #2 

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute group covariance matrix from thickness TSV files.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    p.add_argument(
        "out_dir",
        help="Directory where output CSV files will be saved.",
    )

    p.add_argument(
        "--input-dir",
        default=".",
        help="Directory containing TSV files. Default: current directory.",
    )

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    data = {}

    for tsv_file in sorted(input_dir.glob("*.tsv")):
        df = pd.read_csv(tsv_file, sep="\t")

        cols = [
            c for c in df.columns
            if c.endswith("_thickness") and "MeanThickness" not in c
        ]

        for _, row in df.iterrows():
            subject = str(row["sample"]).replace("_fs", "").replace("_ses-baseline", "")
            data[subject] = pd.to_numeric(row[cols], errors="coerce").to_dict()

    X_df = pd.DataFrame.from_dict(data, orient="index").dropna(axis=1)
    X_df.index.name = "subject"
    X_df.to_csv(out_dir / "subject_region_thickness_matrix.csv")

    X = X_df.to_numpy(dtype=float)

    X_z = (X - np.mean(X, axis=0)) / np.std(X, axis=0)
    group_covariance = np.corrcoef(X_z, rowvar=False)

    regions = list(X_df.columns)
    pd.DataFrame(
        group_covariance,
        index=regions,
        columns=regions,
    ).to_csv(out_dir / "group_covariance.csv")

    print(f"[OK] Saved outputs in: {out_dir}")
    print(f"Subjects: {X_df.shape[0]}")
    print(f"Regions: {X_df.shape[1]}")


if __name__ == "__main__":
    main()