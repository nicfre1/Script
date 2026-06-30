#Code pour script covariance 

#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def build_arg_parser():
    p = argparse.ArgumentParser(
        description="Compute individual and group covariance matrices from thickness TSV files."
    )

    p.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing thickness TSV files.",
    )

    p.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where covariance matrices will be saved.",
    )

    return p


def clean_subject_id(sample):
    subject = str(sample).strip()
    subject = subject.replace("_fs", "")
    subject = subject.replace("_ses-baseline", "")
    return subject


def load_thickness_tsv(input_dir):
    tsv_files = sorted(input_dir.glob("*.tsv"))

    if not tsv_files:
        raise FileNotFoundError(f"No .tsv files found in: {input_dir}")

    subject_data = {}

    for tsv_file in tsv_files:
        df = pd.read_csv(tsv_file, sep="\t")

        thickness_columns = [
            col for col in df.columns
            if col.endswith("_thickness") and "MeanThickness" not in col
        ]

        for _, row in df.iterrows():
            subject_id = clean_subject_id(row["sample"])
            subject_data[subject_id] = row[thickness_columns].astype(float).to_dict()

    thickness_df = pd.DataFrame.from_dict(subject_data, orient="index")
    thickness_df.index.name = "subject"
    thickness_df = thickness_df.sort_index()

    return thickness_df


def compute_individual_covariance(thickness_values):
    thickness_values = np.asarray(thickness_values, dtype=float)

    z = (thickness_values - np.mean(thickness_values)) / np.std(thickness_values)
    diff = z[:, None] - z[None, :]

    return np.exp(-(diff ** 2))


def compute_group_covariance(thickness_matrix):
    thickness_matrix = np.asarray(thickness_matrix, dtype=float)

    x_z = (
        thickness_matrix - np.mean(thickness_matrix, axis=0)
    ) / np.std(thickness_matrix, axis=0)

    return np.corrcoef(x_z, rowvar=False)


def save_matrix(matrix, regions, output_file):
    output_file.parent.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(matrix, index=regions, columns=regions)
    df.to_csv(output_file)


def main():
    args = build_arg_parser().parse_args()

    output_dir = args.output_dir
    individual_dir = output_dir / "individual_covariance"
    output_dir.mkdir(parents=True, exist_ok=True)

    thickness_df = load_thickness_tsv(args.input_dir)
    regions = list(thickness_df.columns)

    thickness_df.to_csv(output_dir / "subject_region_thickness_matrix.csv")

    for subject_id, row in thickness_df.iterrows():
        matrix = compute_individual_covariance(row.values)
        output_file = individual_dir / f"{subject_id}_individual_covariance.csv"
        save_matrix(matrix, regions, output_file)
        print(f"[OK] individual covariance: {output_file}")

    group_matrix = compute_group_covariance(thickness_df.values)
    save_matrix(group_matrix, regions, output_dir / "group_covariance.csv")

    print(f"[OK] group covariance: {output_dir / 'group_covariance.csv'}")
    print(f"Subjects: {thickness_df.shape[0]}")
    print(f"Regions: {thickness_df.shape[1]}")


if __name__ == "__main__":
    main()