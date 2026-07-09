#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute group, individual covariance matrix, group and individual gradients,
and lambdas from left/right thickness TSV files.

Example:
compute_covariance.py lh_stats.tsv rh_stats.tsv
"""

import argparse
from pathlib import Path
import logging

import numpy as np
import pandas as pd
from scipy.stats import zscore
from brainspace.gradient import GradientMaps


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    p.add_argument("i_stats", nargs=2, help="Left and right stats - thickness atlas.")

    p.add_argument(
        "--out_dir",
        default="results",
        help="Directory where output CSV files will be saved.",
    )

    p.add_argument(
        "-v",
        default="WARNING",
        const="INFO",
        nargs="?",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        dest="verbose",
        help="Produces verbose output depending on "
        "the provided level. \nDefault level is warning, "
        "default when using -v is info.",
    )

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    logging.getLogger().setLevel(logging.getLevelName(args.verbose))

    i_stats = args.i_stats

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    i_left = i_stats[0]
    i_right = i_stats[1]

    logging.info("Read TSV files")
    df_left = pd.read_csv(i_left, sep="\t")
    df_right = pd.read_csv(i_right, sep="\t")

    logging.info("Filter out column not thickness")
    left = df_left.loc[
        :,
        df_left.columns.str.contains("thickness") | df_left.columns.str.match("sample"),
    ]

    right = df_right.loc[
        :,
        df_right.columns.str.contains("thickness")
        | df_right.columns.str.match("sample"),
    ]

    thickness_df = pd.merge(left, right, on="sample")

    logging.debug("Compute zscore")
    df_zscore = thickness_df.copy()
    numeric_cols = thickness_df.columns[thickness_df.columns != "sample"]

    df_zscore[numeric_cols] = pd.DataFrame(
        zscore(thickness_df[numeric_cols], axis=1),
        columns=numeric_cols,
        index=thickness_df.index,
    )

    logging.debug("Average left and right zscores")
    zscore_regions = df_zscore[numeric_cols].copy()

    zscore_regions.columns = (
        zscore_regions.columns.str.replace("^lh_", "", regex=True)
        .str.replace("^rh_", "", regex=True)
        .str.replace("_L_", "_", regex=False)
        .str.replace("_R_", "_", regex=False)
    )

    zscore_regions = zscore_regions.T.groupby(level=0).mean().T

    df_zscore = pd.concat(
        [thickness_df["sample"], zscore_regions],
        axis=1,
    )

    numeric_cols = zscore_regions.columns

    logging.debug("Create empty object to store covariance matrices")
    index_multi = pd.MultiIndex.from_product(
        [thickness_df["sample"], numeric_cols],
        names=["subjects", "covariance"],
    )

    df_3d = pd.DataFrame(index=index_multi, columns=numeric_cols, dtype=float)

    for subject in thickness_df["sample"]:

        z = df_zscore.loc[
            thickness_df["sample"] == subject, numeric_cols
        ].values.flatten()

        diff = z[:, None] - z[None, :]
        matrix = np.exp(-(diff**2))

        df_3d.loc[subject] = matrix

    logging.info("Individual Matrices")
    # print(df_3d)

    logging.debug("Group Matrix")
    group_matrix = df_3d.groupby(level="covariance").mean()
    group_matrix = group_matrix.loc[numeric_cols, numeric_cols]

    # print(group_matrix)

    logging.debug("Compute cortical Gradients")
    gradient_model = GradientMaps(
        n_components=10,
        approach="dm",
        kernel="normalized_angle",
        random_state=0,
    )

    gradient_model.fit(group_matrix.to_numpy(dtype=float))

    gradients = pd.DataFrame(
        gradient_model.gradients_,
        index=group_matrix.index,
        columns=[
            f"gradient_{i + 1}" for i in range(gradient_model.gradients_.shape[1])
        ],
    )

    lambdas = pd.DataFrame(
        gradient_model.lambdas_,
        index=[f"gradient_{i + 1}" for i in range(len(gradient_model.lambdas_))],
        columns=["lambda"],
    )

    logging.debug("Cortical Gradients")
    # print(gradients)

    logging.debug("Eigenvalues")
    # print(lambdas)

    logging.debug("Compute individual cortical gradients")

    subjects = thickness_df["sample"].tolist()

    individual_covariances = [
        df_3d.loc[subject].to_numpy(dtype=float) for subject in subjects
    ]

    gm_indiv = GradientMaps(
        n_components=10,
        approach="dm",
        kernel="normalized_angle",
        random_state=0,
        alignment="procrustes",
    )

    gm_indiv.fit(individual_covariances, reference=gradient_model.gradients_)

    individual_gradients = pd.concat(
        {
            subject: pd.DataFrame(
                aligned_gradients,
                index=numeric_cols,
                columns=[
                    f"gradient_{i + 1}" for i in range(aligned_gradients.shape[1])
                ],
            )
            for subject, aligned_gradients in zip(subjects, gm_indiv.aligned_)
        },
        names=["subjects", "region"],
    )

    individual_lambdas = pd.DataFrame(
        gm_indiv.lambdas_,
        index=subjects,
        columns=[f"lambda_{i + 1}" for i in range(len(gm_indiv.lambdas_[0]))],
    )

    individual_lambdas.index.name = "subjects"

    logging.debug("Individual Cortical Gradients")
    # print(individual_gradients)

    logging.debug("Individual Eigenvalues")
    # print(individual_lambdas)

    logging.info("Save CSV files")

    output_prefix = "lh_rh_average"

    df_3d.to_csv(out_dir / f"{output_prefix}_individual_matrices.csv")

    group_matrix.to_csv(out_dir / f"{output_prefix}_group_matrix.csv")

    gradients.to_csv(out_dir / f"{output_prefix}_group_gradients.csv")

    lambdas.to_csv(out_dir / f"{output_prefix}_group_lambdas.csv")

    individual_gradients.to_csv(out_dir / f"{output_prefix}_individual_gradients.csv")

    individual_lambdas.to_csv(out_dir / f"{output_prefix}_individual_lambdas.csv")


if __name__ == "__main__":
    main()
