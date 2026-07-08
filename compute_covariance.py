#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute group and individual covariance matrix from thickness TSV files.

Example: compute_covariance.py lh_stats.tsv rh_stats.tsv

"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import zscore

def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    p.add_argument("i_stats",
        nargs=2,
        help="Left and right stats - thickness atlas."
    )

    p.add_argument(
        "--out_dir",
        default='results',
        help="Directory where output CSV files will be saved.",
    )

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    i_stats = args.i_stats
   

    for current_thickness in i_stats:

        print("Read TSV file")
        df = pd.read_csv(current_thickness, sep="\t")

        print("Filter out column not thickness")
        toto = df.iloc[:, df.columns.str.contains('thickness') | df.columns.str.match('sample')] 
        
        print("Compute zscore")
        df_zscore = toto.copy()
        numeric_cols = toto.columns[toto.columns != "sample"]
        df_zscore[numeric_cols] = pd.DataFrame(zscore(toto[numeric_cols], axis=1),
                                                      columns=numeric_cols,
                                                      index=toto.index)
        


        print("Create empty object to store covariance matrices")
        index_multi = pd.MultiIndex.from_product(
                        [toto["sample"], numeric_cols], 
                        names=['subjects', 'covariance'])
        df_3d = pd.DataFrame(index=index_multi, columns=numeric_cols, dtype=float)
    
        for subject in toto["sample"]:

            z = df_zscore.loc[toto["sample"] == subject, numeric_cols].values.flatten()

            diff = z[:, None] - z[None, :]
            matrix = np.exp(-(diff ** 2))

            df_3d.loc[subject] = matrix
        print("Individual Matrices") 
        print(df_3d)

        print("Group Matrix")
        group_matrix= df_3d.groupby(level="covariance").mean()
        group_matrix = group_matrix.loc[numeric_cols,  numeric_cols]
        
        print(group_matrix)


if __name__ == "__main__":
    main()

    

            







    # out_dir.mkdir(parents=True, exist_ok=True)

    # data = {}

    # for tsv_file in sorted(input_dir.glob("*.tsv")):
    #     df = pd.read_csv(tsv_file, sep="\t")

    #     cols = [
    #         c for c in df.columns
    #         if c.endswith("_thickness") and "MeanThickness" not in c
    #     ]

    #     for _, row in df.iterrows():
    #         subject = str(row["sample"]).replace("_fs", "").replace("_ses-baseline", "")
    #         data[subject] = pd.to_numeric(row[cols], errors="coerce").to_dict()

    # X_df = pd.DataFrame.from_dict(data, orient="index").dropna(axis=1)
    # X_df.index.name = "subject"
    # X_df.to_csv(out_dir / "subject_region_thickness_matrix.csv")

    # X = X_df.to_numpy(dtype=float)

    # X_z = (X - np.mean(X, axis=0)) / np.std(X, axis=0)
    # group_covariance = np.corrcoef(X_z, rowvar=False)

    # regions = list(X_df.columns)
    # pd.DataFrame(
    #     group_covariance,
    #     index=regions,
    #     columns=regions,
    # ).to_csv(out_dir / "group_covariance.csv")

    # print(f"[OK] Saved outputs in: {out_dir}")
    # print(f"Subjects: {X_df.shape[0]}")
    # print(f"Regions: {X_df.shape[1]}")


