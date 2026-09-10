#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate synthetic DrawEM32 thickness TSV files to test the pipeline without
real data.

See others/generate_example_thickness.md for details.

Example:
generate_example_thickness.py --out_dir example_data --n_subjects 25
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

N_PARCELS = 40


def _latent_axes(n_parcels):
    position = np.linspace(0.0, 1.0, n_parcels)
    axis_1 = np.cos(np.pi * position)
    axis_2 = np.sin(2.0 * np.pi * position)
    return axis_1, axis_2


def _hemisphere_frame(rng, subjects, id_column, hemi, suffix,
                      baseline, scale_1, scale_2, noise_sd,
                      score_1, score_2, subject_offset):
    axis_1, axis_2 = _latent_axes(N_PARCELS)
    n_subjects = len(subjects)

    data = {id_column: subjects}
    for idx in range(N_PARCELS):
        label = f"dwm_{idx + 1:02d}"
        values = (
            baseline
            + subject_offset
            + scale_1 * axis_1[idx] * score_1
            + scale_2 * axis_2[idx] * score_2
            + rng.normal(0.0, noise_sd, size=n_subjects)
        )
        data[f"{hemi}_{label}{suffix}"] = np.round(values, 4)
    return pd.DataFrame(data)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--out_dir", default="example_data",
                   help="Directory where the synthetic TSV files are written.")
    p.add_argument("--n_subjects", type=int, default=25,
                   help="Number of synthetic subjects.")
    p.add_argument("--id_column", default="sample",
                   help="Subject identifier column.")
    p.add_argument("--seed", type=int, default=0,
                   help="Random seed for reproducibility.")
    p.add_argument("--with_area", action="store_true",
                   help="Also write matching surface-area files.")
    return p


def main():
    args = _build_arg_parser().parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    subjects = [f"sub-{i + 1:03d}" for i in range(args.n_subjects)]
    n = args.n_subjects

    score_1 = rng.normal(1.0, 0.25, size=n)
    score_2 = rng.normal(1.0, 0.35, size=n)
    subject_offset = rng.normal(0.0, 0.12, size=n)

    for hemi in ("lh", "rh"):
        thickness = _hemisphere_frame(
            rng, subjects, args.id_column, hemi, "_thickness",
            baseline=2.5, scale_1=0.55, scale_2=0.22, noise_sd=0.08,
            score_1=score_1, score_2=score_2, subject_offset=subject_offset,
        )
        thickness_path = out_dir / f"drawem32_thickness_{hemi}_stats.tsv"
        thickness.to_csv(thickness_path, sep="\t", index=False)
        print(f"Wrote {thickness_path}  ({thickness.shape[0]} subjects, "
              f"{thickness.shape[1] - 1} parcels)")

        if args.with_area:
            area = _hemisphere_frame(
                rng, subjects, args.id_column, hemi, "_area",
                baseline=900.0, scale_1=260.0, scale_2=110.0, noise_sd=40.0,
                score_1=score_1, score_2=score_2,
                subject_offset=subject_offset * 200.0,
            )
            area_path = out_dir / f"drawem32_area_{hemi}_stats.tsv"
            area.to_csv(area_path, sep="\t", index=False)
            print(f"Wrote {area_path}")


if __name__ == "__main__":
    main()
