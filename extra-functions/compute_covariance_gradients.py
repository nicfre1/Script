#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute individual and group covariance matrices, group and individual
gradients, and lambdas from MCRIBS thickness files (68 Desikan-Killiany regions).

See others/compute_covariance_gradients.md for details.

Example:
compute_covariance_gradients.py lh_stats.tsv rh_stats.tsv
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from brainspace.gradient import GradientMaps

THICKNESS_SUFFIX = "_thickness"
HEMI_PREFIXES = ("lh_", "rh_")

DK_REGIONS = [
    "bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", "cuneus",
    "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal",
    "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", "lingual",
    "medialorbitofrontal", "middletemporal", "parahippocampal", "paracentral",
    "parsopercularis", "parsorbitalis", "parstriangularis", "pericalcarine",
    "postcentral", "posteriorcingulate", "precentral", "precuneus",
    "rostralanteriorcingulate", "rostralmiddlefrontal", "superiorfrontal",
    "superiorparietal", "superiortemporal", "supramarginal", "frontalpole",
    "temporalpole", "transversetemporal", "insula",
]


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument(
        "i_stats",
        nargs="*",
        metavar="STATS",
        help="Left and right MCRIBS thickness TSV files. Omit with --thickness_dir.",
    )
    p.add_argument(
        "--thickness_dir",
        type=Path,
        default=None,
        help="Directory of per-subject <subject>_thickness.csv (Region, ThickAvg).",
    )
    p.add_argument(
        "--out_dir",
        default="results_mcribs",
        help="Directory where CSV files and figures are written.",
    )
    p.add_argument(
        "--id_column",
        default="sample",
        help="Subject identifier column in the wide TSV files.",
    )
    p.add_argument(
        "--n_components",
        type=int,
        default=10,
        help="Number of gradient components.",
    )
    p.add_argument(
        "--no_brain_plot",
        action="store_true",
        help="Skip the cortical-surface render of the group gradients.",
    )
    p.add_argument(
        "-v",
        default="WARNING",
        const="INFO",
        nargs="?",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        dest="verbose",
        help="Verbosity level. Default: WARNING, or INFO when -v is given.",
    )
    return p


def _region_basename(column):
    name = column
    for prefix in HEMI_PREFIXES:
        if name.startswith(prefix):
            name = name[len(prefix):]
            break
    if name.endswith(THICKNESS_SUFFIX):
        name = name[: -len(THICKNESS_SUFFIX)]
    return name


def load_from_wide_tsv(lh_path, rh_path, id_column):
    lh = pd.read_csv(lh_path, sep="\t")
    rh = pd.read_csv(rh_path, sep="\t")
    for name, df in ((lh_path, lh), (rh_path, rh)):
        if id_column not in df.columns:
            raise ValueError(f"Column '{id_column}' not found in '{name}'.")

    merged = pd.merge(lh, rh, on=id_column)
    if merged.empty:
        raise ValueError(
            f"No subject in common between '{lh_path}' and '{rh_path}' "
            f"(merged on '{id_column}')."
        )

    lh_cols = [c for c in lh.columns if c.endswith(THICKNESS_SUFFIX)]
    rh_cols = [c for c in rh.columns if c.endswith(THICKNESS_SUFFIX)]
    lh_map = {_region_basename(c): c for c in lh_cols}
    rh_map = {_region_basename(c): c for c in rh_cols}

    ordered_lh = [r for r in DK_REGIONS if r in lh_map] + \
                 [r for r in lh_map if r not in DK_REGIONS]
    ordered_rh = [r for r in DK_REGIONS if r in rh_map] + \
                 [r for r in rh_map if r not in DK_REGIONS]

    columns = [lh_map[r] for r in ordered_lh] + [rh_map[r] for r in ordered_rh]
    regions = [f"ctx-lh-{r}" for r in ordered_lh] + \
              [f"ctx-rh-{r}" for r in ordered_rh]

    X = merged[columns].to_numpy(dtype=float)
    subjects = merged[id_column].astype(str).tolist()
    return X, regions, subjects


def load_from_thickness_dir(thickness_dir):
    files = sorted(
        p for p in thickness_dir.glob("*_thickness.csv")
        if not p.name.startswith("._")
    )
    if not files:
        raise FileNotFoundError(f"No *_thickness.csv in {thickness_dir}")

    regions = None
    rows, subjects = [], []
    for path in files:
        frame = pd.read_csv(path)
        current = frame["Region"].astype(str).tolist()
        if regions is None:
            regions = current
        elif current != regions:
            raise ValueError(f"{path.name}: region order differs from {files[0].name}")
        rows.append(frame["ThickAvg"].to_numpy(dtype=float))
        subjects.append(path.name.removesuffix("_thickness.csv"))

    return np.vstack(rows), regions, subjects


def individualized_covariance_from_z(z_subject):
    diff = z_subject[:, None] - z_subject[None, :]
    return 1 / np.exp(diff**2)


def rescale_to_minus1_1(x):
    return 2 * (x - np.min(x)) / (np.max(x) - np.min(x)) - 1


def save_individual_matrix_figure(matrix, subject_id, n_subjects, path):
    figure, axis = plt.subplots(figsize=(10, 10))
    image = axis.imshow(
        matrix, cmap="viridis", vmin=0, vmax=1,
        interpolation="nearest", aspect="equal",
    )
    figure.colorbar(image, ax=axis, label="Individualized covariance")
    axis.set_title(f"{subject_id} individualized covariance (n_cohort={n_subjects})")
    axis.set_xticks([])
    axis.set_yticks([])
    figure.tight_layout()
    figure.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def render_group_gradients_brain(rescaled_csv, out_png):
    script = Path(__file__).with_name("plot_group_gradients_brain.py")
    if script.exists():
        try:
            completed = subprocess.run(
                [sys.executable, str(script),
                 "--gradient-csv", str(rescaled_csv),
                 "--output", str(out_png)],
                capture_output=True, text=True, timeout=600,
            )
            if completed.returncode == 0 and out_png.exists():
                logging.warning("Saved cortical-surface render to %s", out_png)
                return
            logging.warning(
                "Cortical-surface render failed (exit %s). Last output: %s",
                completed.returncode,
                (completed.stderr or completed.stdout or "").strip().splitlines()[-1:]
                or "<none>",
            )
        except (subprocess.TimeoutExpired, OSError) as exc:
            logging.warning("Cortical-surface render could not run (%s).", exc)

    rescaled_df = pd.read_csv(rescaled_csv, index_col=0)
    columns = [c for c in rescaled_df.columns if c.startswith("Gradient_")][:4]
    fig, axes = plt.subplots(1, len(columns), figsize=(14, 9), sharey=True)
    if len(columns) == 1:
        axes = [axes]
    for ax, column in zip(axes, columns):
        im = ax.imshow(
            rescaled_df[[column]].to_numpy(dtype=float),
            cmap="viridis_r", vmin=-1, vmax=1,
            aspect="auto", interpolation="nearest",
        )
        ax.set_title(column.replace("_", " "))
        ax.set_xticks([])
    axes[0].set_ylabel("Region (lh then rh, Desikan-Killiany order)")
    fig.colorbar(im, ax=axes, label="Gradient value (rescaled to [-1, 1])",
                 shrink=0.6)
    fig.suptitle("Group cortical gradients")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logging.warning("Saved fallback group-gradient figure to %s", out_png)


def save_scree_figure(lambdas, path):
    x = range(1, len(lambdas) + 1)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(x, lambdas)
    ax.plot(x, lambdas)
    ax.set_xlabel("Component")
    ax.set_ylabel("Eigenvalue")
    ax.set_title("Group Gradient Eigenvalues")
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    args = _build_arg_parser().parse_args()
    logging.getLogger().setLevel(logging.getLevelName(args.verbose))

    if args.thickness_dir is not None:
        X, regions, subjects = load_from_thickness_dir(args.thickness_dir)
    elif len(args.i_stats) == 2:
        X, regions, subjects = load_from_wide_tsv(
            args.i_stats[0], args.i_stats[1], args.id_column
        )
    else:
        raise SystemExit(
            "Provide either two wide TSV files (LH RH) or --thickness_dir DIR."
        )

    n_subjects, n_regions = X.shape
    if n_subjects < 2:
        raise SystemExit(
            "Cohort z-scoring needs at least 2 subjects "
            f"(got {n_subjects})."
        )
    logging.warning("Cohort: %d subjects, %d regions", n_subjects, n_regions)

    out_dir = Path(args.out_dir)
    indiv_cov_csv_dir = out_dir / "individualized_covariance" / "csv"
    indiv_cov_png_dir = out_dir / "individualized_covariance" / "png"
    group_dir = out_dir / "group_covariance"
    indiv_grad_dir = out_dir / "individual_gradients4"
    for directory in (
        indiv_cov_csv_dir, indiv_cov_png_dir, group_dir, indiv_grad_dir
    ):
        directory.mkdir(parents=True, exist_ok=True)

    X_z = (X - np.mean(X, axis=0)) / np.std(X, axis=0)
    logging.info(
        "z-score check: region means ~0 (%s), region stds ~1 (%s)",
        np.allclose(X_z.mean(axis=0), 0, atol=1e-10),
        np.allclose(X_z.std(axis=0), 1, atol=1e-10),
    )

    individual_covariances = []
    for index, subject_id in enumerate(subjects):
        matrix = individualized_covariance_from_z(X_z[index])
        individual_covariances.append(matrix)
        pd.DataFrame(matrix, index=regions, columns=regions).to_csv(
            indiv_cov_csv_dir / f"{subject_id}_individualized_covariance.csv"
        )
        save_individual_matrix_figure(
            matrix, subject_id, n_subjects,
            indiv_cov_png_dir / f"{subject_id}_individualized_covariance.png",
        )
    logging.warning(
        "Wrote %d individualized covariance matrices (csv + png)", n_subjects
    )

    covariance_group = np.corrcoef(X_z, rowvar=False)
    pd.DataFrame(covariance_group, index=regions, columns=regions).to_csv(
        group_dir / "group_covariance_matrix.csv"
    )

    gm_group = GradientMaps(n_components=args.n_components, random_state=0)
    gm_group.fit(covariance_group)

    grad_columns = [
        f"Gradient_{i + 1}" for i in range(gm_group.gradients_.shape[1])
    ]
    group_grad_df = pd.DataFrame(
        gm_group.gradients_, index=regions, columns=grad_columns
    )
    group_grad_df.to_csv(group_dir / "group_gradients.csv")

    group_rescaled = pd.DataFrame(
        {
            f"Gradient_{i + 1}": rescale_to_minus1_1(gm_group.gradients_[:, i])
            for i in range(min(4, gm_group.gradients_.shape[1]))
        },
        index=regions,
    )
    group_rescaled_csv = group_dir / "group_gradients_rescaled_minus1_1.csv"
    group_rescaled.to_csv(group_rescaled_csv)

    pd.DataFrame(
        gm_group.lambdas_,
        index=[f"Gradient_{i + 1}" for i in range(len(gm_group.lambdas_))],
        columns=["lambda"],
    ).to_csv(group_dir / "group_lambdas.csv")
    save_scree_figure(gm_group.lambdas_, group_dir / "group_lambdas.png")

    gm_indiv = GradientMaps(
        n_components=args.n_components, random_state=0, alignment="procrustes"
    )
    gm_indiv.fit(individual_covariances, reference=gm_group.gradients_)

    for index, subject_id in enumerate(subjects):
        aligned = gm_indiv.aligned_[index]
        columns = [f"Gradient_{j + 1}" for j in range(aligned.shape[1])]
        pd.DataFrame(aligned, index=regions, columns=columns).to_csv(
            indiv_grad_dir / f"{subject_id}_gradients_aligned.csv"
        )
        pd.DataFrame(
            {
                f"Gradient_{j + 1}": rescale_to_minus1_1(aligned[:, j])
                for j in range(min(4, aligned.shape[1]))
            },
            index=regions,
        ).to_csv(
            indiv_grad_dir / f"{subject_id}_gradients_aligned_rescaled_minus1_1.csv"
        )
    logging.warning("Wrote aligned individual gradients (raw + rescaled)")

    if not args.no_brain_plot:
        render_group_gradients_brain(
            group_rescaled_csv, group_dir / "group_gradients_brain.png"
        )

    logging.warning("Done. Outputs written to %s", out_dir)


if __name__ == "__main__":
    main()
