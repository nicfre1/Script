#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Render the first four group gradients on the cortical surface (aparc_fsa5).

Needs a VTK stack (see others/plot_group_gradients_brain.md and
requirements-brain.txt); run with a VTK-capable interpreter.

Example:
plot_group_gradients_brain.py --gradient-csv group_gradients_rescaled_minus1_1.csv --output group_gradients_brain.png
"""

import argparse
import os
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from brainspace.mesh.mesh_creation import build_polydata
from brainspace.plotting import plot_hemispheres
from enigmatoolbox.utils.parcellation import parcel_to_surface
import enigmatoolbox


def load_fsa5_surfaces():
    surf_dir = os.path.join(
        os.path.dirname(enigmatoolbox.__file__), "datasets", "surfaces"
    )
    lh = nib.load(os.path.join(surf_dir, "fsa5_lh.surf.gii"))
    rh = nib.load(os.path.join(surf_dir, "fsa5_rh.surf.gii"))
    surf_lh = build_polydata(
        np.asarray(lh.darrays[0].data, dtype=float),
        np.asarray(lh.darrays[1].data, dtype=np.int32),
    )
    surf_rh = build_polydata(
        np.asarray(rh.darrays[0].data, dtype=float),
        np.asarray(rh.darrays[1].data, dtype=np.int32),
    )
    return surf_lh, surf_rh


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--gradient-csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--gradient-columns", nargs="+",
        default=["Gradient_1", "Gradient_2", "Gradient_3", "Gradient_4"],
    )
    parser.add_argument(
        "--label-text", nargs="+",
        default=["Grad1", "Grad2", "Grad3", "Grad4"],
    )
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.gradient_csv, index_col=0)

    array_list = [
        parcel_to_surface(frame[column].to_numpy(dtype=float), "aparc_fsa5")
        for column in args.gradient_columns
    ]

    surf_lh, surf_rh = load_fsa5_surfaces()
    plot_hemispheres(
        surf_lh,
        surf_rh,
        array_name=array_list,
        size=(1200, 400),
        cmap="viridis_r",
        color_bar=True,
        label_text=args.label_text,
        zoom=1.55,
        interactive=False,
        offscreen=True,
        screenshot=True,
        filename=str(args.output),
        transparent_bg=False,
        background=(1, 1, 1),
    )
    print("Saved brain plot to:", args.output)


if __name__ == "__main__":
    main()
