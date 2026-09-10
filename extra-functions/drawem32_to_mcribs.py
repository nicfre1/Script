#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert DrawEM32 left/right thickness TSV files to the MCRIBS Desikan-Killiany
parcellation (68 regions), using an external region mapping file.

See others/drawem32_to_mcribs.md for details.

Example:
drawem32_to_mcribs.py lh_stats.tsv rh_stats.tsv --mapping mapping_drawem32_to_mcribs.json
"""

import argparse
import json
import logging
from pathlib import Path

import pandas as pd


THICKNESS_SUFFIX = "_thickness"
AREA_SUFFIX = "_area"
HEMI_PREFIXES = ("lh_", "rh_")


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    p.add_argument(
        "i_stats",
        nargs=2,
        metavar=("LH_STATS", "RH_STATS"),
        help="Left and right hemisphere DrawEM32 thickness TSV files.",
    )

    p.add_argument(
        "--mapping",
        default=str(Path(__file__).with_name("mapping_drawem32_to_mcribs.json")),
        help="JSON file mapping each MCRIBS region to its DrawEM32 source regions.",
    )

    p.add_argument(
        "--out_dir",
        default="mcribs_thickness_68",
        help="Directory where the converted TSV files are written.",
    )

    p.add_argument(
        "--id_column",
        default="sample",
        help="Name of the subject identifier column.",
    )

    p.add_argument(
        "--aggregation",
        choices=["mean", "area_weighted"],
        default=None,
        help="How to combine several DrawEM32 sources; overrides the mapping file.",
    )

    p.add_argument(
        "--area_stats",
        nargs=2,
        metavar=("LH_AREA", "RH_AREA"),
        help="Left/right DrawEM32 surface area TSV files, for area_weighted.",
    )

    p.add_argument(
        "--strict",
        action="store_true",
        help="Fail if a MCRIBS region has no available DrawEM32 source.",
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


def region_basename(column):
    name = column
    for prefix in HEMI_PREFIXES:
        if name.startswith(prefix):
            name = name[len(prefix):]
            break
    for suffix in (THICKNESS_SUFFIX, AREA_SUFFIX):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def load_mapping(path):
    with open(path, "r", encoding="utf-8") as handle:
        mapping = json.load(handle)

    if "regions" not in mapping or not isinstance(mapping["regions"], dict):
        raise ValueError(
            f"Mapping file '{path}' must contain a 'regions' object "
            "(MCRIBS region -> list of DrawEM32 source labels)."
        )

    regions = {
        target: [src for src in sources if not str(src).startswith("_")]
        for target, sources in mapping["regions"].items()
        if not str(target).startswith("_")
    }

    empty = sorted(t for t, s in regions.items() if not s)
    if empty:
        logging.warning(
            "%d MCRIBS region(s) have no source label in the mapping file: %s",
            len(empty),
            ", ".join(empty),
        )

    mapping["regions"] = regions
    mapping.setdefault("aggregation", "mean")
    return mapping


def _read_stats(path, id_column, suffix):
    df = pd.read_csv(path, sep="\t")
    if id_column not in df.columns:
        raise ValueError(
            f"Identifier column '{id_column}' not found in '{path}'. "
            f"Available columns start with: {list(df.columns[:5])} ..."
        )
    value_cols = [c for c in df.columns if c.endswith(suffix)]
    if not value_cols:
        raise ValueError(f"No '*{suffix}' column found in '{path}'.")
    return df, value_cols


def convert_hemisphere(
    df,
    thickness_cols,
    id_column,
    hemi,
    mapping,
    aggregation,
    area_df=None,
    area_cols=None,
):
    available = {region_basename(c): c for c in thickness_cols}
    area_available = {}
    if area_df is not None and area_cols is not None:
        area_available = {region_basename(c): c for c in area_cols}

    out = pd.DataFrame({id_column: df[id_column].to_numpy()})
    report_rows = []
    used_columns = set()

    for target, sources in mapping["regions"].items():
        matched = [s for s in sources if s in available]
        missing = [s for s in sources if s not in available]

        if not matched:
            report_rows.append(
                {
                    "mcribs_region": target,
                    "status": "NO_SOURCE",
                    "sources_used": "",
                    "sources_missing": ";".join(missing),
                }
            )
            continue

        cols = [available[s] for s in matched]
        used_columns.update(cols)
        values = df[cols].astype(float)

        if aggregation == "area_weighted":
            weight_cols = [area_available.get(s) for s in matched]
            if any(w is None for w in weight_cols):
                raise ValueError(
                    f"area_weighted aggregation requested but surface area is "
                    f"missing for some sources of '{hemi}_{target}': {matched}"
                )
            weights = area_df[weight_cols].astype(float).to_numpy()
            weighted = (values.to_numpy() * weights).sum(axis=1)
            total = weights.sum(axis=1)
            aggregated = weighted / total
        else:
            aggregated = values.mean(axis=1).to_numpy()

        out[f"{hemi}_{target}{THICKNESS_SUFFIX}"] = aggregated
        report_rows.append(
            {
                "mcribs_region": target,
                "status": "OK",
                "sources_used": ";".join(matched),
                "sources_missing": ";".join(missing),
            }
        )

    unused = sorted(c for c in thickness_cols if c not in used_columns)
    for col in unused:
        report_rows.append(
            {
                "mcribs_region": "",
                "status": "UNUSED_INPUT",
                "sources_used": col,
                "sources_missing": "",
            }
        )

    report = pd.DataFrame(
        report_rows,
        columns=[
            "mcribs_region",
            "status",
            "sources_used",
            "sources_missing",
        ],
    )
    return out, report


def main():
    args = _build_arg_parser().parse_args()
    logging.getLogger().setLevel(logging.getLevelName(args.verbose))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mapping = load_mapping(args.mapping)
    aggregation = args.aggregation or mapping.get("aggregation", "mean")
    logging.info("Aggregation method: %s", aggregation)

    area_frames = {"lh": (None, None), "rh": (None, None)}
    if aggregation == "area_weighted":
        if not args.area_stats:
            raise SystemExit(
                "area_weighted aggregation requires --area_stats LH_AREA RH_AREA"
            )
        for hemi, path in zip(("lh", "rh"), args.area_stats):
            area_df, area_cols = _read_stats(path, args.id_column, AREA_SUFFIX)
            area_frames[hemi] = (area_df, area_cols)

    n_regions = 0
    for hemi, stats_path in zip(("lh", "rh"), args.i_stats):
        logging.info("Convert %s hemisphere from %s", hemi, stats_path)
        df, thickness_cols = _read_stats(stats_path, args.id_column, THICKNESS_SUFFIX)
        area_df, area_cols = area_frames[hemi]

        converted, report = convert_hemisphere(
            df=df,
            thickness_cols=thickness_cols,
            id_column=args.id_column,
            hemi=hemi,
            mapping=mapping,
            aggregation=aggregation,
            area_df=area_df,
            area_cols=area_cols,
        )

        missing = report.loc[report["status"] == "NO_SOURCE", "mcribs_region"].tolist()
        if missing:
            message = (
                f"{hemi} hemisphere: {len(missing)} MCRIBS region(s) without any "
                f"available DrawEM32 source: {', '.join(missing)}"
            )
            if args.strict:
                raise SystemExit(message)
            logging.warning(message)

        out_path = out_dir / f"thickness_68_{hemi}_stats.tsv"
        report_path = out_dir / f"conversion_report_{hemi}.tsv"
        converted.to_csv(out_path, sep="\t", index=False)
        report.to_csv(report_path, sep="\t", index=False)

        produced = converted.shape[1] - 1
        n_regions += produced
        logging.info("Wrote %s (%d regions)", out_path, produced)
        logging.info("Wrote %s", report_path)

    logging.warning(
        "Done. %d cortical regions written across both hemispheres "
        "(expected 68 for a complete Desikan-Killiany parcellation).",
        n_regions,
    )


if __name__ == "__main__":
    main()
