# drawem32_to_mcribs.py

Convert DrawEM32 cortical thickness statistics into the MCRIBS (M-CRIB-S)
Desikan-Killiany parcellation with 68 cortical regions (34 per hemisphere).

The conversion is driven entirely by an external mapping file, so it contains no
study-specific data and can be reused by anyone with their own thickness files.
For every MCRIBS/Desikan-Killiany target region the mapping file lists the
DrawEM32 source regions that overlap it. The thickness of the target region is
obtained by aggregating (mean, or surface-area weighted mean) the thickness of
its DrawEM32 sources, subject by subject.

## Input files

Two tab-separated files (left and right hemisphere), each with:

* one identifier column (default name: `sample`)
* one column per DrawEM32 region, named `<hemi>_<label>_thickness`
  (for instance `lh_dwm_07_thickness`).

## Output files (written in `--out_dir`)

* `thickness_68_lh_stats.tsv` / `thickness_68_rh_stats.tsv` — thickness files in
  the MCRIBS 68-region convention, directly usable as input to
  `compute_covariance.py` or `compute_covariance_gradients.py`.
* `conversion_report_lh.tsv` / `conversion_report_rh.tsv` — per-region log of
  which sources were used and which input columns were left unused.

## Options

| Option | Meaning |
| --- | --- |
| `--mapping` | JSON file mapping each MCRIBS/DK region to its DrawEM32 source regions. Defaults to `mapping_drawem32_to_mcribs.json` next to the script. |
| `--out_dir` | Directory where the converted TSV files are written. |
| `--id_column` | Name of the subject identifier column (default `sample`). |
| `--aggregation` | `mean` or `area_weighted`; overrides the value stored in the mapping file. |
| `--area_stats LH_AREA RH_AREA` | Left/right TSV files holding DrawEM32 surface areas (columns `<hemi>_<label>_area`). Required for `area_weighted`. |
| `--strict` | Fail if a MCRIBS region has no available DrawEM32 source (default: warn and skip). |
| `-v` | Verbosity level. Default `WARNING`, or `INFO` when `-v` is given. |

## The mapping file

Edit `mapping_drawem32_to_mcribs.json`. For every MCRIBS/Desikan-Killiany region,
list the DrawEM32 source labels that overlap it. A source label is the thickness
column name without the hemisphere prefix (`lh_`/`rh_`) and without the
`_thickness` suffix. The same list is applied to both hemispheres. Several
DrawEM32 sources may feed one MCRIBS region; they are combined with the method in
`aggregation`. Keys starting with `_` are ignored (used for comments).

Validate this correspondence against your own atlas documentation before using
the results.

## Example

```bash
python drawem32_to_mcribs.py \
    lh_drawem32_thickness.tsv rh_drawem32_thickness.tsv \
    --mapping mapping_drawem32_to_mcribs.json \
    --out_dir mcribs_thickness_68
```
