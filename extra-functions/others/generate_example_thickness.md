# generate_example_thickness.py

Generate synthetic DrawEM32 cortical thickness files, so the DrawEM32 -> MCRIBS
-> covariance/gradients chain can be run end to end without any real subject
data.

The synthetic atlas has 40 cortical parcels per hemisphere (`dwm_01` ...
`dwm_40`), which is what `mapping_example.json` expects.

To make the downstream covariance/gradient results interpretable, parcel
thickness is built from two smooth latent spatial axes (a dominant one and a
weaker one). Each subject gets one score per axis, shared by both hemispheres,
so that regions covary across subjects along those axes: the cohort covariance
matrix then has real structure and the first group gradient recovers the
dominant axis.

## Options

| Option | Meaning |
| --- | --- |
| `--out_dir` | Directory where the synthetic TSV files are written. |
| `--n_subjects` | Number of synthetic subjects (default 25). |
| `--id_column` | Subject identifier column (default `sample`). |
| `--seed` | Random seed for reproducibility (default 0). |
| `--with_area` | Also write matching surface-area files (for `area_weighted`). |

## Example

```bash
python generate_example_thickness.py --out_dir example_data --n_subjects 30
```
