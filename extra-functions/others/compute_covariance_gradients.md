# compute_covariance_gradients.py

Individualized (cohort-referenced) covariance matrices, group covariance matrix,
group gradients and Procrustes-aligned individual gradients, from MCRIBS
thickness in the Desikan-Killiany parcellation (68 cortical regions, 34 per
hemisphere, no left/right averaging).

This is a data-free port of the validated reference scripts
(`individualized_covariance.py` + `group_and_individual_gradients.py`). Every
equation and BrainSpace call is kept:

```
cohort z-score per region:   X_z = (X - X.mean(0)) / X.std(0)
individual matrix:           M[i, j] = 1 / exp((z_i - z_j) ** 2)
group covariance matrix:     np.corrcoef(X_z, rowvar=False)
group gradients:             GradientMaps(n_components=10, random_state=0)
individual gradients:        GradientMaps(..., alignment="procrustes")
                             .fit(individual_covariances,
                                  reference=gm_group.gradients_)
gradient rescaling:          2 * (x - min) / (max - min) - 1   (first 4)
```

Run it on your own data: nothing here is study specific.

## Inputs

Either

```bash
compute_covariance_gradients.py LH_STATS.tsv RH_STATS.tsv
```

where each file is the wide output of `drawem32_to_mcribs.py` (`sample` column +
one `<hemi>_<dkname>_thickness` column per region), or

```bash
compute_covariance_gradients.py --thickness_dir DIR
```

where `DIR` holds one `<subject>_thickness.csv` per subject with columns
`Region` and `ThickAvg` (68 Desikan-Killiany rows, lh then rh).

## Options

| Option | Meaning |
| --- | --- |
| `--out_dir` | Directory where CSV files and figures are written. |
| `--id_column` | Subject identifier column in the wide TSV files (default `sample`). |
| `--n_components` | Number of gradient components (default 10). |
| `--no_brain_plot` | Skip the cortical-surface render of the group gradients. |
| `-v` | Verbosity level. Default `WARNING`, or `INFO` when `-v` is given. |

## Outputs (in `--out_dir`)

```
individualized_covariance/csv/<subject>_individualized_covariance.csv
individualized_covariance/png/<subject>_individualized_covariance.png
group_covariance/group_covariance_matrix.csv
group_covariance/group_gradients.csv
group_covariance/group_gradients_rescaled_minus1_1.csv
group_covariance/group_lambdas.csv
group_covariance/group_lambdas.png
group_covariance/group_gradients_brain.png
individual_gradients4/<subject>_gradients_aligned.csv
individual_gradients4/<subject>_gradients_aligned_rescaled_minus1_1.csv
```

`group_gradients_brain.png` is the cortical-surface render when a VTK stack is
available (see `plot_group_gradients_brain.md`), otherwise a matplotlib fallback.

## Region order

Regions are ordered lh first, then rh, each in the canonical FreeSurfer
Desikan-Killiany sequence (`ctx-lh-bankssts` ... `ctx-lh-insula`,
`ctx-rh-bankssts` ... `ctx-rh-insula`), which is the order ENIGMA's `aparc_fsa5`
surface mapping expects.

## Example

```bash
python compute_covariance_gradients.py \
    mcribs_thickness_68/thickness_68_lh_stats.tsv \
    mcribs_thickness_68/thickness_68_rh_stats.tsv \
    --out_dir mcribs_thickness_68/results
```
