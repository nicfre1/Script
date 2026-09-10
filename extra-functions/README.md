# extra-functions

Helper scripts to move cortical thickness statistics from the **DrawEM32**
parcellation to the **MCRIBS** (M-CRIB-S) Desikan-Killiany parcellation with
**68 cortical regions** (34 per hemisphere), then run the validated structural
covariance and gradient analysis on them.

No study data is shipped here. Every script is data-agnostic: anyone can run the
chain on their own thickness files and obtain their own results. Synthetic files
are provided only so the pipeline can be exercised end to end.

Tested with Python 3.11 and the packages in `../requirements.txt`. The cortical
surface render of the group gradients additionally needs a working VTK stack
(`enigmatoolbox`, `brainspace[plotting]`, `nibabel`, `vtk`); without it the
pipeline still completes and writes a matplotlib fallback figure.

## Contents

| File | Purpose |
| --- | --- |
| `drawem32_to_mcribs.py` | Convert DrawEM32 thickness TSVs to the MCRIBS 68-region convention, driven by an external mapping file. |
| `mapping_drawem32_to_mcribs.json` | Mapping **template**: the 34 Desikan-Killiany target regions with empty source lists to fill in. |
| `compute_covariance_gradients.py` | Individualized covariance matrices, group covariance matrix, group + Procrustes-aligned individual gradients, figures. Data-free port of the validated reference scripts. |
| `plot_group_gradients_brain.py` | Standalone cortical-surface render of the 4 group gradients (aparc_fsa5). Run with a VTK-capable interpreter. |
| `generate_example_thickness.py` | Write synthetic DrawEM32 thickness files so the chain can be run without real data. |
| `mapping_example.json` | Example mapping matching the synthetic files (demo only, not anatomically valid). |
| `others/` | Detailed notes for each script (kept out of the code headers). |

The method (cohort z-score, `np.corrcoef` group matrix, `GradientMaps` +
Procrustes, rescale to `[-1, 1]`) is summarised in [`others/README.md`](others/README.md)
and unchanged from the validated reference scripts.

## 1. Provide the DrawEM32 -> MCRIBS correspondence

Edit `mapping_drawem32_to_mcribs.json`. For every MCRIBS/Desikan-Killiany region,
list the DrawEM32 source labels that overlap it (column name minus the
`lh_`/`rh_` prefix and the `_thickness` suffix). Same list for both hemispheres.
Validate it against your atlas documentation before trusting the results.

## 2. Convert

```bash
python drawem32_to_mcribs.py \
    lh_drawem32_thickness.tsv rh_drawem32_thickness.tsv \
    --mapping mapping_drawem32_to_mcribs.json \
    --out_dir mcribs_thickness_68
```

Produces `thickness_68_lh_stats.tsv`, `thickness_68_rh_stats.tsv` and a
`conversion_report_{lh,rh}.tsv`.

## 3. Covariance and gradients

```bash
python compute_covariance_gradients.py \
    mcribs_thickness_68/thickness_68_lh_stats.tsv \
    mcribs_thickness_68/thickness_68_rh_stats.tsv \
    --out_dir mcribs_thickness_68/results
```

Or, if you already have one `<subject>_thickness.csv` per subject (columns
`Region`, `ThickAvg`, 68 Desikan-Killiany rows):

```bash
python compute_covariance_gradients.py --thickness_dir dk68_stats --out_dir results
```

Outputs:

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
available, otherwise a matplotlib fallback. To produce (or re-produce) the real
render with a dedicated interpreter:

```bash
/path/to/vtk-python plot_group_gradients_brain.py \
    --gradient-csv results/group_covariance/group_gradients_rescaled_minus1_1.csv \
    --output results/group_covariance/group_gradients_brain.png
```

## Run the full chain on synthetic data

```bash
python generate_example_thickness.py --out_dir example_data --n_subjects 30
python drawem32_to_mcribs.py \
    example_data/drawem32_thickness_lh_stats.tsv \
    example_data/drawem32_thickness_rh_stats.tsv \
    --mapping mapping_example.json --out_dir example_out
python compute_covariance_gradients.py \
    example_out/thickness_68_lh_stats.tsv \
    example_out/thickness_68_rh_stats.tsv \
    --out_dir example_out/results
```
