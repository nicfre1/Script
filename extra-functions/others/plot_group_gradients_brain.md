# plot_group_gradients_brain.py

Render the first four group gradients on the cortical surface (aparc_fsa5).

Data-free port of the validated `render_dk68_myrna_style.py`: each of the 68
Desikan-Killiany gradient values is painted on the fsaverage5 surface with
ENIGMA Toolbox's `parcel_to_surface(..., "aparc_fsa5")` and BrainSpace's
`plot_hemispheres` (colormap `viridis_r`, colorbar, four rows `Grad1` .. `Grad4`).

## Requirements

Needs a working VTK stack: `enigmatoolbox`, `brainspace[plotting]`, `nibabel`,
`vtk` (see `../requirements-brain.txt`). On some interpreters `import vtk` crashes
the process, so this lives in its own script and
`compute_covariance_gradients.py` calls it as a subprocess. A dedicated
environment (Python 3.10, VTK 9.3.x) is recommended.

## Input

A gradient CSV with columns `Gradient_1` .. `Gradient_4` in the canonical
Desikan-Killiany order (lh 34 regions, then rh 34 regions), ideally already
rescaled to `[-1, 1]` (`group_gradients_rescaled_minus1_1.csv`).

## Options

| Option | Meaning |
| --- | --- |
| `--gradient-csv` | Path to the gradient CSV (required). |
| `--output` | Path of the PNG to write (required). |
| `--gradient-columns` | Columns to render (default `Gradient_1 Gradient_2 Gradient_3 Gradient_4`). |
| `--label-text` | Row labels (default `Grad1 Grad2 Grad3 Grad4`). |

## Example

```bash
/path/to/vtk-python plot_group_gradients_brain.py \
    --gradient-csv results/group_covariance/group_gradients_rescaled_minus1_1.csv \
    --output results/group_covariance/group_gradients_brain.png
```
