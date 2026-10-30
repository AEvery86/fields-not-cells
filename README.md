# Fields, not cells

Code, derived data and the **fovplan** tool accompanying the manuscript
*Fields, not cells: design effects, false discoveries and field-number planning for spatial cell-composition studies*.

Imaging-based spatial studies profile a few fields of view (FOVs) per tissue. Cells within a field are correlated, so
the precision of a cell-type proportion is set by the number of fields, not the number of cells. This repository
quantifies the resulting design effects, shows the false discoveries produced by cell-level tests, and provides
`fovplan`, which turns pilot fields into the number of fields needed for a target precision.

## fovplan in one minute

```bash
pip install .                      # installs the `fovplan` command (numpy, pandas, scipy)
fovplan examples/fovplan_example_pilot.tsv --target 0.25 --relative --assurance 0.9
```

Input: one row per pilot field, one column per cell type (counts); non-numeric columns are ignored.
Output per cell type: pilot proportion, cells per field, design effect `D`, implied intraclass correlation, the
number of fields required (point plan and conservative plan), and the number implied if cells were independent.
We recommend `--assurance 0.9` for strongly compartmentalized tissue (see manuscript, Fig. 8d).

```python
import pandas as pd, fovplan
counts = pd.read_csv('pilot.tsv', sep='\t').select_dtypes('number')
print(fovplan.plan(counts, target=0.25, relative=True, assurance=0.9))
```

## Reproducing the analyses

```bash
pip install -r requirements.txt
bash run_all.sh            # runs every analysis and figure script in order (about 1-2 h on 2 cores)
pytest tests
```

| Script | Content | Main outputs |
|---|---|---|
| `a1_design_effect_atlas.py` | Design-effect atlas, exact validation | Fig. 2 |
| `a2_fovplan_validation.py`, `a2b_fields_required.py` | Held-out validation; fields required | Fig. 5 |
| `a3_false_positive_da.py` | False-positive rates, real contrasts | Fig. 3 |
| `a4_model_interval.py` | Interval procedures incl. quasi-binomial field model | Fig. 4, Table 2 |
| `a6_independent_check.py` | Independent recomputation of key results | Supplementary Methods S5 |
| `b1_simulation.py`, `b2_simulation_da.py` | Simulation with known truth | Fig. 6 |
| `c3_mouse_whole_section.py` | Whole-section truth, mouse kidney | Fig. 7 |
| `c0_hickey_assemble.py`, `c4_hickey.py`, `c5_pilot_sensitivity.py` | Human intestine whole-region truth, donor-level nulls | Fig. 8 |
| `a5_figures.py`, `b3_figures.py`, `c6_figure8.py`, `p3_concept.py`, `p3_supp_figures.py` | Figures | `figures/` |

Seeds are fixed in each script (20261027-20261062). Results shipped in `analysis/results/` let every figure be
redrawn without rerunning the simulations.

## Data

`analysis/data/` contains per-field cell-type counts and other derived tables computed from public sources; see
`DATA_SOURCES.md` for provenance, licences and how to reacquire the raw data. Derived lupus (CC BY 4.0) and mouse
(CC BY 4.0) tables are attributed to their original authors; the intestine extract derives from a public-domain
Dryad deposit. No restricted-access human data are included.

The intestine extract (`analysis/data/hickey_extract/`, about 22 MB) is not stored in the repository tree. Download
`fields-not-cells_v1.0.0_hickey_extract.zip` from the v1.0.0 release (also archived on Zenodo) and unzip it in the
parent directory of the repository before running `c0_hickey_assemble.py`, `c4_hickey.py` or `c5_pilot_sensitivity.py`.

## Citation

Please cite the manuscript and the archived release (see `CITATION.cff`; DOI assigned by Zenodo on release).

## Licence

Code: MIT (see `LICENSE`). Derived data keep the licences of their sources (`DATA_SOURCES.md`).
