# Fields, not cells

Code, derived data and the **fovplan** tool accompanying the manuscript
*Fields, not cells, determine precision and replication in spatial cell-composition studies*
(Xiao Lv, Yuanfang Zou, Liuyi Yang).

Imaging-based spatial studies profile a few fields of view (FOVs) per tissue. Cells within a field are correlated, so
the precision of a cell-type proportion is set by the number of fields, not the number of cells. This repository
quantifies the resulting design effects across platforms and organs, shows the false discoveries produced by
cell-level tests, evaluates practical design choices, and provides `fovplan`, which turns pilot fields into the number
of fields needed for a target precision.

## fovplan in one minute

```bash
pip install .                      # installs the `fovplan` command (numpy, pandas, scipy); fovplan 0.2.0
fovplan examples/fovplan_example_pilot.tsv --target 0.25 --relative --assurance 0.9
```

Input: one row per pilot field, one column per cell type (counts). Keep identifiers and areas out of the count table:
every numeric column is treated as a cell-type count. Output per cell type: pilot proportion, cells per field,
design effect `D`, implied intraclass correlation, the number of fields required (point plan and assurance plan),
and the number implied if cells were independent. A worked, verified example is described in `TUTORIAL.md`.

```python
import pandas as pd, fovplan
counts = pd.read_csv('pilot.tsv', sep='\t').select_dtypes('number')
print(fovplan.plan(counts, target=0.25, relative=True, assurance=0.9))
```

## Reproducing the analyses

```bash
pip install -r requirements.txt
bash run_all.sh            # original analyses and figures (about 1-2 h on 2 cores)
pytest tests
```

Scripts in `analysis/scripts/` are run from that directory.

| Script | Content |
|---|---|
| `a1_design_effect_atlas.py` | Design-effect atlas, exact validation |
| `a2_fovplan_validation.py`, `a2b_fields_required.py` | Held-out validation; fields required |
| `a3_false_positive_da.py` | False-positive rates, real contrasts |
| `a4_model_interval.py` | Interval procedures incl. quasi-binomial field model |
| `a6_independent_check.py` | Independent recomputation of key results |
| `a7_glmm_da.py`, `da_models.py`, `a7_summary.py` | Mixed models, beta-binomial, negative-binomial and propeller under the same null splits (v1.0.1) |
| `a7_export_for_R.py`, `../R_validation/` | Export of identical inputs and re-analysis with lme4, glmmTMB, MASS and speckle (v1.0.1) |
| `a8_planner_benchmark.py`, `a8_summary.py` | Benchmark of cell-count, bootstrap and fovplan planners (v1.0.1) |
| `a9_preferential.py`, `a10_planner_sensitivity.py`, `a11_field_size.py`, `a12_annotation_noise.py`, `a9_12_summary.py` | Design choices: field selection, precision target and pilot size, field size, annotation error (v1.0.1) |
| `b1_simulation.py`, `b2_simulation_da.py` | Simulation with known truth |
| `c3_mouse_whole_section.py` | Whole-section references, mouse kidney |
| `c0_hickey_assemble.py`, `c4_hickey.py`, `c5_pilot_sensitivity.py` | Human intestine whole-region references, donor-level nulls |
| `a5_figures.py`, `b3_figures.py`, `c6_figure8.py`, `p3_concept.py`, `p3_supp_figures.py` | Figures |

Seeds are fixed in each script (20261027-20261130). Summary tables of all analyses are in
`analysis/results/` and `analysis/results/supplementary_tables/`; the full long-format outputs of the v1.0.1 analyses
(A7-A12) and of the R re-analysis are attached to the v1.0.1 GitHub release as
`fields-not-cells_v1.0.1_results.zip`.

## Data

`analysis/data/` contains per-field cell-type counts and other derived tables computed from public sources; see
`DATA_SOURCES.md` for provenance, licences and how to reacquire the raw data. Derived lupus (CC BY 4.0) and mouse
(CC BY 4.0) tables are attributed to their original authors; the intestine extract derives from a CC0 Dryad deposit.
No restricted-access human data are included.

The intestine extract (`analysis/data/hickey_extract/`, about 22 MB) is not stored in the repository tree. Download
`fields-not-cells_v1.0.0_hickey_extract.zip` from the assets of the v1.0.0 GitHub release and unzip it in the parent
directory of the repository before running `c0_hickey_assemble.py`, which writes `analysis/data/hickey_7sites.pkl`
used by the intestinal analyses. (The Zenodo archive contains the repository source only, not this release asset.)

## Citation

Please cite the manuscript and the archived release (see `CITATION.cff`; Zenodo concept DOI
https://doi.org/10.5281/zenodo.23797220).

## Licence

Code: MIT (see `LICENSE`). Derived data keep the licences of their sources (`DATA_SOURCES.md`).
