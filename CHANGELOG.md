# Changelog

## v1.0.1 (2026-11-02)

- Authors: Yuanfang Zou added; affiliations added to `.zenodo.json` and `CITATION.cff`.
- fovplan 0.2.0: the smallest number of fields is found by exponential bracketing and bisection. In 3,000 random
  parameter settings it returned exactly the same plans as the exhaustive scan of 0.1.0. A fifth unit test checks
  this agreement. `analysis/scripts/fovplan.py` is the same code.
- `TUTORIAL.md` and `examples/run_immune_example.py`: a verified five-field pilot example (GSE325587, seed 20261071)
  with expected outputs.
- New analyses: mixed models and differential-abundance tools (A7) with an R reference-package re-analysis
  (`analysis/R_validation/`); planner benchmark (A8); design-choice analyses (A9-A12).
- `analysis/results/supplementary_tables/`: diagnostic tables S5-S7 and summary tables S10-S14, A7 and A8.
- README: corrected the statement on the intestine extract (it is a GitHub release asset; the Zenodo archive holds
  the source only) and updated the script table.

## v1.0.0 (2026-10-30)

- Initial release.
