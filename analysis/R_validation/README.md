# R reference-package re-analysis of the A7 comparison

The mixed-model and differential-abundance comparison (A7) was implemented in Python (`analysis/scripts/da_models.py`).
These scripts re-fit the same eight procedures on identical inputs with the reference R implementations:
`lme4::glmer` (individual random intercept, nAGQ = 7; individual and field random intercepts, Laplace),
`glmmTMB` beta-binomial regression, `MASS::glm.nb` and `speckle::propeller.ttest` with `limma`.

1. Export the inputs: `cd analysis/scripts && python a7_export_for_R.py ../R_validation/data`
   (or unzip `data/` from the v1.0.1 release asset), and place the Python A7 results in `python_results/`.
2. Install R (4.6.2 was used) and run `Rscript 00_install_packages.R`.
3. Run `Rscript 01_validate_A7.R quick` and then `Rscript 01_validate_A7.R full [cores]`
   (Windows: `run_quick.bat`, `run_full.bat`). Results are cached in `output/cache`, so runs can be resumed.

The report `output/R_vs_Python_report_full.txt` gives, for each procedure, Python and R rejection rates, decision
agreement, Spearman correlation and differences of log10 P values. In our run (R 4.6.2; lme4 2.1-0, glmmTMB 1.1.15.2,
MASS 7.3-66, speckle 1.12.0, limma 3.68.5), decisions agreed in 99.5-100% of 30,578 comparisons; all
negative-binomial discrepancies arose from non-converged `glm.nb` fits (Supplementary Table S14).
