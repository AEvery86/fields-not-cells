# fovplan tutorial

## Run the verified example

From the repository root:

```text
python examples/run_immune_example.py
```

This uses the bundled `src` directory, validates the counts and checks the computed plans against the expected
outputs supplied in `examples/`. To use the command-line tool:

```text
python -m pip install .
fovplan examples/immune_pilot_G24.tsv --target 0.25 --relative --assurance 0.9
fovplan examples/immune_pilot_G24.tsv --target 0.02 --assurance 0.9
python -m pytest tests
```

The five fields in `immune_pilot_G24.tsv` were sampled from deposited GSE325587 counts with seed 20261071 after
choosing the tissue by identifier; the tissue and field identifiers are in `immune_pilot_G24_metadata.tsv`. The
count table contains immune cells and all other cells. The pilot proportion is 0.0306 and D is 1.89. A relative
target of ±25% (±0.0077 in proportion units) gives a point plan of 7 fields and a 90%-assurance plan of 17; an
absolute target of ±0.02 gives 3 and 5 fields. This retrospective demonstration is not prospective validation.

## Inputs

- Analyze one tissue sampling frame at a time. Supply counts of mutually exclusive, exhaustive categories, or a focal
  category and its complement; do not supply proportions.
- Keep patient, tissue and field identifiers and field areas in a separate metadata file. The command-line tool
  treats every numeric column as a count, so numeric identifiers, areas, totals or overlapping lineages would corrupt
  the denominator.
- Pilot and planned fields should share geometry, eligibility rules and the selection design. Field area is not a
  model input, and the calculation does not correct preferential selection or unequal inclusion probabilities.

## Interpreting plans

- A relative target is a fraction of the pilot proportion; an absolute target is a difference in proportion units.
  An absolute ±0.02 is not a relative ±2%.
- `--fields-available M` describes a real finite sampling frame, not an imaging budget. If a plan exceeds the budget,
  report it as infeasible rather than capping it.
- Two pilot fields are the computational minimum, not an adequate pilot. Boundary proportions (0 or 1) and zero
  residual variance give no identifiable plan; expand the pilot. Rare categories can give unstable plans.
- The assurance plan uses an approximate upper bound of the between-field variance; it is not a guarantee. Taking the
  maximum over cell types does not give simultaneous coverage.
- fovplan plans fields within a tissue. It does not determine the number of individuals or the power of group
  comparisons: fields within an individual are not biological replicates.
