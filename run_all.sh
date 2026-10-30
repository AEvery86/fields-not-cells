#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/analysis/scripts"
python a1_design_effect_atlas.py
python a2_fovplan_validation.py
python a2b_fields_required.py
python a3_false_positive_da.py
python a4_model_interval.py
python a6_independent_check.py
python b1_simulation.py
python b2_simulation_da.py
python c3_mouse_whole_section.py
python c0_hickey_assemble.py
python c4_hickey.py ../data/hickey_7sites.pkl ../results/C4_hickey
python c5_pilot_sensitivity.py
python a5_figures.py && python b3_figures.py && python c6_figure8.py && python p3_concept.py && python p3_supp_figures.py
