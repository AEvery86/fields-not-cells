"""G20-A2b. Fields required per cell type and lineage (complete field collection used as pilot).
Targets: +/-25% relative half-width and +/-0.02 absolute half-width (95%, no finite-population correction)."""
import numpy as np, pandas as pd
from common import load_all, wide
import fovplan
data = load_all()
# fields required for +/-25% relative precision and +/-0.02 absolute, full field collection as pilot
req = []
for (src, tis), g in data.groupby(['source', 'tissue']):
    lin = g.groupby('label').lineage.first()
    (labels, X, fields), = wide(g).values()
    lineages = sorted(set(lin))
    L = np.stack([X[:, [labels.index(l) for l in labels if lin[l] == ln]].sum(1) for ln in lineages], 1)
    for level, labs, Y in [('label', labels, X), ('lineage', lineages, L)]:
        est = fovplan.plan(pd.DataFrame(Y, columns=labs), 0.25, relative=True)
        est2 = fovplan.plan(pd.DataFrame(Y, columns=labs), 0.02, relative=False)
        est['fields_abs002'] = est2.fields_required; est['naive_abs002'] = est2.fields_if_cells_independent
        est['source'] = src; est['tissue'] = tis; est['level'] = level; est['M'] = len(Y)
        req.append(est)
req = pd.concat(req, ignore_index=True)
req['le_M'] = req.fields_required <= req.M
req.to_csv('../results/A2_fields_required.tsv', sep='\t', index=False)
q = req[(req.p >= 0.01)]
tab = q.groupby(['source', 'level']).agg(combos=('p', 'size'), rel25_median=('fields_required', 'median'),
                                         rel25_naive_median=('fields_if_cells_independent', 'median'),
                                         abs002_median=('fields_abs002', 'median'), abs002_naive_median=('naive_abs002', 'median'),
                                         share_rel25_le_M=('le_M', 'mean')).reset_index()
tab.to_csv('../results/A2_fields_required_summary.tsv', sep='\t', index=False)
print(tab.round(2).to_string())
lin = req[(req.level == 'lineage')].groupby(['source', 'cell_type']).agg(p=('p', 'median'), rel25=('fields_required', 'median'),
                                                                         rel25_naive=('fields_if_cells_independent', 'median'),
                                                                         abs002=('fields_abs002', 'median')).reset_index()
lin.to_csv('../results/A2_fields_required_lineage.tsv', sep='\t', index=False)
print(lin.round(3).to_string())
