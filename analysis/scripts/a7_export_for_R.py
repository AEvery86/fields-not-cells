"""Export the exact inputs of the A7 comparison (kidney, intestine, simulation) for re-analysis in R
with the reference packages (lme4, glmmTMB, MASS, speckle/limma).

Reproduces the field tables, field subsets, categories and pseudo-group splits of a7_glmm_da.py
(same seeds and code paths) and writes them as plain CSV files.
usage: a7_export_for_R.py OUTDIR
"""
import sys, os, itertools
import numpy as np, pandas as pd
sys.path.insert(0, '.')
import a7_glmm_da as A7

out = sys.argv[1]; os.makedirs(out, exist_ok=True)
w = lambda df, name: df.to_csv(os.path.join(out, name), index=False)

# ------------------------------------------------------------------ kidney (as run_kidney)
from common import load_primary, load_lupus
rng = np.random.default_rng(20261080)
for name, loader, disease, n1 in [('GSE325587', load_primary, 'DKD', 3), ('Lupus', load_lupus, 'SLE', 4)]:
    long = loader()
    lin = long.groupby('label').lineage.first()
    W = long.pivot_table(index=['patient', 'tissue', 'field'], columns='label', values='count', aggfunc='sum', fill_value=0)
    labcols = list(W.columns)
    overall = W.sum(0) / W.to_numpy().sum()
    cats = list(overall[overall >= 0.005].index) + ['[lineage] ' + l for l in sorted(set(lin))]
    F = W.reset_index()
    grp_of = long.groupby('patient').group.first()
    dis = sorted(grp_of[grp_of == disease].index)
    Fd = F[F.patient.isin(dis)]
    sub3 = pd.concat([g.iloc[rng.choice(len(g), min(3, len(g)), replace=False)] for _, g in Fd.groupby('tissue')])
    F = F.copy(); F.insert(0, 'row', np.arange(len(F)))
    F.insert(1, 'in_3_fields', F.index.isin(sub3.index).astype(int))
    F.insert(2, 'disease', F.patient.isin(dis).astype(int))
    F.insert(3, 'group', F.patient.map(grp_of))
    lab = pd.DataFrame(dict(column=[f'L{j + 1}' for j in range(len(labcols))], label=labcols, lineage=lin.reindex(labcols).to_numpy()))
    F = F.rename(columns=dict(zip(labcols, lab.column)))
    w(F, f'kidney_{name}_fields.csv'); w(lab, f'kidney_{name}_labels.csv')
    w(pd.DataFrame(dict(category=cats)), f'kidney_{name}_categories.csv')
    sp = [dict(split=s_i, member=m) for s_i, s in enumerate(itertools.combinations(dis, n1)) for m in s]
    w(pd.DataFrame(sp), f'kidney_{name}_splits.csv')
    w(pd.DataFrame(dict(individual=dis)), f'kidney_{name}_disease_individuals.csv')
    print(name, F.shape, len(cats), 'splits', sp[-1]['split'] + 1, flush=True)

# ------------------------------------------------------------------ intestine (as run_hickey)
import c2_whole_section as c2
PX = 0.37744
d = pd.read_pickle(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'hickey_7sites.pkl')).dropna(subset=['cell_type', 'unique_region', 'x', 'y'])
sites = sorted(d.Tissue_location.unique())
for si, site in enumerate(sites):
    dd = d[d.Tissue_location == site]
    cells = pd.DataFrame(dict(section=dd.unique_region, donor=dd.donor, x_um=dd.x * PX, y_um=dd.y * PX, label=dd.cell_type))
    cells['field'] = pd.concat([c2.tile(g, 500) for _, g in cells.groupby('section')]).reindex(cells.index)
    W = cells.pivot_table(index=['donor', 'section', 'field'], columns='label', values='x_um', aggfunc='size', fill_value=0)
    W = W[W.sum(1) >= 100]
    labcols = list(W.columns)
    prop = W.sum(0) / W.to_numpy().sum()
    cats = list(prop[prop >= 0.01].index)
    donors = sorted(W.index.get_level_values('donor').unique())
    W = W.reset_index()
    r = np.random.default_rng(20261081 + si)
    sub10 = W.iloc[r.permutation(len(W))].groupby('section').head(10)
    F = W.copy(); F.insert(0, 'row', np.arange(len(F))); F.insert(1, 'in_10_fields', F.index.isin(sub10.index).astype(int))
    lab = pd.DataFrame(dict(column=[f'L{j + 1}' for j in range(len(labcols))], label=labcols))
    F = F.rename(columns=dict(zip(labcols, lab.column)))
    w(F, f'intestine_site{si}_fields.csv'); w(lab, f'intestine_site{si}_labels.csv')
    w(pd.DataFrame(dict(category=cats)), f'intestine_site{si}_categories.csv')
    sp = [dict(split=s_i, member=m) for s_i, s in enumerate(itertools.combinations(donors, len(donors) // 2)) if s[0] == donors[0] for m in s]
    w(pd.DataFrame(sp), f'intestine_site{si}_splits.csv')
    print(site, F.shape, len(cats), flush=True)
w(pd.DataFrame(dict(site_index=range(len(sites)), site=sites)), 'intestine_sites.csv')

# ------------------------------------------------------------------ simulation (as sim_job)
rows = []
for rep in range(200):
    rg = np.random.default_rng(20261082 + rep)
    for delta in (0.0, 0.4):
        for g in (0, 1):
            for i in range(12):
                u = rg.normal(0, 0.3) + (delta if g == 1 else 0)
                C = A7.section_all(1.0, 200, 0.1, 6, 6, rg, u)
                C = C[rg.permutation(len(C))]
                for f in range(len(C)):
                    rows.append((rep, delta, g, i, f + 1, *C[f].astype(int)))
    if rep % 50 == 49:
        print('sim rep', rep + 1, flush=True)
S = pd.DataFrame(rows, columns=['rep', 'delta', 'group', 'individual', 'field_order'] + [f'T{t + 1}' for t in range(6)])
S.to_csv(os.path.join(out, 'simulation_fields.csv.gz'), index=False, compression='gzip')
print('simulation', S.shape)
