"""G22-C5. Sensitivity of fovplan to pilot size (5, 8, 10 fields) and assurance level (0.8, 0.9) against
whole-section truth in the Hickey intestine regions and the two mouse kidney sections (seed 20261062)."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '.')
import c2_whole_section as c2, fovplan
rng = np.random.default_rng(20261062)
d = pd.read_pickle('../data/hickey_7sites.pkl')
cells = [pd.DataFrame(dict(section=d.unique_region, organ=d.tissue.map({'SB': 'human small intestine', 'CL': 'human colon'}),
                           x_um=d.x * 0.37744, y_um=d.y * 0.37744, label=d.cell_type))]
for s in ['ShamR', 'Day14R']:
    m = pd.read_csv(f'../data/mouse_example_{s}.tsv.gz', sep='\t')
    cells.append(pd.DataFrame(dict(section=s, organ='mouse kidney', x_um=m.x, y_um=m.y, label=m.celltype)))
cells = pd.concat(cells, ignore_index=True)
cells['field'] = pd.concat([c2.tile(g, 500) for _, g in cells.groupby('section')]).reindex(cells.index)
rows = []
for sec, g in cells.groupby('section'):
    truth = g.label.value_counts(normalize=True)
    W = g.pivot_table(index='field', columns='label', values='x_um', aggfunc='size', fill_value=0)
    W = W[W.sum(axis=1) >= 100]; X = W.to_numpy(float); n = X.sum(1); M = len(n)
    if M < 20: continue
    for j, lab in enumerate(W.columns):
        P = truth[lab]
        if P < 0.01: continue
        y = X[:, j]
        for m in [5, 8, 10]:
            for rep in range(10):
                perm = rng.permutation(M); pil, pool = perm[:m], perm[m:]
                yp, npl = y[pil], n[pil]; pp = yp.sum() / npl.sum()
                if not 0 < pp < 1: continue
                Sp = (yp - pp * npl).var(ddof=1)
                for q in [None, 0.8, 0.9]:
                    S = Sp if q is None else fovplan.se2_upper(Sp, m, q)
                    k = fovplan.required_fields(pp, S, npl.mean(), .25, True, M=len(pool))
                    if not np.isfinite(k) or k > len(pool):
                        rows.append(dict(organ=g.organ.iloc[0], section=sec, label=lab, pilot=m, q=q or 0, k=np.inf, assur=np.nan)); continue
                    ss = np.array([rng.choice(len(pool), int(k), replace=False) for _ in range(150)])
                    est = y[pool][ss].sum(1) / n[pool][ss].sum(1)
                    rows.append(dict(organ=g.organ.iloc[0], section=sec, label=lab, pilot=m, q=q or 0, k=k, assur=np.mean(np.abs(est - P) <= .25 * P)))
r = pd.DataFrame(rows); r.to_csv('../results/C5_pilot_sensitivity.tsv', sep='\t', index=False)
print(r.groupby(['organ', 'pilot', 'q']).assur.mean().unstack().round(3))
print(r.groupby(['pilot', 'q']).assur.mean().unstack().round(3))
print(r.groupby(['pilot', 'q']).k.apply(lambda s: np.median(s[np.isfinite(s)])).unstack())
print('infeasible share', r.groupby(['pilot', 'q']).k.apply(lambda s: np.mean(~np.isfinite(s))).unstack().round(3))
