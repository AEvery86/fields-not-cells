"""Summaries for A8 planner benchmark and fovplan runtime (fovplan 0.2.0)."""
import glob, time, numpy as np, pandas as pd
import fovplan
R = '../results/a8/'
d = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(R + '*.tsv.gz'))])
order = ['cellcount', 'bootstrap', 'fovplan', 'fovplan80', 'fovplan90']
src = ['simulated (sigma = 0.5)', 'simulated (sigma = 1.0)', 'simulated (sigma = 1.5)', 'mouse kidney', 'human small intestine',
       'human colon', 'human kidney (GSE325587)', 'human kidney (GSE211785)', 'human kidney (Lupus)']
fin = lambda s: np.nanmedian(np.where(np.isfinite(s), s, np.nan))
g = d.groupby(['source', 'method']).agg(share=('share', 'mean'), k=('k', fin), infeasible=('infeasible', 'mean')).reset_index()
g.to_csv(R + 'A8_planner_long.tsv', sep='\t', index=False)
tab = pd.DataFrame(index=src)
for m in order:
    sub = g[g.method == m].set_index('source').reindex(src)
    tab[m] = [f'{100 * a:.1f} ({b:.0f})' for a, b in zip(sub.share, sub.k)]
u = d[d.method == 'fovplan'].groupby('source').agg(units=('unit', 'nunique'), plans=('rep', 'size')).reindex(src)
tab = tab.join(u); tab.to_csv(R + 'A8_planner_table.tsv', sep='\t')
print(tab.to_string())
print(d.groupby('method').share.mean().reindex(order).round(3))
print('min source share', g.groupby('method').share.min().reindex(order).round(3))
rt = (d.groupby('method').seconds.median().reindex(order) * 1e3); print('median ms per plan (A8 run, fovplan 0.1.0 grid search)\n', rt.round(3))
rng = np.random.default_rng(1); res = []
for K, m in [(30, 5), (100, 10), (1000, 50)]:
    X = rng.poisson(rng.gamma(1, 20, size=K), size=(m, K))
    reps = 20 if K < 1000 else 3
    t = time.perf_counter()
    for _ in range(reps):
        fovplan.plan(pd.DataFrame(X), 0.25, relative=True, assurance=0.9)
    res.append((K, m, (time.perf_counter() - t) / reps))
rs = pd.DataFrame(res, columns=['cell_types', 'pilot_fields', 'seconds']); rs.to_csv(R + 'A8_runtime_scale.tsv', sep='\t', index=False)
print(rs)
yp = rng.poisson(30, 5).astype(float); np_ = rng.poisson(1000, 5).astype(float)
import a8_planner_benchmark as A
t = time.perf_counter()
for _ in range(20):
    A.bootstrap_k(yp, np_, 0.25 * yp.sum() / np_.sum(), 200, rng)
print('bootstrap planner ms per plan', 1e3 * (time.perf_counter() - t) / 20)
t = time.perf_counter()
for _ in range(200):
    fovplan.required_fields(yp.sum() / np_.sum(), (yp - yp.sum() / np_.sum() * np_).var(ddof=1), np_.mean(), 0.25, True, M=200)
print('fovplan 0.2.0 ms per plan', 1e3 * (time.perf_counter() - t) / 200)
