"""G20-A4. Field-level model-based interval as an additional comparator for the immune proportion.

Quasi-binomial logistic model fitted to the k sampled fields (intercept only):
  p_hat = sum(y)/sum(n); phi = Pearson X^2/(k-1); Var(logit p_hat) = (1-k/M) * phi / (sum(n) p_hat (1-p_hat))
  95% interval: expit(logit p_hat +/- t_{k-1,0.975} * SE). If p_hat is 0 or 1, 0.5 is added to y and 1 to n.
It is compared on the same draws with cell-level Wilson, region-ratio t and Fieller sets
(implementations copied unchanged from the G18 package: region_methods.py, interval_methods.py).
Three annotated human sources; k = 3, 5, 10 where k < M; 500 draws per tissue and k. Seed 20261043.
"""
import sys, os
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'data'))
from interval_methods import wilson, region_interval
from region_methods import fieller_coefficients, polynomial_set
from common import load_all, wide

rng = np.random.default_rng(20261043)
DRAWS = 500


def quasibin(y, n, M):
    k = len(y); Y, N = y.sum(), n.sum()
    if Y == 0 or Y == N:
        Yc, Nc = Y + .5, N + 1
    else:
        Yc, Nc = Y, N
    p = Yc / Nc
    x2 = np.sum((y - n * p) ** 2 / (n * p * (1 - p)))
    phi = x2 / (k - 1)
    se = np.sqrt((1 - k / M) * phi / (Nc * p * (1 - p)))
    tq = stats.t.ppf(.975, k - 1)
    lg = np.log(p / (1 - p))
    lo, hi = 1 / (1 + np.exp(-(lg - tq * se))), 1 / (1 + np.exp(-(lg + tq * se)))
    return lo, hi, phi


rows = []
data = load_all()
for (src, tis), g in data.groupby(['source', 'tissue']):
    (labels, X, fields), = wide(g).values()
    lin = g.groupby('label').lineage.first()
    y_all = X[:, [i for i, l in enumerate(labels) if lin[l] == 'Immune']].sum(1)
    n_all = X.sum(1); M = len(n_all); truth = y_all.sum() / n_all.sum()
    arr = np.stack([y_all, n_all - y_all], 1)
    for k in [3, 5, 10]:
        if k >= M:
            continue
        for d in range(DRAWS):
            ids = rng.choice(M, k, replace=False)
            y, n = y_all[ids], n_all[ids]; phat = y.sum() / n.sum()
            wl, wh = wilson(float(y.sum()), float(n.sum()))
            rl, rh = region_interval(arr[ids], M)
            fs = polynomial_set(fieller_coefficients(y, n, M))
            ql, qh, phi = quasibin(y, n, M)
            for meth, segs in [('cell_Wilson', [(wl, wh)]), ('region_ratio_t', [(rl[0], rh[0])]), ('Fieller', fs), ('quasibinomial', [(ql, qh)])]:
                width = sum(h - l for l, h in segs)
                rows.append(dict(source=src, tissue=tis, k=k, draw=d, method=meth, width=width,
                                 full=width >= 1 - 1e-12, covered=any(l - 1e-12 <= truth <= h + 1e-12 for l, h in segs),
                                 phi=phi if meth == 'quasibinomial' else np.nan))
r = pd.DataFrame(rows)
r.to_csv('../results/A4_model_interval_draws.tsv.gz', sep='\t', index=False, compression='gzip')
t = r.groupby(['source', 'tissue', 'k', 'method']).agg(coverage=('covered', 'mean'), width=('width', 'mean'), full=('full', 'mean')).reset_index()
t.to_csv('../results/A4_model_interval_by_tissue.tsv', sep='\t', index=False)
s = t.groupby(['source', 'k', 'method']).agg(units=('tissue', 'nunique'), coverage=('coverage', 'median'), cov_min=('coverage', 'min'),
                                             cov_max=('coverage', 'max'), width=('width', 'median'), full=('full', 'median')).reset_index()
s.to_csv('../results/A4_model_interval_summary.tsv', sep='\t', index=False)
print(s.round(3).to_string())
print('median phi at k=3:', r[(r.method == 'quasibinomial') & (r.k == 3)].groupby('source').phi.median().round(1).to_dict())
