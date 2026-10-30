"""G20-A2. Held-out validation of fovplan.

For every tissue with at least 15 observed fields, and every native cell category or lineage with
reference proportion >= 1%:
  repeat R = 50 random splits: m = 5 pilot fields / remaining fields = validation pool (size P)
    for k in {3, 5, 10} with k <= P - 2:
      predicted 95% half-width from the pilot (FPC with M = P):
         naive  : cells treated as independent (design effect 1)
         fovplan: pilot design effect (point estimate)
         fovplan-80: 80% upper confidence bound of S_e^2
      realised: 200 random k-field samples from the pool; assurance = fraction with
                |p_hat - p_pool| <= predicted half-width (target 0.95)
Seed 20261041. Also tabulates fields required for a +/-25% relative precision target
using each tissue's complete field collection.
"""
import numpy as np, pandas as pd
from scipy import stats
from common import load_all, wide
import fovplan

rng = np.random.default_rng(20261041)
R, m, DRAWS = 50, 5, 200
data = load_all()
rows = []
for (src, tis), g in data.groupby(['source', 'tissue']):
    lin = g.groupby('label').lineage.first()
    (labels, X, fields), = wide(g).values()
    lineages = sorted(set(lin))
    L = np.stack([X[:, [labels.index(l) for l in labels if lin[l] == ln]].sum(1) for ln in lineages], 1)
    Y = np.concatenate([X, L], 1); names = labels + ['[lineage] ' + l for l in lineages]
    n = X.sum(1); M = len(n)
    if M < 15:
        continue
    pref = Y.sum(0) / n.sum()
    keep = np.nonzero(pref >= 0.01)[0]
    for rep in range(R):
        perm = rng.permutation(M); pil, pool = perm[:m], perm[m:]
        P = len(pool); npil = n[pil]; nbar_p = npil.mean()
        npool = n[pool]; ppool = Y[pool].sum(0) / npool.sum()
        for k in [3, 5, 10]:
            if k > P - 2:
                continue
            subs = np.array([rng.choice(P, k, replace=False) for _ in range(DRAWS)])
            ns = npool[subs].sum(1)
            PH = Y[pool][subs].sum(1) / ns[:, None]
            for i in keep:
                y = Y[pil, i]; p = y.sum() / npil.sum()
                if not 0 < p < 1:
                    rows.append(dict(source=src, tissue=tis, label=names[i], rep=rep, k=k, pilot_zero=True))
                    continue
                Se2 = (y - p * npil).var(ddof=1)
                h_f = float(fovplan.half_width(k, p, Se2, nbar_p, P))
                h_c = float(fovplan.half_width(k, p, fovplan.se2_upper(Se2, m, 0.8), nbar_p, P))
                h_n = 1.959964 * np.sqrt((1 - k / P) * p * (1 - p) / (k * nbar_p))
                phat = PH[:, i]
                err = np.abs(phat - ppool[i])
                rows.append(dict(source=src, tissue=tis, label=names[i], level='lineage' if names[i].startswith('[') else 'label',
                                 rep=rep, k=k, pilot_zero=False, p_pool=ppool[i],
                                 h_naive=h_n, h_fovplan=h_f, h_fovplan80=h_c, emp_rmse=np.sqrt(np.mean((phat - ppool[i]) ** 2)),
                                 assur_naive=np.mean(err <= h_n), assur_fovplan=np.mean(err <= h_f), assur_fovplan80=np.mean(err <= h_c)))
v = pd.DataFrame(rows)
v.to_csv('../results/A2_fovplan_validation_draws.tsv.gz', sep='\t', index=False, compression='gzip')
vv = v[~v.pilot_zero]
summ = vv.groupby(['source', 'k']).agg(tissues=('tissue', 'nunique'), combos=('label', 'size'),
                                         naive=('assur_naive', 'median'), fovplan=('assur_fovplan', 'median'),
                                         fovplan80=('assur_fovplan80', 'median'),
                                         naive_mean=('assur_naive', 'mean'), fovplan_mean=('assur_fovplan', 'mean'),
                                         fovplan80_mean=('assur_fovplan80', 'mean')).reset_index()
allk = vv.groupby('k').agg(combos=('label', 'size'), naive=('assur_naive', 'median'), fovplan=('assur_fovplan', 'median'),
                           fovplan80=('assur_fovplan80', 'median'), naive_mean=('assur_naive', 'mean'),
                           fovplan_mean=('assur_fovplan', 'mean'), fovplan80_mean=('assur_fovplan80', 'mean')).reset_index()
allk.insert(0, 'source', 'All')
summ = pd.concat([summ, allk])
summ.to_csv('../results/A2_fovplan_validation_summary.tsv', sep='\t', index=False)
print(summ.round(3).to_string())
print('pilot zero fraction', round(v.pilot_zero.mean(), 4), 'tissues', v.tissue.nunique(), sorted(v.tissue.unique()))
# ratio of predicted to empirical SE
vv = vv.assign(ratio_f=vv.h_fovplan / (1.96 * vv.emp_rmse), ratio_n=vv.h_naive / (1.96 * vv.emp_rmse))
print('median predicted/empirical half-width: fovplan', round(vv.ratio_f.median(), 3), 'naive', round(vv.ratio_n.median(), 3))

