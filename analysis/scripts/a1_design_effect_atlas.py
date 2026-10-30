"""G20-A1. Design-effect atlas for cell-type proportions in three annotated human CosMx sources.

For each tissue and native cell category (and for five coarse lineages):
  linearised variance of the ratio estimator under SRS of k of M fields:
      V_F(k) = (1 - k/M) * S_e^2 / (k * nbar^2),  e_j = y_j - p * n_j
  variance under SRS of n cells without replacement (hypergeometric):
      V_C(n) = (1 - n/N) * p(1-p) * N/(N-1) / n
  intrinsic design effect (no finite-population terms):  D = S_e^2 / (nbar * p(1-p))
  implied intraclass correlation:                          rho = (D - 1)/(nbar - 1)
Validation: exact design effect for k = 3, i.e. the mean squared error of p_hat over all
3-field subsets (or 20,000 random subsets when M > 40) divided by the mean V_C(n_s) of
cell-count-matched random-cell samples.
"""
import itertools, sys
import numpy as np, pandas as pd
from scipy import stats
from common import load_all, wide

rng = np.random.default_rng(20261040)
K = 3
rows = []
data = load_all()
for (src, tis), g in data.groupby(['source', 'tissue']):
    group = g.group.iloc[0]
    lin = g.groupby('label').lineage.first()
    (labels, X, fields), = wide(g).values()
    # add lineage aggregates
    lineages = sorted(set(lin))
    L = np.stack([X[:, [labels.index(l) for l in labels if lin[l] == ln]].sum(1) for ln in lineages], 1)
    n = X.sum(1); M = len(n); N = n.sum(); nbar = n.mean()
    if M > 40:
        subs = np.array([rng.choice(M, K, replace=False) for _ in range(20000)])
    else:
        subs = np.array(list(itertools.combinations(range(M), K)))
    ns = n[subs].sum(1)
    for level, labs, Y in [('label', labels, X), ('lineage', lineages, L)]:
        for i, lab in enumerate(labs):
            y = Y[:, i]; p = y.sum() / N
            if p <= 0 or p >= 1:
                continue
            e = y - p * n; Se2 = e.var(ddof=1)
            D = Se2 / (nbar * p * (1 - p))
            vf = (1 - K / M) * Se2 / (K * nbar ** 2)
            vc = (1 - K * nbar / N) * p * (1 - p) * N / (N - 1) / (K * nbar)
            phat = Y[subs, i].sum(1) / ns
            mse_f = np.mean((phat - p) ** 2)
            vc_s = np.mean((1 - ns / N) * p * (1 - p) * N / (N - 1) / ns)
            rows.append(dict(source=src, tissue=tis, group=group, level=level, label=lab,
                             lineage=(lin[lab] if level == 'label' else lab), p=p, M=M, N=int(N), nbar=nbar,
                             cv_field_size=n.std(ddof=1) / nbar, D=D, rho=(D - 1) / (nbar - 1),
                             deff3_lin=vf / vc, deff3_exact=mse_f / vc_s, neff3=K * nbar / (mse_f / vc_s),
                             subsets=len(subs), exact_enumeration=M <= 40))
r = pd.DataFrame(rows)
r.to_csv('../results/A1_design_effect_atlas.tsv', sep='\t', index=False)

# --- summaries ---
lab = r[r.level == 'label']
ok = lab[lab.p >= 0.01]
summ = []
for src, g in ok.groupby('source'):
    summ.append(dict(source=src, tissues=g.tissue.nunique(), labels=g.label.nunique(), combos=len(g),
                     D_median=g.D.median(), D_q25=g.D.quantile(.25), D_q75=g.D.quantile(.75),
                     D_min=g.D.min(), D_max=g.D.max(), share_D_gt10=(g.D > 10).mean(), share_D_gt2=(g.D > 2).mean(),
                     rho_median=g.rho.median(), nbar_median=g.nbar.median(),
                     neff3_median=g.neff3.median(), cells3_median=(3 * g.nbar).median()))
summ = pd.DataFrame(summ)
summ.to_csv('../results/A1_summary_by_source.tsv', sep='\t', index=False)
print(summ.round(3).to_string())
allok = ok
print('ALL p>=1%: combos', len(allok), 'D median', round(allok.D.median(), 1), 'IQR', allok.D.quantile([.25, .75]).round(1).tolist(),
      'share>10', round((allok.D > 10).mean(), 3), 'share>2', round((allok.D > 2).mean(), 3), 'rho median', round(allok.rho.median(), 4))

# validation: linearised vs exact
v = lab.copy()
lr = np.log(v.deff3_lin / v.deff3_exact)
print('validation combos', len(v), 'median |log ratio|', round(np.median(np.abs(lr)), 4),
      'spearman', round(stats.spearmanr(v.deff3_lin, v.deff3_exact)[0], 4),
      'within 20%:', round((np.abs(lr) < np.log(1.2)).mean(), 3), 'p>=1% within 20%:', round((np.abs(lr[v.p >= .01]) < np.log(1.2)).mean(), 3))
v[['source', 'tissue', 'label', 'p', 'M', 'deff3_lin', 'deff3_exact', 'exact_enumeration']].to_csv('../results/A1_validation_lin_vs_exact.tsv', sep='\t', index=False)

# lineage summary
lin = r[r.level == 'lineage']
ls = lin.groupby(['source', 'label']).agg(tissues=('tissue', 'nunique'), p_median=('p', 'median'), D_median=('D', 'median'),
                                          D_min=('D', 'min'), D_max=('D', 'max'), rho_median=('rho', 'median')).reset_index()
ls.to_csv('../results/A1_lineage_summary.tsv', sep='\t', index=False)
print(ls.round(3).to_string())

# per-label medians (primary) and correlations
for src, g in lab.groupby('source'):
    rs = stats.spearmanr(g.D, g.p)
    print(src, 'spearman D~p', round(rs[0], 3), 'p=', '%.2g' % rs[1], 'n', len(g))

# determinants: OLS of log10 D on log10 p, log10 nbar, lineage and source (descriptive)
d = ok.copy()
d['lp'] = np.log10(d.p); d['ln'] = np.log10(d.nbar); d['lD'] = np.log10(d.D)
Xd = pd.get_dummies(d[['lp', 'ln', 'source', 'lineage']], columns=['source', 'lineage'], drop_first=True, dtype=float)
Xd.insert(0, 'const', 1.0)
beta, *_ = np.linalg.lstsq(Xd.to_numpy(), d.lD.to_numpy(), rcond=None)
pred = Xd.to_numpy() @ beta
r2 = 1 - np.sum((d.lD - pred) ** 2) / np.sum((d.lD - d.lD.mean()) ** 2)
coef = pd.Series(beta, index=Xd.columns)
print('OLS log10 D, R2 =', round(r2, 3)); print(coef.round(3).to_string())
coef.to_frame('beta').assign(R2=r2).to_csv('../results/A1_determinants_ols.tsv', sep='\t')
# variance decomposition: how much of log D varies between cell types vs between tissues (primary)
for src, g in ok.groupby('source'):
    if g.tissue.nunique() < 3: continue
    lD = np.log10(g.D)
    between_label = g.assign(lD=lD).groupby('label').lD.transform('mean')
    between_tis = g.assign(lD=lD).groupby('tissue').lD.transform('mean')
    print(src, 'share of var(log D) explained by cell type', round(np.var(between_label) / np.var(lD), 3),
          'by tissue', round(np.var(between_tis) / np.var(lD), 3))
