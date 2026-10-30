"""G20-A6. Independent recomputation of key G20 results with code that does not reuse common.py
loaders, fovplan or the analysis functions. Writes ../results/G20_independent_check.json."""
import json, itertools
import numpy as np, pandas as pd
from math import comb

D_ = '../data/'
out = {}
# ---------- 1. design effects from raw per-field files ----------
p = pd.read_csv(D_ + 'fov_celltype_counts.tsv', sep='\t')
g = pd.read_csv(D_ + 'GSE211785_fov_celltype_counts.tsv', sep='\t')
l = pd.read_csv(D_ + 'Danaher_fov_label_counts.tsv', sep='\t')
lw = l.pivot_table(index=['tissue', 'fovkey'], columns='label', values='count', aggfunc='sum', fill_value=0).reset_index()
frames = {'GSE325587': (p.drop(columns=['slide', 'fov']).rename(columns={'sample': 'unit'})),
          'GSE211785': (g.drop(columns=['slide', 'fov']).rename(columns={'sample': 'unit'})),
          'Lupus': (lw.drop(columns=['fovkey']).rename(columns={'tissue': 'unit'}))}
rec = []
for src, f in frames.items():
    for u, h in f.groupby('unit'):
        Y = h.drop(columns=['unit']).astype(float)
        n = Y.sum(axis=1)
        for c in Y.columns:
            pp = Y[c].sum() / n.sum()
            if 0 < pp < 1:
                resid = Y[c] - pp * n
                rec.append((src, u, c, pp, (resid ** 2).sum() / (len(n) - 1) / (n.mean() * pp * (1 - pp))))
ind = pd.DataFrame(rec, columns=['source', 'tissue', 'label', 'p', 'D_ind'])
a1 = pd.read_csv('../results/A1_design_effect_atlas.tsv', sep='\t')
a1 = a1[a1.level == 'label'].merge(ind, on=['source', 'tissue', 'label'])
out['design_effect'] = dict(combos=len(a1), max_rel_diff=float(np.max(np.abs(a1.D / a1.D_ind - 1))),
                            median_D_p1pct=float(ind[ind.p >= 0.01].D_ind.median()))

# ---------- 2. false-positive rates: cell-level z test and individual-level exact rank-sum ----------
def ranksum_exact_p(x, y):
    allv = np.concatenate([x, y]); r = pd.Series(allv).rank().to_numpy(); n1 = len(x); N = len(allv)
    W = r[:n1].sum()
    dist = np.array([r[list(c)].sum() for c in itertools.combinations(range(N), n1)])
    mu = dist.mean()
    return np.mean(np.abs(dist - mu) >= abs(W - mu) - 1e-9)

def ztest_p(a, n1, c, n0):
    from scipy.stats import norm
    p1, p0, pp = a / n1, c / n0, (a + c) / (n1 + n0)
    se = np.sqrt(pp * (1 - pp) * (1 / n1 + 1 / n0))
    return 2 * norm.sf(abs(p1 - p0) / se) if se > 0 else 1.0

fpr = {}
st = pd.read_csv(D_ + 'sample_mapping.tsv', sep='\t').set_index('sample').stage
dkd = p[p['sample'].map(st) != 'Control']
cats = dkd.drop(columns=['sample', 'slide', 'fov']).sum()
tot = cats.sum(); keep = list(cats[cats / p.drop(columns=['sample', 'slide', 'fov']).sum().sum() >= 0.005].index)
indiv = dkd.groupby('sample')[keep].sum(); ni = dkd.drop(columns=['slide', 'fov']).groupby('sample').sum().sum(axis=1)
units = list(indiv.index)
cell_rej = []; ind_rej = []
for s in itertools.combinations(units, 3):
    m = indiv.index.isin(s)
    for c in keep:
        cell_rej.append(ztest_p(indiv[c][m].sum(), ni[m].sum(), indiv[c][~m].sum(), ni[~m].sum()) < 0.05)
for c in keep[:3]:
    prop = (indiv[c] / ni).to_numpy()
    for s in itertools.combinations(range(len(units)), 3):
        m = np.isin(np.arange(len(units)), s)
        ind_rej.append(ranksum_exact_p(prop[m], prop[~m]) < 0.05)
a3 = pd.read_csv('../results/A3_null_pvalues.tsv.gz', sep='\t')
a3 = a3[(a3.source == 'GSE325587') & (a3.scenario == 'all_fields') & (~a3.category.str.startswith('[lineage]'))]
out['false_positive_DKD_labels'] = dict(cell_level_independent=float(np.mean(cell_rej)), cell_level_A3=float((a3.cell_chisq < 0.05).mean()),
                                        indiv_wilcoxon_independent_3cats=float(np.mean(ind_rej)),
                                        indiv_wilcoxon_A3=float((a3.indiv_wilcoxon < 0.05).mean()))

# ---------- 3. interval coverage, immune proportion, k=3, independent code and seed ----------
from scipy.stats import t as tdist, norm
rng = np.random.default_rng(99)
lin = pd.read_csv(D_ + 'harmonized_category_mapping.tsv', sep='\t')
imm_p = list(lin[lin.coarse == 'Immune'].original_label)
imm_l = sorted(l[l.is_immune.astype(str).str.upper() == 'TRUE'].label.unique())
sets = {'GSE325587': (frames['GSE325587'], imm_p), 'GSE211785': (frames['GSE211785'], ['Immune']), 'Lupus': (frames['Lupus'], imm_l)}
res = {}
for src, (f, im) in sets.items():
    cov = {'wilson': [], 'ratio_t': [], 'quasibin': []}
    for u, h in f.groupby('unit'):
        y = h[im].sum(axis=1).to_numpy(float); n = h.drop(columns=['unit']).sum(axis=1).to_numpy(float)
        M = len(n); P = y.sum() / n.sum(); k = 3
        hits = {kk: 0 for kk in cov}
        for _ in range(400):
            ix = rng.choice(M, k, replace=False); yy, nn = y[ix], n[ix]; ph = yy.sum() / nn.sum(); N = nn.sum()
            z = 1.959964; cen = (ph + z * z / (2 * N)) / (1 + z * z / N); hw = z * np.sqrt(ph * (1 - ph) / N + z * z / (4 * N * N)) / (1 + z * z / N)
            hits['wilson'] += cen - hw <= P <= cen + hw
            e = yy - ph * nn; v = (1 - k / M) * e.var(ddof=1) / (k * nn.mean() ** 2); q = tdist.ppf(.975, k - 1)
            hits['ratio_t'] += max(ph - q * np.sqrt(v), 0) <= P <= min(ph + q * np.sqrt(v), 1)
            pc = ph if 0 < ph < 1 else (yy.sum() + .5) / (N + 1)
            phi = np.sum((yy - nn * pc) ** 2 / (nn * pc * (1 - pc))) / (k - 1)
            se = np.sqrt((1 - k / M) * phi / ((N if 0 < ph < 1 else N + 1) * pc * (1 - pc)))
            lo, hi = 1 / (1 + np.exp(-(np.log(pc / (1 - pc)) - q * se))), 1 / (1 + np.exp(-(np.log(pc / (1 - pc)) + q * se)))
            hits['quasibin'] += lo <= P <= hi
        for kk in cov:
            cov[kk].append(hits[kk] / 400)
    res[src] = {kk: float(np.median(v)) for kk, v in cov.items()}
a4 = pd.read_csv('../results/A4_model_interval_summary.tsv', sep='\t'); a4 = a4[a4.k == 3]
cmp = {}
for src in res:
    s = a4[a4.source == src].set_index('method').coverage
    cmp[src] = dict(independent=res[src], A4=dict(wilson=float(s['cell_Wilson']), ratio_t=float(s['region_ratio_t']), quasibin=float(s['quasibinomial'])))
out['interval_coverage_k3'] = cmp
json.dump(out, open('../results/G20_independent_check.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
