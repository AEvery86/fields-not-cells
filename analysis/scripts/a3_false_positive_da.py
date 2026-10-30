"""G20-A3. False-positive rates of differential-abundance tests that treat cells or fields as replicates.

Data: GSE325587 (3 control vs 14 DKD tissues; one tissue per individual) and the pediatric lupus
source (4 control vs 8 SLE individuals; the three SLE8 biopsies are pooled into one individual).
Cell categories: native labels with overall proportion >= 0.5% plus five coarse lineages.

Tests (two-sided; a = 0.05):
  cell_chisq     cells as replicates: Pearson chi-square (2x2) on pooled counts
  field_wilcoxon fields as replicates: Wilcoxon rank-sum on per-field proportions
  field_quasibin fields as replicates: quasi-binomial logistic model (Pearson dispersion, t on F-2 df)
  field_cluster  field-level logistic model with individual-clustered sandwich variance (t on G-2 df)
  indiv_wilcoxon individuals as replicates: exact Wilcoxon rank-sum on per-individual proportions
  indiv_quasibin individuals as replicates: quasi-binomial logistic model on individual totals (t on G-2 df)
Null: all random splits of the disease individuals into two pseudo-groups of the real sizes
(GSE325587: 3 vs 11 of 14 DKD; lupus: 4 vs 4 of 8 SLE individuals), (a) using all fields and
(b) using 3 random fields per tissue. Real contrast: disease vs control, BH-adjusted.
Seed 20261042.
"""
import itertools
import numpy as np, pandas as pd
from scipy import stats
from common import load_primary, load_lupus

rng = np.random.default_rng(20261042)
MINP = 0.005


def prepare(long):
    """Return fields table: patient, tissue, field, counts per category (labels + lineages)."""
    lin = long.groupby('label').lineage.first()
    w = long.pivot_table(index=['patient', 'tissue', 'field'], columns='label', values='count', aggfunc='sum', fill_value=0)
    tot = w.sum(1)
    overall = w.sum(0) / tot.sum()
    cats = list(overall[overall >= MINP].index)
    W = w[cats].copy()
    for ln in sorted(set(lin)):
        W['[lineage] ' + ln] = w[[c for c in w.columns if lin[c] == ln]].sum(1)
    W['__n'] = tot
    return W.reset_index()


def logit_diff_quasi(y1, n1, y0, n0, df_resid):
    """Quasi-binomial two-group logistic model on observations (vectors). Returns t, p."""
    Y1, N1, Y0, N0 = y1.sum(), n1.sum(), y0.sum(), n0.sum()
    if Y1 == 0 or Y0 == 0 or Y1 == N1 or Y0 == N0:
        Y1, N1, Y0, N0 = Y1 + .5, N1 + 1, Y0 + .5, N0 + 1
    p1, p0 = Y1 / N1, Y0 / N0
    x2 = np.sum((y1 - n1 * p1) ** 2 / (n1 * p1 * (1 - p1))) + np.sum((y0 - n0 * p0) ** 2 / (n0 * p0 * (1 - p0)))
    phi = x2 / df_resid
    se = np.sqrt(phi * (1 / (N1 * p1 * (1 - p1)) + 1 / (N0 * p0 * (1 - p0))))
    t = (np.log(p1 / (1 - p1)) - np.log(p0 / (1 - p0))) / se
    return t, 2 * stats.t.sf(abs(t), df_resid)


def cluster_logit(yf, nf, cl, grp):
    """Field-level logistic (independence working model) with cluster-robust variance.
    yf, nf: per-field counts; cl: cluster id per field; grp: 1/0 per field."""
    out = {}
    for g in (1, 0):
        m = grp == g
        Y, N = yf[m].sum(), nf[m].sum()
        if Y == 0 or Y == N:
            Y, N = Y + .5, N + 1
        p = Y / N
        U = pd.Series(yf[m] - nf[m] * p).groupby(cl[m]).sum().to_numpy()
        G = len(U)
        var = (G / (G - 1)) * np.sum(U ** 2) / (N * p * (1 - p)) ** 2
        out[g] = (np.log(p / (1 - p)), var, G)
    G = out[1][2] + out[0][2]
    t = (out[1][0] - out[0][0]) / np.sqrt(out[1][1] + out[0][1])
    return t, 2 * stats.t.sf(abs(t), G - 2)


def run_tests(F, cats, g1_patients):
    """F: fields table (subset); g1_patients: set of patients in group 1. Returns DataFrame of p-values."""
    grp = F.patient.isin(g1_patients).to_numpy().astype(int)
    nf = F['__n'].to_numpy(float); cl = F.patient.to_numpy()
    P = F.groupby('patient')[cats + ['__n']].sum()
    pg = P.index.isin(g1_patients).astype(int)
    res = []
    for c in cats:
        yf = F[c].to_numpy(float)
        # cells
        a, b = yf[grp == 1].sum(), nf[grp == 1].sum(); cc, d = yf[grp == 0].sum(), nf[grp == 0].sum()
        tab = np.array([[a, b - a], [cc, d - cc]])
        p_cell = stats.chi2_contingency(tab, correction=False)[1] if tab.min() >= 0 and tab.sum(0).min() > 0 else 1.0
        # fields
        pf = yf / nf
        p_fw = stats.mannwhitneyu(pf[grp == 1], pf[grp == 0], alternative='two-sided').pvalue
        p_fq = logit_diff_quasi(yf[grp == 1], nf[grp == 1], yf[grp == 0], nf[grp == 0], len(yf) - 2)[1]
        p_fc = cluster_logit(yf, nf, cl, grp)[1]
        # individuals
        yi, ni = P[c].to_numpy(float), P['__n'].to_numpy(float)
        pi = yi / ni
        p_iw = stats.mannwhitneyu(pi[pg == 1], pi[pg == 0], alternative='two-sided', method='exact').pvalue
        p_iq = logit_diff_quasi(yi[pg == 1], ni[pg == 1], yi[pg == 0], ni[pg == 0], len(yi) - 2)[1]
        res.append(dict(category=c, cell_chisq=p_cell, field_wilcoxon=p_fw, field_quasibin=p_fq,
                        field_cluster=p_fc, indiv_wilcoxon=p_iw, indiv_quasibin=p_iq))
    return pd.DataFrame(res)


TESTS = ['cell_chisq', 'field_wilcoxon', 'field_quasibin', 'field_cluster', 'indiv_wilcoxon', 'indiv_quasibin']


def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def subsample(F, k):
    parts = []
    for t, g in F.groupby('tissue'):
        parts.append(g.iloc[rng.choice(len(g), min(k, len(g)), replace=False)])
    return pd.concat(parts)


null_rows, real_rows = [], []
for name, loader, disease, n1 in [('GSE325587', load_primary, 'DKD', 3), ('Lupus', load_lupus, 'SLE', 4)]:
    long = loader()
    F = prepare(long)
    cats = [c for c in F.columns if c not in ('patient', 'tissue', 'field', '__n')]
    grp_of = long.groupby('patient').group.first()
    dis = sorted(grp_of[grp_of == disease].index); ctl = sorted(grp_of[grp_of == 'Control'].index)
    # real contrast (all fields)
    r = run_tests(F, cats, set(dis)); r['source'] = name; r['contrast'] = f'{disease} vs control'
    for t in TESTS:
        r[t + '_q'] = bh(r[t])
    real_rows.append(r)
    # null splits within disease individuals
    Fd = F[F.patient.isin(dis)]
    splits = list(itertools.combinations(dis, n1))
    for s_i, s in enumerate(splits):
        for scen, FF in [('all_fields', Fd), ('3_fields', subsample(Fd, 3))]:
            r = run_tests(FF, cats, set(s)); r['source'] = name; r['split'] = s_i; r['scenario'] = scen
            null_rows.append(r)
    print(name, 'categories', len(cats), 'splits', len(splits), 'disease', len(dis), 'control', len(ctl))

null = pd.concat(null_rows, ignore_index=True)
null.to_csv('../results/A3_null_pvalues.tsv.gz', sep='\t', index=False, compression='gzip')
real = pd.concat(real_rows, ignore_index=True)
real.to_csv('../results/A3_real_contrast.tsv', sep='\t', index=False)

null['level'] = np.where(null.category.str.startswith('[lineage]'), 'lineage', 'label')
fpr = null.groupby(['source', 'scenario'])[TESTS].apply(lambda d: (d < 0.05).mean()).reset_index()
fpr.to_csv('../results/A3_false_positive_rates.tsv', sep='\t', index=False)
print('Per-test false-positive rate (nominal 0.05):'); print(fpr.round(3).to_string())
# family-level: share of splits with >=1 BH-significant category
fam = []
for (src, scen, sp), d in null.groupby(['source', 'scenario', 'split']):
    row = dict(source=src, scenario=scen, split=sp)
    for t in TESTS:
        row[t] = (bh(d[t]) < 0.05).any()
        row[t + '_n'] = (bh(d[t]) < 0.05).sum()
    fam.append(row)
fam = pd.DataFrame(fam)
famsum = fam.groupby(['source', 'scenario'])[TESTS].mean().reset_index()
famsum.to_csv('../results/A3_familywise_any_discovery.tsv', sep='\t', index=False)
print('Share of null splits with >=1 BH discovery:'); print(famsum.round(3).to_string())
nsum = fam.groupby(['source', 'scenario'])[[t + '_n' for t in TESTS]].mean().reset_index()
print('Mean BH discoveries per null split:'); print(nsum.round(2).to_string())
nsum.to_csv('../results/A3_mean_null_discoveries.tsv', sep='\t', index=False)
# MC uncertainty is not binomial (splits share individuals); report range across categories
bycat = null[null.scenario == 'all_fields'].groupby(['source', 'category'])[TESTS].apply(lambda d: (d < 0.05).mean()).reset_index()
bycat.to_csv('../results/A3_fpr_by_category.tsv', sep='\t', index=False)
print(bycat.groupby('source')[TESTS].agg(['min', 'median', 'max']).round(3).T.to_string())
# real contrast discoveries
disc = real.groupby('source')[[t + '_q' for t in TESTS]].apply(lambda d: (d < 0.05).sum()).reset_index()
disc['categories'] = real.groupby('source').size().values
disc.to_csv('../results/A3_real_discoveries.tsv', sep='\t', index=False)
print('Real contrast BH discoveries:'); print(disc.to_string())
