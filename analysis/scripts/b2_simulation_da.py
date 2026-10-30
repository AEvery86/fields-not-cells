"""G21-B2. Differential-abundance simulation with known truth.

Each individual has a 3 x 3 mm section (36 fields of 500 um) generated as in B1 (sigma = 1,
ell = 200 um, focal base proportion 0.1). Individual-level biological variation: the focal
log-intensity is shifted by u_i ~ N(0, tau^2), tau = 0.3. Group 2 receives an additional shift
delta (0 = null, 0.4 = true effect, about a 45% relative increase). k fields are sampled per individual.
Configurations: n individuals per group in {3, 6, 12}; k in {3, 10, 36}; delta in {0, 0.4};
300 replicates each (seed 20261051). Six tests as in A3.
"""
import numpy as np, pandas as pd
from scipy import stats
from sim_core import section

rng = np.random.default_rng(20261051)
TESTS = ['cell_chisq', 'field_wilcoxon', 'field_quasibin', 'field_cluster', 'indiv_wilcoxon', 'indiv_quasibin']


def quasi2(y1, n1, y0, n0, df):
    Y1, N1, Y0, N0 = y1.sum(), n1.sum(), y0.sum(), n0.sum()
    if min(Y1, Y0) == 0 or Y1 == N1 or Y0 == N0:
        Y1, N1, Y0, N0 = Y1 + .5, N1 + 1, Y0 + .5, N0 + 1
    p1, p0 = Y1 / N1, Y0 / N0
    x2 = np.sum((y1 - n1 * p1) ** 2 / (n1 * p1 * (1 - p1))) + np.sum((y0 - n0 * p0) ** 2 / (n0 * p0 * (1 - p0)))
    se = np.sqrt(x2 / df * (1 / (N1 * p1 * (1 - p1)) + 1 / (N0 * p0 * (1 - p0))))
    return 2 * stats.t.sf(abs((np.log(p1 / (1 - p1)) - np.log(p0 / (1 - p0))) / se), df)


def cluster(yf, nf, cl, grp):
    out = {}
    for g in (1, 0):
        m = grp == g; Y, N = yf[m].sum(), nf[m].sum(); p = Y / N
        U = pd.Series(yf[m] - nf[m] * p).groupby(cl[m]).sum().to_numpy(); G = len(U)
        out[g] = (np.log(p / (1 - p)), (G / (G - 1)) * np.sum(U ** 2) / (N * p * (1 - p)) ** 2, G)
    t = (out[1][0] - out[0][0]) / np.sqrt(out[1][1] + out[0][1])
    return 2 * stats.t.sf(abs(t), out[1][2] + out[0][2] - 2)


def tests(yf, nf, cl, grp):
    a, b, c, d = yf[grp == 1].sum(), nf[grp == 1].sum(), yf[grp == 0].sum(), nf[grp == 0].sum()
    p_cell = stats.chi2_contingency(np.array([[a, b - a], [c, d - c]]), correction=False)[1]
    pf = yf / nf
    p_fw = stats.mannwhitneyu(pf[grp == 1], pf[grp == 0]).pvalue
    p_fq = quasi2(yf[grp == 1], nf[grp == 1], yf[grp == 0], nf[grp == 0], len(yf) - 2)
    p_fc = cluster(yf, nf, cl, grp)
    I = pd.DataFrame(dict(y=yf, n=nf, cl=cl, g=grp)).groupby('cl').agg(y=('y', 'sum'), n=('n', 'sum'), g=('g', 'first'))
    pi = (I.y / I.n).to_numpy(); gi = I.g.to_numpy()
    p_iw = stats.mannwhitneyu(pi[gi == 1], pi[gi == 0], method='exact').pvalue
    p_iq = quasi2(I.y.to_numpy()[gi == 1], I.n.to_numpy()[gi == 1], I.y.to_numpy()[gi == 0], I.n.to_numpy()[gi == 0], len(gi) - 2)
    return [p_cell, p_fw, p_fq, p_fc, p_iw, p_iq]


if __name__ == '__main__':
    rows = []
    for delta in [0.0, 0.4]:
        for nper in [3, 6, 12]:
            for k in [3, 10, 36]:
                for rep in range(300):
                    yf, nf, cl, grp = [], [], [], []
                    for g in (0, 1):
                        for i in range(nper):
                            u = rng.normal(0, 0.3) + (delta if g == 1 else 0)
                            y, n = section(1.0, 200, 0.1, nx=6, ny=6, rng=rng, shift=u)
                            ix = rng.choice(len(n), k, replace=False)
                            yf += list(y[ix]); nf += list(n[ix]); cl += [f'{g}_{i}'] * k; grp += [g] * k
                    ps = tests(np.array(yf), np.array(nf), np.array(cl), np.array(grp))
                    rows.append(dict(delta=delta, n_per_group=nper, k=k, rep=rep, **dict(zip(TESTS, ps))))
                print(delta, nper, k, flush=True)
    r = pd.DataFrame(rows)
    r.to_csv('../results/B2_sim_da_pvalues.tsv.gz', sep='\t', index=False, compression='gzip')
    s = r.groupby(['delta', 'n_per_group', 'k'])[TESTS].apply(lambda d: (d < 0.05).mean()).reset_index()
    s.to_csv('../results/B2_sim_da_rejection.tsv', sep='\t', index=False)
    print(s.round(3).to_string())
