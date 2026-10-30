"""G21-C2. Generic whole-section analysis (applies to any cell table with coordinates).

Input: TSV/CSV/parquet with columns  section, donor, organ, x_um, y_um, label  (one row per cell).
For each section, the tissue is tiled into square virtual fields of width F um (default 500);
fields with fewer than MIN cells are not sampled. The TRUE proportion is computed from all cells
in the section (whole-section truth), so error caused by restricting sampling to well-populated
fields is included.

Outputs (prefix given by --out):
  _atlas.tsv      section x label: p_true, D (all retained fields), rho, effective cells of 3 fields
  _sampling.tsv   section x k: median TV of k random fields and of matched random cells vs truth
  _intervals.tsv  section x label(lineage or focal) x k x method: coverage of the TRUE proportion
  _planning.tsv   section x label: fovplan (5 pilot fields) assurance vs truth, cell-count plan for comparison
  _da.tsv         organ: false-positive rates under random donor splits (when >= 6 donors)
Seed 20261060.
"""
import argparse, itertools, sys, os
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
import fovplan

rng = np.random.default_rng(20261060)


def tile(df, F):
    fx = np.floor((df.x_um - df.x_um.min()) / F).astype(int)
    fy = np.floor((df.y_um - df.y_um.min()) / F).astype(int)
    return fx.astype(str) + '_' + fy.astype(str)


def wilson(y, n, z=1.959964):
    p = y / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def ratio_t(y, n, M):
    k = len(y); p = y.sum() / n.sum(); e = y - p * n
    v = (1 - k / M) * e.var(ddof=1) / (k * n.mean() ** 2); h = stats.t.ppf(.975, k - 1) * np.sqrt(v)
    return max(p - h, 0), min(p + h, 1)


def quasibin(y, n, M):
    k = len(y); Y, N = y.sum(), n.sum()
    Yc, Nc = (Y + .5, N + 1) if (Y == 0 or Y == N) else (Y, N)
    p = Yc / Nc; phi = np.sum((y - n * p) ** 2 / (n * p * (1 - p))) / (k - 1)
    se = np.sqrt((1 - k / M) * phi / (Nc * p * (1 - p))); q = stats.t.ppf(.975, k - 1); lg = np.log(p / (1 - p))
    return 1 / (1 + np.exp(-(lg - q * se))), 1 / (1 + np.exp(-(lg + q * se)))


def analyse(cells, F=500, MIN=100, out='whole', focal=None, draws=500):
    cells = cells.reset_index(drop=True); cells['field'] = pd.concat([tile(d, F) for _, d in cells.groupby('section')]).reindex(cells.index)
    atlas, samp, ints, plan = [], [], [], []
    for sec, d in cells.groupby('section'):
        organ, donor = d.organ.iloc[0], d.donor.iloc[0]
        truth = d.label.value_counts(normalize=True)
        W = d.pivot_table(index='field', columns='label', values='x_um', aggfunc='size', fill_value=0)
        W = W[W.sum(1) >= MIN]
        labels = list(truth.index); W = W.reindex(columns=labels, fill_value=0)
        X = W.to_numpy(float); n = X.sum(1); M = len(n)
        if M < 15:
            continue
        nbar = n.mean(); cover = n.sum() / len(d)
        for i, lab in enumerate(labels):
            P = truth[lab]
            if not 0 < P < 1: continue
            y = X[:, i]; pf = y.sum() / n.sum()
            Se2 = (y - pf * n).var(ddof=1); D = Se2 / (nbar * pf * (1 - pf)) if 0 < pf < 1 else np.nan
            atlas.append(dict(section=sec, organ=organ, donor=donor, label=lab, p_true=P, p_fields=pf, M=M, nbar=nbar,
                              coverage_of_cells=cover, D=D, rho=(D - 1) / (nbar - 1), neff3=3 * nbar / D if D else np.nan))
        # composition error vs whole-section truth
        q = truth.reindex(labels).to_numpy()
        for k in [1, 3, 5, 10]:
            if k >= M: continue
            tvf, tvc = [], []
            for _ in range(draws):
                ix = rng.choice(M, k, replace=False); s = X[ix].sum(0)
                tvf.append(0.5 * np.abs(s / s.sum() - q).sum())
                c = rng.multivariate_hypergeometric(np.round(truth.reindex(labels).to_numpy() * len(d)).astype(int), int(s.sum()))
                tvc.append(0.5 * np.abs(c / c.sum() - q).sum())
            samp.append(dict(section=sec, organ=organ, donor=donor, k=k, M=M, tv_fields=np.median(tvf), tv_cells=np.median(tvc)))
        # intervals and planning for labels with p >= 1%
        for i, lab in enumerate(labels):
            P = truth[lab]
            if P < 0.01 or (focal and lab not in focal): continue
            y = X[:, i]
            for k in [3, 5, 10]:
                if k >= M: continue
                hit = {'Wilson': 0, 'ratio_t': 0, 'quasibinomial': 0}
                for _ in range(200):
                    ix = rng.choice(M, k, replace=False); yy, nn = y[ix], n[ix]
                    for m, (l, u) in {'Wilson': wilson(yy.sum(), nn.sum()), 'ratio_t': ratio_t(yy, nn, M), 'quasibinomial': quasibin(yy, nn, M)}.items():
                        hit[m] += l <= P <= u
                for m in hit:
                    ints.append(dict(section=sec, organ=organ, label=lab, k=k, method=m, coverage=hit[m] / 200))
            for pr in range(20):
                perm = rng.permutation(M); pil, pool = perm[:5], perm[5:]
                yp, npl = y[pil], n[pil]; pp = yp.sum() / npl.sum()
                if not 0 < pp < 1: continue
                Sp = (yp - pp * npl).var(ddof=1); h = 0.25 * P
                row = dict(section=sec, organ=organ, label=lab, p_true=P, pilot=pr)
                for name, kk in {'cells': fovplan.naive_required_fields(pp, npl.mean(), .25, True),
                                 'fovplan': fovplan.required_fields(pp, Sp, npl.mean(), .25, True, M=len(pool)),
                                 'fovplan80': fovplan.required_fields(pp, fovplan.se2_upper(Sp, 5, .8), npl.mean(), .25, True, M=len(pool))}.items():
                    if not np.isfinite(kk) or kk > len(pool):
                        row[name + '_k'] = np.inf; row[name + '_assur'] = np.nan; continue
                    ss = np.array([rng.choice(len(pool), int(kk), replace=False) for _ in range(200)])
                    est = y[pool][ss].sum(1) / n[pool][ss].sum(1)
                    row[name + '_k'] = kk; row[name + '_assur'] = np.mean(np.abs(est - P) <= h)
                plan.append(row)
    for name, rows in [('atlas', atlas), ('sampling', samp), ('intervals', ints), ('planning', plan)]:
        pd.DataFrame(rows).to_csv(f'{out}_{name}.tsv', sep='\t', index=False)
    return pd.DataFrame(atlas), pd.DataFrame(samp), pd.DataFrame(ints), pd.DataFrame(plan)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cells'); ap.add_argument('--field', type=float, default=500); ap.add_argument('--min-cells', type=int, default=100)
    ap.add_argument('--out', default='whole')
    a = ap.parse_args()
    f = a.cells
    cells = pd.read_parquet(f) if f.endswith('.parquet') else pd.read_csv(f, sep='\t' if f.endswith(('.tsv', '.tsv.gz')) else ',')
    at, sa, it, pl = analyse(cells, a.field, a.min_cells, a.out)
    print(at.groupby('organ').D.median(), sa.groupby(['organ', 'k'])[['tv_fields', 'tv_cells']].median(), sep='\n')
