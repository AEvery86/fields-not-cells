"""A9 (exploratory). Non-random (preferential) field selection against whole-tissue truth.

Schemes for choosing k fields among the eligible 500-um virtual fields of a section or region:
  random       simple random sampling (reference)
  systematic   every (M/k)-th field along a serpentine path with a random start
  interior     simple random sampling restricted to fields whose four neighbours are all eligible fields
  dense        probability proportional to cell count (cellular areas preferred), without replacement
  enriched     probability proportional to exp(2 z), z = standardized proportion of the region's most
               compartmentalized cell type (largest design effect among types >= 1%): lesion- or
               structure-enriched selection
Outcomes (k = 3, 5, 10; 300 draws): median total variation distance from the whole-tissue composition;
relative bias and quasi-binomial interval coverage for cell types with proportion >= 1%; and fovplan
90%-assurance planning (5 pilot fields and held-out fields both chosen by the scheme; 10 pilot draws,
100 held-out samples; +/-25% of the whole-tissue proportion). Seed 20261100.
"""
import sys, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, '.')
import fovplan, wt_common as W
from c2_whole_section import quasibin

SCHEMES = ['random', 'systematic', 'interior', 'dense', 'enriched']
KS = [3, 5, 10]


def serpentine(fx, fy):
    return np.lexsort((np.where(fy % 2 == 0, fx, -fx), fy))


def chooser(scheme, fx, fy, n, X, focal_j, rng):
    M = len(n)
    if scheme == 'interior':
        S = set(zip(fx, fy))
        ok = np.array([all((a + dx, b + dy) in S for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]) for a, b in zip(fx, fy)])
        elig = np.nonzero(ok)[0]
        return lambda k, pool=None: (rng.choice(np.intersect1d(elig, pool) if pool is not None else elig, k, replace=False)
                                     if len(np.intersect1d(elig, pool) if pool is not None else elig) >= k else None)
    if scheme == 'systematic':
        order = serpentine(fx, fy)
        def f(k, pool=None):
            o = order if pool is None else order[np.isin(order, pool)]
            if len(o) < k:
                return None
            step = len(o) / k; s = rng.uniform(0, step)
            return o[np.floor(s + step * np.arange(k)).astype(int)]
        return f
    if scheme == 'random':
        w = np.ones(M)
    elif scheme == 'dense':
        w = n.astype(float)
    else:
        pf = X[:, focal_j] / n; z = (pf - pf.mean()) / (pf.std() + 1e-12); w = np.exp(2 * z)
    def f(k, pool=None):
        idx = np.arange(M) if pool is None else np.asarray(pool)
        if len(idx) < k:
            return None
        pw = w[idx] / w[idx].sum()
        return rng.choice(idx, k, replace=False, p=pw)
    return f


def analyse(args):
    sec, organ, truth, Wt, fx, fy, seed = args
    rng = np.random.default_rng(seed)
    X = Wt.to_numpy(float); n = X.sum(1); M = len(n); labs = list(Wt.columns)
    q = truth.reindex(labs).to_numpy()
    keep = [j for j, l in enumerate(labs) if truth[l] >= 0.01]
    D = {}
    for j in keep:
        p = X[:, j].sum() / n.sum(); D[j] = (X[:, j] - p * n).var(ddof=1) / (n.mean() * p * (1 - p))
    focal_j = max(D, key=D.get)
    comp, cov, plan = [], [], []
    for sc in SCHEMES:
        pick = chooser(sc, fx, fy, n, X, focal_j, rng)
        for k in KS:
            if k >= M - 5:
                continue
            tv, est, hit, ok = [], {j: [] for j in keep}, {j: 0 for j in keep}, 0
            for _ in range(300):
                ix = pick(k)
                if ix is None:
                    break
                ok += 1
                s = X[ix].sum(0); tv.append(0.5 * np.abs(s / s.sum() - q).sum())
                for j in keep:
                    est[j].append(s[j] / s.sum())
                    lo, hi = quasibin(X[ix, j], n[ix], M); hit[j] += lo <= truth[labs[j]] <= hi
            if ok == 0:
                continue
            comp.append(dict(section=sec, organ=organ, scheme=sc, k=k, M=M, tv=np.median(tv)))
            for j in keep:
                P = truth[labs[j]]
                cov.append(dict(section=sec, organ=organ, scheme=sc, k=k, label=labs[j], focal=j == focal_j, p=P,
                                rel_bias=(np.mean(est[j]) - P) / P, coverage=hit[j] / ok))
        # planning: pilot and held-out fields chosen by the same scheme
        for j in keep:
            P = truth[labs[j]]
            for rep in range(10):
                pil = pick(5)
                if pil is None:
                    break
                pool = np.setdiff1d(np.arange(M), pil)
                yp, npl = X[pil, j], n[pil]; pp = yp.sum() / npl.sum()
                if not 0 < pp < 1:
                    continue
                Se2 = (yp - pp * npl).var(ddof=1)
                kk = fovplan.required_fields(pp, fovplan.se2_upper(Se2, 5, 0.9), npl.mean(), 0.25, True, M=len(pool))
                if not np.isfinite(kk) or kk > len(pool):
                    plan.append(dict(section=sec, organ=organ, scheme=sc, label=labs[j], focal=j == focal_j, k=np.inf, share=np.nan)); continue
                hits, tot = 0, 0
                for _ in range(100):
                    ix = pick(int(kk), pool)
                    if ix is None:
                        break
                    e = X[ix, j].sum() / n[ix].sum(); hits += abs(e - P) <= 0.25 * P; tot += 1
                plan.append(dict(section=sec, organ=organ, scheme=sc, label=labs[j], focal=j == focal_j, k=kk, share=hits / tot if tot else np.nan))
    return comp, cov, plan


if __name__ == '__main__':
    out = sys.argv[1]
    cells = W.load_cells()
    tasks = [(sec, organ, truth, Wt, fx, fy, 20261100 + i) for i, (sec, organ, truth, Wt, fx, fy) in enumerate(W.field_tables(cells))]
    with Pool(2) as pool:
        res = pool.map(analyse, tasks, chunksize=2)
    for k, name in enumerate(['composition', 'coverage', 'planning']):
        pd.DataFrame([r for rr in res for r in rr[k]]).to_csv(f'{out}_{name}.tsv.gz', sep='\t', index=False, compression='gzip')
    print('sections', len(tasks))
