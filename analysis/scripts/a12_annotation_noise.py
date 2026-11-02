"""A12 (exploratory). Sensitivity to cell-type annotation error.

Perturbations of the field-by-cell-type count tables:
  none           original labels
  random 5%      each cell relabelled with probability 0.05 to a label drawn from the tissue's composition
  random 10%     as above with probability 0.10
  field-varying  each field has its own relabelling probability drawn from Uniform(0, 0.20)
                 (heterogeneous segmentation or staining quality between fields)
Outcomes: design effects (kidney sources and intestine; cell types >= 1%); rejection rates of the six
A3 tests under the kidney null splits (DKD 364, lupus 70; all fields); and the share of held-out samples
within +/-25% for fovplan point and 90%-assurance plans in the 16 kidney tissues with >= 15 fields
(five pilot fields, 10 pilot draws, 100 held-out samples; reference = perturbed complete collection).
Seed 20261130.
"""
import sys, itertools, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, '.')
import fovplan, wt_common as W
from b2_simulation_da import tests as six_tests, TESTS as SIX

PERT = [('none', None), ('random 5%', 0.05), ('random 10%', 0.10), ('field-varying 0-20%', 'field')]


def perturb(X, eps, rng):
    """X: fields x labels integer counts. Returns perturbed counts with the same field totals."""
    if eps is None:
        return X.copy()
    X = X.astype(int); comp = X.sum(0) / X.sum()
    e = rng.uniform(0, 0.2, size=len(X)) if eps == 'field' else np.full(len(X), eps)
    out = np.zeros_like(X)
    for j in range(len(X)):
        flip = rng.binomial(X[j], e[j])
        out[j] = X[j] - flip + rng.multinomial(flip.sum(), comp)
    return out


def design_effects(X, keep):
    n = X.sum(1); res = []
    for j in keep:
        p = X[:, j].sum() / n.sum()
        if 0 < p < 1:
            res.append((X[:, j] - p * n).var(ddof=1) / (n.mean() * p * (1 - p)))
    return res


def kidney_tables():
    from common import load_primary, load_lupus, load_gse211785
    out = []
    for name, loader, disease, n1 in [('GSE325587', load_primary, 'DKD', 3), ('Lupus', load_lupus, 'SLE', 4), ('GSE211785', load_gse211785, None, None)]:
        long = loader()
        Wt = long.pivot_table(index=['patient', 'tissue', 'field'], columns='label', values='count', aggfunc='sum', fill_value=0)
        grp = long.groupby('patient').group.first()
        out.append((name, Wt, grp, disease, n1))
    return out


def job_kidney(args):
    name, Wt, grp, disease, n1, pert, eps, seed = args
    rng = np.random.default_rng(seed)
    labs = list(Wt.columns); idx = Wt.index.to_frame(index=False)
    X = np.zeros(Wt.shape, int)
    for t, rows in idx.groupby('tissue').groups.items():       # perturb within tissue
        X[rows] = perturb(Wt.to_numpy()[rows], eps, rng)
    rows_D, rows_fp, rows_pl = [], [], []
    for t, rr in idx.groupby('tissue').groups.items():
        Xt = X[rr]; n = Xt.sum(1); comp = Xt.sum(0) / n.sum()
        keep = [j for j in range(len(labs)) if comp[j] >= 0.01]
        for d in design_effects(Xt, keep):
            rows_D.append(dict(source=name, perturbation=pert, D=d))
        if len(n) >= 15:
            for j in keep:
                y = Xt[:, j].astype(float); P = comp[j]
                for rep in range(10):
                    perm = rng.permutation(len(n)); pil, pool = perm[:5], perm[5:]
                    pp = y[pil].sum() / n[pil].sum()
                    if not 0 < pp < 1:
                        continue
                    Se2 = (y[pil] - pp * n[pil]).var(ddof=1)
                    for meth, S in [('fovplan', Se2), ('fovplan90', fovplan.se2_upper(Se2, 5, 0.9))]:
                        k = fovplan.required_fields(pp, S, n[pil].mean(), 0.25, True, M=len(pool))
                        if not np.isfinite(k) or k > len(pool):
                            continue
                        ss = np.array([rng.choice(len(pool), int(k), replace=False) for _ in range(100)])
                        est = y[pool][ss].sum(1) / n[pool][ss].sum(1)
                        rows_pl.append(dict(source=name, perturbation=pert, method=meth, share=np.mean(np.abs(est - P) <= 0.25 * P)))
    if disease is not None:
        dis = sorted(grp[grp == disease].index)
        m = idx.patient.isin(dis).to_numpy()
        Xd = X[m]; ind = idx.patient.to_numpy()[m]; nf = Xd.sum(1).astype(float)
        comp = Xd.sum(0) / nf.sum(); keep = [j for j in range(len(labs)) if comp[j] >= 0.005]
        for s_i, s in enumerate(itertools.combinations(dis, n1)):
            g = np.isin(ind, s).astype(int)
            for j in keep:
                ps = six_tests(Xd[:, j].astype(float), nf, ind, g)
                rows_fp.append(dict(source=name, perturbation=pert, split=s_i, label=labs[j], **dict(zip(SIX, ps))))
    return rows_D, rows_fp, rows_pl


def job_intestine(args):
    sec, organ, Wt, pert, eps, seed = args
    rng = np.random.default_rng(seed)
    X = perturb(Wt.to_numpy(), eps, rng); comp = X.sum(0) / X.sum()
    return [dict(source=organ, perturbation=pert, D=d) for d in design_effects(X, [j for j in range(X.shape[1]) if comp[j] >= 0.01])]


if __name__ == '__main__':
    out = sys.argv[1]
    kt = kidney_tables()
    jobs = [(name, Wt, grp, dis, n1, p, e, 20261130 + 10 * i + k) for i, (name, Wt, grp, dis, n1) in enumerate(kt) for k, (p, e) in enumerate(PERT)]
    cells = W.load_cells(); cells = cells[cells.organ != 'mouse kidney']
    ij = [(sec, organ, Wt, p, e, 20262130 + 10 * i + k) for i, (sec, organ, truth, Wt, fx, fy) in enumerate(W.field_tables(cells)) for k, (p, e) in enumerate(PERT)]
    with Pool(2) as pool:
        rk = pool.map(job_kidney, jobs, chunksize=1)
        ri = pool.map(job_intestine, ij, chunksize=8)
    D = pd.DataFrame([r for rr in rk for r in rr[0]] + [r for rr in ri for r in rr])
    FP = pd.DataFrame([r for rr in rk for r in rr[1]]); PL = pd.DataFrame([r for rr in rk for r in rr[2]])
    D.to_csv(out + '_design.tsv.gz', sep='\t', index=False, compression='gzip')
    FP.to_csv(out + '_fp.tsv.gz', sep='\t', index=False, compression='gzip')
    PL.to_csv(out + '_planning.tsv.gz', sep='\t', index=False, compression='gzip')
    print(D.groupby(['source', 'perturbation']).D.median().unstack().round(1))
    print((100 * FP.groupby(['source', 'perturbation'])[SIX].apply(lambda d: (d < 0.05).mean())).round(1))
    print((100 * PL.groupby(['perturbation', 'method']).share.mean()).unstack().round(1))
