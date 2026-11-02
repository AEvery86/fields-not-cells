"""A10 (exploratory). Sensitivity of the planner benchmark (A8) to the precision target and pilot size.

Settings: relative targets +/-10% and +/-50% and an absolute target +/-0.02 with five pilot fields, and
pilot sizes of 3, 8 and 10 fields with the +/-25% relative target (the A8 setting, +/-25% with five
pilot fields, is not repeated). Same references, planners and evaluation as A8, with 10 pilot draws and
100 held-out samples per plan. Seed 20261110.
usage: a10_planner_sensitivity.py OUT.tsv.gz
"""
import sys, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, '.')
import fovplan
import a8_planner_benchmark as A8

SETTINGS = [('rel', 0.10, 5), ('rel', 0.50, 5), ('abs', 0.02, 5), ('rel', 0.25, 3), ('rel', 0.25, 8), ('rel', 0.25, 10)]
METHODS = ['cellcount', 'bootstrap', 'fovplan', 'fovplan80', 'fovplan90']


def plans(yp, np_, kmax, kind, target, m, rng):
    pp = yp.sum() / np_.sum()
    if not 0 < pp < 1:
        return None
    rel = kind == 'rel'; h = target * pp if rel else target
    Se2 = (yp - pp * np_).var(ddof=1); nbar = np_.mean()
    out = dict(cellcount=fovplan.naive_required_fields(pp, nbar, target, rel),
               fovplan=fovplan.required_fields(pp, Se2, nbar, target, rel, M=kmax),
               fovplan80=fovplan.required_fields(pp, fovplan.se2_upper(Se2, m, 0.8), nbar, target, rel, M=kmax),
               fovplan90=fovplan.required_fields(pp, fovplan.se2_upper(Se2, m, 0.9), nbar, target, rel, M=kmax),
               bootstrap=A8.bootstrap_k(yp, np_, h, kmax, rng))
    return out, h


def job(task):
    unit, y, n, P, seed = task
    rng = np.random.default_rng(seed + 7)
    M = len(n); rows = []
    for kind, target, m in SETTINGS:
        if M - m < 5:
            continue
        for rep in range(10):
            perm = rng.permutation(M); pil, pool = perm[:m], perm[m:]
            r = plans(y[pil], n[pil], len(pool), kind, target, m, rng)
            if r is None:
                continue
            pl, h_pilot = r
            hP = target * P if kind == 'rel' else target
            for meth in METHODS:
                k = pl[meth]
                row = dict(unit, setting=f'{kind} {target} / {m} pilot', method=meth, k=k, pool=len(pool))
                if not np.isfinite(k) or k > len(pool):
                    row.update(infeasible=True, share=np.nan)
                else:
                    k = int(max(k, 1))
                    ss = np.array([rng.choice(len(pool), k, replace=False) for _ in range(100)])
                    est = y[pool][ss].sum(1) / n[pool][ss].sum(1)
                    row.update(infeasible=False, share=np.mean(np.abs(est - P) <= hP))
                rows.append(row)
    return rows


if __name__ == '__main__':
    tasks = A8.tasks_sim() + A8.tasks_mouse() + A8.tasks_kidney() + A8.tasks_intestine()
    with Pool(2) as pool:
        res = pool.map(job, tasks, chunksize=4)
    df = pd.DataFrame([r for rr in res for r in rr])
    df.to_csv(sys.argv[1], sep='\t', index=False, compression='gzip')
    print(df.groupby(['setting', 'method']).agg(share=('share', 'mean'), infeasible=('infeasible', 'mean')).unstack().round(3))
