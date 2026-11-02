"""A11 (exploratory). Field size and imaging cost against whole-tissue truth.

Each section or region is tiled into square virtual fields of 250, 500 or 1,000 um (mouse: stored
coordinate units); fields with fewer than 25, 100 or 400 cells (area-proportional thresholds) are not
sampled. For cell types with a whole-tissue proportion >= 1%:
  design effect D and mean cells per field for each field size;
  area-matched precision: for imaging budgets equivalent to 1, 2 and 4 fields of 1,000 um (1, 2 and 4 mm2),
  300 random samples of 1/2/4 large, 4/8/16 medium or 16/32/64 small fields; share of estimates within
  +/-25% of the whole-tissue proportion and relative root-mean-square error;
  area required by the fovplan point plan (+/-25%), computed with all eligible fields as the pilot.
Seed 20261120.
"""
import sys, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, '.')
import fovplan, wt_common as W

SIZES = [(250.0, 25, 16), (500.0, 100, 4), (1000.0, 400, 1)]   # (F, min cells, fields per 1 mm2-equivalent)
BUDGETS = [1, 2, 4]


def job(args):
    sec, g, seed = args
    rng = np.random.default_rng(seed)
    rows, arows = [], []
    for F, mn, per in SIZES:
        tabs = list(W.field_tables(g, F, mn, min_fields=4))
        if not tabs:
            continue
        _, organ, truth, Wt, fx, fy = tabs[0]
        X = Wt.to_numpy(float); n = X.sum(1); M = len(n)
        for j, lab in enumerate(Wt.columns):
            P = truth[lab]
            if P < 0.01:
                continue
            y = X[:, j]; pf = y.sum() / n.sum(); Se2 = (y - pf * n).var(ddof=1)
            D = Se2 / (n.mean() * pf * (1 - pf)) if 0 < pf < 1 else np.nan
            k = fovplan.required_fields(pf, Se2, n.mean(), 0.25, True, M=M) if 0 < pf < 1 else np.nan
            area = k * (F / 1000) ** 2 if np.isfinite(k) else np.nan
            arows.append(dict(section=sec, organ=organ, F=F, label=lab, p=P, M=M, nbar=n.mean(), D=D, k_req=k, area_req_mm2=area))
            for b in BUDGETS:
                kk = b * per
                if kk > M:
                    continue
                ss = np.array([rng.choice(M, kk, replace=False) for _ in range(300)])
                est = y[ss].sum(1) / n[ss].sum(1)
                rows.append(dict(section=sec, organ=organ, F=F, label=lab, p=P, budget_mm2=b, k=kk, M=M,
                                 share=np.mean(np.abs(est - P) <= 0.25 * P), rrmse=np.sqrt(np.mean((est - P) ** 2)) / P))
    return rows, arows


if __name__ == '__main__':
    out = sys.argv[1]
    cells = W.load_cells()
    tasks = [(sec, g, 20261120 + i) for i, (sec, g) in enumerate(cells.groupby('section'))]
    with Pool(2) as pool:
        res = pool.map(job, tasks, chunksize=2)
    a = pd.DataFrame([r for rr in res for r in rr[0]]); b = pd.DataFrame([r for rr in res for r in rr[1]])
    a.to_csv(out + '_precision.tsv.gz', sep='\t', index=False, compression='gzip')
    b.to_csv(out + '_design.tsv.gz', sep='\t', index=False, compression='gzip')
    print(b.groupby(['organ', 'F']).agg(D=('D', 'median'), nbar=('nbar', 'median'), area=('area_req_mm2', 'median'), M=('M', 'median')).round(2))
    print(a.groupby(['organ', 'budget_mm2', 'F']).share.mean().unstack().round(3))
