"""A8 (exploratory, post-G26; for the CMPB version). Benchmark of field-number planners.

Planners (all from the same 5 pilot fields, +/-25% relative precision target for the tissue proportion):
  cellcount   cells treated as independent (design effect 1)
  fovplan     point plan (pilot S_e^2)
  fovplan80   80%-assurance plan (chi-square upper bound of S_e^2)
  fovplan90   90%-assurance plan
  bootstrap   resampling planner: smallest k such that the 95th percentile of |p* - p_pilot| over
              B = 400 bootstrap samples of k fields drawn with replacement from the pilot, multiplied by the
              finite-population factor sqrt(1 - k/M), is <= h
Evaluation: 150 random samples of the planned number of fields from the remaining fields; share whose
estimate lies within +/-25% of the reference proportion (target 95%); planned k; infeasible plans
(k > fields remaining) are recorded and excluded from the share. Wall-clock time per plan is recorded.
References: whole-section mouse kidney (2 sections) and human intestine (57 regions), kidney CosMx
tissues with >= 15 fields (observed-field reference), simulated sections with known truth.
10 pilot draws per tissue and cell type (cell types with reference proportion >= 1%). Seed 20261090.
usage: a8_planner_benchmark.py OUT.tsv.gz [part]   part in {kidney, mouse, intestine, sim}
"""
import sys, os, time, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, '.')
import fovplan
import c2_whole_section as c2

TARGET, M_PILOT, REPS, DRAWS, B = 0.25, 5, 10, 150, 400
METHODS = ['cellcount', 'fovplan', 'fovplan80', 'fovplan90', 'bootstrap']


def bootstrap_k(yp, np_, h, kmax, rng):
    """Finite-population-corrected percentile bootstrap planner (M = kmax fields available)."""
    M = kmax
    p0 = yp.sum() / np_.sum()
    def q95(k):
        ix = rng.integers(0, len(yp), size=(B, k))
        est = yp[ix].sum(1) / np_[ix].sum(1)
        return np.quantile(np.abs(est - p0), 0.95) * np.sqrt(max(1 - k / M, 0))
    if q95(1) <= h:
        return 1
    lo, hi = 1, 2
    while hi <= kmax and q95(hi) > h:
        lo, hi = hi, hi * 2
    if hi > kmax:
        if q95(kmax) > h:
            return np.inf
        hi = kmax
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if q95(mid) <= h:
            hi = mid
        else:
            lo = mid
    return hi


def plan_all(yp, np_, kmax, rng):
    pp = yp.sum() / np_.sum()
    if not 0 < pp < 1:
        return None
    Se2 = (yp - pp * np_).var(ddof=1); nbar = np_.mean(); h = TARGET * pp
    out, tm = {}, {}
    t = time.perf_counter(); out['cellcount'] = fovplan.naive_required_fields(pp, nbar, TARGET, True); tm['cellcount'] = time.perf_counter() - t
    t = time.perf_counter(); out['fovplan'] = fovplan.required_fields(pp, Se2, nbar, TARGET, True, M=kmax); tm['fovplan'] = time.perf_counter() - t
    for q, name in [(0.8, 'fovplan80'), (0.9, 'fovplan90')]:
        t = time.perf_counter(); out[name] = fovplan.required_fields(pp, fovplan.se2_upper(Se2, M_PILOT, q), nbar, TARGET, True, M=kmax); tm[name] = time.perf_counter() - t
    t = time.perf_counter(); out['bootstrap'] = bootstrap_k(yp, np_, h, kmax, rng); tm['bootstrap'] = time.perf_counter() - t
    return out, tm


def evaluate(unit, y, n, P, seed):
    """y, n: per-field counts of one category and totals (eligible fields); P: reference proportion."""
    rng = np.random.default_rng(seed)
    M = len(n); rows = []
    for rep in range(REPS):
        perm = rng.permutation(M); pil, pool = perm[:M_PILOT], perm[M_PILOT:]
        r = plan_all(y[pil], n[pil], len(pool), rng)
        if r is None:
            continue
        plans, tm = r
        for meth in METHODS:
            k = plans[meth]
            row = dict(unit, rep=rep, method=meth, k=k, seconds=tm[meth])
            if not np.isfinite(k) or k > len(pool):
                row.update(infeasible=True, share=np.nan)
            else:
                k = int(max(k, 1))
                ss = np.array([rng.choice(len(pool), k, replace=False) for _ in range(DRAWS)])
                est = y[pool][ss].sum(1) / n[pool][ss].sum(1)
                row.update(infeasible=False, share=np.mean(np.abs(est - P) <= TARGET * P))
            rows.append(row)
    return rows


def job(a):
    return evaluate(*a)


def tasks_wholetissue(cells, organ_of, seed0):
    cells = cells.copy()
    cells['field'] = pd.concat([c2.tile(g, 500) for _, g in cells.groupby('section')]).reindex(cells.index)
    tasks = []
    for si, (sec, g) in enumerate(cells.groupby('section')):
        truth = g.label.value_counts(normalize=True)
        W = g.pivot_table(index='field', columns='label', values='x_um', aggfunc='size', fill_value=0)
        W = W[W.sum(axis=1) >= 100]; X = W.to_numpy(float); n = X.sum(1)
        if len(n) < 20:
            continue
        for j, lab in enumerate(W.columns):
            if truth[lab] < 0.01:
                continue
            tasks.append((dict(source=organ_of(g), unit=sec, label=lab), X[:, j], n, truth[lab], seed0 + 1000 * si + j))
    return tasks


def tasks_mouse():
    cells = []
    for s in ['ShamR', 'Day14R']:
        m = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', f'mouse_example_{s}.tsv.gz'), sep='\t')
        cells.append(pd.DataFrame(dict(section=s, x_um=m.x, y_um=m.y, label=m.celltype)))
    return tasks_wholetissue(pd.concat(cells, ignore_index=True), lambda g: 'mouse kidney', 20261090)


def tasks_intestine():
    d = pd.read_pickle(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'hickey_7sites.pkl')).dropna(subset=['cell_type', 'unique_region', 'x', 'y'])
    cells = pd.DataFrame(dict(section=d.unique_region, organ=d.tissue.map({'SB': 'human small intestine', 'CL': 'human colon'}),
                              x_um=d.x * 0.37744, y_um=d.y * 0.37744, label=d.cell_type))
    return tasks_wholetissue(cells, lambda g: g.organ.iloc[0], 20262090)


def tasks_kidney():
    from common import load_all, wide
    data = load_all(); tasks = []
    for ti, ((src, tis), g) in enumerate(data.groupby(['source', 'tissue'])):
        lin = g.groupby('label').lineage.first()
        (labels, X, fields), = wide(g).values()
        n = X.sum(1)
        if len(n) < 15:
            continue
        lineages = sorted(set(lin))
        L = np.stack([X[:, [labels.index(l) for l in labels if lin[l] == ln]].sum(1) for ln in lineages], 1)
        Y = np.concatenate([X, L], 1); names = labels + ['[lineage] ' + l for l in lineages]
        pref = Y.sum(0) / n.sum()
        for j in np.nonzero(pref >= 0.01)[0]:
            tasks.append((dict(source='human kidney (' + src + ')', unit=tis, label=names[j]), Y[:, j].astype(float), n.astype(float), pref[j], 20263090 + 1000 * ti + j))
    return tasks


def tasks_sim():
    from sim_core import section
    rng = np.random.default_rng(20264090); tasks = []
    for sigma in (0.5, 1.0, 1.5):
        for s in range(10):
            y, n = section(sigma, 200, 0.1, rng=rng)          # 96 fields of ~1,000 cells
            tasks.append((dict(source=f'simulated (sigma = {sigma})', unit=f'{sigma}_{s}', label='focal'), y, n, y.sum() / n.sum(),
                          20265090 + int(sigma * 100) + s))
    return tasks


if __name__ == '__main__':
    out = sys.argv[1]; part = sys.argv[2]
    tasks = dict(kidney=tasks_kidney, mouse=tasks_mouse, intestine=tasks_intestine, sim=tasks_sim)[part]()
    with Pool(2) as pool:
        res = pool.map(job, tasks, chunksize=4)
    df = pd.DataFrame([r for rr in res for r in rr]); df.to_csv(out, sep='\t', index=False, compression='gzip')
    print(part, 'tasks', len(tasks), 'rows', len(df))
    print(df.groupby(['source', 'method']).agg(share=('share', 'mean'), k_median=('k', 'median'), infeasible=('infeasible', 'mean')).round(3))
