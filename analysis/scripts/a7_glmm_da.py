"""A7 (exploratory, post-G25). Mixed models and commonly used differential-abundance tools under the
same null comparisons as A3/C4 and in the B2-type case-control simulation.

New comparators (see da_models.py): glmm_ri (LRT, Wald), glmm_olre (LRT, Wald), betabin, nb_glm,
propeller (logit, asin). The six A3 tests are recomputed on identical data for reference.

usage:
  a7_glmm_da.py kidney  OUT.tsv.gz                 # DKD 364 + lupus 70 null splits, real contrasts
  a7_glmm_da.py hickey  OUT.tsv.gz SITE_INDEX      # one intestinal site (0-6), 35 donor splits
  a7_glmm_da.py sim     OUT.tsv.gz REP_START REP_END
Seeds: kidney 20261080, hickey 20261081, sim 20261082 + rep.
"""
import sys, os, itertools, warnings
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, '.')
import da_models as M
from b2_simulation_da import tests as six_tests, TESTS as SIX

warnings.filterwarnings('ignore')
NEW = ['glmm_ri_lrt', 'glmm_ri_wald', 'glmm_olre_lrt', 'glmm_olre_wald', 'betabin', 'nb_glm', 'propeller_logit', 'propeller_asin']
ALL = SIX + NEW


def model_tests(yf, nf, ind, grp):
    """Field-level inputs for one category. Returns dict of the 6 A3 tests and the field-input GLMMs."""
    out = dict(zip(SIX, six_tests(yf, nf, ind, grp)))
    a, b, _ = M.glmm_ri(yf, nf, ind, grp); out['glmm_ri_lrt'], out['glmm_ri_wald'] = a, b
    a, b, _ = M.glmm_olre(yf, nf, ind, grp); out['glmm_olre_lrt'], out['glmm_olre_wald'] = a, b
    I = pd.DataFrame(dict(y=yf, n=nf, i=ind, g=grp)).groupby('i').agg(y=('y', 'sum'), n=('n', 'sum'), g=('g', 'first'))
    out['betabin'] = M.betabin(I.y.to_numpy(float), I.n.to_numpy(float), I.g.to_numpy(float))
    out['nb_glm'] = M.nb_glm(I.y.to_numpy(float), I.n.to_numpy(float), I.g.to_numpy(float))
    return out


def one_contrast(F, labcols, cats, lineage_of, g1, ind_col):
    """F: field table with all native label columns; cats: tested columns (natives and '[lineage] X').
    g1: set of individuals in group 1. Returns list of dicts (one per tested category)."""
    ind = F[ind_col].to_numpy(); grp = np.isin(ind, list(g1)).astype(int)
    nf = F[labcols].sum(1).to_numpy(float)
    # propeller inputs: per-individual counts of all mutually exclusive clusters
    P = F.groupby(ind_col)[labcols].sum()
    gi = P.index.isin(list(g1)).astype(float)
    prop = {}
    for tr in ('logit', 'asin'):
        pv = M.propeller(P.to_numpy(), gi, tr); prop.update({(tr, c): p for c, p in zip(labcols, pv)})
        if lineage_of is not None:
            L = P.T.groupby(lineage_of.reindex(P.columns).to_numpy()).sum().T
            pv = M.propeller(L.to_numpy(), gi, tr); prop.update({(tr, '[lineage] ' + c): p for c, p in zip(L.columns, pv)})
    rows = []
    for c in cats:
        if c.startswith('[lineage] '):
            yf = F[[x for x in labcols if lineage_of[x] == c[10:]]].sum(1).to_numpy(float)
        else:
            yf = F[c].to_numpy(float)
        r = model_tests(yf, nf, ind, grp)
        r['propeller_logit'], r['propeller_asin'] = prop[('logit', c)], prop[('asin', c)]
        r['category'] = c
        rows.append(r)
    return rows


# ------------------------------------------------------------------------------------- kidney
def kidney_job(args):
    name, scen, s_i, s, F, labcols, cats, lineage_of = args
    rows = one_contrast(F, labcols, cats, lineage_of, set(s), 'patient')
    for r in rows:
        r.update(source=name, scenario=scen, split=s_i)
    return rows


def run_kidney(out):
    from common import load_primary, load_lupus
    rng = np.random.default_rng(20261080)
    jobs, real = [], []
    for name, loader, disease, n1 in [('GSE325587', load_primary, 'DKD', 3), ('Lupus', load_lupus, 'SLE', 4)]:
        long = loader()
        lin = long.groupby('label').lineage.first()
        W = long.pivot_table(index=['patient', 'tissue', 'field'], columns='label', values='count', aggfunc='sum', fill_value=0)
        labcols = list(W.columns)
        overall = W.sum(0) / W.to_numpy().sum()
        cats = list(overall[overall >= 0.005].index) + ['[lineage] ' + l for l in sorted(set(lin))]
        F = W.reset_index()
        grp_of = long.groupby('patient').group.first()
        dis = sorted(grp_of[grp_of == disease].index)
        # real contrast
        rows = one_contrast(F, labcols, cats, lin, set(dis), 'patient')
        for r in rows:
            r.update(source=name, scenario='real_contrast', split=-1)
        real += rows
        Fd = F[F.patient.isin(dis)]
        sub3 = pd.concat([g.iloc[rng.choice(len(g), min(3, len(g)), replace=False)] for _, g in Fd.groupby('tissue')])
        for s_i, s in enumerate(itertools.combinations(dis, n1)):
            for scen, FF in [('all_fields', Fd), ('3_fields', sub3)]:
                jobs.append((name, scen, s_i, s, FF, labcols, cats, lin))
        print(name, 'categories', len(cats), 'disease individuals', len(dis), flush=True)
    with Pool(2) as pool:
        res = pool.map(kidney_job, jobs, chunksize=8)
    df = pd.DataFrame([r for rr in res for r in rr] + real)
    df.to_csv(out, sep='\t', index=False, compression='gzip')


# ------------------------------------------------------------------------------------- intestine
def hickey_job(args):
    site, scen, s_i, s, F, labcols, cats = args
    rows = one_contrast(F, labcols, cats, None, set(s), 'donor')
    for r in rows:
        r.update(site=site, scenario=scen, split=s_i)
    return rows


def run_hickey(out, site_index):
    import c2_whole_section as c2
    PX = 0.37744
    d = pd.read_pickle(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'hickey_7sites.pkl')).dropna(subset=['cell_type', 'unique_region', 'x', 'y'])
    sites = sorted(d.Tissue_location.unique())
    site = sites[site_index]
    d = d[d.Tissue_location == site]
    cells = pd.DataFrame(dict(section=d.unique_region, donor=d.donor, x_um=d.x * PX, y_um=d.y * PX, label=d.cell_type))
    cells['field'] = pd.concat([c2.tile(g, 500) for _, g in cells.groupby('section')]).reindex(cells.index)
    W = cells.pivot_table(index=['donor', 'section', 'field'], columns='label', values='x_um', aggfunc='size', fill_value=0)
    W = W[W.sum(1) >= 100]
    labcols = list(W.columns)
    prop = W.sum(0) / W.to_numpy().sum()
    cats = list(prop[prop >= 0.01].index)
    donors = sorted(W.index.get_level_values('donor').unique())
    W = W.reset_index()
    rng = np.random.default_rng(20261081 + site_index)
    sub10 = W.iloc[rng.permutation(len(W))].groupby('section').head(10)
    jobs = []
    for scen, F in [('all_fields', W), ('10_fields', sub10)]:
        for s_i, s in enumerate(itertools.combinations(donors, len(donors) // 2)):
            if s[0] != donors[0]:
                continue
            jobs.append((site, scen, s_i, s, F, labcols, cats))
    with Pool(2) as pool:
        res = pool.map(hickey_job, jobs, chunksize=2)
    pd.DataFrame([r for rr in res for r in rr]).to_csv(out, sep='\t', index=False, compression='gzip')
    print(site, 'donors', len(donors), 'categories', len(cats), 'jobs', len(jobs), flush=True)


# ------------------------------------------------------------------------------------- simulation
def section_all(sigma, ell, p0, nx, ny, rng, shift, K=6):
    """As sim_core.section, but returns counts of all K types per field (fields x K)."""
    from sim_core import grf, FIELD, LAM
    H, Wd = ny * FIELD, nx * FIELD
    b = np.r_[p0, np.full(K - 1, (1 - p0) / (K - 1))]
    eta = np.stack([np.log(b[t]) + sigma * grf((H, Wd), ell, rng) for t in range(K)])
    eta[0] += shift
    eta -= eta.max(0)
    pi = np.exp(eta); pi /= pi.sum(0)
    tot = rng.poisson(LAM, (H, Wd))
    # multinomial per pixel via sequential binomials
    counts = np.zeros((K, H, Wd)); rem = tot.copy(); rest = np.ones((H, Wd))
    for t in range(K - 1):
        q = np.clip(pi[t] / np.maximum(rest, 1e-12), 0, 1)
        c = rng.binomial(rem, q); counts[t] = c; rem = rem - c; rest = rest - pi[t]
    counts[K - 1] = rem
    return counts.reshape(K, ny, FIELD, nx, FIELD).sum((2, 4)).reshape(K, -1).T


def sim_job(rep):
    rng = np.random.default_rng(20261082 + rep)
    rows = []
    for delta in (0.0, 0.4):
        secs = {}
        for g in (0, 1):
            for i in range(12):
                u = rng.normal(0, 0.3) + (delta if g == 1 else 0)
                C = section_all(1.0, 200, 0.1, 6, 6, rng, u)
                secs[(g, i)] = C[rng.permutation(len(C))]
        for nper in (3, 6, 12):
            for k in (3, 10, 36):
                Cs, ind, grp = [], [], []
                for g in (0, 1):
                    for i in range(nper):
                        Cs.append(secs[(g, i)][:k]); ind += [f'{g}_{i}'] * k; grp += [g] * k
                C = np.vstack(Cs); ind = np.array(ind); grp = np.array(grp)
                r = model_tests(C[:, 0].astype(float), C.sum(1).astype(float), ind, grp)
                P = pd.DataFrame(C).groupby(ind).sum()
                gi = np.array([1.0 if x.startswith('1_') else 0.0 for x in P.index])
                r['propeller_logit'] = M.propeller(P.to_numpy(), gi, 'logit')[0]
                r['propeller_asin'] = M.propeller(P.to_numpy(), gi, 'asin')[0]
                r.update(delta=delta, n_per_group=nper, k=k, rep=rep)
                rows.append(r)
    return rows


def run_sim(out, a, b):
    with Pool(2) as pool:
        res = pool.map(sim_job, range(a, b), chunksize=2)
    pd.DataFrame([r for rr in res for r in rr]).to_csv(out, sep='\t', index=False, compression='gzip')


if __name__ == '__main__':
    mode, out = sys.argv[1], sys.argv[2]
    if mode == 'kidney':
        run_kidney(out)
    elif mode == 'hickey':
        run_hickey(out, int(sys.argv[3]))
    elif mode == 'sim':
        run_sim(out, int(sys.argv[3]), int(sys.argv[4]))
