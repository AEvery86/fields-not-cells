"""G21-B1. Simulation with known truth.

Synthetic tissue sections are generated on a 25-um pixel grid. Cell density is 4,000 cells/mm^2
(2.5 expected cells per pixel). Six cell types have log-intensities  eta_t = log b_t + sigma * Z_t,
where Z_t are independent unit-variance Gaussian random fields (white noise smoothed with a Gaussian
kernel of SD ell) and type probabilities are softmax(eta). Cell counts per pixel and type are Poisson
(thinning), so field counts are sums over pixels. Fields are 500 x 500 um (20 x 20 pixels, about
1,000 cells); a 6 x 4 mm section is tiled by 96 fields, and the true proportion is the whole-section
proportion. The focal type has target mean proportion p0.

Experiments (seed 20261050):
  S1  design effect: true D from all fields vs empirical variance ratio of 3-field samples;
      pilot (5 fields) estimate vs true D.
  S2  fovplan validation against the true section proportion: plan k for +/-25% relative precision
      from 5 random pilot fields (point and 80%-assurance), then draw k new fields 200 times;
      assurance = share within the target half-width. Cell-count planning shown for comparison.
  S3  interval coverage of the true section proportion at k = 3, 5, 10 (Wilson, ratio t, Fieller,
      quasi-binomial).
Scenario grid: sigma in {0, 0.5, 1, 1.5}; ell in {50, 200, 800} um; p0 in {0.02, 0.1, 0.3};
20 sections per scenario.
"""
import sys, os, itertools
import numpy as np, pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'data'))
from interval_methods import wilson
from region_methods import fieller_coefficients, polynomial_set
import fovplan

PIX = 25.0; LAM = 2.5; FIELD = 20  # pixels per field side
rng = np.random.default_rng(20261050)


def grf(shape, ell_um, rng):
    z = rng.standard_normal(shape)
    if ell_um > 0:
        z = gaussian_filter(z, ell_um / PIX, mode='wrap')
    return (z - z.mean()) / z.std()


def section(sigma, ell, p0, nx=12, ny=8, rng=rng, shift=0.0, K=6):
    """Return field counts array (fields x 2: focal, total) and the true focal proportion."""
    H, W = ny * FIELD, nx * FIELD
    b = np.r_[p0, np.full(K - 1, (1 - p0) / (K - 1))]
    eta = np.stack([np.log(b[t]) + sigma * grf((H, W), ell, rng) for t in range(K)])
    eta[0] += shift
    eta -= eta.max(0)
    pi = np.exp(eta); pi /= pi.sum(0)
    tot = rng.poisson(LAM, (H, W))
    foc = rng.binomial(tot, pi[0])
    f = foc.reshape(ny, FIELD, nx, FIELD).sum((1, 3)).ravel()
    n = tot.reshape(ny, FIELD, nx, FIELD).sum((1, 3)).ravel()
    return f.astype(float), n.astype(float)


def quasibin(y, n, M):
    k = len(y); Y, N = y.sum(), n.sum()
    Yc, Nc = (Y + .5, N + 1) if (Y == 0 or Y == N) else (Y, N)
    p = Yc / Nc
    phi = np.sum((y - n * p) ** 2 / (n * p * (1 - p))) / (k - 1)
    se = np.sqrt((1 - k / M) * phi / (Nc * p * (1 - p)))
    q = stats.t.ppf(.975, k - 1); lg = np.log(p / (1 - p))
    return 1 / (1 + np.exp(-(lg - q * se))), 1 / (1 + np.exp(-(lg + q * se)))


def ratio_t(y, n, M):
    k = len(y); p = y.sum() / n.sum(); e = y - p * n
    v = (1 - k / M) * e.var(ddof=1) / (k * n.mean() ** 2); h = stats.t.ppf(.975, k - 1) * np.sqrt(v)
    return max(p - h, 0), min(p + h, 1)


S1, S2, S3 = [], [], []
grid = list(itertools.product([0, 0.5, 1.0, 1.5], [50, 200, 800], [0.02, 0.1, 0.3]))
for si, (sigma, ell, p0) in enumerate(grid):
    for rep in range(20):
        y, n = section(sigma, ell, p0)
        M = len(n); P = y.sum() / n.sum(); N = n.sum(); nbar = n.mean()
        Se2 = (y - P * n).var(ddof=1); D_true = Se2 / (nbar * P * (1 - P))
        # empirical design effect for k = 3 against the true proportion
        subs = np.array([rng.choice(M, 3, replace=False) for _ in range(2000)])
        ph = y[subs].sum(1) / n[subs].sum(1); ns = n[subs].sum(1)
        vc = np.mean((1 - ns / N) * P * (1 - P) * N / (N - 1) / ns)
        D_emp3 = np.mean((ph - P) ** 2) / vc
        # pilot estimates and planning
        for pr in range(20):
            perm = rng.permutation(M); pil, pool = perm[:5], perm[5:]
            yp, npl = y[pil], n[pil]; pp = yp.sum() / npl.sum()
            if not 0 < pp < 1:
                S2.append(dict(sigma=sigma, ell=ell, p0=p0, rep=rep, pilot=pr, pilot_zero=True)); continue
            Sp = (yp - pp * npl).var(ddof=1); D_pil = Sp / (npl.mean() * pp * (1 - pp))
            row = dict(sigma=sigma, ell=ell, p0=p0, rep=rep, pilot=pr, pilot_zero=False, P=P, D_true=D_true, D_pilot=D_pil)
            h = 0.25 * P
            plans = {'cells': fovplan.naive_required_fields(pp, npl.mean(), 0.25, relative=True),
                     'fovplan': fovplan.required_fields(pp, Sp, npl.mean(), 0.25, relative=True, M=len(pool)),
                     'fovplan80': fovplan.required_fields(pp, fovplan.se2_upper(Sp, 5, 0.8), npl.mean(), 0.25, relative=True, M=len(pool))}
            for name, k in plans.items():
                if not np.isfinite(k) or k > len(pool):
                    row[name + '_k'] = np.inf; row[name + '_assur'] = np.nan; continue
                k = int(k)
                ss = np.array([rng.choice(len(pool), k, replace=False) for _ in range(200)])
                est = y[pool][ss].sum(1) / n[pool][ss].sum(1)
                row[name + '_k'] = k; row[name + '_assur'] = np.mean(np.abs(est - P) <= h)
            S2.append(row)
        S1.append(dict(sigma=sigma, ell=ell, p0=p0, rep=rep, P=P, nbar=nbar, D_true=D_true, D_emp3=D_emp3,
                       rho=(D_true - 1) / (nbar - 1)))
        # intervals
        for k in [3, 5, 10]:
            hits = {'Wilson': 0, 'ratio_t': 0, 'Fieller': 0, 'quasibinomial': 0}; wid = {m: 0.0 for m in hits}; full = 0
            for d in range(100):
                ix = rng.choice(M, k, replace=False); yy, nn = y[ix], n[ix]
                ints = {'Wilson': [wilson(yy.sum(), nn.sum())], 'ratio_t': [ratio_t(yy, nn, M)],
                        'Fieller': polynomial_set(fieller_coefficients(yy, nn, M)), 'quasibinomial': [quasibin(yy, nn, M)]}
                for m, segs in ints.items():
                    hits[m] += any(l - 1e-12 <= P <= u + 1e-12 for l, u in segs); wid[m] += sum(u - l for l, u in segs)
                full += sum(u - l for l, u in ints['Fieller']) >= 1 - 1e-12
            for m in hits:
                S3.append(dict(sigma=sigma, ell=ell, p0=p0, rep=rep, k=k, method=m, coverage=hits[m] / 100, width=wid[m] / 100,
                               full=(full / 100 if m == 'Fieller' else 0.0)))
    print(si + 1, '/', len(grid), sigma, ell, p0, flush=True)

S1 = pd.DataFrame(S1); S2 = pd.DataFrame(S2); S3 = pd.DataFrame(S3)
S1.to_csv('../results/B1_sim_design_effect.tsv', sep='\t', index=False)
S2.to_csv('../results/B1_sim_planning.tsv.gz', sep='\t', index=False, compression='gzip')
S3.to_csv('../results/B1_sim_intervals.tsv', sep='\t', index=False)
