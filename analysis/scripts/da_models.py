"""Additional differential-abundance comparators (exploratory, added after G25).

All tests are two-sided for a two-group comparison of one cell category.

glmm_ri      binomial GLMM on field counts, logit(p_ij) = b0 + b1*g_i + u_i, u_i ~ N(0, s_u^2);
             adaptive Gauss-Hermite quadrature (7 nodes; cf. lme4::glmer nAGQ = 7).
glmm_olre    as glmm_ri plus a field-level (observation-level) random effect v_ij ~ N(0, s_v^2);
             Laplace approximation (cf. glmer(... + (1|individual) + (1|field)), the default nAGQ = 1).
             Both GLMMs return a likelihood-ratio P value (chi-square, 1 df) and a Wald z P value.
betabin      beta-binomial regression on per-individual totals (cf. glmmTMB/aod betabin); LRT.
nb_glm       negative-binomial GLM on per-individual counts with log(total) offset and a
             category-specific dispersion (edgeR/DESeq2-style count model); LRT.
propeller    speckle::propeller re-implementation: per-individual proportions of all mutually exclusive
             clusters, logit (pseudo-count 0.5) or arcsine-square-root transform, linear model and
             limma empirical-Bayes moderated t (eBayes, robust = FALSE).
"""
import numpy as np
from scipy import optimize, special, stats

GH_X, GH_W = np.polynomial.hermite.hermgauss(7)   # physicists' nodes: integral exp(-x^2) f(x)


def _bincount(idx, w, n):
    return np.bincount(idx, weights=w, minlength=n)


# ---------------------------------------------------------------- binomial GLMM, random intercept
_WARM = {}


def _glmm_ri_negll(theta, y, n, idx, gi, G, fixed_b1=None):
    if fixed_b1 is None:
        b0, b1, ls = theta
    else:
        (b0, ls), b1 = theta, fixed_b1
    su = np.exp(ls)
    off = b0 + b1 * gi                                  # per individual
    # mode of h(u) = sum_j loglik(y | off+u) - u^2/(2 su^2) via Newton (vectorized over individuals),
    # warm-started from the previous evaluation
    key = ('ri', id(y))
    u = _WARM.get(key, np.zeros(G)).copy()
    if len(u) != G:
        u = np.zeros(G)
    for _ in range(100):
        p = special.expit(off[idx] + u[idx])
        g = _bincount(idx, y - n * p, G) - u / su ** 2
        h = _bincount(idx, n * p * (1 - p), G) + 1 / su ** 2
        step = np.clip(g / h, -2, 2)
        u += step
        if np.max(np.abs(step)) < 1e-10:
            break
    if np.all(np.isfinite(u)):
        _WARM[key] = u.copy()
    else:
        _WARM.pop(key, None)
    p = special.expit(off[idx] + u[idx])
    h = _bincount(idx, n * p * (1 - p), G) + 1 / su ** 2
    sd = 1 / np.sqrt(h)
    # adaptive GHQ: integral exp(h(u)) N(u;0,su) du, nodes u = mode + sqrt(2) sd x
    ll = np.zeros(G)
    terms = []
    for x, w in zip(GH_X, GH_W):
        uu = u + np.sqrt(2) * sd * x
        eta = off[idx] + uu[idx]
        lb = _bincount(idx, y * eta - n * np.logaddexp(0, eta), G)
        terms.append(lb - uu ** 2 / (2 * su ** 2) - np.log(su) - 0.5 * np.log(2 * np.pi)
                     + np.log(w) + x ** 2 + np.log(np.sqrt(2) * sd))
    ll = special.logsumexp(np.vstack(terms), axis=0)
    return -ll.sum()


def glmm_ri(yf, nf, ind, grp):
    """yf, nf: field counts; ind: individual label per field; grp: 0/1 per field."""
    _WARM.clear()
    yf = np.array(yf, float)
    labs, idx = np.unique(ind, return_inverse=True)
    G = len(labs)
    gi = np.zeros(G); gi[idx] = grp
    p0 = (yf.sum() + .5) / (nf.sum() + 1)
    start = np.array([np.log(p0 / (1 - p0)), 0.0, np.log(0.5)])
    bnd = [(-30, 30), (-20, 20), (-8, 3)]
    f1 = optimize.minimize(_glmm_ri_negll, start, args=(yf, nf, idx, gi, G), method='L-BFGS-B', bounds=bnd)
    f0 = optimize.minimize(_glmm_ri_negll, start[[0, 2]], args=(yf, nf, idx, gi, G, 0.0), method='L-BFGS-B',
                           bounds=[bnd[0], bnd[2]])
    lr = max(2 * (f0.fun - f1.fun), 0.0)
    p_lrt = stats.chi2.sf(lr, 1)
    # Wald: numerical Hessian of the full model at the optimum
    H = _num_hess(lambda t: _glmm_ri_negll(t, yf, nf, idx, gi, G), f1.x)
    try:
        se = np.sqrt(np.linalg.inv(H)[1, 1])
        p_wald = 2 * stats.norm.sf(abs(f1.x[1] / se))
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        p_wald = np.nan
    return p_lrt, p_wald if np.isfinite(p_wald) else np.nan, f1.x


def _num_hess(f, x, eps=1e-4):
    k = len(x); H = np.zeros((k, k))
    for i in range(k):
        for j in range(i, k):
            e_i = np.zeros(k); e_i[i] = eps; e_j = np.zeros(k); e_j[j] = eps
            H[i, j] = H[j, i] = (f(x + e_i + e_j) - f(x + e_i - e_j) - f(x - e_i + e_j) + f(x - e_i - e_j)) / (4 * eps ** 2)
    return H


# ---------------------------------------------- binomial GLMM, individual + field random effects (Laplace)
def _glmm_olre_negll(theta, y, n, idx, gi, G, fixed_b1=None):
    if fixed_b1 is None:
        b0, b1, lu, lv = theta
    else:
        (b0, lu, lv), b1 = theta, fixed_b1
    su2, sv2 = np.exp(2 * lu), np.exp(2 * lv)
    off = b0 + b1 * gi
    key = ('olre', id(y))
    u, v = _WARM.get(key, (np.zeros(G), np.zeros(len(y))))
    u, v = u.copy(), v.copy()
    if len(u) != G or len(v) != len(y):
        u, v = np.zeros(G), np.zeros(len(y))
    for _ in range(200):
        eta = off[idx] + u[idx] + v
        p = special.expit(eta); w = n * p * (1 - p)
        gv = (y - n * p) - v / sv2
        gu = _bincount(idx, y - n * p, G) - u / su2
        hvv = w + 1 / sv2
        huu = _bincount(idx, w, G) + 1 / su2
        S = huu - _bincount(idx, w ** 2 / hvv, G)          # Schur complement (arrowhead Hessian)
        du = (gu - _bincount(idx, w * gv / hvv, G)) / S
        dv = (gv - w * du[idx]) / hvv
        du = np.clip(du, -2, 2); dv = np.clip(dv, -2, 2)
        u += du; v += dv
        if max(np.max(np.abs(du)), np.max(np.abs(dv))) < 1e-9:
            break
    if np.all(np.isfinite(u)) and np.all(np.isfinite(v)):
        _WARM[key] = (u.copy(), v.copy())
    else:
        _WARM.pop(key, None)
    eta = off[idx] + u[idx] + v
    p = special.expit(eta); w = n * p * (1 - p)
    hvv = w + 1 / sv2
    huu = _bincount(idx, w, G) + 1 / su2
    S = huu - _bincount(idx, w ** 2 / hvv, G)
    logdet = np.log(S) + _bincount(idx, np.log(hvv), G)
    J = np.bincount(idx, minlength=G)
    lb = _bincount(idx, y * eta - n * np.logaddexp(0, eta) - v ** 2 / (2 * sv2), G)
    ll = lb - u ** 2 / (2 * su2) - 0.5 * np.log(su2) - 0.5 * J * np.log(sv2) - 0.5 * logdet
    return -ll.sum()


def glmm_olre(yf, nf, ind, grp):
    _WARM.clear()
    yf = np.array(yf, float)
    labs, idx = np.unique(ind, return_inverse=True)
    G = len(labs)
    gi = np.zeros(G); gi[idx] = grp
    p0 = (yf.sum() + .5) / (nf.sum() + 1)
    start = np.array([np.log(p0 / (1 - p0)), 0.0, np.log(0.3), np.log(0.3)])
    bnd = [(-30, 30), (-20, 20), (-8, 3), (-8, 3)]
    f1 = optimize.minimize(_glmm_olre_negll, start, args=(yf, nf, idx, gi, G), method='L-BFGS-B', bounds=bnd)
    f0 = optimize.minimize(_glmm_olre_negll, start[[0, 2, 3]], args=(yf, nf, idx, gi, G, 0.0), method='L-BFGS-B',
                           bounds=[bnd[0], bnd[2], bnd[3]])
    lr = max(2 * (f0.fun - f1.fun), 0.0)
    p_lrt = stats.chi2.sf(lr, 1)
    H = _num_hess(lambda t: _glmm_olre_negll(t, yf, nf, idx, gi, G), f1.x)
    try:
        se = np.sqrt(np.linalg.inv(H)[1, 1])
        p_wald = 2 * stats.norm.sf(abs(f1.x[1] / se))
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        p_wald = np.nan
    return p_lrt, p_wald if np.isfinite(p_wald) else np.nan, f1.x


# ----------------------------------------------------------------- individual-level count models
def _bb_negll(theta, y, n, g, fixed_b1=None):
    if fixed_b1 is None:
        b0, b1, lphi = theta
    else:
        (b0, lphi), b1 = theta, fixed_b1
    mu = special.expit(b0 + b1 * g); phi = np.exp(lphi)
    a, b = mu * phi, (1 - mu) * phi
    return -np.sum(special.betaln(y + a, n - y + b) - special.betaln(a, b))


def betabin(yi, ni, gi):
    p0 = (yi.sum() + .5) / (ni.sum() + 1)
    start = np.array([np.log(p0 / (1 - p0)), 0.0, np.log(50.0)])
    bnd = [(-30, 30), (-20, 20), (-5, 20)]
    f1 = optimize.minimize(_bb_negll, start, args=(yi, ni, gi), method='L-BFGS-B', bounds=bnd)
    f0 = optimize.minimize(_bb_negll, start[[0, 2]], args=(yi, ni, gi, 0.0), method='L-BFGS-B', bounds=[bnd[0], bnd[2]])
    return stats.chi2.sf(max(2 * (f0.fun - f1.fun), 0), 1)


def _nb_negll(theta, y, off, g, fixed_b1=None):
    if fixed_b1 is None:
        b0, b1, lk = theta
    else:
        (b0, lk), b1 = theta, fixed_b1
    mu = np.exp(off + b0 + b1 * g); r = np.exp(-lk)          # dispersion k = 1/r, var = mu + k mu^2
    return -np.sum(special.gammaln(y + r) - special.gammaln(r) - special.gammaln(y + 1)
                   + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu)))


def nb_glm(yi, ni, gi):
    off = np.log(ni)
    b0 = np.log((yi.sum() + .5) / ni.sum())
    start = np.array([b0, 0.0, np.log(0.1)])
    bnd = [(-40, 10), (-20, 20), (-15, 5)]
    f1 = optimize.minimize(_nb_negll, start, args=(yi, off, gi), method='L-BFGS-B', bounds=bnd)
    f0 = optimize.minimize(_nb_negll, start[[0, 2]], args=(yi, off, gi, 0.0), method='L-BFGS-B', bounds=[bnd[0], bnd[2]])
    return stats.chi2.sf(max(2 * (f0.fun - f1.fun), 0), 1)


# ---------------------------------------------------------------------------------- propeller
def _trigamma_inverse(x):
    """limma::trigammaInverse (Newton iteration)."""
    x = np.asarray(x, float)
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.max(-dif / y) < 1e-8:
            break
    return y


def _fit_fdist(s2, df):
    """limma::fitFDist moment estimates (equal residual df, no covariate)."""
    s2 = np.maximum(s2, 1e-5 * np.median(s2[s2 > 0]) if np.any(s2 > 0) else 1e-12)
    z = np.log(s2)
    e = z - special.digamma(df / 2) + np.log(df / 2)
    emean = e.mean()
    evar = np.sum((e - emean) ** 2) / (len(e) - 1) - special.polygamma(1, df / 2)
    if evar > 0:
        d0 = 2 * _trigamma_inverse(evar)
        s0 = np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0, s0 = np.inf, np.exp(emean)
    return s0, d0


def propeller(counts, gi, transform='logit'):
    """counts: individuals x clusters (all mutually exclusive clusters); gi: 0/1.
    Returns moderated-t P values for every cluster (column)."""
    counts = np.asarray(counts, float)
    if transform == 'logit':
        c = counts + 0.5
        props = c / c.sum(1, keepdims=True)
        Y = np.log(props / (1 - props))
    else:
        props = counts / counts.sum(1, keepdims=True)
        Y = np.arcsin(np.sqrt(props))
    X = np.column_stack([np.ones(len(gi)), gi])
    XtXi = np.linalg.inv(X.T @ X)
    B = XtXi @ X.T @ Y
    R = Y - X @ B
    df = len(gi) - 2
    s2 = (R ** 2).sum(0) / df
    s0, d0 = _fit_fdist(s2, df)
    if np.isfinite(d0):
        s2p = (d0 * s0 + df * s2) / (d0 + df)
        dft = df + d0
    else:
        s2p, dft = np.full_like(s2, s0), np.inf
    t = B[1] / np.sqrt(s2p * XtXi[1, 1])
    return 2 * (stats.norm.sf(np.abs(t)) if not np.isfinite(dft) else stats.t.sf(np.abs(t), dft))
