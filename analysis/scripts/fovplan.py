#!/usr/bin/env python3
"""fovplan: plan the number of fields of view (FOVs) needed to estimate cell-type proportions
from spatial data, using a pilot set of fields.

Model
-----
For k fields drawn at random from a tissue, the ratio estimator p_hat = sum(y)/sum(n) has
approximate variance
    Var(p_hat) = (1 - k/M) * D * p(1-p) / (k * nbar)
where nbar is the mean number of cells per field, M the number of fields available in the
section (omit for an effectively unlimited section) and D = S_e^2 / (nbar p(1-p)) the design
effect, with e_j = y_j - p n_j. D is estimated from pilot fields. A conservative plan replaces
S_e^2 by an upper confidence bound  (m-1) S_e^2 / chi2_{1-q, m-1}  (m pilot fields).

The required number of fields is the smallest k with  t_{k-1, 1-alpha/2} * sqrt(Var) <= h,
where h is the target half-width (absolute, or relative to p when --relative is used).

Input: a CSV/TSV with one row per pilot field and one column per cell type (counts).
Non-numeric columns are ignored.

Example
-------
    python fovplan.py pilot_counts.tsv --target 0.25 --relative --assurance 0.9   # 0.9 recommended for compartmentalized tissue

Dependencies: numpy, pandas, scipy.
"""
import argparse, sys
import numpy as np, pandas as pd
from scipy import stats

__version__ = '0.1.0'


def estimate(counts):
    """counts: (fields x types) array-like or DataFrame. Returns DataFrame of pilot estimates."""
    df = pd.DataFrame(counts)
    X = df.to_numpy(float)
    n = X.sum(1); m = len(n); N = n.sum(); nbar = n.mean()
    if m < 2:
        raise ValueError('at least two pilot fields are required')
    out = []
    for i, c in enumerate(df.columns):
        y = X[:, i]; p = y.sum() / N
        e = y - p * n; Se2 = e.var(ddof=1) if m > 1 else np.nan
        D = Se2 / (nbar * p * (1 - p)) if 0 < p < 1 else np.nan
        out.append(dict(cell_type=c, p=p, pilot_fields=m, nbar=nbar, Se2=Se2, D=D,
                        rho=(D - 1) / (nbar - 1) if np.isfinite(D) else np.nan))
    return pd.DataFrame(out)


def se2_upper(Se2, m, assurance):
    """Upper `assurance` confidence bound for S_e^2 from m pilot fields (chi-square)."""
    return (m - 1) * Se2 / stats.chi2.ppf(1 - assurance, m - 1)


def half_width(k, p, Se2, nbar, M=None, conf=0.95):
    """Expected half-width of a region-level t interval with k fields."""
    k = np.asarray(k, float)
    fpc = 1.0 if M is None else np.clip(1 - k / M, 0, 1)
    tq = stats.t.ppf(1 - (1 - conf) / 2, np.maximum(k - 1, 1))
    return tq * np.sqrt(fpc * Se2 / (k * nbar ** 2))


def required_fields(p, Se2, nbar, target, relative=False, M=None, conf=0.95, kmax=100000):
    """Smallest k >= 2 with expected half-width <= target (absolute or relative to p)."""
    if not (np.isfinite(Se2) and 0 < p < 1):
        return np.nan
    h = target * p if relative else target
    ks = np.arange(2, kmax + 1)
    if M is not None:
        ks = ks[ks <= M]
    hw = half_width(ks, p, Se2, nbar, M, conf)
    ok = np.nonzero(hw <= h)[0]
    return int(ks[ok[0]]) if len(ok) else np.inf


def naive_required_fields(p, nbar, target, relative=False, conf=0.95):
    """Fields implied by treating cells as independent (design effect = 1)."""
    h = target * p if relative else target
    z = stats.norm.ppf(1 - (1 - conf) / 2)
    cells = z * z * p * (1 - p) / h ** 2
    return max(1, int(np.ceil(cells / nbar)))


def plan(counts, target, relative=False, M=None, conf=0.95, assurance=None):
    est = estimate(counts)
    rec, rec_c, naive = [], [], []
    for _, r in est.iterrows():
        rec.append(required_fields(r.p, r.Se2, r.nbar, target, relative, M, conf))
        naive.append(naive_required_fields(r.p, r.nbar, target, relative, conf) if 0 < r.p < 1 else np.nan)
        if assurance:
            rec_c.append(required_fields(r.p, se2_upper(r.Se2, r.pilot_fields, assurance), r.nbar, target, relative, M, conf))
    est['fields_required'] = rec
    if assurance:
        est[f'fields_required_assurance{int(assurance * 100)}'] = rec_c
    est['fields_if_cells_independent'] = naive
    return est


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('counts', help='CSV/TSV, one row per pilot field, one column per cell type')
    ap.add_argument('--target', type=float, required=True, help='target half-width (proportion, or fraction of p with --relative)')
    ap.add_argument('--relative', action='store_true')
    ap.add_argument('--fields-available', type=int, default=None, help='M, number of fields in the section (optional)')
    ap.add_argument('--conf', type=float, default=0.95)
    ap.add_argument('--assurance', type=float, default=None, help='e.g. 0.8 for a conservative plan')
    ap.add_argument('--out', default=None)
    a = ap.parse_args(argv)
    sep = '\t' if a.counts.endswith(('.tsv', '.txt')) else ','
    df = pd.read_csv(a.counts, sep=sep)
    df = df.select_dtypes('number')
    res = plan(df, a.target, a.relative, a.fields_available, a.conf, a.assurance)
    if a.out:
        res.to_csv(a.out, sep='\t', index=False)
    with pd.option_context('display.width', 200, 'display.max_columns', 20):
        print(res.round(4).to_string(index=False))


if __name__ == '__main__':
    main()
