import numpy as np, pandas as pd
import fovplan


def test_design_effect_near_one_for_unclustered_counts():
    rng = np.random.default_rng(1)
    n = rng.poisson(1000, 200)
    y = rng.binomial(n, 0.1)
    est = fovplan.estimate(pd.DataFrame({'a': y, 'b': n - y}))
    assert 0.7 < est.loc[0, 'D'] < 1.3


def test_required_fields_increase_with_clustering_and_precision():
    p, nbar = 0.1, 1000
    k1 = fovplan.required_fields(p, 1.0 * nbar * p * (1 - p), nbar, 0.25, relative=True)
    k10 = fovplan.required_fields(p, 10.0 * nbar * p * (1 - p), nbar, 0.25, relative=True)
    k10_tight = fovplan.required_fields(p, 10.0 * nbar * p * (1 - p), nbar, 0.10, relative=True)
    assert k1 <= k10 < k10_tight


def test_assurance_plan_is_more_conservative():
    rng = np.random.default_rng(2)
    n = rng.poisson(1000, 6); y = rng.binomial(n, rng.beta(5, 45, 6))
    df = pd.DataFrame({'t': y, 'rest': n - y})
    pl = fovplan.plan(df, 0.25, relative=True, assurance=0.9)
    assert pl.loc[0, 'fields_required_assurance90'] >= pl.loc[0, 'fields_required']


def test_naive_plan_matches_binomial_formula():
    k = fovplan.naive_required_fields(0.1, 1000, 0.02)
    assert k == int(np.ceil(1.96 ** 2 * 0.09 / 0.02 ** 2 / 1000))


def test_bisection_search_matches_exhaustive_scan():
    from scipy import stats
    rng = np.random.default_rng(3)
    for _ in range(300):
        p = rng.uniform(0.005, 0.5); nbar = rng.uniform(100, 2000); D = np.exp(rng.uniform(0, 5))
        Se2 = D * nbar * p * (1 - p); M = None if rng.random() < 0.5 else int(rng.integers(3, 300))
        ks = np.arange(2, 5001) if M is None else np.arange(2, M + 1)
        hw = fovplan.half_width(ks, p, Se2, nbar, M)
        ok = np.nonzero(hw <= 0.25 * p)[0]
        expected = int(ks[ok[0]]) if len(ok) else np.inf
        got = fovplan.required_fields(p, Se2, nbar, 0.25, relative=True, M=M, kmax=5000)
        assert got == expected
