"""Shared synthetic-tissue generator for G21 simulations (see b1_simulation.py)."""
import numpy as np
from scipy.ndimage import gaussian_filter

PIX = 25.0; LAM = 2.5; FIELD = 20  # pixels per field side


def grf(shape, ell_um, rng):
    z = rng.standard_normal(shape)
    if ell_um > 0:
        z = gaussian_filter(z, ell_um / PIX, mode='wrap')
    return (z - z.mean()) / z.std()


def section(sigma, ell, p0, nx=12, ny=8, rng=None, shift=0.0, K=6):
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


