"""
Gaussian prior on the nodal values of a(x) (no structure on the increments):

    a ~ N(A_MEAN, C),   C = GAMMA^2 I                                      (KERNEL = "iid")
                        C_ij = GAMMA^2 exp(-(x_i - x_j)^2 / (2 ELL^2))     (KERNEL = "se", Gaussian process prior)

    log_prior(a), grad_log_prior(a)     up to a constant; for MAP, MALA and HMC
    from_whitened(w)                    a = A_MEAN + R w with C = R R^T; under the prior w ~ N(0, I) (pCN works in w)
    sample_prior(size, rng)             draws of shape (size, M_P)

    python gaussian_prior.py            plots draws under both kernels against the true a(x)
"""
import numpy as np
import matplotlib.pyplot as plt

from bayes.python.fault import x, M_P, B, a_true, plot_draws, summarize
from priors import se_covariance

KERNEL = "iid"
A_MEAN = B                  # prior mean on the VW/VS boundary: neither regime is favoured
GAMMA = 0.002               # pointwise standard deviation
ELL = 2.0                   # correlation length (km), KERNEL = "se" only


def covariance(kernel=KERNEL, gamma=GAMMA, ell=ELL):
    if kernel == "iid":
        return gamma**2 * np.eye(M_P)
    if kernel == "se":
        return se_covariance(x, gamma, ell)
    raise ValueError("kernel must be 'iid' or 'se'")


def factors(C):
    """Cholesky factor R of C and the precision C^{-1} (a small ridge keeps the SE covariance factorizable)."""
    R = np.linalg.cholesky(C + 1e-10 * np.trace(C) / len(C) * np.eye(len(C)))
    R_inv = np.linalg.inv(R)
    return R, R_inv.T @ R_inv


R, P = factors(covariance())


def log_prior(a):
    d = a - A_MEAN
    return -0.5 * d @ P @ d


def grad_log_prior(a):
    return -P @ (a - A_MEAN)


def from_whitened(w):
    return A_MEAN + w @ R.T


def sample_prior(size=1, rng=None):
    rng = np.random.default_rng() if rng is None else rng
    return from_whitened(rng.standard_normal((size, M_P)))


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.5), sharey=True, constrained_layout=True)
    for ax, kernel in zip(axs, ["iid", "se"]):
        R, P = factors(covariance(kernel))
        draws = sample_prior(10000, rng)
        title = f"$a_i \\sim \\mathcal{{N}}({A_MEAN}, \\gamma^2)$, $\\gamma = {GAMMA}$" if kernel == "iid" else \
            f"GP, squared exponential, $\\gamma = {GAMMA}$, $\\ell = {ELL}$ km"
        plot_draws(ax, draws, title)
        summarize(kernel, draws)
        print(f"{kernel}: squared Mahalanobis distance of the truth {-2 * log_prior(a_true(x)):.1f} (draws: about {M_P})")
    plt.show()
