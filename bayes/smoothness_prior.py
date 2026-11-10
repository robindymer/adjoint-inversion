"""
Gaussian smoothness prior on the increments of a(x):

    a_1 ~ N(A_MEAN, GAMMA0^2),   (a_{i+1} - a_i) ~ N(0, GAMMA^2)          (TRANSFORM = "differences")
    a_1, a_2 ~ N(A_MEAN, GAMMA0^2),   (a_{i+2} - 2 a_{i+1} + a_i) ~ N(0, GAMMA^2)   (TRANSFORM = "second-differences")

i.e. p(a) ~ exp(-(a - A_MEAN)^T P (a - A_MEAN) / 2) with the precision P = L^T L / GAMMA^2 + (level terms) / GAMMA0^2
from the course code. GAMMA is per node spacing (0.8 km); the level a_1 needs its own prior because L does not see
constant shifts. Under "differences" the draws are random walks started at a_1, so their spread grows along the fault.

    log_prior(a), grad_log_prior(a)     up to a constant; for MAP, MALA and HMC
    from_whitened(w)                    a = A_MEAN + R^{-T} w with P = R R^T; under the prior w ~ N(0, I) (pCN works in w)
    sample_prior(size, rng)             draws of shape (size, M_P)

    python smoothness_prior.py          plots draws for three values of GAMMA against the true a(x)
"""
import numpy as np
import matplotlib.pyplot as plt

from fault import x, M_P, B, a_true, plot_draws, summarize
from priors import gaussian_precision

TRANSFORM = "differences"
A_MEAN = B                  # mean of the level a_1 (and so of every a_i)
GAMMA = 1e-3                # standard deviation of one increment
GAMMA0 = 0.002              # standard deviation of the level a_1


def factors(gamma=GAMMA, gamma0=GAMMA0, transform=TRANSFORM):
    """Precision P and its Cholesky factor R (P = R R^T)."""
    P = gaussian_precision(M_P, gamma, transform, gamma0)
    return P, np.linalg.cholesky(P)


P, R = factors()


def log_prior(a):
    d = a - A_MEAN
    return -0.5 * d @ P @ d


def grad_log_prior(a):
    return -P @ (a - A_MEAN)


def from_whitened(w):
    """A_MEAN + R^{-T} w for w of shape (M_P,) or (size, M_P); then Cov(a) = R^{-T} R^{-1} = P^{-1}."""
    return A_MEAN + np.linalg.solve(R.T, np.transpose(w)).T


def sample_prior(size=1, rng=None):
    rng = np.random.default_rng() if rng is None else rng
    return from_whitened(rng.standard_normal((size, M_P)))


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    gammas = [2.5e-4, 1e-3, 2.5e-3]
    fig, axs = plt.subplots(1, len(gammas), figsize=(14, 3.5), sharey=True, constrained_layout=True)
    for ax, gamma in zip(axs, gammas):
        P, R = factors(gamma)            # the functions above read the module-level P, R
        draws = sample_prior(10000, rng)
        plot_draws(ax, draws, f"increments $\\sim \\mathcal{{N}}(0, \\gamma^2)$, $\\gamma = {gamma}$, "
                              f"$\\gamma_0 = {GAMMA0}$")
        std = np.sqrt(np.diag(np.linalg.inv(P)))
        ax.fill_between(x, A_MEAN - 2 * std, A_MEAN + 2 * std, color="0.85", zorder=0)   # +- 2 prior std
        summarize(f"gamma = {gamma}", draws)
        print(f"gamma = {gamma}: squared Mahalanobis distance of the truth {-2 * log_prior(a_true(x)):.1f} (draws: about {M_P})")
    plt.show()
