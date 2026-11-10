"""
Hierarchical (GSBL) prior on the increments z = L a of a(x), with an inverse-gamma hyper-prior:

    a_1 ~ N(A_MEAN, GAMMA0^2),   z_i | theta_i ~ N(0, theta_i),   theta_i ~ InvGamma(BETA, VARTHETA),  i = 1..M_P-1.

Most theta_i stay near the mode VARTHETA / (BETA + 1), so most increments are tiny, but the heavy tail lets a few
theta_i be large: the draws are close to piecewise constant, like the true a(x). Marginally z_i is Student-t with
2 BETA degrees of freedom and scale sqrt(VARTHETA / BETA). With BETA = 1/2, VARTHETA = 1e-8 a draw has on average
about 1.5 increments larger than 0.0015 and a median |z_i| of about 1.4e-4 (the true a(x) has two jumps of 0.004).

    log_prior(a, theta)                 joint log density up to a constant
    grad_log_prior(a, theta)            gradient in a for fixed theta (Gibbs-type and IAS steps)
    ias_theta_update(a)                 closed-form theta-step of the iterative alternating scheme for the joint MAP
    to_normalized(a, theta), from_normalized(w, tau)
                                        prior-normalized coordinates: under the prior (w, tau) ~ N(0, I) with
                                        w_0 = (a_1 - A_MEAN) / GAMMA0, w_i = z_i / sqrt(theta_i), tau_i = Phi^{-1}(F(theta_i)),
                                        F the InvGamma cdf; pCN then targets likelihood x standard normal
    sample_prior(size, rng, return_theta)   draws of shape (size, M_P) (and theta of shape (size, M_P - 1))

    python hierarchical_prior.py        plots draws of a and their theta against the true a(x)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import invgamma, norm

from bayes.python.fault import x, M_P, B, a_true, plot_draws, summarize
from priors import sample_gsbl, neg_log_invgamma
from solvers import difference_matrix

A_MEAN = B                  # mean of the level a_1
GAMMA0 = 0.002              # standard deviation of the level a_1
BETA = 0.5                  # hyper-prior shape: smaller is heavier-tailed (fewer, larger jumps)
VARTHETA = 1e-8             # hyper-prior scale: sets the size of the small increments

L = difference_matrix(M_P)
hyper = invgamma(BETA, scale=VARTHETA)


def log_prior(a, theta):
    if np.any(theta <= 0):
        return -np.inf
    z = L @ a
    return (-0.5 * (a[0] - A_MEAN) ** 2 / GAMMA0**2 + np.sum(-0.5 * np.log(theta) - 0.5 * z**2 / theta)
            - np.sum(neg_log_invgamma(theta, BETA, VARTHETA)))


def grad_log_prior(a, theta):
    g = -L.T @ (L @ a / theta)
    g[0] -= (a[0] - A_MEAN) / GAMMA0**2
    return g


def ias_theta_update(a):
    """argmax over theta of log_prior(a, theta): theta_i = (VARTHETA + z_i^2 / 2) / (BETA + 3/2)."""
    return (VARTHETA + 0.5 * (L @ a) ** 2) / (BETA + 1.5)


def theta_from_tau(tau):
    """theta = F^{-1}(Phi(tau)), accurate in both tails (upper tail via survival functions)."""
    tau = np.asarray(tau, dtype=float)
    return np.where(tau < 0, hyper.ppf(norm.cdf(tau)), hyper.isf(norm.sf(tau)))


def tau_from_theta(theta):
    c = hyper.cdf(theta)
    return np.where(c < 0.5, norm.ppf(c), norm.isf(hyper.sf(theta)))


def to_normalized(a, theta):
    w = np.concatenate([[(a[0] - A_MEAN) / GAMMA0], L @ a / np.sqrt(theta)])
    return w, tau_from_theta(theta)


def from_normalized(w, tau):
    theta = theta_from_tau(tau)
    a = A_MEAN + GAMMA0 * w[0] + np.concatenate([[0.0], np.cumsum(np.sqrt(theta) * w[1:])])
    return a, theta


def sample_prior(size=1, rng=None, return_theta=False):
    rng = np.random.default_rng() if rng is None else rng
    a, theta = sample_gsbl(M_P, BETA, VARTHETA, "differences", size, rng, return_theta=True)
    a += A_MEAN + GAMMA0 * rng.standard_normal((size, 1))
    return (a, theta) if return_theta else a


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    draws, theta = sample_prior(10000, rng, return_theta=True)
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.5), constrained_layout=True)
    plot_draws(axs[0], draws, f"GSBL on increments, $\\beta = {BETA}$, $\\vartheta = {VARTHETA:g}$")
    xm = 0.5 * (x[:-1] + x[1:])
    for k in range(5):
        axs[1].semilogy(xm, theta[k], ".-", lw=0.8, ms=3)
    axs[1].axhline(VARTHETA / (BETA + 1), color="0.4", ls=":", label="hyper-prior mode")
    axs[1].set_xlabel("$x$ (km)"); axs[1].set_title("$\\theta_i$ of the same draws")
    axs[1].legend(frameon=False, fontsize=7)
    summarize("hierarchical", draws)
    jumps = np.sum(np.abs(np.diff(draws, axis=1)) > 1.5e-3, axis=1)
    print(f"hierarchical: mean number of jumps > 0.0015 per draw {jumps.mean():.2f}, P(no jump) {np.mean(jumps == 0):.2f}")
    plt.show()
