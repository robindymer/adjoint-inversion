"""Posterior mean and 95 % bands of a(x) against the true a, one panel per (MAP file, chain glob):
MCMC mean and 2.5-97.5 percentile band (shaded) and the Laplace band at the MAP (dashed).

    python plot_posterior_bands.py OUT.png MAPFILE CHAINGLOB [MAPFILE CHAINGLOB ...]
"""
import sys
from glob import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.stats import norm


def panel(ax, map_file, pattern):
    M = loadmat(map_file, squeeze_me=True, struct_as_record=False)
    a_map, H = M["aMap"], M["H"]
    m_p = len(a_map)
    x = np.linspace(-10, 10, m_p)
    std = np.sqrt(np.diag(np.linalg.inv(H)))
    q = norm.ppf(0.975)

    runs = [loadmat(f, squeeze_me=True, struct_as_record=False)["chain"] for f in sorted(glob(pattern))]
    burn = int(getattr(runs[0].opts, "nAdapt", 0))
    S = np.vstack([r.a[:, burn:].T for r in runs])
    lo, hi = np.percentile(S, [2.5, 97.5], axis=0)

    xt = np.linspace(-10, 10, 1001)   # true a on the grid of the data, a step at x = -5 and 6
    at = np.where((xt >= -5) & (xt <= 6), 0.009, 0.013)
    ax.fill_between(x, lo, hi, color="0.8", label="MCMC 95 %")
    ax.plot(x, S.mean(axis=0), "k-", lw=1.5, label="MCMC mean")
    ax.plot(x, a_map - q*std, "r--", lw=1, label="Laplace 95 %")
    ax.plot(x, a_map + q*std, "r--", lw=1)
    ax.plot(xt, at, color="C0", lw=2, label="true")
    ax.set_title(f"m = {int(M['m'])}, m_p = {m_p}, {len(runs)} chains x {S.shape[0]//len(runs)} samples")
    ax.set_ylabel("a")
    ax.grid(alpha=0.3)


def main(out, args):
    pairs = list(zip(args[::2], args[1::2]))
    fig, axes = plt.subplots(len(pairs), 1, figsize=(9, 4*len(pairs)), sharex=True, squeeze=False)
    for ax, (map_file, pattern) in zip(axes[:, 0], pairs):
        panel(ax, map_file, pattern)
    axes[0, 0].legend(fontsize=9, loc="upper right")
    axes[-1, 0].set_xlabel("x (km)")
    fig.tight_layout()
    fig.savefig(out, dpi=200)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2:])
