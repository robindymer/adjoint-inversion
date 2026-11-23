"""
Diagnostics of several chains of pcn.m (run_pcn.m, run_hpcn.m), as in Sec. 4.6 and Remark 4.21 of the script:

    rank-normalized split-Rhat and ESS (samplers.rank_rhat, samplers.ess) of U and of every node of a,
    forward solves per effective sample, and the posterior mean and sd against the Laplace approximation.

    python diagnostics.py <glob of chain files> [burn-in]
    e.g. python diagnostics.py "../results/hpcn_m36_eta0.02_s1_laplace_seed*.mat"

The burn-in defaults to max(nAdapt, 10 % of the chain). The figure is saved next to the chains, in figures/.
"""
import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

import fault  # also puts the course code on the path
from samplers import ess, rank_rhat


def load_chains(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"No files match {pattern}")
    runs = [loadmat(f, squeeze_me=True, struct_as_record=False) for f in files]
    n = min(r["chain"].a.shape[1] for r in runs)   # chains of unequal length are cut to the shortest
    A = np.stack([r["chain"].a[:, :n].T for r in runs])   # (chains, iterations, m_p)
    U = np.stack([r["chain"].U[:n] for r in runs])        # (chains, iterations)
    acc = np.stack([r["chain"].accepted[:n] for r in runs])
    n_adapt = int(getattr(runs[0]["chain"].opts, "nAdapt", 0))
    return files, runs, A, U, acc, n_adapt


def main(pattern, burn_in=None):
    files, runs, A, U, acc, n_adapt = load_chains(pattern)
    C, N, m_p = A.shape
    if burn_in is None:
        burn_in = max(n_adapt, N // 10)
    A, U, acc = A[:, burn_in:], U[:, burn_in:], acc[:, burn_in:]
    a_true = runs[0]["aTrue"]
    print(f"{C} chains of {N} iterations, burn-in {burn_in}, kept {C * (N - burn_in)} samples")
    print(f"acceptance rate per chain: {np.round(acc.mean(axis=1), 2)}")

    # Rhat and ESS of U and of every node
    rhat = np.array([rank_rhat(A[:, :, i]) for i in range(m_p)])
    ess_a = np.array([ess(A[:, :, i]) for i in range(m_p)])
    rhat_U, ess_U = rank_rhat(U), ess(U)
    n_solves = C * N   # one forward solve per iteration, burn-in included
    print(f"U:     Rhat {rhat_U:.3f}, ESS {ess_U:.0f}")
    print(f"nodes: Rhat max {rhat.max():.3f} (node {rhat.argmax() + 1}), "
          f"ESS min {ess_a.min():.0f} (node {ess_a.argmin() + 1}), median {np.median(ess_a):.0f}")
    print(f"forward solves per effective sample (worst node): {n_solves / ess_a.min():.1f}")
    ok = rhat.max() < 1.01 and rhat_U < 1.01 and ess_a.min() > 400 and ess_U > 400
    print(f"Rhat < 1.01 and ESS > 400 for U and all nodes: {'yes' if ok else 'NO'}")

    # Posterior summaries, and the Laplace approximation if the chains were run with it
    samples = A.reshape(-1, m_p)
    mean, sd = samples.mean(axis=0), samples.std(axis=0, ddof=1)
    lo, hi = np.percentile(samples, [2.5, 97.5], axis=0)
    mcse = sd / np.sqrt(ess_a)   # Monte Carlo error of the posterior mean
    print(f"|mean - a_true|/sd: max {np.max(np.abs(mean - a_true) / sd):.2f}; "
          f"a_true inside the 95 % band at {np.sum((a_true >= lo) & (a_true <= hi))} of {m_p} nodes")
    print(f"Monte Carlo error of the mean / sd: max {np.max(mcse / sd):.3f}")
    laplace = None
    if "H" in runs[0]:
        H, a_map = runs[0]["H"], runs[0]["aMap"]
        sd_lap = np.sqrt(np.diag(np.linalg.inv(H)))
        laplace = (a_map, sd_lap)
        print(f"vs Laplace: |mean - a_MAP|/sd max {np.max(np.abs(mean - a_map) / sd):.2f}, "
              f"sd/sd_Laplace in [{np.min(sd / sd_lap):.3f}, {np.max(sd / sd_lap):.3f}]")

    # Figure: traces of U and of the slowest node, the posterior band, and the posterior sd
    x = fault.x
    worst = ess_a.argmin()
    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    it = np.arange(burn_in, N)
    for c in range(C):
        ax[0, 0].plot(it, U[c], lw=0.6)
        ax[0, 1].plot(it, A[c, :, worst], lw=0.6)
    ax[0, 0].set(xlabel="iteration", ylabel="$U$", title=f"Traces of $U$ ($\\hat R$ = {rhat_U:.3f})")
    ax[0, 1].set(xlabel="iteration", ylabel=f"$a$ at node {worst + 1}",
                 title=f"Slowest node (ESS = {ess_a[worst]:.0f})")
    ax[1, 0].fill_between(x, lo - a_true, hi - a_true, alpha=0.3, label="95 % band")
    ax[1, 0].plot(x, mean - a_true, ".-", label="mean")
    ax[1, 0].axhline(0, color="k", lw=1, label="$a$, true")
    ax[1, 0].set(xlabel="$x$ (km)", ylabel="$a - a_{true}$", title="Posterior relative to the truth")
    ax[1, 0].legend(frameon=False, fontsize=8)
    ax[1, 1].semilogy(x, sd, ".-", label="samples")
    if laplace is not None:
        ax[1, 1].semilogy(x, laplace[1], "--", label="Laplace")
    ax[1, 1].set(xlabel="$x$ (km)", ylabel="posterior sd of $a$", title="Posterior standard deviation")
    ax[1, 1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    out_dir = os.path.join(os.path.dirname(files[0]), "figures")
    os.makedirs(out_dir, exist_ok=True)
    name = os.path.basename(files[0]).replace("_seed1.mat", "").replace(".mat", "")
    out = os.path.join(out_dir, f"diagnostics_{name}.png")
    fig.savefig(out, dpi=150)
    print(f"Saved {out}")

    # Figure: posterior of a with 50 % and 95 % credible bands, on top of the true a(x)
    lo50, hi50 = np.percentile(samples, [25, 75], axis=0)
    xx = np.linspace(fault.X_MIN, fault.X_MAX, 1001)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(xx, fault.a_true(xx), color="0.75", lw=4, label="$a$, true", zorder=0)
    ax.fill_between(x, lo, hi, color="C0", alpha=0.2, lw=0, label="95 % band")
    ax.fill_between(x, lo50, hi50, color="C0", alpha=0.4, lw=0, label="50 % band")
    ax.plot(x, mean, ".-", color="C0", lw=1.2, label="posterior mean")
    ax.axhline(fault.B, color="0.4", ls=":", label="$b$ (VW below)")
    ax.set(xlabel="$x$ (km)", ylabel="$a$", title=f"Posterior of $a$ ({C} chains, {C * (N - burn_in)} samples)")
    ax.legend(frameon=False, fontsize=8, loc="upper center", ncol=3)
    fig.tight_layout()
    out = os.path.join(out_dir, f"posterior_a_{name}.png")
    fig.savefig(out, dpi=150)
    print(f"Saved {out}")

    plt.show()


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else None)
