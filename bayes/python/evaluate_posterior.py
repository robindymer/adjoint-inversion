"""
Calibration of a high-res posterior against the true a(x): is a_true inside the 95 % band?

    python evaluate_posterior.py MAPFILE [CHAINGLOB ...]

MAPFILE    results/highres_map_*.mat: the Laplace approximation N(a_MAP, H^-1)
CHAINGLOB  e.g. "results/mala_highres_m126_mp126_eta0.02_seed*.mat": the pooled chains after burn-in
           (the adaptation phase nAdapt)

For each posterior the nodes are split into the VW middle (-5 < x < 6 km), the VS edges, and the two
nodes on the steps (where a_true is ambiguous). Reported: the fraction of nodes with a_true in the 95 %
band, the relative L2 error of the mean (MAP for Laplace) and |mean - a_true|/std. a_true is only used
here, for evaluation.
"""
from glob import glob
import os
import sys

import numpy as np
from scipy.io import loadmat
from scipy.stats import norm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "context", "bayes",
                                "Computational-and-Bayesian-inverse-problems-code", "course_code"))
from samplers import ess, rank_rhat

VW = (-5.0, 6.0)


def regions(m_p):
    x = np.linspace(-10, 10, m_p)
    step = np.isclose(x, VW[0]) | np.isclose(x, VW[1])
    mid = (x > VW[0]) & (x < VW[1]) & ~step
    edge = ~mid & ~step
    return x, {"middle": mid, "edges": edge, "steps": step, "all": ~step}


def calibration(name, center, lo, hi, std, a_true):
    _, reg = regions(len(a_true))
    inside = (a_true >= lo) & (a_true <= hi)
    z = np.abs(center - a_true) / std
    parts = [f"{r} {inside[k].sum()}/{k.sum()}" for r, k in reg.items()]
    err = np.linalg.norm(center - a_true) / np.linalg.norm(a_true)
    mid = reg["middle"]
    print(f"{name:10s} 95 % coverage: {', '.join(parts)} | rel. err {err:.3e} | "
          f"|err|/std middle median {np.median(z[mid]):.1f} max {z[mid].max():.1f} | "
          f"std middle median {np.median(std[mid]):.1e}")
    return inside


def main(map_file, chain_globs):
    M = loadmat(map_file, squeeze_me=True, struct_as_record=False)
    a_true, a_map, H, noise = M["aTrue"], M["aMap"], M["H"], M["noise"]
    m_p = len(a_true)
    mtrue = getattr(noise, "Mtrue", np.nan)
    print(f"{os.path.basename(map_file)}: m_p = {m_p}, kappa = {M['inflation']:.3f}, "
          f"M(a_MAP)/M_noise = {np.atleast_2d(M['infHist'])[-1][1] / noise.Mnoise if np.size(M['infHist']) else np.nan:.3f}, "
          f"M(a_true)/M_noise = {mtrue / noise.Mnoise:.3f}")
    std_lap = np.sqrt(np.diag(np.linalg.inv(H)))
    q = norm.ppf(0.975)
    calibration("Laplace", a_map, a_map - q * std_lap, a_map + q * std_lap, std_lap, a_true)

    for pattern in chain_globs:
        files = sorted(glob(pattern))
        if not files:
            print(f"No files match {pattern}")
            continue
        runs = [loadmat(f, squeeze_me=True, struct_as_record=False)["chain"] for f in files]
        n = min(r.a.shape[1] for r in runs)
        burn = int(getattr(runs[0].opts, "nAdapt", 0))
        A = np.stack([r.a[:, burn:n].T for r in runs])
        U = np.stack([r.U[burn:n] for r in runs])
        rhat = max(rank_rhat(A[:, :, i]) for i in range(m_p))
        ess_min = min(ess(A[:, :, i]) for i in range(m_p))
        acc = np.mean([r.accepted[burn:n].mean() for r in runs])
        print(f"{os.path.basename(pattern)}: {len(runs)} chains, {sum(r.a.shape[1] - burn for r in runs)} samples kept, "
              f"acceptance {acc:.2f}, "
              f"Rhat max {rhat:.3f} (U {rank_rhat(U):.3f}), ESS min {ess_min:.0f}")
        # Pooled samples use the full length of each chain; R-hat and ESS above use the common length
        S = np.vstack([r.a[:, burn:].T for r in runs])
        lo, hi = np.percentile(S, [2.5, 97.5], axis=0)
        calibration("MCMC", S.mean(axis=0), lo, hi, S.std(axis=0, ddof=1), a_true)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2:])
