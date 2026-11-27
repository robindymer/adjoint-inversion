"""Summary figures of the posterior calibration runs (context/bayes/auto_improve_results.md), saved to
adjoint_rate_state_inversion/assets/bayes/final_*.png:

    final_maps.png         MAP estimates of a(x) for the baseline and the m = m_p = 201 priors
    final_bands.png        95 % bands against the true a: baseline MCMC, m201 MCMC, m201 GSBL Laplace
    final_zscores.png      |posterior mean - a_true|/std along x, the 95 % level marked
    final_coverage.png     95 % coverage per region

    python plot_final_results.py
"""
import os
from glob import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
OUT = os.path.join(HERE, "..", "..", "adjoint_rate_state_inversion", "assets", "bayes")
Q = norm.ppf(0.975)


def true_a(x):
    return np.where((x >= -5) & (x <= 6), 0.009, 0.013)


def load_map(name):
    M = loadmat(os.path.join(RES, name), squeeze_me=True, struct_as_record=False)
    a = M["aMap"]
    std = np.sqrt(np.diag(np.linalg.inv(M["H"]))) if "H" in M else None
    return np.linspace(-10, 10, len(a)), a, std


def load_chains(pattern):
    runs = [loadmat(f, squeeze_me=True, struct_as_record=False)["chain"]
            for f in sorted(glob(os.path.join(RES, pattern)))]
    burn = int(getattr(runs[0].opts, "nAdapt", 0))
    S = np.vstack([r.a[:, burn:].T for r in runs])
    lo, hi = np.percentile(S, [2.5, 97.5], axis=0)
    return S.mean(axis=0), S.std(axis=0, ddof=1), lo, hi


def posteriors():
    """(label, x, center, std, lo, hi) of each posterior that is compared."""
    out = []
    x, _, _ = load_map("highres_map_m151_mp26_eta0.02.mat")
    m, s, lo, hi = load_chains("mala_highres_m151_mp26_eta0.02_seed*.mat")
    out.append(("baseline m151/mp26, H-MALA", x, m, s, lo, hi))
    x, a, s = load_map("highres_map_m201_mp201_eta0.02.mat")
    out.append(("m201, Gaussian prior, Laplace", x, a, s, a - Q*s, a + Q*s))
    m, s, lo, hi = load_chains("mala_highres_m201_mp201_eta0.02_seed*.mat")
    out.append(("m201, Gaussian prior, H-MALA", x, m, s, lo, hi))
    x, a, s = load_map("highres_map_m201_mp201_eta0.02_gsbl.mat")
    out.append(("m201, GSBL prior (θ fixed), Laplace", x, a, s, a - Q*s, a + Q*s))
    return out


def regions(x):
    step = np.isclose(x, -5) | np.isclose(x, 6)
    mid = (x > -5) & (x < 6) & ~step
    return {"middle (VW)": mid, "edges (VS)": ~mid & ~step}


def plot_maps():
    xt = np.linspace(-10, 10, 1001)
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(xt, true_a(xt), color="C0", lw=2.5, label="true")
    for name, label, style in [
            ("highres_map_m151_mp26_eta0.02.mat", "baseline m151/mp26 (Gaussian)", dict(color="0.6", ls=":", marker="o", ms=3)),
            ("highres_map_m201_mp201_eta0.02.mat", "m201, Gaussian", dict(color="k", ls="--")),
            ("highres_map_m201_mp201_eta0.02_tv.mat", "m201, TV", dict(color="C2")),
            ("highres_map_m201_mp201_eta0.02_gsbl.mat", "m201, GSBL", dict(color="C3"))]:
        x, a, _ = load_map(name)
        mid = regions(x)["middle (VW)"]
        rms = np.sqrt(np.mean((a - true_a(x))[mid]**2))
        ax.plot(x, a, lw=1.3, label=f"{label}: RMS err middle {rms:.1e}", **style)
    ax.set(xlabel="x (km)", ylabel="a", title="MAP estimates")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc="upper center", ncol=2)
    ax.set_ylim(0.0060, 0.0195)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "final_maps.png"), dpi=200)


def plot_bands(posts):
    xt = np.linspace(-10, 10, 1001)
    sel = [posts[0], posts[2], posts[3]]
    fig, axes = plt.subplots(len(sel), 1, figsize=(10, 3.6*len(sel)), sharex=True, sharey=True)
    for ax, (label, x, c, s, lo, hi) in zip(axes, sel):
        at = true_a(x)
        out = (at < lo) | (at > hi)
        ax.fill_between(x, lo, hi, color="0.75", label="95 % band")
        ax.plot(x, c, "k-", lw=1.2, label="posterior mean / MAP")
        ax.plot(xt, true_a(xt), color="C0", lw=2, label="true")
        ax.plot(x[out], at[out], "rx", ms=5, label="true a outside band")
        ax.set_title(f"{label}: {np.sum(~out)}/{len(x)} nodes covered")
        ax.set_ylabel("a")
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=9, loc="upper right")
    axes[0].set_ylim(0.0060, 0.0190)
    axes[-1].set_xlabel("x (km)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "final_bands.png"), dpi=200)


def plot_zscores(posts):
    fig, ax = plt.subplots(figsize=(10, 4.5))
    styles = [dict(color="0.6", marker="o", ms=3, ls=":"), dict(color="k", ls="--"), dict(color="C1"), dict(color="C3")]
    for (label, x, c, s, lo, hi), st in zip(posts, styles):
        ax.semilogy(x, np.abs(c - true_a(x))/s, lw=1.2, label=label, **st)
    ax.axhline(Q, color="C0", lw=2, label="1.96 (95 % level)")
    ax.axvspan(-5, 6, color="C0", alpha=0.06)
    ax.set(xlabel="x (km)", ylabel="|mean - true| / std", title="Calibration along the fault (below the blue line = covered)")
    ax.set_ylim(1e-2, 1e2)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=9, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "final_zscores.png"), dpi=200)


def plot_coverage(posts):
    fig, ax = plt.subplots(figsize=(8, 4))
    names = list(regions(posts[0][1]))
    w = 0.8/len(posts)
    colors = ["0.6", "k", "C1", "C3"]
    for k, (label, x, c, s, lo, hi) in enumerate(posts):
        at = true_a(x)
        cov = [np.mean((at >= lo) & (at <= hi) & r)/np.mean(r) for r in regions(x).values()]
        bars = ax.bar(np.arange(len(names)) + (k - (len(posts) - 1)/2)*w, 100*np.array(cov), w,
                      color=colors[k], label=label)
        ax.bar_label(bars, fmt="%.0f", fontsize=8)
    ax.axhline(95, color="C0", ls="--", label="nominal 95 %")
    ax.set_xticks(range(len(names)), names)
    ax.set(ylabel="nodes covered (%)", ylim=(0, 105), title="95 % coverage of the true a")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "final_coverage.png"), dpi=200)


if __name__ == "__main__":
    posts = posteriors()
    plot_maps()
    plot_bands(posts)
    plot_zscores(posts)
    plot_coverage(posts)
