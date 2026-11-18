"""
Figures for diagnostics_explained.md, made from the H-pCN and pCN chains in bayes/results/.

    python make_figures.py      (from this directory) writes figures/*.png
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import fault  # noqa: E402  (puts the course code on the path)
from diagnostic_reference import load_chains  # noqa: E402
from samplers import autocorrelation, iat, ess  # noqa: E402

RES = os.path.join(HERE, "..", "..", "results")
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

_, runs_h, A_h, U_h, _, _ = load_chains(os.path.join(RES, "hpcn_m36_eta0.02_s1_laplace_seed*.mat"))
_, runs_p, A_p, U_p, _, nad_p = load_chains(os.path.join(RES, "pcn_m36_eta0.02_prior_seed*.mat"))
BURN_H, BURN_P = 100, nad_p
x = fault.x


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name), dpi=130)
    plt.close(fig)
    print("saved", name)


# 1) Burn-in: where the chains start and when they have forgotten it
fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
for c in range(U_h.shape[0]):
    ax[0].plot(np.arange(1, 201), U_h[c, :200], lw=0.8)
ax[0].axvline(BURN_H, color="k", ls="--", lw=1)
ax[0].text(BURN_H + 3, ax[0].get_ylim()[1] * 0.99999, "burn-in ends", va="top", fontsize=8)
ax[0].set(xlabel="iteration", ylabel="$U$", title="H-pCN: dispersed starts, at equilibrium after a few steps")
for c in range(U_p.shape[0]):
    ax[1].semilogy(np.arange(1, U_p.shape[1] + 1), U_p[c], lw=0.8)
ax[1].axvline(BURN_P, color="k", ls="--", lw=1)
ax[1].axhline(np.mean(U_h[:, BURN_H:]), color="0.5", ls=":", lw=1)
ax[1].text(50, np.mean(U_h[:, BURN_H:]) * 1.3, "posterior level of $U$", fontsize=8, color="0.4")
ax[1].set(xlabel="iteration", ylabel="$U$ (log scale)", title="pCN: never reaches the posterior")
save(fig, "01_burn_in.png")

# 2) Autocorrelation and the integrated autocorrelation time (Geyer's initial positive sequence)
fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
node_h = 5   # node 6, the slowest for H-pCN
node_p = 17  # node 18, the slowest for pCN
xh, xp = A_h[0, BURN_H:, node_h], A_p[0, BURN_P:, node_p]
rh, rp = autocorrelation(xh, 40), autocorrelation(xp, 1000)
ax[0].bar(np.arange(41), rh, width=0.8, color="C0")
ax[0].axhline(0, color="k", lw=0.6)
ax[0].set(xlabel="lag $j$", ylabel=r"$\hat\rho(j)$", title=f"H-pCN, node {node_h + 1}: forgets after ~1 step")
ax[1].plot(np.arange(1001), rp, color="C3")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set(xlabel="lag $j$", ylabel=r"$\hat\rho(j)$", title=f"pCN, node {node_p + 1}: remembers for 1000+ steps")
# running tau = 1 + 2 sum rho over pairs, and where Geyer stops
rho = autocorrelation(xh)
pairs = rho[0:-1:2] + rho[1::2]
tau_run = -1 + 2 * np.cumsum(pairs)
stop = np.argmax(pairs <= 0) if np.any(pairs <= 0) else len(pairs)
ax[2].plot(2 * np.arange(len(tau_run[:150])), tau_run[:150], color="C0")
ax[2].axvline(2 * stop, color="k", ls="--", lw=1)
ax[2].annotate(f"first negative pair: stop here,\n$\\tau$ = {iat(xh):.2f}", xy=(2 * stop, tau_run[stop - 1]),
               xytext=(80, 1.1), fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.8))
ax[2].set(xlabel="lag up to which pairs are summed", ylabel=r"running estimate of $\tau$",
          title=f"H-pCN, node {node_h + 1}: summing all lags adds noise")
save(fig, "02_autocorrelation.png")

# 3) Rhat: do the chains agree? Histograms of each chain at one node
fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
for c in range(A_h.shape[0]):
    ax[0].hist(A_h[c, BURN_H:, 16], bins=40, histtype="step", lw=1.2, density=True, label=f"chain {c + 1}")
ax[0].set(xlabel="$a$ at node 17", ylabel="density", title=r"H-pCN: chains overlap, $\hat R \approx 1$")
ax[0].legend(frameon=False, fontsize=7)
for c in range(A_p.shape[0]):
    ax[1].hist(A_p[c, BURN_P:, node_p], bins=40, histtype="stepfilled", alpha=0.7, density=True, label=f"chain {c + 1}")
ax[1].set(xlabel=f"$a$ at node {node_p + 1}", ylabel="density",
          title=r"pCN: between-chain spread $\gg$ within, $\hat R = 4.5$")
ax[1].legend(frameon=False, fontsize=7)
save(fig, "03_rhat.png")

# 4) Rank normalization: values -> ranks -> normal scores, robust to heavy tails
rng = np.random.default_rng(0)
raw = np.exp(rng.standard_normal(2000) * 1.2)   # heavy right tail
r = rankdata(raw)
z = norm.ppf((r - 3 / 8) / (raw.size + 1 / 4))
fig, ax = plt.subplots(1, 3, figsize=(11, 3))
ax[0].hist(raw, bins=60, color="0.5")
ax[0].set(title="raw draws (heavy tail)", xlabel="value")
ax[1].hist(r, bins=40, color="0.5")
ax[1].set(title="ranks 1..S (uniform)", xlabel="rank")
ax[2].hist(z, bins=40, color="C0")
ax[2].set(title=r"normal scores $\Phi^{-1}$(rank) (Gaussian)", xlabel="z")
save(fig, "04_rank_normalization.png")

# 5) From samples at one node to a credible band, and the Laplace approximation on top
samples = A_h[:, BURN_H:].reshape(-1, A_h.shape[2])
H, a_map = runs_h[0]["H"], runs_h[0]["aMap"]
sd_lap = np.sqrt(np.diag(np.linalg.inv(H)))
fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
for k, i in enumerate([2, 12]):
    s = samples[:, i]
    ax[k].hist(s, bins=50, density=True, color="C0", alpha=0.35, label="samples")
    g = np.linspace(s.min(), s.max(), 300)
    ax[k].plot(g, norm.pdf(g, a_map[i], sd_lap[i]), "C1--", lw=1.5, label="Laplace $N(a_{MAP}, H^{-1})$")
    p = np.percentile(s, [2.5, 25, 75, 97.5])
    for q, lab in zip(p, ["2.5 %", "25 %", "75 %", "97.5 %"]):
        ax[k].axvline(q, color="C0", lw=1, ls=":")
        ax[k].text(q, ax[k].get_ylim()[1] * 0.97, lab, rotation=90, va="top", ha="right", fontsize=7)
    ax[k].axvline(s.mean(), color="C0", lw=1.5, label="sample mean")
    ax[k].axvline(runs_h[0]["aTrue"][i], color="k", lw=1.5, label="$a_{true}$")
    ax[k].set(xlabel=f"$a$ at node {i + 1} ($x$ = {x[i]:.1f} km)", ylabel="density")
    ax[k].ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useOffset=False)
ax[0].set_title("VS end: prior-dominated, wide")
ax[1].set_title("VW interior: data-dominated, narrow")
fig.legend(*ax[1].get_legend_handles_labels(), loc="lower center", ncol=4, frameon=False, fontsize=8)
fig.tight_layout(rect=(0, 0.07, 1, 1))
fig.savefig(os.path.join(OUT, "05_band_and_laplace.png"), dpi=130)
plt.close(fig)
print("saved 05_band_and_laplace.png")

# 6) ESS per chain and the summed ESS (what samplers.ess returns for shape (C, N))
fig, ax = plt.subplots(figsize=(5.5, 3.2))
per_chain = [ess(A_h[c, BURN_H:, node_h]) for c in range(A_h.shape[0])]
ax.bar(np.arange(1, 5), per_chain, color="C0")
ax.axhline(A_h.shape[1] - BURN_H, color="k", ls="--", lw=1)
ax.text(0.6, A_h.shape[1] - BURN_H + 15, "kept iterations per chain", fontsize=8)
ax.set(xlabel="chain", ylabel="ESS", ylim=(0, 1100), xticks=[1, 2, 3, 4],
       title=f"H-pCN node {node_h + 1}: ESS summed over chains = {sum(per_chain):.0f}")
save(fig, "06_ess_per_chain.png")
