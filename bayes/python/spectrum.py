"""Spectrum of the posterior Hessian at a_MAP in the coordinates each sampler works in.

python spectrum.py [hessian .mat file]
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

here = os.path.dirname(os.path.abspath(__file__))
f = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "results", "hessian_m36_eta0.02.mat")
R = loadmat(f, squeeze_me=True)
H = (R["H"] + R["H"].T) / 2

spectra = {
    # eig(∇²U(a_MAP))
    "none ($\\mathbf{a}$)": np.sort(np.linalg.eigvalsh(H))[::-1],
    # eig(R_p^-T ∇²U R_p^-1) = 1 + eig(R_p^-T (∇²U - P) R_p^-1), P = R_p^T R_p the prior precision
    "prior": np.sort(1 + np.real(R["lambda"]))[::-1],
    # eig(R_H^-T ∇²U R_H^-1) = 1, R_H the Cholesky factor of ∇²U
    "Laplace": np.ones(H.shape[0]),
}

print(f"{'preconditioner':16s} {'mu_max':>9s} {'mu_min':>9s} {'cond':>9s} {'eps_HMC':>9s} {'h_MALA':>9s} "
      f"{'HMC steps':>10s} {'MALA steps':>10s}")
for name, mu in spectra.items():
    cond = mu[0] / mu[-1]
    # leapfrog stable for eps < 2/sqrt(mu_max), MALA drift for h < 4/mu_max; steps to cross the widest direction
    print(f"{name.split()[0]:16s} {mu[0]:9.2e} {mu[-1]:9.2e} {cond:9.2e} {2/np.sqrt(mu[0]):9.2e} {4/mu[0]:9.2e} "
          f"{np.sqrt(cond):10.1e} {cond:10.1e}")

fig, ax = plt.subplots(figsize=(5, 3.2))
for (name, mu), mk in zip(spectra.items(), "osd"):
    ax.semilogy(np.arange(1, mu.size + 1), mu, mk + "-", ms=4, label=name)
ax.set_xlabel("index $i$")
ax.set_ylabel("eigenvalue of the Hessian of $U$")
ax.legend(title="preconditioner", frameon=False)
fig.tight_layout()
out = os.path.join(here, "..", "..", "bayes_report", "figures", "spectrum_preconditioning.png" if len(sys.argv) == 1
                   else f"spectrum_preconditioning_{os.path.basename(f)[:-4]}.png")
fig.savefig(out, dpi=300)
print("saved", os.path.normpath(out))
