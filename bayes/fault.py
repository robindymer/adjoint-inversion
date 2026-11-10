"""
Fault setup shared by the prior scripts (Table 6.1 of Stiernström, Almquist & Dunham, JCP 2024):

    x            parameter grid on the fault, m_p = 26 nodes on [-10, 10] km
    a_true(x)    direct-effect parameter: 0.009 on the velocity-weakening (VW) part [-5, 6] km, 0.013 elsewhere
    B            state-evolution parameter b = 0.011 (a < b is VW, a > b is velocity-strengthening)
    plot_draws, summarize    plotting and prior predictive checks for draws of shape (size, m_p)

Importing this module also puts the course code (priors.py, solvers.py) on the path.
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "Computational-and-Bayesian-inverse-problems-code", "course_code"))

X_MIN, X_MAX, M_P = -10.0, 10.0, 26
VW = (-5.0, 6.0)
A_VW, A_VS, B = 0.009, 0.013, 0.011

x = np.linspace(X_MIN, X_MAX, M_P)


def a_true(x):
    return np.where((x >= VW[0]) & (x <= VW[1]), A_VW, A_VS)


def plot_draws(ax, draws, title, n_show=5):
    """The first n_show draws on top of the true a(x) and the VW/VS boundary a = b."""
    xx = np.linspace(X_MIN, X_MAX, 1001)
    for k in range(min(n_show, len(draws))):
        ax.plot(x, draws[k], ".-", lw=0.8, ms=3)
    ax.plot(xx, a_true(xx), "k", lw=1.5, label="$a$, true")
    ax.axhline(B, color="0.4", ls=":", label="$b$")
    ax.set_xlabel("$x$ (km)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=7, loc="upper right")


def summarize(name, draws):
    """Prior predictive checks: unphysical draws (a <= 0 somewhere) and how much of the fault is VW (a < b)."""
    vw = draws < B
    print(f"{name}: P(a <= 0 somewhere) = {np.mean(np.any(draws <= 0, axis=1)):.3f}, "
          f"P(some VW node) = {np.mean(np.any(vw, axis=1)):.3f}, "
          f"mean VW fraction = {np.mean(vw):.3f} (true {np.mean(a_true(x) < B):.3f})")
