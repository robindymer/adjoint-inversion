
"""
Minimal iid Gaussian prior on the nodal values of a(x):  a_i ~ N(MEAN, GAMMA^2), independently.
Change MEAN, GAMMA or N_SHOW and rerun to see how the draws change.
"""
import numpy as np
import matplotlib.pyplot as plt
 
MEAN = 0.011     # prior mean (on the VW/VS boundary b)
GAMMA = 0.002    # pointwise standard deviation
N_SHOW = 3       # number of draws to plot
 
x = np.linspace(-10, 10, 26)                                 # 26 nodes on the fault (km)
a_true = np.where((x >= -5) & (x <= 6), 0.009, 0.013)        # true a(x): VW patch on [-5, 6] km
 
rng = np.random.default_rng(0)
draws = MEAN + GAMMA * rng.standard_normal((N_SHOW, len(x)))  # each row is one draw of a
 
for d in draws:
    plt.plot(x, d, ".-", lw=0.8)
plt.step(x, a_true, "k", where="mid", lw=1.5, label="a, true")
plt.axhline(0.011, color="0.4", ls=":", label="b")
plt.xlabel("x (km)")
plt.title(f"a_i ~ N({MEAN}, {GAMMA}²)")
plt.legend()
plt.show()
