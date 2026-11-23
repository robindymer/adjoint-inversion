""" Diagnostic and verification, based on 4.6 in Glaubitz course litterature.
Want to answer
1) Has the chain forgotten its start (diagnostic)
2) Value of ESS, i.e. how many independent samples is it worth (diagnostic)
3) Does it sample the right posterior (verification)
"A sampler never reports its own failure, and every result needs several chains,
R hat, effective sample sizes, and a verification against a known answer" - chapter 4 summary, Glaubitz course litterature.
See Figure 4.18 for inspiration.
"""
from glob import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

# Such that we cant import from samplers.py
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "context", "bayes",
                                "Computational-and-Bayesian-inverse-problems-code", "course_code"))
from samplers import ess, rank_rhat

# First, plot QOI (a(x)) against iterations for each chain
def load_chains(pattern): # e.g. pattern is hpcn_m36_eta0.02_s1_laplace_seed*.mat
    files = sorted(glob(pattern))
    if not files:
        sys.exit(f"No files match {pattern}")
    # Load matlab files
    # runs has one element for each chain
    runs = [loadmat(f, squeeze_me=True, struct_as_record=False) for f in files]
    n = min(r["chain"].a.shape[1] for r in runs) # cut to shortest chain length
    # E.g. access chain 1, iteration 100 at fault point 10 with A[0,99,9]
    A = np.stack([r["chain"].a[:, :n].T for r in runs]) # a(x) values, size (chains, iterations, m_p)
    U = np.stack([r["chain"].U[:n] for r in runs]) # value of potential, size (chains, iterations)
    acc = np.stack([r["chain"].accepted[:n] for r in runs])
    n_adapt = int(getattr(runs[0]["chain"].opts, "nAdapt", 0))
    return files, runs, A, U, acc, n_adapt

def main(pattern, burn_in=None):
    files, runs, A, U, acc, n_adapt = load_chains(pattern)
    # print(files)
    # print(runs[3]["chain"].alpha)
    print(A.size, A)
    print(A[0,10,5])
    # print(U)
    # print(acc)
    # print(n_adapt)

    """ Trace plots and burn-in """
    # With Laplace preconditoning, no burn-in is seen. Need to have step-size s high enough for burn-in
    plt.figure()
    for i, r in enumerate(runs):
        print(f"Chain {i+1}: {r['chain'].a.shape[1]} iterations, acceptance rate {np.round(r['chain'].accepted.mean(), 2)}")

        sample_size = r['chain'].a.shape[1]
        a_index = 13
        a_values = A[i,:,a_index]
        plt.plot([s for s in range(a_values.size)], a_values)
        # for j in range(r['chain'].a.shape[1]):
        #     print(A[i, j, 10])
    plt.figure()
    for i, r in enumerate(runs):
        sample_size = r['chain'].a.shape[1]
        U_values = U[i,:]
        plt.plot([s for s in range(U_values.size)], U_values)

    plt.show()

    # NOTE: If burn_in is seen, we need to add something like
    # burn_in = max(n_adapt, N // 10) or manually

    """ Autocorrelation time (can be derived from ESS) """
    # ESS, autocorrelation and Rhat at every node AND for U
    C, N, m_p = A.shape
    print(f"C: {C}, N: {N}, m_p: {m_p}")
    rhat_a = np.array([rank_rhat(A[:,:,i]) for i in range(m_p)])
    ess_a = np.array([ess(A[:,:,i]) for i in range(m_p)])
    rhat_U, ess_U = rank_rhat(U), ess(U)
    n_solves = C * N
    print(f"For U: Rhat {rhat_U}, ESS {ess_U}")
    print(f"For a(x): Rhat max {rhat_a.max()} at node {rhat_a.argmax() + 1}")
    print(f"ESS min {ess_a.min()} at node {ess_a.argmin() + 1}, median {np.median(ess_a)}")
    print(f"Forward solves per effective sample at worst node: {n_solves / ess_a.min()}")
    # Run at least four chains and to requirê R < 1.01 and Neff > 400 for every reported quantity 
    ok = rhat_a.max() < 1.01 and rhat_U < 1.01 and ess_a.min() > 400 and ess_U > 400
    print(f"Rhat < 1.01 and ESS > 400 for U and all nodes: {'yes' if ok else 'NO'}")

    # Posterior summaries
    a_true = runs[0]["aTrue"]
    samples = A.reshape(-1, m_p)
    print(samples.size) # 4 * 1000 * 26
    mean, std = samples.mean(axis=0), samples.std(axis=0, ddof=1)
    print(f"Mean: {mean}, std: {std}")
    lo, hi = np.percentile(samples, [2.5, 97.5], axis=0)
    mcse = std / np.sqrt(ess_a) # monte carlo error of the posterior mean
    print(f"|mean - a_true|/std: max {np.max(np.abs(mean - a_true) / std):.2f}; "
          f"a_true inside the 95 % band at {np.sum((a_true >= lo) & (a_true <= hi))} of {m_p} nodes")
    print(f"Monte Carlo error of the mean / std: max {np.max(mcse / std):.3f}")

    # How well Laplace approximates posterior
    if "H" in runs[0]:
        H, a_map = runs[0]["H"], runs[0]["aMap"]
        std_lap = np.sqrt(np.diag(np.linalg.inv(H)))
        laplace = (a_map, std_lap)
        print(f"vs Laplace: |mean - a_MAP|/std max {np.max(np.abs(mean - a_map) / std):.2f}, "
                f"std/std_Laplace in [{np.min(std / std_lap):.3f}, {np.max(std / std_lap):.3f}]")

if __name__ == "__main__":
    try:
        pattern = sys.argv[1]
    except IndexError:
        print("ERROR: Must provide pattern, e.g. \"../results/hpcn_m36_eta0.02_s1_laplace_seed*.mat\"")
    burn_in = int(sys.argv[2]) if len(sys.argv) > 2 else None
    main(pattern, burn_in)

""" Structure of runs """
# run = loadmat(f)                      # dict, one per chain/seed
# ├── m        int 36                   # mesh / grid resolution
# ├── m_p      int 26                   # number of parameters in a(x)
# ├── s        1                        # pCN step size
# ├── seed     1
# ├── start    'laplace'                # chain started at the Laplace approx. (MAP)
# ├── aTrue    (26,)                    # true a(x)
# ├── aMap     (26,)                    # MAP estimate
# ├── H        (26, 26)                 # Hessian at the MAP
# ├── noise    struct
# │   ├── eta    0.02                   # relative noise level
# │   ├── sigma  0.0325                 # absolute noise std
# │   ├── seed   1
# │   ├── dt_d   0.003                  # data time step
# │   └── N      176088                 # number of data points
# └── chain    struct
#     ├── a         (26, 1000)          # samples: parameters × iterations  ← the QoI
#     ├── U         (1000,)             # neg. log posterior per iteration
#     ├── phi       (1000,)             # data misfit (neg. log-likelihood)
#     ├── alpha     (1000,)             # acceptance probability min(1, ·)
#     ├── accepted  (1000,) uint8       # 0/1 per iteration
#     ├── s         (1000,)             # step size used per iteration
#     ├── nFwd      1000                # forward solves
#     ├── time      13071 s
#     └── opts      struct
#         ├── ref        struct         # reference Gaussian for hpCN
#         │   ├── mu  (26,)             #   mean (≈ aMap)
#         │   └── R   (26, 26)          #   Cholesky / precision factor
#         ├── nAdapt     0              # adaptation iterations (0 = s fixed)
#         ├── target     0.25           # target acceptance rate
#         ├── saveFile, saveEvery, printEvery
