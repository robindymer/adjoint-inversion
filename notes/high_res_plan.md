# Plan: high-resolution data inversion and regularization

**Goal:** check whether Newton-CG still beats L-BFGS when the data are not an inverse crime.
- Data: high-res data generated locally at m = m_p = 251, k = 0.0005 (`matlab -batch "data_m=251; fractal_fault_generate_data"`), so that no cluster is needed. m = 76 differs from m = 251 by 20% (`optimizer_comparison.md` §7), which is what breaks the inverse crime. The m = 1001 data (cluster only) can be used with the same scripts by leaving out `data_m`.
- Inversion: m = 76, m_p = 26, `a` only, initial guess a = 0.0135, 88 receivers.
- Solvers: `lbfgs` (fmincon, unbounded) vs `newton-cg-ls`.
- Also: test whether TV works without smoothing (eps = 0).

**Why the outcome is open:**
- In the inverse crime, Newton-CG won in the late, convex phase (f < 1e-2).
- With high-res data, the m = 76 discretization error alone is worth a misfit of about 1.5 (`optimizer_comparison.md` §7), against an initial misfit of 29.8.
- Convergence may therefore stop before that phase, and the two solvers could look about equal.

**Inverse-crime baseline at m = 76** (saved runs; iterations / wall time to reach each level):

| | f < 1 | f < 1e-2 | f < 1e-5 | final f |
|---|---|---|---|---|
| `newton-cg-ls` | 11 / 23 min | 17 / 35 min | 21 / 46 min | 8.7e-7 |
| `CG` (fmincon ip) | 15 / 71 min | 21 / 89 min | 25 / 107 min | 8.7e-7 |
| `lbfgs` (unbounded) | – | – | – | 1.36 (stalled) |

## Metrics (all steps)
- **Objective against wall time**, one plot with both solvers (`history.t`). For each run, also the wall time to reach (f − f\*)/(f0 − f\*) = 1e-1, 1e-2 and 1e-3, where f\* is the best objective any run reached.
- **Relative error in `a`**, both over all nodes and over the resolved nodes 6–25.
  - Over all nodes, the error is dominated by the unresolved edge nodes 1–5 and 26.
  - At m = 251 the global error was 0.0215 against 0.0266, but the resolved-node error was 2.3e-4 against 1.9e-2.
- **Final `a`**: look for oscillations (compare Fig. 1f in the report).
- **Speedup**: L-BFGS time / Newton-CG time to reach the same objective level.

## 0. Preparation (local)
1. Let `high_res_data_inversion.m` take `m` from outside (it is fixed at 251).
2. Support eps = 0 in `+adjopt/regularization.m`.
   - At present the initial guess is constant, so Dp = 0, and Dp./|Dp| = 0/0 gives NaN.
   - With eps = 0, use the subgradient `sign(Dp)` and a zero Hessian term.
3. **Misfit floor:** compute f(a_true sampled on m_p = 26) at m = 76 against the high-res data. This is one forward solve.
   - The floor includes the representation error of the step at x = −5, which is not a grid node.
4. Write a Slurm array script for the sweeps.

## 1. Regularization in the inverse-crime setting
- The true `a` is constant (0.013) out to both edges. TV or tikhonov-grad should therefore pull the unresolved edge nodes toward their neighbours and remove most of the 2% error.
- Plain Tikhonov keeps them at the 0.0135 prior instead.
- Runs: TV and tikhonov-grad, a few values of β each, `newton-cg-ls`.
- This separates two effects: fixing the null space, and damping the fit to the discretization error (step 3).

## 2. High-res data, no regularization
- `lbfgs` and `newton-cg-ls`.
- **Expected:** both stall near the floor from 0.3, at similar cost, possibly with oscillations in `a`.

## 3. β sweep, both solvers
- **TV**, eps = 1e-4 (2.5% of the 0.004 jump): β ∈ {1, 10, 100, 1000}. For scale, TV(a_true) = 0.008.
- **tikhonov-grad:** β ∈ {1e3, …, 1e6}. For scale, Σ(Da_true)² ≈ 3e-5.
- **Choice of β:** discrepancy principle (final misfit ≈ the floor) or an L-curve.
- **Main result:** the speedup as a function of β.
  - Regularization adds curvature exactly in the poorly resolved, oscillatory modes where L-BFGS fails.
  - It may therefore help L-BFGS more than Newton-CG and narrow the gap.

## 4. eps study (TV at the chosen β, both solvers)
- eps ∈ {0, 1e-5, 1e-4, 1e-3}, plus optionally a continuation run (1e-3 → 1e-5, warm-started).
- **Compare all runs on the exact objective F + β·TV**, not on each run's own smoothed objective.
- **Stopping:** at eps = 0 the gradient test never triggers, so runs stop at maxIter. Compare them by time to reach a given objective level.
- **Expected:**
  - L-BFGS often copes with nonsmooth objectives (Lewis & Overton 2013).
  - `newton-cg-ls` gets no curvature from R and may zig-zag around the steps.
- **"Got away with it"** means: the exact objective is no worse than with smoothing, the step is sharp, and there are no NaNs or stalls.

## 5. If Newton-CG still does not stand out
In order of effort:
1. **m_p = 51:** more unknowns and worse conditioning.
2. **m = 151:** the misfit floor is about 3× lower, so more of the convex phase is reached.
3. **Joint inversion** (`a` + `tau0` or `D_c`): needs parameter scaling first (`optimizer_comparison.md` §6).

## Cost
- Each m = 76 run takes about 1–3 h on 8 cores.
- Steps 1–4 need about 25 jobs. The jobs within a step are independent, and steps 1 and 2 can be submitted together.

## For the report (new §3.4)
- **Setup table:** Table 1 describes the Hessian test (m = 26, order 6, T = 4), not the inversions. The new section needs its own table: data at m = 1001 with k = 0.0005, inversion at m = 76, order 8, T = 6, regularization type, β and eps.
- **Convergence plot:** objective against wall time, both solvers.
- **§3.3:** state which error the comparison refers to. The global error hides the 80× difference in the resolved-node error.

## Running
**Cluster** (steps 0.3 and 2–4, with the m = 1001 data already there). Submit from `scripts/` after pulling the code:
```bash
sbatch high-res_floor.sh
sbatch --array=1-$(grep -cv '^#' high-res_runs_step23.txt) high-res_sweep_a.sh high-res_runs_step23.txt
# after choosing beta in step 3 (edit beta in the run file first):
sbatch --array=1-$(grep -cv '^#' high-res_runs_step4.txt) high-res_sweep_a.sh high-res_runs_step4.txt
```
Add `251` after the run file to use the m = 251 data instead (it must then be generated or copied there first).

**Local** (step 1, inverse crime): `scripts/run_local_queue.sh direct_effect_inverse_crime scripts/inverse-crime_runs_step1.txt 5 3`.

**Summary:** `matlab -batch "pattern='rough_*_m76_mp26_velocity_highres_*'; summarize_inversions"` from `adjoint_rate_state_inversion/`.

## Findings so far
- Misfit of the true `a` sampled on the parameter grid, at m = 76 against the m = 251 data: 9.07 with m_p = 26, 3.85 with m_p = 76. The time step (k = 0.003 or 0.0005) makes no difference. The floor is dominated by how the step at x = −5 (between the nodes at −5.2 and −4.4) is represented with m_p = 26. It is only an upper bound on the reachable misfit; use the unregularized minimum from step 2 for the discrepancy principle.
