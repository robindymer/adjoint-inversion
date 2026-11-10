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

## 5. Making Newton-CG stand out
With high-res data the unregularized minimizer fits the model error, so compare the solvers on regularized problems only, by the time to reach the regularized minimizer (J*). Procedure:
1. **Step 4 (eps study):** if eps = 1e-3 makes Newton-CG faster on TV, the 1/eps stiffness is the cause. Use the best eps from then on. Expected: with eps = 0 Newton-CG loses, since R adds no curvature.
2. **Step 5a, m_p = 51** (`high-res_runs_step5.txt`, 4 runs): doubles the unknowns and worsens the conditioning, where L-BFGS should suffer more. Same β as the best in step 3 (tikhonov-grad 1e4, TV 10): for the step in the true `a`, Σ(Da)² and Σ|Da| do not change with the node spacing. Add the best eps from step 4 to the TV lines before submitting.
3. **Step 5b, m = 151** (`high-res_runs_step5b.txt`, 8 runs, m_p = 26): step 5a was not convincing (speedup about 2×, not growing with m_p). The model error drops from 1.5 to 0.47 in misfit, so more of the convergence is in the late, badly conditioned phase. About 2–3× the cost per run.
4. **Step 5c, hybrid L-BFGS → Newton-CG** (`x0_iter=<k>` in the drivers, `high-res_runs_hybrid.txt` on the cluster, `inverse-crime_runs_hybrid.txt` locally): `newton-cg-ls` starts from iteration k of the `lbfgs` run with the same settings, which must be in `mat/inversion/`. The output is `rough_lbfgs<k>-newton-cg-ls_...`, and its history contains the L-BFGS iterations 0..k, so `history.t` is the total wall time. k = 20 and 40 (J − J* about 1e-1 and 2e-2) for TV 10, TV 100 (both eps = 1e-4) and tikhonov-grad 1e4; k = 30 and 60 in the inverse crime. Original idea: L-BFGS through the nonconvex phase, then `newton-cg-ls` from that iterate, switching when the relative decrease per L-BFGS iteration stays below about 1e-3. It also removes the large-β stall, which comes from the 1/eps stiffness at the constant initial guess. Cheap test without new L-BFGS runs: add `x0_from=<file>; x0_iter=<k>` to the drivers and start `newton-cg-ls` from saved `history.x` iterates (inverse crime at f ≈ 2; high-res TV β = 100). Expected: much faster than L-BFGS in the inverse crime but only slightly faster than Newton-CG alone; removes the 6–7 h stall at large β.
5. **Joint inversion** (`a` + `tau0` or `D_c`): needs parameter scaling first (`optimizer_comparison.md` §6).

**Claim for the report:** Newton-CG reaches the regularized minimizer about 2× faster at the same accuracy in `a` (m = 76, tikhonov-grad 1e4: 0.31 h against 0.63 h; m = 251, TV β = 1: 3.7 h against 6.2 h, to rel. level 1e-3). The inverse crime remains the case where L-BFGS stalls completely.

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
# step 5a, m_p = 51 (third argument; "" keeps the m = 1001 data):
sbatch --array=1-4 high-res_sweep_a.sh high-res_runs_step5.txt "" 51
# step 5b, m = 151 (fourth argument), m_p = 26:
sbatch --array=1-8 high-res_sweep_a.sh high-res_runs_step5b.txt "" "" 151
# step 5c, hybrid (needs the step 3 lbfgs results in mat/inversion/):
sbatch --array=1-6 high-res_sweep_a.sh high-res_runs_hybrid.txt
```
Hybrid, inverse crime (local): `scripts/run_local_queue.sh direct_effect_inverse_crime scripts/inverse-crime_runs_hybrid.txt 2 4`.
Add `251` after the run file to use the m = 251 data instead (it must then be generated or copied there first). The m = 251 high-res results from before this plan (`mat/inversion/inversion/*_m251_mp26_*`) give the second data point for the claim in step 5. `summarize_inversions` assumes m_p = 26 for the true `a` and the resolved nodes, so it must be adapted for m_p = 51.

**Local** (step 1, inverse crime): `scripts/run_local_queue.sh direct_effect_inverse_crime scripts/inverse-crime_runs_step1.txt 5 3`.

**Summary:** `matlab -batch "pattern='rough_*_m76_mp26_velocity_highres_*'; summarize_inversions"` from `adjoint_rate_state_inversion/`.

## Findings so far
- Misfit of the true `a` sampled on the parameter grid, at m = 76 against the m = 251 data: 9.07 with m_p = 26, 3.85 with m_p = 76. The time step (k = 0.003 or 0.0005) makes no difference. The floor is dominated by how the step at x = −5 (between the nodes at −5.2 and −4.4) is represented with m_p = 26. It is only an upper bound on the reachable misfit; use the unregularized minimum from step 2 for the discrepancy principle.
- With the m = 1001 data (cluster, job 7536888) the same misfit is 8.95, 11.09 and 11.68 at m = 76, 151 and 251, against 9.07, 11.23 and 11.85 with the m = 251 data. The m = 251 data is therefore a good local stand-in.
- Unregularized L-BFGS reached f = 1.55 after 18 iterations, below the 8.95 of the true `a` sampled on the m_p = 26 grid, so the inversion fits better than the truth on this grid.
- **Step 1** (inverse crime, `newton-cg-ls`, local, so the wall times are inflated by sharing the laptop): regularization fixes the null space as expected. The global error drops from 0.0200 to 0.0009 (tikhonov-grad β = 1) and 0.0002 (TV β = 1e-2), and the resolved-node error from 0.0012 to 0.0003 and 0.0002. Larger β biases the result (TV β = 1: 0.018). It costs more iterations (24–118 against 24).
- **Step 2** (high-res, no regularization): both solvers reach F = 0.73, far below the 8.95 of the sampled truth, so that floor is useless for the discrepancy principle. `newton-cg-ls` converges in 19 iterations (0.86 h), `lbfgs` runs to 200 iterations (1.53 h), but the time to rel. level 1e-3 is about equal (0.50 h against 0.61 h). Newton-CG fits the model error better and gets the worse `a`: resolved-node error 0.29 against 0.12 (implicit regularization of L-BFGS).
- **Step 3** (β sweep): the smallest resolved-node errors are tikhonov-grad β = 1e4 (`lbfgs` 0.073, `newton-cg-ls` 0.086) and TV β = 10 (0.085 and 0.103). They are 100× worse than in the inverse crime. At equal β both solvers reach the same objective; `newton-cg-ls` is faster for tikhonov-grad 1e3–1e4 (0.31–0.38 h against 0.63 h to level 1e-3), slower for TV.
- **Newton-CG stalls at large β** (TV β ≥ 100, tikhonov-grad 1e6): about 50 iterations with J unchanged and R alternating between 1e-5 and 1e-12 before it suddenly converges (6–7 h; TV β = 1000 still stalled after 65 iterations, 10 h). The regularization Hessian at the constant initial guess is huge (β/eps·D'D for TV, β·D'D for tikhonov-grad), and the step zig-zags in the stiff modes. The CG flags in the cluster logs should show whether CG stops on negative curvature.
- **Step 4** ran at TV β = 100 (the cluster had the old run file), which is the stiff case; four runs were still unfinished at download (about 1.8 h in). Exact objective F + 100·TV and resolved-node error:
  - eps = 1e-5: `newton-cg-ls` makes no progress at all (J = 32.86 after 15 iterations; curvature β/eps = 1e7). `lbfgs` is fine (1.56, error 0.084).
  - eps = 1e-3: `newton-cg-ls` converges in 19 iterations (1.42 h, J = 1.52) without the stall seen at eps = 1e-4 (64 iterations, 7 h). `lbfgs` was at 1.53 after 186 iterations. This confirms that the 1/eps stiffness causes the stall. The smoothing biases the minimizer: J = 1.52 against 1.38 with eps = 1e-4.
  - eps = 0: `lbfgs` works (J = 1.65 after 184 iterations, not yet converged) and gives the best global error (0.074, resolved 0.087). `newton-cg-ls` crawls (J = 4.89 after 10 iterations, about 10 min each, probably from line-search backtracking).
  - Verdict: eps = 0 is fine for L-BFGS, not for Newton-CG. Newton-CG needs eps large enough that β/eps stays moderate (eps = 1e-3 at β = 100).
  - Still to run: β = 10 (the run file now says β = 10).

- **Step 5a** (m_p = 51, all 4 runs finished). Error = relative L2 error of the linearly interpolated `a` against the step on x ∈ [−6, 9.2], the resolved region; it can be compared across m_p. Speedup = L-BFGS time for 200 iterations / Newton-CG time to the same J (L-BFGS never reaches level 1e-4: it ends at about 2e-4).

  | m_p | case | error L-BFGS / NCG | NCG iter. | L-BFGS 200 it. | NCG to L-BFGS's J | speedup |
  |---|---|---|---|---|---|---|
  | 26 | tikhonov-grad 1e4 | 0.092 / 0.099 | 13 | 1.53 h | 0.34 h | 4.5 |
  | 51 | tikhonov-grad 1e4 | 0.096 / 0.099 | 28 | 1.66 h | 0.72 h | 2.3 |
  | 26 | TV 10, eps 1e-4 | 0.113 / 0.118 | 27 | 1.52 h | 1.80 h | 0.8 |
  | 26 | TV 100, eps 1e-3 | — | 19 | 1.88 h | 1.08 h | 1.7 |
  | 51 | TV 10, eps 1e-3 | 0.105 / 0.104 | 30 | 1.38 h | 0.76 h | 1.8 |

  Newton-CG wins both m_p = 51 cases by about 2×, but the gap did not grow with m_p: L-BFGS gets no slower (one gradient per iteration either way), while Newton-CG needs twice the Newton iterations for tikhonov-grad. Finer m_p does not improve `a` (the error is set by the model error). So step 5a is not the convincing case; next is the hybrid (L-BFGS → Newton-CG) or step 5b (m = 151).
- **Step 5c, inverse crime** (local, unfinished at 16:13 but clear): Newton-CG started from L-BFGS iteration 30 (f = 2.28) or 60 (f = 1.38) stops at f = 1.3616 with ||g|| ≈ 1–5 (against 2e3 at the start), the same place as L-BFGS (1.3624 after 200 iterations), with error 0.19–0.20 in `a`. Plain Newton-CG from the constant guess reaches f = 9e-7 and error 0.02 in 24 iterations (0.6 h). So L-BFGS's plateau at 1.36 is a **local minimum (or a very flat saddle) that L-BFGS has already entered by iteration 30**, not slow convergence in an ill-conditioned valley. Newton-CG wins in the inverse crime because its early steps take a different path, not because it converges faster near the end. This corrects the "late, convex valley" explanation above. The hybrid cannot help here; the high-res runs (regularized) are still a fair test.

## Jobs running after 3 november
[robind@pelle1 scripts]$ squeue -u robind
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
         7540683_4     pelle high-res   robind  R    8:49:38      1 p62
         7540683_2     pelle high-res   robind  R    8:58:21      1 p59
        7536889_10     pelle high-res   robind  R   20:44:36      1 p105
         7546731_1     pelle high-res   robind  R      39:05      1 p60
         7546731_2     pelle high-res   robind  R      39:05      1 p114
         7546731_3     pelle high-res   robind  R      39:05      1 p114
         7546731_4     pelle high-res   robind  R      39:05      1 p114
         7546731_5     pelle high-res   robind  R      39:05      1 p114
         7546731_6     pelle high-res   robind  R      39:05      1 p52
         7546723_1     pelle high-res   robind  R      39:36      1 p68
         7546723_2     pelle high-res   robind  R      39:36      1 p54
         7546723_3     pelle high-res   robind  R      39:36      1 p54
         7546723_4     pelle high-res   robind  R      39:36      1 p54
         7546723_5     pelle high-res   robind  R      39:36      1 p54
         7546723_6     pelle high-res   robind  R      39:36      1 p60
         7546723_7     pelle high-res   robind  R      39:36      1 p60
         7546723_8     pelle high-res   robind  R      39:36      1 p60
[robind@pelle1 scripts]$ - **Step 5b** (m = 151, m_p = 26, all 8 finished). Error as in step 5a (interpolated, x ∈ [−6, 9.2]); speedup = L-BFGS time for 200 iterations / Newton-CG time to the same J.

  | case | error L-BFGS / NCG | NCG iter. | L-BFGS 200 it. | NCG to L-BFGS's J | speedup | (m = 76) |
  |---|---|---|---|---|---|---|
  | tikhonov-grad 1e3 | 0.082 / 0.177 | 16 | 5.36 h | never | NCG fails | 2.3 |
  | tikhonov-grad 1e4 | 0.087 / 0.094 | 13 | 4.47 h | 1.82 h | 2.5 | 4.5 |
  | TV 1, eps 1e-3 | 0.088 / 0.088 | 22 | 4.09 h | 2.34 h | 1.7 | — |
  | TV 10, eps 1e-3 | 0.084 / 0.100 | 16 | 4.05 h | 1.61 h | 2.5 | — |

  Less model error helps `a` (L-BFGS 0.082–0.088 against 0.092–0.113 at m = 76), but the speedup does not grow: about 2× as at m = 76. Newton-CG at tikhonov-grad 1e3 never leaves the constant initial guess (J 32.86 → 32.84 in 16 iterations, 3.5 h, flag 0). β·D'D is mild there, so this stall comes from the misfit's curvature at the start, not only from the regularization stiffness.
- **Step 5c, hybrid, high-res m = 76** (Newton-CG from L-BFGS iteration k; times are totals, including the L-BFGS part):

  | case | L-BFGS (200 it.) | Newton-CG | hybrid k = 20 | hybrid k = 40 |
  |---|---|---|---|---|
  | tikhonov-grad 1e4, t to 1e-4 | > 1.53 h | **0.35 h** | 0.57 h | 0.92 h |
  | TV 10 (eps 1e-4), t to 1e-4 | > 1.52 h | 1.89 h | **0.98 h** | 1.75 h |
  | TV 100 (eps 1e-4), t to 1e-4 | > 1.46 h | 6.47 h (stall) | **1.39 h** | 1.39 h |

  All runs end at the same J* and `a` as Newton-CG. The hybrid with k = 20 removes the Newton-CG stall (4.7× faster at TV 100) and beats L-BFGS where L-BFGS is slow to finish (TV 10: 0.82 h against 1.52 h to L-BFGS's final J), but at TV 100 it is only about as fast as L-BFGS (1.35 h against 1.46 h to L-BFGS's final J; L-BFGS reaches 1e-3 first). When Newton-CG alone does not stall, starting it directly is fastest. Earlier switching is better: k = 40 spends too long in L-BFGS.
- **TV and edges** (`assets/high_res_tv_a_eps.pdf`): φ_ε(x) ≈ x²/(2ε) for |x| ≪ ε, and the node differences of `a` are about 1e-4 (median |Da| 1.4e-4 at ε = 1e-3, 4.6e-5 at 1e-4), so smoothed TV with ε ≥ 1e-4 acts as tikhonov-grad with weight β/ε and gives ramps. Only exact TV (ε = 0, L-BFGS, β = 100) is piecewise constant (median |Da| 2.7e-6), with the best error (0.086 against 0.104 and 0.120). Newton-CG with ε = 0 fails: zero R Hessian, gradient jumps at the kinks, J = 4.43 after 50 iterations (14 h, stopped) against 1.65 for L-BFGS. Local runs with exact TV and β = 3, 10, 30 (`high-res_runs_tv_exact.txt`, m = 251 data) test whether a smaller β removes the remaining contrast loss.
