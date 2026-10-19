<!-- the cost of a gradient versus an Hv, and the Hv correctness checks;
why the start is non-convex;
the benchmark tables at m=36 and m=76;
the null-space finding at the left edge of the fault;
the resolution study;
how trustRegionNewtonCG.m works and how to call it;
the next steps. -->

# Optimizer comparison for the rate-and-state inversion

Comparison of optimization routines for inverting the rate-and-state parameter `a` on the rough (fractal) fault, `direct_effect_inverse_crime.m` setup, using the gradient and the Hessian-vector product (Hv).

**Problem and settings**
- Synthetic (inverse-crime) data and velocity misfit.
- `optPars = {'a'}`, constant initial guess `a = 0.0135`, `m_p = 26` parameter grid points.
- Grid sizes tested: `m = 36`, a confirmation run at `m = 76`, and a single `newton-tr` run at `m = 251` (cluster).

**Recommendation:** use the new trust-region Newton-CG option, `opt_routine_str = 'newton-tr'` (`+adjopt/trustRegionNewtonCG.m`). It was the only method that converged to the misfit floor at both `m = 36` and `m = 76`. With unchanged settings it also converged at `m = 251` (section 4b).

## 1. Cost of the building blocks

Measured at `m = 36` with 6 threads:

| operation | time | in gradients |
|---|---|---|
| gradient (forward + adjoint) | 5.6 s | 1 |
| misfit only (forward) | 4.3 s | 0.77 |
| Hv, first-order solutions reused | 4.1 s | 0.73 |
| Hv, first-order solutions recomputed | 9.0 s | 1.6 |

An Hv costs less than a gradient when the forward and adjoint solutions at the current point are reused, which `computeHessianVector` already does via `firstOrderParsVec`.

**Hv correctness**
- Agreement with a central finite difference of the gradient: relative error 7e-6.
- Symmetry of the full 26×26 Hessian: relative asymmetry 1e-15.

Below, *cost* means equivalent gradient evaluations, `(4.3·nFwd + 1.3·nAdj + 4.1·nHv)/5.6`. This counts PDE solves, so it is independent of machine load.

## 2. The landscape is strongly non-convex at the initial guess

The full Hessian at `a = 0.0135` has **all 26 eigenvalues negative**, ranging from −165 to −5e-6, with |condition number| ≈ 3e7. With `a > b` everywhere the fault is velocity strengthening and does not rupture. Lowering `a` switches rupture on, so the misfit is locally concave.

After the first steps the misfit enters a narrow, stiff valley around f ≈ 9:
- |g| jumps from 0.37 to about 1e4.
- Every method spends several iterations there.

Consequences for the methods:
- A line-search Newton-CG falls back to the steepest-descent direction when it detects negative curvature. Examples are `matlab_impl_constr/optimization/newton_cg.m` and fmincon interior-point with `'cg'`.
- A trust-region (Steihaug) CG instead follows the negative-curvature direction to the trust-region boundary and lets the actual reduction decide the step length.
- A loose CG tolerance turns Newton-CG into scaled steepest descent in the valley. The forcing term is therefore capped at η ≤ 0.1.

## 3. Benchmark at `m = 36`

All methods had a budget of 250 gradient-equivalents.

| method | final f | final err | cost to f < 1 | f < 1e-1 | f < 1e-2 | f < 1e-3 |
|---|---|---|---|---|---|---|
| **TR Newton-CG + L-BFGS preconditioner** | 8.8e-8 | 0.018* | **41** | **74** | **81** | **98** |
| TR Newton-CG, no preconditioner | 1.9e-7 | 0.023 | 73 | 89 | 102 | 102 |
| L-BFGS (own implementation, backtracking line search) | 4.9e-5 | 0.026 | 45 | 88 | 121 | 154 |
| fmincon interior-point + `'cg'` + Hv (`'CG'` in commit 7164513, `lb = 0`) | 1.3e-3 | 4.03 | 120 | 143 | 147 | 171 |
| fmincon trust-region-reflective + Hv | 0.34 | 1.06 | 232 | – | – | – |
| fmincon interior-point L-BFGS, no lower bound (current `'lbfgs'`) | 4.3e-4 | 0.030 | 45 | 93 | 147 | 219 |
| fmincon interior-point L-BFGS, `lb = 0` (`'lbfgs'` in commit 7164513) | 2.7 | 0.20 | – | – | – | – |
| fminunc quasi-Newton (BFGS) | 2.7 | 0.17 | – | – | – | – |

`err` is the relative L2 error in `a` on the parameter grid; see section 5 for why it bottoms out around 0.02.

\* This error is from the variant that switches to L-BFGS at the end. The plain preconditioned run reached the same f but then took a large step in the null space (err 0.26), which led to the fixes in section 6.

**Observations**
- **The lower bound `lb = 0` is what stalls fmincon L-BFGS; it is not path dependence.**
  - Your commit 7164513 added `lb = zeros` to the `'lbfgs'` and `'CG'` options. Interior-point then adds a log-barrier term μ·log(a), and with a ≈ 0.01 that barrier distorts the steps.
  - With the bound, the run never got below f = 2.7.
  - Without it, the run reproduces the saved pre-7164513 run (`mat/inversion/rough_lbfgs_maxiter200_m36_..._1.mat`) iterate for iterate, and reaches f = 4.3e-4.
  - The bound is not needed in practice: `a` never went below 0.0074 in unbounded runs.
  - The unbounded variant is therefore the fair L-BFGS baseline. Even so, `newton-tr` reaches f < 1e-3 at about half its cost (98 against 219 here).
- **The `'CG'` option (with `lb = 0`, as tested) converges in misfit but drifts far off in `a`** (err 4.0), in directions the data cannot see.
- **The L-BFGS preconditioner built from CG pairs (Morales & Nocedal) helps most in the non-convex early phase.** The cost to reach f < 1 drops from 73 to 41.
- **Switching from Newton-CG to L-BFGS once the Hessian looks positive definite gave no speed-up.** The switch only triggered after f had already converged.

## 4. Confirmation at `m = 76`

The three most relevant methods were run in parallel, 3 threads each. Cost weights are those measured at `m = 36`; wall time is also given.

| method | final f | final err | cost to f < 1 / 1e-1 / 1e-2 / 1e-3 / 1e-5 | minutes to the same |
|---|---|---|---|---|
| **`newton-tr` (repo version)** | **8.7e-7** | **0.022** | 49 / 92 / 117 / 122 / 144 | 22 / 42 / 54 / 56 / 66 |
| fmincon interior-point L-BFGS, `lb = 0` (`'lbfgs'` in commit 7164513) | 1.0e-3 (budget hit) | 0.081 | 39 / 61 / 93 / – / – | 18 / 29 / 44 / – / – |
| L-BFGS (own implementation) | 1.36 (stalled) | 0.126 | – | – |

`newton-tr` stopped on its own (relative gradient tolerance) at cost 172, about 79 minutes.

The fmincon run here used the lower bound. The unbounded variant (section 3) was not run at `m = 76`.

The pattern holds at both resolutions:
- **Early, non-convex phase:** L-BFGS is comparable or slightly cheaper down to f ≈ 1e-2.
- **Late phase, once the Hessian is positive definite:** Newton-CG converges much faster and reaches the misfit floor, while L-BFGS crawls.

## 4b. `m = 251` (cluster): `newton-tr`, Hv check and L-BFGS comparison

One run on the cluster (`scripts/newton_tr_m251.sh`, 8 cores, Slurm job 7213723). It used the same harness and the same settings as at `m = 76`. Results are in `optimizer_benchmark/results/res_newton_tr_repo_m251.mat`.

| | `m = 36` | `m = 76` | `m = 251` |
|---|---|---|---|
| initial f | 29.8 | 29.8 | 36.5 |
| final f | 8.8e-8 | 8.7e-7 | 2.4e-5 |
| final err | 0.018 | 0.022 | 0.021 |
| outer iterations | 54 | 52 | 49 |
| cost (gradient-equivalents) | 201* | 172 | 166 |
| wall time | 21 min (3 threads) | 79 min (3 threads) | 5.9 h (8 cores) |

\* Prototype with the switch to L-BFGS (`trncg_pc_switch`); the repo version was not run to convergence at `m = 36`.

- The run **stopped on the relative gradient tolerance** after 49 iterations (‖g‖ = 4e-5 against a peak of about 3e4). There was no null-space blow-up and no retuning was needed.
- **The convergence history is nearly identical across resolutions:**
  - the same stiff valley around f ≈ 8.7, with |g| ≈ 1e4;
  - negative curvature in roughly the first 15 iterations;
  - fast convergence once CG converges without hitting negative curvature.
- **Time to reach** f < 1, 1e-1, 1e-2, 1e-3, 1e-4: iterations 16, 27, 32, 41, 44, at 1.2, 2.3, 3.0, 4.1 and 4.4 h.
- **Solve costs at this resolution:** a gradient takes about 110 s and an Hv about 97 s. At `m = 36` an Hv was 0.73 of a gradient; here it is 0.9, so the cost column underestimates the Hv share slightly. Wall time is the reliable number here.
- **The misfit floor is higher** (2.4e-5 against 8.7e-7 at `m = 76`). The interior parameters are recovered equally well at both resolutions (section 5), so the difference is in f, not in the model.
- **Hv at `m = 251` agrees with central finite differences of the gradient** (`scripts/hv_check_m251.sh`, job 7213722). The relative error is 2.9e-3, 2.9e-5 and 2.7e-7 for steps of 1e-3, 1e-4 and 1e-5. One Hv took 132 s and one gradient 135 s on 8 cores.

### Comparison with L-BFGS at `m = 251`

The `'lbfgs'` option (fmincon interior-point with L-BFGS, no lower bound) was run with the same harness, on 8 cores, with a budget of 400 gradient-equivalents (`scripts/lbfgs_m251.sh`, job 7215641).

Wall time to reach each misfit level:

| | f < 1 | f < 1e-1 | f < 1e-2 | f < 1e-3 | f < 2e-4 | f < 1e-4 | f < 2.4e-5 |
|---|---|---|---|---|---|---|---|
| **`newton-tr`** | 1.2 h | 2.3 h | 3.0 h | **4.1 h** | **4.1 h** | 4.4 h | 5.5 h |
| `'lbfgs'` | 1.7 h | 2.9 h | 4.4 h | 7.3 h | 12.5 h | – | – |

Final state of each run:

| | final f | stop reason | wall time | iterations | solves |
|---|---|---|---|---|---|
| **`newton-tr`** | 2.4e-5 | converged (relative gradient tolerance) | 5.9 h | 49 | 58 forward + adjoint, 148 Hv |
| `'lbfgs'` | 1.5e-4 | budget (401 gradients) | 13.3 h | 271 | 401 forward + adjoint |

Recovery of the resolved nodes 6–25:

| | relative error | max abs deviation |
|---|---|---|
| **`newton-tr`** | 2.3e-4 | 1e-5 |
| `'lbfgs'` | 1.9e-2 | 5e-4 |

**Verdict:** at `m = 251`, `newton-tr` is faster than L-BFGS at every misfit level.
- The gap is moderate early on: about 1.3–1.5× to f < 1e-2.
- It grows in the convex phase: 1.8× to f < 1e-3 and 3× to f < 2e-4.
- `newton-tr` converged in less than half the wall time L-BFGS was given, while L-BFGS was still crawling at 1.5e-4 when its budget ran out.
- The model is about 80× more accurate on the resolved nodes.
- Both leave the unresolved nodes 1–5 and 26 at the initial guess (section 5).

This confirms the `m = 36` and `m = 76` conclusions at a resolution where the forward problem is close to converged (section 7).

## 5. Null space: the edges of the fault are not constrained by the data

At convergence, some edge nodes stay at the initial guess:
- at `m = 36`: nodes 1–4 (x ≈ −10 to −8.4 km);
- at `m = 76` and `m = 251`: nodes 1–5 (x ≈ −10 to −6.8 km) and node 26 (x = 10 km), all exactly at the initial 0.013500.

These nodes lie outside the rupture. At `m = 36`, one Newton-CG step after convergence moved `a` at node 4 from 0.0138 to 0.0268 with **no change in the misfit to 5 digits**.

| node | true | converged | one step later |
|---|---|---|---|
| 1 | 0.0130 | 0.0135 | 0.0150 |
| 2 | 0.0130 | 0.0135 | 0.0164 |
| 3 | 0.0130 | 0.0135 | 0.0159 |
| 4 | 0.0130 | 0.0135 | 0.0268 |
| 5–26 | | recovered to about 1e-4 | unchanged |

At `m = 251` the resolved nodes 6–25 match the true values to about 1e-5.

- The remaining ~2% model error comes entirely from these edge nodes. The interior is recovered essentially exactly at all three resolutions.
- No optimizer can fix this. It needs regularization (Tikhonov or smoothness), a prior, or fixing those nodes.
- The model error should therefore be reported over the resolved region, or the null space handled explicitly.

## 6. Implementation (`+adjopt/trustRegionNewtonCG.m`)

**Algorithm**
- Trust-region Newton-CG (Steihaug–Toint), using the exact Hv from `adjopt.hessianVector`.
- CG stops at negative curvature (follows the direction to the boundary), at the trust-region boundary, or at the forcing tolerance η·‖g‖ with η = min(0.1, √(‖g‖/max‖g‖)).
- Preconditioner: L-BFGS from the positive-curvature CG pairs (d, Hd) of the previous accepted iteration. The trust region is measured in the preconditioner norm; the needed norms are updated by recurrences, so the preconditioner matrix is never formed.

**Cost savings**
- Trial points are evaluated with gradient and misfit together (`vectorGradient`). Most trial points are accepted, and an accepted point then needs no extra forward solve before the Hv products.
- Hv products are memoized per iterate. After a rejected step, the re-solve with a smaller radius replays the same CG iterations at no cost.

**Safeguards**
- Positivity of `a` via a fraction-to-boundary rule.
- A NaN trial point counts as a rejected step.
- Euclidean step cap (`maxStep = 0.25‖x0‖`) against huge steps along near-null directions.
- Stops when ‖g‖ ≤ 1e-8 · max‖g‖, or when the trust region collapses.

**Usage**
```matlab
optimize_antiplaneshear(maxIter, m, m_p, misfitType, optPars, val, scalingStudy, [], [], [], 'newton-tr');
```
`optimize_antiplaneshear.m` selects `AntiplaneShear2DRSFrictionOptHessian` for `'newton-tr'`, just as for `'CG'`.

### Tuning: which settings are standard and which may need revisiting

The algorithm itself is standard:
- Steihaug–Toint truncated CG (Nocedal & Wright, *Numerical Optimization*, Alg. 7.2).
- Preconditioned trust region (Conn, Gould & Toint, *Trust-Region Methods*, 2000).
- L-BFGS preconditioner from CG pairs (Morales & Nocedal, SIAM J. Optim. 10(4), 2000).

The following settings are **ad hoc choices, tuned only on this problem at `m = 36`** (the `a`-only inversion, `m_p = 26`, inverse-crime data):

| setting | value | why it was chosen |
|---|---|---|
| `etaMax` (forcing-term cap) | 0.1 | With η = 0.5, CG stopped after one iteration and the method stalled in the stiff valley |
| forcing term | min(etaMax, √(‖g‖/max‖g‖)) | Ad hoc. Relative to \|g0\| it failed, because \|g\| grows by about 1e5 after the first step |
| radius expansion | ×4 when ρ > 0.9, ×2 when ρ > 0.75 | Faster recovery after the radius had been cut |
| `maxStep` (Euclidean step cap) | 0.25‖x0‖ | Added after a step along near-null directions blew up the model error |
| `gtolRel` | 1e-8 | Relative gradient stop, added for the same reason |
| `maxCG`, `memL` | 20, 20 | Generous for 26 unknowns |
| `Delta0`, `DeltaMax` | 0.1‖x0‖, 10‖x0‖ | Rough guesses; the radius adapts anyway |

The 0.25/0.75 thresholds, the 1e-4 acceptance threshold and the fraction to boundary (0.9) are standard values.

The **grid size m should not require retuning**:
- The optimizer only sees the `m_p` parameter values.
- The length-scale settings are relative to ‖x0‖.
- The stopping test is scale-free.
- The same settings converged cleanly at `m = 76` and `m = 251` without changes (section 4b).

**Revisit the settings when:**
1. **Inverting for several parameters with different magnitudes** (e.g. `a` ≈ 0.01 together with `D_c` ≈ 0.2–1). ‖x‖, and with it `Delta0`, `DeltaMax` and `maxStep`, would be dominated by the largest parameter, so the steps in the smallest one could be far too large.
   - Prefer **parameter scaling** over retuning: optimize in x/x0 or log x, so that all components are O(1).
   - Scaling enters the gradient and Hv by the chain rule. For log x the Hessian also gets a diagonal term from the gradient.
   - Not implemented yet.
2. **A much finer parameter grid** (`m_p` ≈ 100+). `maxCG = 20` may then truncate CG too early, and `memL = 20` pairs cover a smaller part of the space.
3. **A different misfit landscape:** non-inverse-crime data, filtered or other misfit types, other initial guesses. `etaMax` and the radius rules were chosen for the stiff valley seen here. A poor choice slows convergence but should not make it wrong.

## 7. Grid resolution: `m = 36` is too coarse

Receiver recordings at the true parameters, compared with `m = 251`. All grid sizes use the same time step (`k = 0.003`, 8000 time levels).

| m | forward time (3 threads) | relative difference vs m = 251 | misfit-equivalent of the difference |
|---|---|---|---|
| 36 | 6 s | 55% | 11.2 |
| 76 | 15 s | 20% | 1.5 |
| 151 | 40 s | 11% | 0.47 |
| 251 | 89 s | – | – |

For reference, the initial misfit is 29.8 and 0.5∫v² of the data is 36.5.
- At `m = 36` the discretization error is about a third of the initial misfit, so it is effectively a different physical problem.
- In the inverse-crime setting the optimization problem is still self-consistent. However, the features that decide the optimizer ranking (non-convexity, stiffness of rupture onset) could depend on resolution. The ranking did hold at `m = 76`.
- The forward cost grows roughly linearly in m at a fixed time step. A `newton-tr` run at `m = 251` took 5.9 h on 8 cores (section 4b).

## 8. Next steps

1. ~~Run `newton-tr` at `m = 251`, check the Hv there with finite differences, and compare with L-BFGS~~: done (section 4b). `newton-tr` converged in 5.9 h; L-BFGS reached only f = 1.5e-4 in 13.3 h. Hv agrees with finite differences to 2.7e-7.
2. Check robustness to the initial guess (for example `a = 0.011` and `0.015`) at `m = 76`, comparing `newton-tr` with fmincon L-BFGS.
3. Invert for more than one parameter (for example `a` and `D_c`). Different scales and coupling are where Newton should gain most over L-BFGS.
   - Check the mixed second-derivative terms of the Hv with finite differences first.
   - Add parameter scaling (x/x0 or log x) before this, since the `newton-tr` settings are relative to ‖x‖ (section 6, *Tuning*).
4. Decide how to treat the null space (section 5): regularization, a prior, or a restricted error metric.
5. Other uses of the Hv:
   - Dominant Hessian eigenpairs at the solution (Lanczos, about 30–50 Hv) for uncertainty and resolution analysis.
   - A hybrid that runs L-BFGS early and switches to `newton-tr` once f has dropped a few orders of magnitude (expected saving about 20–30%).
   - Tests on non-inverse-crime data (`loadData`).

## Reproducing

The harness, the saved iteration histories and the run logs are in `adjoint_rate_state_inversion/optimizer_benchmark/`; see its `README.md`. To regenerate the tables in sections 3 and 4, run this from `adjoint_rate_state_inversion/`:

```matlab
run('optimizer_benchmark/summarize.m')
```
