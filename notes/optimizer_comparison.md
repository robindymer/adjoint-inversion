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
- Grid sizes tested: `m = 36`, with a confirmation run at `m = 76`.

**Recommendation:** use the new trust-region Newton-CG option, `opt_routine_str = 'newton-tr'` (`+adjopt/trustRegionNewtonCG.m`). It was the only method that converged to the misfit floor at both `m = 36` and `m = 76`.

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
| fmincon interior-point + `'cg'` + Hv (current `'CG'`) | 1.3e-3 | 4.03 | 120 | 143 | 147 | 171 |
| fmincon trust-region-reflective + Hv | 0.34 | 1.06 | 232 | – | – | – |
| fmincon interior-point L-BFGS (current `'lbfgs'`) | 2.7 | 0.20 | – | – | – | – |
| fminunc quasi-Newton (BFGS) | 2.7 | 0.17 | – | – | – | – |

`err` is the relative L2 error in `a` on the parameter grid; see section 5 for why it bottoms out around 0.02.

\* This error is from the variant that switches to L-BFGS at the end. The plain preconditioned run reached the same f but then took a large step in the null space (err 0.26), which led to the fixes in section 6.

**Observations**
- **fmincon interior-point L-BFGS is very path-dependent.** In this run it never got below f = 2.7. An earlier saved run with identical settings (`mat/inversion/rough_lbfgs_maxiter200_m36_..._1.mat`) reached f ≈ 4e-4 after about 300 gradients.
- **The current `'CG'` option converges in misfit but drifts far off in `a`** (err 4.0), in directions the data cannot see.
- **The L-BFGS preconditioner built from CG pairs (Morales & Nocedal) helps most in the non-convex early phase.** The cost to reach f < 1 drops from 73 to 41.
- **Switching from Newton-CG to L-BFGS once the Hessian looks positive definite gave no speed-up.** The switch only triggered after f had already converged.

## 4. Confirmation at `m = 76`

The three most relevant methods were run in parallel, 3 threads each. Cost weights are those measured at `m = 36`; wall time is also given.

| method | final f | final err | cost to f < 1 / 1e-1 / 1e-2 / 1e-3 / 1e-5 | minutes to the same |
|---|---|---|---|---|
| **`newton-tr` (repo version)** | **8.7e-7** | **0.022** | 49 / 92 / 117 / 122 / 144 | 22 / 42 / 54 / 56 / 66 |
| fmincon interior-point L-BFGS (current `'lbfgs'`) | 1.0e-3 (budget hit) | 0.081 | 39 / 61 / 93 / – / – | 18 / 29 / 44 / – / – |
| L-BFGS (own implementation) | 1.36 (stalled) | 0.126 | – | – |

`newton-tr` stopped on its own (relative gradient tolerance) at cost 172, about 79 minutes.

The pattern holds at both resolutions:
- **Early, non-convex phase:** L-BFGS is comparable or slightly cheaper down to f ≈ 1e-2.
- **Late phase, once the Hessian is positive definite:** Newton-CG converges much faster and reaches the misfit floor, while L-BFGS crawls.

## 5. Null space: the left edge of the fault is not constrained by the data

At convergence, parameter nodes 1–4 (x ≈ −10 to −8.4 km, outside the rupture) stay at the initial guess. One Newton-CG step moved `a` there from 0.0138 to 0.0268 with **no change in the misfit to 5 digits**.

| node | true | converged | one step later |
|---|---|---|---|
| 1 | 0.0130 | 0.0135 | 0.0150 |
| 2 | 0.0130 | 0.0135 | 0.0164 |
| 3 | 0.0130 | 0.0135 | 0.0159 |
| 4 | 0.0130 | 0.0135 | 0.0268 |
| 5–26 | | recovered to about 1e-4 | unchanged |

- The remaining ~2% model error comes entirely from these nodes. The interior is recovered essentially exactly.
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
- The same settings converged cleanly at `m = 76` without changes. That is one data point; the `m = 251` run would confirm it.

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
- The forward cost grows roughly linearly in m at a fixed time step. A `newton-tr` run at `m = 251` is estimated at 5–7 h.

## 8. Next steps

1. Run `newton-tr` once at `m = 251`, after a quick Hv finite-difference check at that resolution.
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
