# High-res inversion of a(x): what we learned

The goal was to recover the rate-and-state parameter `a(x)` on the rough fault with **L-BFGS** and with **Newton-CG** (line search, Hessian-vector products from second-order adjoints). We mainly wanted to show where Newton-CG helps.

- **Data:** "high-res" data come from a forward run on a finer grid (m = 1001, or m = 251 for local runs), so the inversion does not commit the inverse crime.
- **True a:** 0.009 on [−5, 6] km and 0.013 outside.
- **Speed:** t₁₀⁻⁴ is the wall time until (J − J\*)/(J₀ − J\*) = 10⁻⁴, where J\* is the best J reached by either solver.
- **Error:** the relative error in a, either on nodes 6–25 (nodal, m_p = 26 only) or as an interpolated L2 error on [−6, 9.2] km (comparable across grids).

## 1. Inverse crime (data and inversion on the same grid)
- Newton-CG converges fully (J ≈ 9e-7, error ≈ 0.02).
- L-BFGS gets stuck in a local minimum at J ≈ 1.36 (error ≈ 0.19).
- Newton-CG started from L-BFGS iteration 30 or 60 gets stuck at the same J. So Newton-CG succeeds because of the path it takes from the initial guess, not because of better local convergence.

## 2. High-res data needs regularization
With high-res data the problem is ill-posed, so we added tikhonov-grad (β/2‖Dp‖²) or TV (βΣφ_ε(Dp), with φ_ε(x) = √(x² + ε²) − ε).

At m = 76 and m_p = 26, Newton-CG reaches 10⁻⁴ faster when the objective is smooth:

| Regularization | Speedup (L-BFGS time / Newton-CG time) |
|---|---|
| tikhonov-grad, β = 10³ | 1.6 |
| tikhonov-grad, β = 10⁴ | > 4.4 |
| TV, β = 100, ε = 10⁻³ | 1.5 |
| TV, ε = 10⁻⁴ | Newton-CG is not faster, and slower at β = 1 and 100 |

Final errors are similar for both solvers, about 0.07–0.12. L-BFGS is often slightly better, even when Newton-CG reaches a lower J.

## 3. Finer grids (m_p = 51, m = 151)
- **m_p = 51:** Newton-CG is about 2× faster (> 1.9 for tikhonov-grad 10⁴, > 1.6 for TV 10 with ε = 10⁻³). The error does not improve over m_p = 26.
- **m = 151:** the interpolated errors drop to 0.08–0.10.
  - Speedups: > 2.5 (tikhonov-grad 10⁴), 1.3 (TV 1), > 2.4 (TV 10).
  - Newton-CG fails for tikhonov-grad 10³: J goes from 32.86 to 32.84 in 3.5 h.
- The speedup does not grow systematically with problem size.

## 4. Hybrid: L-BFGS for k iterations, then Newton-CG (m = 76)
- **It helps where Newton-CG is slow from the initial guess:**
  - TV 10: 0.98 h, against 1.89 h for Newton-CG alone.
  - TV 100: 1.39 h, against 6.47 h.
- **It costs time where Newton-CG is already fast:** tikhonov-grad 10⁴ takes 0.57 h (k = 20) and 0.92 h (k = 40), against 0.35 h.
- The early iterations decide most of the outcome: they set the basin, and they decide how stiff the problem is where Newton-CG starts.

## 5. TV, ε and edges
- **With ε ≥ 10⁻⁴, TV behaves like tikhonov-grad.** φ_ε is quadratic for |Da| ≪ ε, and the final |Da| is about 5e-5. Edges become ramps.
- **Sharp, piecewise constant a needs ε ≤ 10⁻⁵ or exact TV (ε = 0).** With L-BFGS at β = 100, the error drops from 0.120 to 0.084 and 0.086.
- **TV loses contrast.** The plateau levels move towards each other as β grows, and at β = 10³ the step is gone.
- **The best β is about 10–30.** Smaller β gives oscillations, larger β gives flattening.

## 6. Where Newton-CG fails
- **Large β/ε.** The initial guess is constant, so Da = 0 and the TV Hessian is β/ε there, which reaches 10⁷. The Newton steps become negligible and J stays at its initial value. This happened for β = 100 with ε = 10⁻⁵, and for β = 10³ with ε = 10⁻⁴.
- **Exact TV.** The Hessian of the regularization is zero and its gradient is discontinuous, so CG returns huge directions (‖s‖ up to 10⁸) that the line search cuts to step lengths of 1e-10. Examples:
  - β = 100: Newton-CG reached J = 4.43 in 14 h, while L-BFGS reached 1.65 in 2 h.
  - β = 10 on the m = 251 data: Newton-CG reached J = 1.52 in 3.5 h, a value L-BFGS reached in 0.45 h.
  - A hybrid started from L-BFGS iteration 20 stalls at J = 2.68.
- **The regime where TV gives the best-looking a(x) is exactly where Newton-CG fails.** There, L-BFGS is the method to use.

## 7. Reconstruction quality
- Exact TV (L-BFGS, m = 251 data) gives errors of 0.082, 0.076 and 0.069 for β = 3, 10 and 30. tikhonov-grad 10⁴ gives about the same (0.070), so TV mainly improves the shape, not the error.
- The best errors are about 0.07 in all cases. The remaining error comes from the edge positions (the right edge is placed about 1 km too far left) and from the poorly constrained outer plateaus, not from the solver.

## Bottom line
Newton-CG is the better method when the objective is smooth on the scale of the iterates: tikhonov-grad, or TV with ε ≥ 10⁻³ and moderate β/ε. There it is 1.5–4× faster to a tight tolerance, and it avoids the inverse-crime local minimum. It breaks down when the regularization is stiff or non-smooth, which is the case needed for blocky reconstructions; L-BFGS is then faster and gives better results. Warm-starting Newton-CG from L-BFGS helps in between, but not for exact TV.

## Files
- **Figures:**
  - `adjoint_rate_state_inversion/high_res_report_figures.m`, with outputs in `adjoint_rate_state_inversion/assets/` (copied to `assets/compare/` in Overleaf).
  - `adjoint_rate_state_inversion/exact_tv_report_figures.m` for `high_res_tv_beta_eps.pdf`, `high_res_exact_tv.pdf` and `high_res_table_tv.tex`.
- **Runs:**
  - `scripts/high-res_sweep_a.sh` (cluster).
  - `scripts/run_local_queue.sh` (local). It now waits for its MATLAB jobs before printing "all runs finished".
  - Hybrids are started with `x0_iter`, and their `history.t` includes the L-BFGS part.
- **Plan and step-by-step findings:** `notes/high_res_plan.md`.
