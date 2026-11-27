# Results: calibrating the high-res posterior (autonomous run, started 2026-11-26)

Plan: `auto_improve_plan.md`. Evaluation: `bayes/python/evaluate_posterior.py MAPFILE [CHAINGLOB ...]`
(95 % coverage of a_true per region; middle = −5 < x < 6 km, edges = the rest, nodes on a step excluded).

## Summary table

| run | m | m_p | prior | κ | posterior | middle covered | edges covered | median \|err\|/std (middle) | max | rel. err |
|---|---|---|---|---|---|---|---|---|---|---|
| baseline | 151 | 26 | Gaussian, γ = 1e-3 | 1.89 | Laplace | 2/13 | 10/12 | 11.7 | 40.4 | 9.1 % |
| baseline | 151 | 26 | Gaussian | 1.89 | H-MALA (4 × 450) | 2/13 | 5/12 | 12.3 | 41.9 | 15.3 % |
| baseline | 151 | 26 | Gaussian | 1.89 | H-HMC (4 × 80) | 2/13 | 5/12 | 10.1 | 38.0 | 14.9 % |
| baseline | 151 | 26 | Gaussian | 1.89 | H-pCN s = 0.5 | 2/13 | 6/12 | 13.5 | 43.3 | 16.0 % |
| phase 1 | 126 | 126 | Gaussian, γ = 4.47e-4 (same continuum prior) | 2.0 | Laplace | 34/68 | 51/57 | 1.9 | 29.7 | 7.8 % |
| phase 2 | 201 | 201 | Gaussian, γ = 3.54e-4 (same continuum prior) | 1.24 | Laplace (low rank, r = 120) | **75/109** | 78/90 | **1.1** | 30.1 | 8.9 % |
| phase 2 | 201 | 201 | Gaussian | 1.24 | H-MALA (3 × 300 + 1 × 170, 870 kept) | **83/109** | 84/90 | **1.1** | 15.9 | 7.0 % |
| phase 3 | 201 | 201 | GSBL (β = 0.5, ϑ = 1.25e-9), θ fixed at θ_MAP | 1.24 | Laplace of a \| θ (low rank) | 38/109 | 38/90 | 2.6 | 13.9 | 19.1 % |

## Findings

1. **Laplace screening is valid.** For the baseline, the Laplace approximation and all three samplers
   agree in the middle: 2/13 covered, errors of about 12 std. In the middle the posterior is Gaussian.
2. **m = m_p = 126 helps a lot but does not fully calibrate.** The median error in the middle falls
   from 11.7 to 1.9 std, and coverage rises from 15 % to 50 %. The remaining failures are:
   - a Gibbs-type overshoot just before the step at x = 6 (a = 0.0149 at x = 5.84; 30 std). The
     Gaussian smoothness prior penalizes the jump, so the MAP rings around it;
   - oscillations of ±1.3e-3 in x ≈ 0–1.5 (5–8 std), from the remaining structured model error.
   M(a_true)/M_noise = 2.23 against the noisy data (≈ 1 noise + 1.2 model error), and κ = 2.0.
   Note: x = −5 is *not* a node at m_p = 126 (spacing 0.16 km); x = 6 is.
3. **A hierarchical noise level would not help.** With N ≈ 5e4 data, σ² | a is so concentrated that
   κ ≈ M(a)/M_noise for every sample. That is the discrepancy choice again, so it cannot widen the
   posterior. The model error is structured (correlated in time and receivers), not extra iid noise.
4. **γ sensitivity runs (γ × 2, × ½) stopped:** in the first iterations they tracked γ × 1 exactly,
   because the data dominate. They are not worth the CPU before the base case is resolved.

5. **Hierarchical (GSBL) prior: no gain at m_p = 126.** After 3 IAS outer iterations (β = 0.5,
   ϑ = 2e-9), its MAP is nearly identical to the Gaussian one: the same middle oscillations, and the spike at
   x = 5.84 grows (0.0160). Rel. err 8.4 % against 7.8 %. So the remaining bias comes from the data and
   model, not from the prior. Stopped. (`bayes/gsbl_map.m`; the Gibbs sampler `bayes/gsbl_gibbs.m` is
   written and tested but not run.)
6. **Where the model error comes from** (M(a_true) against the *clean* m = 1001 data, / M_noise;
   `discretization_error.m`):

   | m | m_p | M_model/M_noise | forward solve (4 threads, desktop) |
   |---|---|---|---|
   | 126 | 26 | 25.2 | |
   | 251 | 126 | **13.6** | 119 s |
   | 101 | 101 | 1.7 | |
   | 126 | 126 | 1.2 | |
   | 151 | 151 | 0.8 | |
   | 201 | 201 | **0.36** | 81 s |
   | 251 | 251 | **0.33** | 119 s |

   Refining only the PDE grid makes it *worse* (m = 251, m_p = 126: 13.6). For m_p < m the code maps
   a to the grid with high-order SBP glue interpolation (`AntiplaneShear2DRSFrictionOpt`, `GlueOpsOP`),
   which rings (Gibbs) at a step. So for a step-like a, **m = m_p (no interpolation) is essential**,
   and the model error then falls below the noise level from m = m_p = 201 on. This also explains
   the baseline: with m_p = 26 the interpolated step is badly represented.

7. **m = m_p = 201: close to calibrated away from x ≈ 3.** Middle coverage is 69 % (50 % at 126), and the
   median error is 1.1 std. κ = 1.24 is the discrepancy fixed point; M(a_true)/M_noise = 1.38. The
   failures are now local:
   - an oscillation at x = 2.6–3.4 (±5e-4, up to 30 std);
   - the smoothed step next to x = −5 and x = 6 (e.g. a = 0.0126 at x = −4.9, 16 std). The Gaussian
     increment prior cannot make a sharp jump.
   The rel. err (8.9 %) is dominated by the step nodes, so it is a poor metric here.
   (`results/highres_map_m201_mp201_eta0.02.mat`)

8. **TV prior (MAP only): no better overall.** `tv_prior.m` uses √(Δa² + ε²)/b with ε = 1e-5 and
   b = γ√(2/π), the same mean |increment| as the Gaussian prior. Warm-started from the Gaussian MAP and
   stopped after 35 TR-Newton iterations (about 8.5 h), when U fell by about 1 per 15 min.
   - **Middle:** RMS error 5.16e-4 against 5.35e-4, and the same x ≈ 3 oscillation.
   - **Steps:** marginally steeper, with a larger overshoot at x = −5.1 (0.0162 against 0.0157).
   - **Unconstrained left edge (x < −5):** becomes flat at about 0.0155, against the Gaussian ramp
     toward 0.0134; the truth is 0.013. Rel. err 13.1 % against 9.3 %.
   TV charges a monotone ramp the same as a jump, so it cannot sharpen a smeared step. The smearing
   and the overshoot come from the data and model. Figure:
   `adjoint_rate_state_inversion/assets/bayes/map_m201_gauss_vs_tv.png`. The file is
   `results/highres_map_m201_mp201_eta0.02_tv.mat`; its MAP was stopped by hand, so infHist is empty.

9. **H-MALA at m = m_p = 201 confirms and slightly improves on Laplace.** The chains cover 76 % of the
   middle nodes (Laplace 69 %, baseline 15 %) and 93 % of the edge nodes. The worst error halves to
   16 std, and the posterior-mean rel. err is 7.0 %. The bands have the same width as Laplace in the
   middle (median std 1.8e-4). At the left edge they are narrower than Laplace: the posterior is
   non-Gaussian there, since only a > 0 and the nonlinear forward model constrain it.
   - **Mixing:** U and most nodes mix well (R-hat(U) = 1.16 over 4 chains). In the oscillation patch at
     x = 2.4–3.3 the chains drift slowly at separate levels, about 1e-4 apart; R-hat there reaches 2.1, and
     the global minimum ESS is 21. The bands in that patch are less reliable, but they miss the truth by
     much more than the chain spread anyway.
   - **Remaining misses:** the nodes next to the steps (x = −4.9: a ≈ 0.0125; x = 5.9: a ≈ 0.0117;
     truth 0.009), where the posterior is confident in a smeared step, and the x ≈ 3 oscillation.
   - Figures: `assets/bayes/posterior_bands_baseline_vs_m201.png` (bands, against the baseline) and
     `assets/bayes/mala_traces_m201.png`. Script: `bayes/python/plot_posterior_bands.py`.

## Numerical notes (MAP at m_p = 126)
- Newton-CG in the original coordinates is useless at m_p = 126 (badly scaled prior precision): 115
  iterations in 4 h only reached M = 1.42. In the whitened prior coordinates (`whiten = true`) the
  first step already reaches M = 1.05.
- The landscape is strongly nonlinear: even with exact Newton steps (CG to 1e-3, ρ ≈ 1), U falls by
  only tens per step near the MAP. Stopped once the decrease per step was ≪ m_p/2 (the spread
  of U over the posterior), so that a_MAP is within the typical set. κ = 2 was set directly from
  M/M_noise = 1.97 of the κ = 1 round, which saves one round.
- Hessian (126 Hv): 3 processes × 2 threads × 42 columns, 42 min (`hessian_columns.m`,
  `combine_hessian.m`). The H-v product takes ≈ 60 s at 2 threads; forward + adjoint ≈ 45 s at 6 threads.

- H-MALA at m_p = 126 (2 chains) was stopped after 70 iterations to free the CPU for m = 201. With
  20 samples per chain after adaptation the result is inconclusive (R-hat 2.5).

- m = 201 MAP: warm start from the m126 MAP, κ0 = 1.4, whitened TR-Newton. It was stopped at
  M/M_noise = 1.244 (last step ΔU ≈ 1.2), and κ was set to that fixed point.
- m = 201 Hessian: low-rank randomized (`lowrank_hessian.m`, r = 120, 3 × 40 Hv at 2 threads, 97 min).
  60 eigenvalues are above 1, and the smallest kept is 0.024, so the rank is sufficient.

10. **Hierarchical (GSBL) prior at m = m_p = 201: best MAP in the middle, but no usable UQ without
    sampling θ.** After 8 IAS outer iterations (about 9 h):
    - **MAP:** the middle RMS error is 4.13e-4, against 5.35e-4 for the Gaussian MAP and 5.16e-4 for TV. The
      large θ sit exactly at the true jumps (x = −4.95, 5.75–5.95), and elsewhere √θ ≈ 3e-5, ten times below the Gaussian γ.
    - **No convergence:** IAS never met |Δa|/|a| < 1e-3. The poorly constrained left-edge level kept
      moving (0.015–0.017; the last update was 1.7 %) and ends at about 0.017 against the true 0.013.
    - **Laplace of a | θ_MAP is badly overconfident:** the median middle std is 3.4e-5, against 1.8e-4
      for the Gaussian. Coverage is 38/109 in the middle and 38/90 at the edges, and the left edge is
      0.0170 ± 3.6e-4 (12 std off). Only 17 eigenvalues of the prior-preconditioned data Hessian are
      above 1, so a | θ is prior-dominated. Fixing θ removes exactly the uncertainty that matters.
    - **Proper UQ:** it would need the Gibbs sampler (`gsbl_gibbs.m`). At about 3 min per iteration
      and with slow joint (a, θ) mixing, that is days of CPU on this machine. Not run.
    - **Files:** `results/highres_map_m201_mp201_eta0.02_gsbl.mat` (aMap, theta, hist, low-rank H).

## Recommendation
- **For uncertainty quantification:** m = m_p = 201, the Gaussian smoothness prior with γ ∝ √Δx (the
  same continuum prior as the baseline), κ by the discrepancy principle (1.24), and H-MALA
  preconditioned by the low-rank Laplace Hessian. This is the best-calibrated posterior: 76 % of the
  middle nodes and 93 % of the edge nodes are covered, against 15 % and 50 % for the baseline.
  Laplace alone is a good cheap substitute (69 %).
- **For a point estimate:** the GSBL MAP is best in the middle but unreliable at the poorly
  constrained left edge. The Gaussian MAP is the safer overall estimate (rel. err 9.3 %).
- **What remains:** the misses come from structured model error. M(a_true)/M_noise is 1.38, against
  1.24 at the MAP: the data prefer the MAP to the truth by about 2 800 in U, which no prior of
  reasonable strength can overcome. The next step would be the Bayesian approximation error
  approach: estimate the mean and covariance of the m = 1001 → 201 model error from prior draws and
  add them to the likelihood.
- **Dead ends:** m_p < m (the SBP interpolation of a rings at the steps), γ × 2 / × ½ (tracked γ × 1),
  a hierarchical noise level (it reduces to the discrepancy principle), TV (it charges a ramp the same as
  a jump), and GSBL with θ fixed (overconfident).
