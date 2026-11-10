# Plan: Bayesian inference of frictional parameters from seismograms

## Context

The goal is to carry out the proposal in `context/plan.md` and the report outlined in `context/Bayesian_adjoint_inversion.pdf`. The work has two halves:

1. **Deterministic:** TV-regularized L-BFGS inversion for a(x), with λ chosen by the discrepancy principle. An earlier report (`context/project_report.pdf`) concluded that λ had to be tuned by hand because the misfit had no noise scale. This project fixes that.
2. **Bayesian:** pCN vs MALA vs HMC under a Gaussian smoothness prior, compared by ESS per PDE solve. Then the Laplace approximation and the hierarchical (GSBL) prior.

Decisions already taken:
- Samplers are written in **MATLAB**. Python is used only for priors (reference implementations), diagnostics and figures.
- MCMC runs on the **m = 36** forward grid, with data from a finer grid.
- Time horizon is about **6 weeks**, to the December discussion. Laplace and the hierarchical prior are therefore in scope.

What exists already (in `code/adjoint-inversion/adjoint_rate_state_inversion/`, its own git repo):
- `+adjopt/regularization.m` already provides `tv`, `tikhonov-grad` (H1), `tikhonov` and `lasso`, each with value, gradient and Hessian-vector product. TV uses β·Σ√((Dp)²+ε²)−ε with default ε = 1e-4.
- `high_res_data_inversion.m` (m = 251, m_p = 26, 88 receivers) together with `scripts/high-res_reg_a.sh` runs regularized inversions on Pelle. A TV β-sweep with L-BFGS has been run, but every run stopped at maxIter = 200. The model-error misfit floor is about 0.22.
- `optimizer_benchmark/InvProblem.m` wraps the solver as `f(x)`, `fg(x)` and `Hv(x,V)`. It counts solves, caches results and returns `inf` for x ≤ 0. Samplers will build on it.
- Misfit: M(a) = ½ Σ_r ∫ (v_r − d_r)² dt, with **no noise and no 1/σ²**. High-res data are stored as spline function handles in `mat/fractal_fault_m1001_mp1001/receiverData.mat`, a file that exists **only on the cluster**.
- Cost per solve:

  | m | forward | gradient | Hv |
  |---|---|---|---|
  | 36 | 4.3 s | 5.6 s | 4.1 s |
  | 251 | – | ~120 s | ~110 s |

- The forward solve can return NaN; such points must be rejected.
- m = 36 has about 55 % discretization error relative to m = 251, and m = 76 about 20 %.
- Your Python priors in `code/` (`fault.py`, `smoothness_prior.py`, `hierarchical_prior.py`) and the course code (`samplers.hmc`, `ess`, `rank_rhat`, `approximations.lowrank_*`) are the reference implementations.

New code goes in `code/bayes/` (MATLAB) and `code/*.py` (Python). Changes inside `code/adjoint-inversion/` are kept minimal and committed in that repo.

---

## Phase 0 — Noisy data, likelihood scale, model-error budget (week 1)

**0.1 Noisy data.** New script `code/bayes/make_noisy_data.m`:
- Load the m = 1001 recordings and sample them at t_n with Δt_d = 0.01 s (600 samples per receiver).
- Add iid noise N(0, σ²) with σ = η·max|v|. Start with η = 2 % and a fixed seed.
- Save the clean and noisy samples as **plain arrays** in `.mat`. This keeps them readable from Python.
- Also save the noisy samples as spline handles in the existing `receiverData` format, so `optimize_general` loads them unchanged.
- Run once on Pelle, where the data file lives.

**0.2 Likelihood scale.** Take −log L(a) = M(a)/(σ²Δt_d).
- Calibrate without any PDE solve: evaluate M on the pure noise interpolant (d_noisy − d_clean).
- Check that 2M/(σ²Δt_d) ≈ N = N_rec·N_t.
- If the spline quadrature biases this, use the empirical factor instead.

**0.3 Model-error budget.**
- One forward solve at a_true on each of m = 36, 76 and 251, against the clean high-res data. This gives M_model(m).
- Compare M_model(m) with the noise level M_noise ≈ N σ² Δt_d / 2.
- **Decision gate for the Bayesian grid:** if M_model(36) ≫ M_noise, inflate the noise to σ_tot² = σ² + σ_model² (a scalar, so the misfit code is unchanged). If that makes the data uninformative, sample on m = 51 or 76 instead.
- Record the result; it goes into the "Assumptions" section of the report.

## Phase 1 — TV + discrepancy principle (weeks 1–2, Pelle)

**Changes in `adjoint-inversion`:**
- `optimize_general.m` / `high_res_data_inversion.m`:
  - Add an option to scale the misfit by 1/(σ²Δt_d), so that β is relative to the likelihood.
  - Add an option to load the noisy data file.
  - Log M(a_true) through one extra forward solve. This is the target discrepancy δ² and includes model error.
- `regularization.m`: run with `reg_eps = 1e-6`. The default 1e-4 is about 1/40 of a jump and smooths the small increments. Do a brief ε-sensitivity check at the chosen β.

**Experiments:**
1. Verify on inverse-crime data: m = 76 locally, data with the same noise model. Runs take 1–3 h each.
   - FD gradient test of the scaled TV objective with the existing `adjopt.compareWithFD` / `gradientTestFriction.m`.
   - TV with β by the discrepancy principle (τ ≈ 1.1). The truth's discrepancy is ≈ M_noise here.
2. High-res noisy data at m = 251:
   - Run an initial log-spaced β bracket as parallel Slurm jobs.
   - Refine with 2–3 secant/bisection steps in log β on M(a_β) = τ·M(a_true).
   - Use continuation (warm start from the neighbouring β) to stay inside the 48 h limit, and resume with `loadOptPath`.
   - Report any run that stops at maxIter as unconverged.
3. Comparison at each method's own β_DP: TV vs H1 (`tikhonov-grad`) vs unregularized. Report:
   - ‖a − a_true‖;
   - positions of the VW edges (where a crosses b);
   - behaviour of the unconstrained nodes 1–5 and 26 (the null space found in `notes/optimizer_comparison.md`);
   - the L-curve or discrepancy-vs-β plot.

## Phase 2 — Bayesian, smoothness prior (weeks 2–4)

Posterior: U(a) = M(a)/(σ²Δt_d) + ½(a−μ)ᵀP(a−μ). P comes from `smoothness_prior.py` (`gaussian_precision`, "differences", with GAMMA and GAMMA0 fixed by the existing prior predictive checks). The MAP of this posterior is H1-Tikhonov with β = 1/γ² plus a level term; the report should state this link to Phase 1.

**New MATLAB code in `code/bayes/`:**
- `prior_smoothness.m`: P, R = chol(P), sampling, log-density and whitening a = μ + R⁻ᵀw.
  - Check it against `smoothness_prior.py` by comparing log-prior values on a few exported vectors.
- `BayesProblem.m`: wraps `InvProblem` to give U, ∇U and Hv_U.
  - Returns `inf` on NaN or a ≤ 0, so those points are rejected.
  - Keeps forward, adjoint and Hv counters.
- `gaussian` regularization type in `+adjopt/regularization.m`: R = ½(p−μ)ᵀP(p−μ), with gradient P(p−μ) and Hv P·V. The existing optimizers then compute the MAP.
- Samplers, all sharing these features:
  - burn-in-only step adaptation (Robbins–Monro toward acceptance 0.25 for pCN, 0.57 for MALA, 0.7 for HMC), then fixed steps;
  - a `.mat` checkpoint every N iterations, with resume;
  - solve counters stored with the chain.

  The three samplers:
  - `pcn.m`: in whitened coordinates, w′ = √(1−s²)w + sξ. Needs only the likelihood, so one **forward** solve per step.
  - `mala.m`: preconditioned with C. One gradient per step.
  - `hmc.m`: a port of course `samplers.hmc`, with jitter and mass matrix M = C⁻¹. L gradients per step.
- `run_chain.m` and `code/bayes/scripts/chain.sh <sampler> <seed>` for Slurm on Pelle.

**MAP and Laplace (cheap at m_p = 26):**
- Compute the MAP with the optimizer and the `gaussian` regularization.
- Compute the full misfit Hessian from 26 Hv products, about 2 min at m = 36.
- Export to Python. In new `code/laplace.py`, take the eigenvalues of the prior-preconditioned misfit Hessian RᵀH_misR, following the formulas of `approximations.lowrank_laplace`, `lowrank_std` and `lowrank_sample`, but starting from H instead of J.
  - Show the eigenvalue decay and the rank r with λ > 1.
  - Use H⁻¹ as the preconditioner C for MALA and as the HMC mass matrix.

**Verification before spending PDE solves:**
- Run all three samplers on the Gaussian surrogate U(a) = ½(a−a_MAP)ᵀH(a−a_MAP), which costs no solves.
- Check that they recover the mean and covariance and that ESS is plausible.
- Then do 300–500-step pilot runs on the real posterior to tune step sizes and L.

**Production on Pelle:** 4 chains per sampler.
- Starts: the MAP plus 3 over-dispersed draws from the Laplace approximation. Avoid the constant 0.0135, where the misfit is concave.
- Each chain is a 48 h job and is resumed as needed. One job gives about 3·10⁴ pCN steps, 2.5·10⁴ MALA steps or 2–3·10³ HMC steps (L ≈ 10).
- Measure the m = 36 timings on Pelle first.

**Diagnostics in new `code/diagnostics.py`** (reads the exported chains and uses course `samplers.ess` and `rank_rhat`). Everything in Remark 4.21 of the script:
- rank-R̂ and ESS for all 26 components and the QoIs;
- trace plots of the slowest quantities;
- cost, in **ESS per gradient-equivalent**: (4.3·nFwd + 1.3·nAdj + 4.1·nHv)/5.6, as in `notes/optimizer_comparison.md`;
- the sampler comparison table.

**Quantities of interest (QoIs):**
- posterior mean and 95 % credible bands, compared with a_true and with the Laplace bands;
- P(a(x_i) < b) per node, i.e. the probability that the fault is velocity-weakening there;
- posterior distributions of the VW edge positions;
- coverage of the truth.

## Phase 3 — Hierarchical (GSBL) prior (week 5)

- Port `hierarchical_prior.py` to `code/bayes/prior_gsbl.m`.
  - θ_i ~ InvGamma(β, ϑ) on the increments.
  - Inverse CDF via `θ = ϑ ./ gammaincinv(p, β, 'upper')`, which is base MATLAB and needs no toolbox.
  - Check against Python as in Phase 2.
- **MAP via IAS:**
  - The a-step is the optimizer with a new `weighted-grad` regularization, ½Σ zᵢ²/θᵢ.
  - The θ-step is the closed form `ias_theta_update`.
  - Compare the result with the TV-MAP from Phase 1, using the course's link between IAS and TV.
- **Sampling:** pCN in the normalized coordinates (w, τ), using `to_normalized`/`from_normalized` from your Python file. It is prior-reversible and avoids the funnel between a and θ.
- **Optional:** Metropolis-within-Gibbs, with an exact draw θ|a ~ InvGamma(β+½, ϑ + zᵢ²/2) and a pCN step for a|θ.
- Compare with the smoothness posterior on the QoIs, in particular the sharpness of the edge positions.

## Phase 4 — Report (week 6)

Map the figures to the sections of `context/Bayesian_adjoint_inversion.pdf`:

| Section | Content |
|---|---|
| §2 Problem formulation | ill-posedness, from the Hessian eigenvalues and the null-space nodes; noise model; model-error budget from Phase 0 |
| §3 Computational method | algorithms; verification on the Gaussian surrogate |
| §4 Results | DP/TV figures; sampler comparison table; credible bands; VW probability; Laplace vs MCMC; GSBL |

Every figure gets a script under `code/` so it can be reproduced.

---

## Critical files

- Modify, in `code/adjoint-inversion/adjoint_rate_state_inversion/`:
  - `+adjopt/regularization.m` (types `gaussian` and `weighted-grad`);
  - `+adjopt/optimize_general.m` (misfit scaling, noisy data);
  - `high_res_data_inversion.m`;
  - `scripts/high-res_reg_a.sh`.
- Reuse: `optimizer_benchmark/InvProblem.m`, `optimizer_benchmark/setup_problem.m`, `+adjopt/vectorGradient.m`, `+adjopt/hessianVector.m`, `+adjopt/compareWithFD.m`.
- New: `code/bayes/{make_noisy_data,prior_smoothness,prior_gsbl,BayesProblem,pcn,mala,hmc,run_chain}.m`, `code/bayes/scripts/chain.sh`, `code/laplace.py`, `code/diagnostics.py`.
- Python references: `code/smoothness_prior.py`, `code/hierarchical_prior.py`, `code/fault.py`, and from the course code `samplers.py` (`ess`, `rank_rhat`, `hmc`) and `approximations.py` (`lowrank_*`).

## Verification

- **Noise calibration:** 2M(noise)/(σ²Δt_d) ≈ N_rec·N_t (Phase 0.2).
- **Derivatives:** FD gradient test and Hv test of the scaled objective with TV, gaussian and weighted-grad, using `compareWithFD` and `hv_check.m`.
- **Priors:** MATLAB and Python give the same log-prior to 1e-12, and the MATLAB prior samples match the Python prior predictive statistics (`fault.summarize`).
- **Discrepancy principle:** inverse-crime runs recover a_true with M(a_β) ≈ τ·M_noise, and the discrepancy is monotone in β on the bracket.
- **Samplers:** recover the mean and covariance of the Gaussian surrogate. On the real posterior, rank-R̂ < 1.01 and ESS > 100 for every reported quantity.
- **Laplace:** Laplace and MCMC bands agree where the posterior is near-Gaussian. Disagreements are reported, not hidden.
