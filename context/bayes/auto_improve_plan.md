# Plan: calibrated posterior for the high-res data (autonomous run)

Goal: a posterior whose 95 % bands contain the true a(x) at most nodes, i.e. not confidently wrong.
Baseline: m = 151, m_p = 26, Gaussian smoothness prior, κ ≈ 1.89 from the discrepancy at the MAP
(truth outside the band at most nodes in −5…6 km). Background: `highres_rerun_mp126.md`.

## Ground rules
- **a_true is used only for evaluation**, never to choose hyperparameters. Any setting picked by
  looking at a_true is reported as "oracle-tuned", separately.
- **Screen cheaply, sample expensively.** Each configuration is first judged by MAP + Laplace
  (≈ 3 h on the laptop at m = 126). Only the 1–2 best get full MCMC.
- New code goes beside the old (new files or optional arguments with the current defaults), so the
  m151/mp26 results stay reproducible. Outputs carry the configuration in the filename.

## Metrics (one script, `bayes/evaluate_posterior.m`)
For a MAP + H file, or a set of chains:
- coverage: fraction of nodes where a_true lies in the 95 % band, in the middle (−5…6 km), at the
  edges, and overall. The nodes exactly on the steps are reported separately (a_true is ambiguous there);
- relative L2 error of the mean (or MAP) vs a_true, and the error at the VW/VS edges;
- κ, M(a_MAP)/M_noise and M(a_true)/M_noise;
- for chains, also R-hat and ESS (`python/diagnostic.py`).
Step 0 runs it on the existing m151/mp26 MAP and chains, which gives the baseline row and
checks that Laplace coverage ≈ MCMC coverage in the middle (justifies the screening).

## Phase 1: resolve the step (plan of `highres_rerun_mp126.md`)
1. `smoothness_prior(m_p, gamma, ...)` with γ scaled as γ = 1e-3·√(25/(m_p−1)) (same continuum
   prior). Thread `gamma` (and later a prior choice) through `run_highres_*.m`.
2. MAP + κ loop + Hessian at m = m_p = 126. Expected: κ near 1.2 and much better coverage.
   This is the new reference that the prior experiments are compared against.

## Phase 2: better priors (screened by MAP + Laplace at m = m_p = 126)
a. **Hierarchical (GSBL) prior on the increments**, port of `python/hierarchical_prior.py`:
   z_i | θ_i ~ N(0, θ_i), θ_i ~ InvGamma(β, ϑ). Suits the piecewise-constant truth: it allows a
   few large jumps without forcing oscillations elsewhere.
   - ϑ rescaled for the finer grid (ϑ ∝ Δx, i.e. 1e-8 → 2e-9 at m_p = 126), β = 0.5; one
     sensitivity variant each in β and ϑ.
   - MAP by IAS: θ-step in closed form, a-step = the existing trust-region Newton-CG with the
     Gaussian prior P(θ) (warm start). Combined with the κ loop.
   - Laplace of a | θ_MAP for screening (underestimates the spread, since θ is fixed; stated as such).
b. **Hierarchical noise level** instead of κ at the MAP: σ² ~ InvGamma, sampled in Gibbs together
   with a (σ² | a is closed form). Removes the ad-hoc κ step; the posterior over κ shows whether the
   MAP-based κ was too small. Only evaluated with MCMC (Phase 3), since the MAP gives back κ_MAP.
c. **Smoothness-prior strength**: γ × {1/2, 2} around the scaled value, plus a choice of γ by an
   honest criterion (e.g. discrepancy/L-curve, or a hyperprior on γ in Gibbs). No tuning on a_true.
d. Optional, only if a and b leave a clear model-error gap: Bayesian approximation error (mean and
   covariance of F_fine(a) − F_126(a) over prior draws, with m = 251 as "fine").

## Phase 3: MCMC for the best 1–2 configurations
- Gaussian priors: H-MALA (the sampler that mixed best), 4 chains, ~1000 solves each, same as before.
- GSBL (and noise hyperprior): Metropolis-within-Gibbs: θ | a (and σ² | a) exact draws, a | θ by one
  H-MALA step preconditioned by c·H_data(MAP) + P(θ) (126 × 126, re-factored each iteration at
  negligible cost).
- Check acceptance, R-hat, ESS; evaluate coverage; mean/band plot and P(a < b) per node.

## Compute (laptop: 12 cores, 15 GB, forward solve ≈ 20 s at m = 126)
- Screening: MAP ≈ 1–2 h + Hessian ≈ 1.5 h per configuration, 2–3 in parallel.
- Chains: ≈ 6 h per chain, 4 in parallel → ≈ 6–12 h per sampler/configuration.
- Phases 1–2 ≈ 1 day; Phase 3 ≈ 1 more day. If cluster jobs are preferred, I prepare the
  `scripts/*.sh` with `M`/`MP`/prior exports and you submit.

## Output
- Code in `bayes/`, results in `bayes/results/` (config in names), figures in `bayes/results/figures/`.
- `context/bayes/auto_improve_results.md`: table of all configurations (metrics above), what worked,
  what did not and why, and the recommended result for the report. Updated as the work progresses.
- Small commits on `bayes`, no push.
