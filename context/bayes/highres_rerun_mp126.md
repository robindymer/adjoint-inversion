# High-res rerun with m = m_p = 126: findings and handoff

## What we found (run m = 151, m_p = 26, κ ≈ 1.89)
- In the data-informed middle (x ≈ −5 … 6 km) the posterior collapses onto the MAP: the mean equals the MAP, and the 95 % bands are very narrow. The MAP and the mean oscillate around the true a = 0.009, by up to ±0.0025, and the truth lies outside the 95 % band at most nodes. So the posterior is **overconfident**, and the samplers are not at fault.
- At the poorly informed edges (x < −5, x > 8) the bands are wide and contain the truth. The mean moves away from the MAP there, so the posterior is non-Gaussian and sampling adds information.
- **Cause: model error from the a parametrization, not from the PDE grid.** From `bayes/results/discretization_error_mp*.mat` (M_model = M(a_true) against the clean m = 1001 data):

  | m | m_p | M_model / M_noise | forward solve |
  |---|---|---|---|
  | 126 | 26 | 25.2 | 22 s |
  | 151 | 26 | 26.3 (current run) | 29 s |
  | 101 | 101 | 1.7 | 19 s |
  | 126 | 126 | 1.2 | 21 s |
  | 151 | 151 | 0.8 | 28 s |

  With m_p = 26 (node spacing 0.8 km) the sharp step at x = −5 cannot be represented. Refining m does not help, but refining m_p does.
- **Why the noise inflation did not catch it.** κ comes from the discrepancy principle at the MAP: κ = M(a_MAP)/M_noise ≈ 1.9. The MAP absorbs most of the model error by oscillating, so its misfit is low, while M(a_true) ≈ 26 M_noise. The likelihood is therefore far too sharp. This weakness of choosing κ at the MAP is worth stating in the report.

## Plan: rerun with m = m_p = 126
Model error then ≈ 1.2 M_noise, the forward solve is cheaper than now, and both step locations (−5 and 6 km) fall on nodes (spacing 0.16 km).

### Steps
1. **Parametrize m and m_p in the job scripts.** `scripts/highres_{map,hpcn,mala,hmc}.sh` hardcode `m=151; m_p=26` in the `matlab -batch` line. Add `M="${M:-151}"` and `MP="${MP:-26}"` and pass `m=$M; m_p=$MP`. The `.m` scripts already take `m` and `m_p` from the workspace and name their output files `*_m<m>_mp<m_p>_*`, so the old results are not overwritten.
2. **Check the prior scaling (important).** `bayes/smoothness_prior.m` puts N(0, γ²) on each node increment, with γ = 1e-3 and no dependence on the spacing. With 125 increments instead of 25, the random-walk variance across the domain becomes 5× larger, which gives a much weaker prior. To keep the same continuum prior, scale γ ∝ √Δx: γ = 1e-3·√(25/125) ≈ 4.47e-4. Decide whether to keep the prior equivalent (recommended, so the comparison is fair), and pass γ through `run_highres_*.m`, since they all call `smoothness_prior(m_p)`.
3. **MAP + Hessian:** `sbatch --export=ALL,M=126,MP=126 highres_map.sh`.
   - The MAP is now 126-dimensional, so Newton-CG needs more iterations. Keep `maxIter=150` and use `RESUME=true` if the job hits the time limit.
   - The Hessian is now `post.Uhv(aMap, eye(126))`, which takes 126 Hessian-vector products (about 5× the current cost).
   - Record the new κ. Also report M(a_true)/M_noise for m_p = 126 (setup prints it) next to κ, to show whether the discrepancy-at-MAP gap closed.
4. **Chains:** submit with the same exports and `--dependency=afterok:<map job>`, or extend `highres_submit.sh` to forward `M`/`MP`.
   - Watch acceptance after adaptation (targets: MALA 0.574, HMC 0.65; H-pCN with `s` adaptive).
   - The higher dimension may need smaller steps. Check that H is positive definite.
5. **Evaluation:** `python bayes/python/diagnostic.py` on the new chains, then the mean/band plot as before.
   - Success criterion: the true a (sampled at the 126 nodes) lies inside the 95 % band at most nodes in the middle region, i.e. the bands are calibrated.
   - Compare with the m_p = 26 run in the report.

## Context
- Setup: `bayes/setup_highres_problem.m` (noisy data from `mat/fractal_fault_m1001_mp1001`, eta = 0.02, κ applied as σ² = κσ0²).
- MAP and κ loop: `bayes/run_highres_map.m`. Earlier todo and history: `context/bayes/highres_todo.md`.
- Current results: `bayes/results/*_highres_m151_mp26_*.mat`, `highres_map_m151_mp26_eta0.02.mat`.
- Fallback if there is no time to rerun: report the m_p = 26 result as a finding (overconfident posterior due to an unresolved step, hidden by setting κ at the MAP), with the m_p table above as evidence.
