# High-res TODO

Posterior for the noisy m = 1001 data, solved on m = 151, m_p = 26, with H-pCN, H-MALA and H-HMC.

## Done
- [x] Code, data and scripts on the cluster; smoke test passed.
- [x] First full run (jobs 7969554–7). The MAP hit the 60-iteration cap in every round, and a bug saved inflation 2.933 while `aMap` and H belonged to 7.877. Result: H-pCN frozen (acceptance about 0%), while H-MALA and H-HMC mixed (acceptance 0.4–0.8). Keep these runs as a robustness experiment (inaccurate Laplace approximation) and state the cause. They are not final results.
- [x] Fixed `run_highres_map.m`: it saves the inflation that `aMap` and H belong to; `maxIter = 150`, `maxRounds = 8`. Added a `MAXITER` override to `highres_map.sh`.

## Now
- [ ] Resumed MAP job 7987365 (16 h). The inflation went 2.933 → 2.188 → 1.889, with the error vs. a_true at 8.3%. Round 3 is nearly flat (M ≈ 0.739), so κ ≈ 1.88 should meet the 5% check, and the job then computes H.
  - If it times out while computing H: `sbatch -t 04:00:00 --export=ALL,RESUME=true,MAXITER=10 highres_map.sh`
  - When done, check the log: inflation converged, and "Eigenvalues of the prior-preconditioned data Hessian" printed, with no "not positive definite" error.
- [ ] Rerun the chains on the final MAP: `./highres_submit.sh --skip-map`. Move the old `*_highres_*seed*.mat` files aside first (e.g. into `results/highres_run1/`), since the new chains overwrite them.
- [ ] Check the chain logs.
  - H-pCN: if `std(phi)` is much larger than 1, rerun with `sbatch --export=ALL,S=0.5 highres_hpcn.sh`.
  - H-MALA and H-HMC: check the acceptance rate after adaptation. The adaptation targets are 0.574 for MALA and 0.65 for HMC.

The budget is about 1000 PDE solves per chain for each sampler (H-pCN 1000 iterations, H-MALA 500, H-HMC 100 × 5 leapfrog steps). Change it with `NITER` (and `NLEAP`) through `--export`.

## Inverse-crime Laplace comparison (independent of the above)
- [ ] `sbatch laplace_comparison.sh`: H-MALA and H-HMC at m = 36, 4 seeds each, matching the existing H-pCN chains. This completes the prior- vs. Laplace-preconditioning table.
- [ ] Fetch the results: `rsync -rvP --include='*_laplace_m36_*seed[1-4].mat' --exclude='*' <cluster>:.../bayes/results/ bayes/results/`

## Back home
- [ ] Copy back `bayes/results/{highres_map,hpcn_highres,mala_highres,hmc_highres}_*.mat`.
- [ ] Plot the MAP with `direct_effect_inverse_crime_plot.m` (`resultsDir = fullfile('..','bayes','results')`, `resultFiles = {'highres_map_m151_mp26_eta0.02.mat'}`).
- [ ] Run `python diagnostic.py` on each set of chains. Compare R-hat and ESS per PDE solve across the three samplers.
- [ ] Figures and numbers for the report:
  - the mean with 95% bands
  - $P(a(x) < b)$ per node
  - the VW edge positions
  - Laplace vs. MCMC
  - the inflation factor κ (fill in the exact value below)
  - the choice of m = 151, from `discretization_error_mp*.mat`. With m_p = 26 the model error is about 26 × M_noise at a_true, about 1.9 × at the MAP.
  - a side-by-side with the inverse-crime results.

<!-- The high-resolution data contain, besides the measurement noise, a model error: the forward
model is solved on a coarser grid ($m = 151$) and represents $a$ on $m_p = 26$ nodes, which
cannot resolve the sharp step. Following the simple approach of the course literature, we
absorb it by inflating the noise variance, $\sigma^2 = \kappa\sigma_0^2$. The factor $\kappa$
is chosen by the discrepancy principle: the misfit at the MAP should equal the noise level,
$M(\vect a_{\mathrm{MAP}}) \approx M_{\mathrm{noise}}$. Since the MAP itself depends on
$\kappa$ through the relative weight of the prior, we alternate between computing
$\vect a_{\mathrm{MAP}}$ for fixed $\kappa$ (trust-region Newton-CG) and updating
$\kappa \leftarrow M(\vect a_{\mathrm{MAP}})/M_{\mathrm{noise}}$, until $\kappa$ changes by
less than $5\,\%$. This gives $\kappa \approx 1.9$, i.e.\ a model error of the same order as
the noise. The Laplace approximation is then computed at the final $\vect a_{\mathrm{MAP}}$. -->
