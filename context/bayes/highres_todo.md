# High-res TODO

Posterior for the noisy m = 1001 data, solved on m = 151, m_p = 26, with H-pCN, H-MALA and H-HMC.

## Before the cluster
- [ ] Commit and push `bayes/run_highres_{mala,hmc}.m`, `scripts/highres_{mala,hmc}.sh` and the updated `scripts/highres_submit.sh`.

## On the cluster
- [ ] `git pull`.
- [ ] Copy `adjoint_rate_state_inversion/mat/fractal_fault_m1001_mp1001/` with `rsync`. It is gitignored and needs `receiverData.mat` and `receiverData_noisy_eta0.02.mat`.
- [ ] Optional smoke test: run `matlab -batch "maxIter=2; run_highres_map"` interactively to catch path or data errors.
- [ ] From `scripts/`, run `./highres_submit.sh`. It queues the MAP job, then 4 chains each of H-pCN, H-MALA and H-HMC once the MAP job succeeds.
- [ ] Check `highres_map_*.out`. The inflation should converge and the Hessian should be positive definite. Otherwise:
  - resubmit with `sbatch --export=ALL,RESUME=true highres_map.sh`
  - cancel the waiting chain jobs, which will never start because their dependency failed
  - run `./highres_submit.sh --skip-map` once the MAP job is done.
- [ ] Check the chain logs.
  - H-pCN: if `std(phi)` is much larger than 1, rerun with `sbatch --export=ALL,S=0.5 highres_hpcn.sh`.
  - H-MALA and H-HMC: check the acceptance rate after adaptation. The adaptation targets are 0.574 for MALA and 0.65 for HMC.

The budget is about 1000 PDE solves per chain for each sampler (H-pCN 1000 iterations, H-MALA 500, H-HMC 100 × 5 leapfrog steps). Change it with `NITER` (and `NLEAP`) through `--export`.

## Back home
- [ ] Copy back `bayes/results/{highres_map,hpcn_highres,mala_highres,hmc_highres}_*.mat`.
- [ ] Run `python diagnostic.py` on each set of chains. Compare R-hat and ESS per PDE solve across the three samplers.
- [ ] Figures and numbers for the report:
  - the mean with 95% bands
  - $P(a(x) < b)$ per node
  - the VW edge positions
  - Laplace vs. MCMC
  - the inflation factor
  - the choice of m = 151, from `discretization_error_mp*.mat`
  - a side-by-side with the inverse-crime results.
