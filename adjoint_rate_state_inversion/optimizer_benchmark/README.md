<!-- NOTE: This folder is to reproduce or extend the method comparison -->
# Optimizer benchmark

Harness behind `notes/optimizer_comparison.md` (repo root). All scripts are run from MATLAB with this folder on the path. They `cd` to `adjoint_rate_state_inversion/` themselves.

| file | purpose |
|---|---|
| `setup_problem.m` | Builds the inverse-crime problem (`a` only, velocity misfit) with the Hessian class, as in `adjopt.optimize_general` |
| `InvProblem.m` | Counted, cached wrapper: `f` (forward), `fg` (+ adjoint), `Hv`, iterate log and cost/time budget |
| `run_method.m` | Runs one method and saves its history to `results/res_<method>_m<m>.mat` |
| `tr_ncg.m`, `lbfgs_ls.m`, `lbfgs_twoloop.m` | Prototype trust-region Newton-CG and line-search L-BFGS. The repo version of the former is `+adjopt/trustRegionNewtonCG.m` (method `newton_tr_repo`) |
| `summarize.m` | Table of all `results/res_*.mat`: final f, error, solve counts, cost to reach f < 1, 1e-1, 1e-2, 1e-3 |
| `time_evals.m` | Timing of gradient, misfit and Hv, an FD check of Hv, and the full Hessian and its eigenvalues at x0 (`results/H0.mat`) |
| `resolution_check.m` | Receiver recordings at the true parameters for m = 36, 76, 151, 251, compared with m = 251 |

## Running

Methods for `run_method`:

| method | what it runs |
|---|---|
| `fmincon_ip_lbfgs` | current `'lbfgs'` option |
| `fmincon_ip_cg` | current `'CG'` option |
| `fmincon_trr` | fmincon trust-region-reflective with Hv |
| `fminunc_qn` | fminunc quasi-Newton (BFGS) |
| `lbfgs` | own L-BFGS |
| `trncg` | prototype TR Newton-CG, no preconditioner |
| `trncg_pc` | prototype TR Newton-CG with L-BFGS preconditioner |
| `trncg_pc_switch` | `trncg_pc` with a switch to L-BFGS |
| `newton_tr_repo` | `+adjopt/trustRegionNewtonCG.m` |

Each run stops after 250 gradient-equivalents or after `maxTime` seconds.

Example: several methods in parallel from the shell, 3 threads each, m = 76.
```sh
cd adjoint_rate_state_inversion
for meth in newton_tr_repo fmincon_ip_lbfgs; do
  matlab -batch "addpath('optimizer_benchmark'); run_method('$meth', 10800, 3, 76)" > optimizer_benchmark/results/log_m76_$meth.txt 2>&1 &
done
```

Then summarize:
```matlab
run('optimizer_benchmark/summarize.m')
```

## Notes

- Cost is `(4.3·nFwd + 1.3·nAdj + 4.1·nHv)/5.6` gradient-equivalents. The weights are solve times measured at m = 36 (`time_evals.m`). Wall time is also logged in `hist.t`.
- `results/resolution_data.mat` (the recordings, about 20 MB) is not kept. `resolution_check.m` regenerates it in a few minutes.
- `results/res_*_m36.mat` without `newton_tr_repo` are from the prototypes. The repo implementation was benchmarked at m = 76.
