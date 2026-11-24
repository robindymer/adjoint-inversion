## Revision of plan.md and claude_plan.md

- [ ] When we do hierarchical prior: Include TV-regularized L-BFGS. If we do this
via the noise model we can choose lambda via the
discrepancy principle, but for now we can just take
something that looks good.
- [x] Include H1 regularized 
- [x] Add noise to misfit
1) Add to data
2) Remember to scale misfit by 1 / (sigma^2 delta t_d)
- LATER ON: When using high-res data, remember to add noise to it and
compare m=36, 51, 76 against high-res data.
If M_model >> M_noise, inflate noise (simple approach). If data
becomes uninformative, increase m. More rigorous approach,
enhanced error model of Kaipio and
Somersalo as stated in course litterature.
- [x] At some point: Include deterministic, noise, inverse crime, H1 and TV regularized.
This is a more fair comparison compare to noise-free inversion.
- [x] Implement smoothness prior
- [x] Judge it by its samples
- [x] Implement pCN, MALA and HMC samplers
- [x] Sample from posterior and compute
1) Mean with 95% bands, compare against a_true
2) P(a(x) < b) per node
3) VW edge positions
- [x] Missing diagnostic.py plots (Remark 4.21):
1) Posterior std vs Laplace std per node (see diagnostic_reference.py)
2) a_MAP and a few samples on the posterior plot
3) Histogram of VW edge positions / VW extent
- [ ] diagnostic.py fixes:
1) Cost from chain.nFwd + chain.nAdj, not C*N (HMC undercounted 10x); print chain.time
2) Actually apply burn_in
3) Check hpCN starts are dispersed (widened Laplace), else Rhat is meaningless
- [ ] Verification (Sec. 4.6):
1) Compare pCN/MALA/HMC/hpCN mean and std within MC errors
2) Fixed-truth coverage also tests the prior (a_true not drawn from prior); full coverage out of scope, say so
- [ ] Report (template):
1) Dependence on noise level eta and prior gamma
2) Table: rel. error, coverage, wall-clock for det., MAP, posterior mean
3) 95 % bands need ESS of a few thousand, not 400
- Hierarchical prior: use prior normalization (Remark 4.15) against the funnel, expect high Rhat at jumps

[x] Implement low rank laplace approximation from Hv products
[x] Center the Laplace approximation on the converged MAP (as in the course, Sec. 5.2): replace the
one-Newton-step a_MAP with the newton-tr H1 solution (polished with a few Newton steps on U, since
tikhonov-grad lacks the a_1 anchor of the prior), recompute H there, before the production runs.
[ ] Compare with computing the hessian by simply using finite differences
[ ] Implement hierarchical prior

This is done in the inverse crime setting
Remember that it can be wise to precondition MALA with
the prior covariance C, and use C inverse as HMC
mass matrix.
Also remember to reject bad points (a <= 0 or solver yields NaN)

- Lastly, do the above for high-resolution data!
#### Misc
- Decide on if we want to keep Robbins-Monroe or just manually tune
- When comparing chains we start from different initial points. For pCN, this starting
point is drawn from the prior, reasonable. For H-pCN, the Laplace approximation is centered
on the posterior, but we want to be a bit further away to verify burn-in. Therefore we widen
the laplace.
- TODO: Read about the Laplace approximation. Seems like a no-brainer for our case where computing
the hessian is feasible and only needs to be done once?
When would H-pCN possibly beat MALA and HMC? Probably where we use high-res data and have less of a
Gaussian posterior. Also, using a hiearchichal prior would make the posterior non-gaussian.
- LBFGS gets stuck, even when regularized. Newton-TR continues
working well, so we go with that for the MAP. The objective is simply
not nice enough with noise, for LBFGS.

#### Maximum likelihood, MAP and sampling
This is three different levels of complexity.
In Vidar's article, we simply minimze the L2 norm of the residual.
This is maximum likelihood. If we add noise to the data, we also need
to regularize so we don't fit noise. Solving this is then our MAP estimate.
Still, we only get one solution. The final step is to sample from our posterior,
and this gives credibility bands.