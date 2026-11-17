## Revision of plan.md and claude_plan.md

- [ ] Include TV-regularized L-BFGS. If we do this
via the noise model we can choose lambda via the
discrepancy principle, but for now we can just take
something that looks good.
- [ ] Include H1 regularized 
- [x] Add noise to misfit
1) Add to data
2) Remember to scale misfit by 1 / (sigma^2 delta t_d)
- LATER ON: When using high-res data, remember to add noise to it and
compare m=36, 51, 76 against high-res data.
If M_model >> M_noise, inflate noise (simple approach). If data
becomes uninformative, increase m. More rigorous approach,
enhanced error model of Kaipio and
Somersalo as stated in course litterature.
- [ ] At some point: Include deterministic, noise, inverse crime, H1 and TV regularized.
This is a more fair comparison compare to noise-free inversion.
- [x] Implement smoothness prior
- [ ] Judge it by its samples
- [ ] Implement pCN, MALA and HMC samplers
- [ ] Sample from posterior and computer
1) Mean with 95% bands, compare against a_true
2) P(a(x) < b) per node
3) VW edge positions
This is done in the inverse crime setting

Remember that it can be wise to precondition MALA with
the prior covariance C, and use C inverse as HMC
mass matrix.
Also remember to reject bad points (a <= 0 or solver yields NaN)

#### If time is left
[ ] Implement low rank laplace approximation from Hv products
[ ] Implement hierarchical prior

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

#### Maximum likelihood, MAP and sampling
This is three different levels of complexity.
In Vidar's article, we simply minimze the L2 norm of the residual.
This is maximum likelihood. If we add noise to the data, we also need
to regularize so we don't fit noise. Solving this is then our MAP estimate.
Still, we only get one solution. The final step is to sample from our posterior,
and this gives credibility bands.