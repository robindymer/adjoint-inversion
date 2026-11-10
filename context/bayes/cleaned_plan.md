## Revision of plan.md and claude_plan.md

- Include TV-regularized L-BFGS. If we do this
via the noise model we can choose lambda via the
discrepancy principle, but for now we can just take
something that looks good.
- Include H1 regularized 
- Add noise to misfit
1) Add to data
2) Remember to scale misfit by 1 / (sigma^2 delta t_d)
- LATER ON: When using high-res data, remember to add noise to it and
compare m=36, 51, 76 against high-res data.
If M_model >> M_noise, inflate noise (simple approach). If data
becomes uninformative, increase m. More rigorous approach,
enhanced error model of Kaipio and
Somersalo as stated in course litterature.
- Implement smoothness prior
- Implement pCN, MALA and HMC samplers
- Sample from posterior and computer
1) Mean with 95% bands, compare against a_true
2) P(a(x) < b) per node
3) VW edge positions
This is done in the inverse crime setting

Remember that it can be wise to precondition MALA with
the prior covariance C, and use C inverse as HMC
mass matrix.
Also remember to reject bad points (a <= 0 or solver yields NaN)

#### If time is left
- Implement low rank laplace approximation from Hv products
- Implement hierarchical prior
