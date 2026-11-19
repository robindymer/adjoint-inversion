- Sampling the posterior -> uncertainty quantification
- UQ - compute posterior quantities (expectations, probs, dists) of QOI and report them with their own error
- Roberts, Gelman, and Gilks [104] showed that the best compromise is a step size s ≈ 2.38/√n,
AR 0.234
- Sampling Gaussian: μ + Rξ ∼ N (μ, Γ), härlett från faktun Bx + c ∼ N (Bμ + c, BΓBT)
- Remark 4.21: What to report about a sampler - viktig
- Tolkar nu som att med Gaussian och Smootness prior kan vi slippa HMC, blir lättare

- "Sampling the bathymetry posterior" har vi nog mycket att lära oss från
- Viktigt att ta med "Remark 4.21" i rapporten
- Key förklaring: "A posterior that is narrow in some directions and wide
in others forces steps small enough for the narrowest direction, and the walker then crawls along
the widest one. The remedy is to shape the steps like the posterior, with C a Cholesky factor of
the posterior covariance, if it is known, or of the inverse Hessian of U at the MAP estimate."

Rhat - potential scale reduction factor,
measures if chains have converged or not,
using their variance and mean values.
Rhat > 1 -> chains have not forgotten their starts
ESS - effective sample size,
ESS = total sample size / autocorrelation time.
ESS > 400 is usually required.
Autocorrelation time - measures the amount of samples
in between independent samples

## Chapter 5
- "The samplers of Chapter 4 remain the reference, and the test problems of that chapter, whose posteriors we
know, are where we check.", good to have in mind
- Explore table 5.1? Or motivate why we should do Laplace approximation?

## Frågor
- Behöver man i vanliga fall oroa sig över att modellfelet ska överskugga vårt gaussian noise?
Behöver vara lite noggrann och försiktig när vi går till high-res data
