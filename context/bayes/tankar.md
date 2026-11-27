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

## Tankar
- Just realized that it would probably have been better to use TV even
for Gaussian prior, we don't have to mimic regularization term and prior.
- Everything that is still unclear at the end, find correspondence in the course literature.

The standard name for this is white-noise spectral density. Discrete noise with variance σ² sampled every Δt has spectral density (intensity) R = σ²Δt. Going the other way, R_k = R/Δt. Your likelihood M/(σ²Δt_d) is just (1/(2R))∫r² dt with that R. This relation shows up mostly in estimation and filtering texts. In seismic waveform inversion papers people rarely state it explicitly.

References, from most to least directly checked:

MIT 2.160 lecture notes (Identification, Estimation and Learning), lecture 7. These derive the continuous Kalman–Bucy filter as the Δt → 0 limit of the discrete one, using R_k = R/Δt. They credit Brown & Hwang, Introduction to Random Signals and Applied Kalman Filtering, Section 7.1, which would be a citable textbook source. I've confirmed the notes; I haven't checked the book itself. (ocw.mit.edu)
Gelb (ed.), Applied Optimal Estimation, MIT Press 1974. This is the classic source for the same discrete-to-continuous step, but I couldn't confirm the exact section.
arXiv 1707.08131. This paper uses R_k = R/Δ explicitly and shows the discrete model converges to the continuous one. It's a good example if you want one in writing. (arxiv.org)
Stuart (2010), "Inverse problems: a Bayesian perspective", Acta Numerica. This is the Bayesian inverse problems version: data observed continuously in time with white noise, dy = G(u)dt + √R dW. The likelihood then contains (1/(2R))∫|G(u)|²dt, which is the rigorous form of your Φ. I'm citing this from memory and the search didn't surface the article, so check the section before citing it.