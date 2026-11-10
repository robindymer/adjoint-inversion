- Sampling the posterior -> uncertainty quantification
- UQ - compute posterior quantities (expectations, probs, dists) of QOI and report them with their own error
- Roberts, Gelman, and Gilks [104] showed that the best compromise is a step size s ≈ 2.38/√n,
AR 0.234
- Sampling Gaussian: μ + Rξ ∼ N (μ, Γ), härlett från faktun Bx + c ∼ N (Bμ + c, BΓBT)
- Remark 4.21: What to report about a sampler - viktig
- Tolkar nu som att med Gaussian och Smootness prior kan vi slippa HMC, blir lättare

## Frågor
- Behöver man i vanliga fall oroa sig över att modellfelet ska överskugga vårt gaussian noise?
Behöver vara lite noggrann och försiktig när vi går till high-res data
