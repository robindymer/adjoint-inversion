## TODO
* As Vidar verified discrete gradient, we aim at verifying discrete Hessian against finite difference approximation
* Do "Inverse crime" since it is simple and demonstrates method
* Lastly we can do something slightly more real

## Notes on Vidar's code-based
gradientConvergenceFriction.m - entrypoint for convergence study of discrete gradient
runFDConvergence.m - set synthetic data with true parameter run, and calls for gradient computations
AntiplaneShear2DRSFrictionOpt.m - optimizer object where the main bulk lives! Has e.g. computeGradient, runForward, etc. 1D fault.
rsFriction2DFractalFaultVerification.m - paremeter generator.

computeGradient() -> runForward() -> updateAdjointDiscr() -> runAdjoint() -> gradientFormula()
computeGradientFD(deltaG) -> updateForwardDiscr() -> runForward() -> computeMisfit()
compareGradients() and gradientNorm() turn those two results into the relative error used in the plot.

## Next
Upgrade from FD second-order to fully analytic second-order FP/AP by generating and wiring second derivative friction functions and replacing the FD increments inside runSecondOrderForward/runSecondOrderAdjoint.

* Update the second order discr objects!!

## Other nice files
* rsFrictionFunctions.m - has derivatives of functions
* AntiplaneShearRSFrictionFwdDiscr.m - forward discretization
* rsFrictionOffFailtTransport.m - param file

## Funderingar
* Vad är det för interpolation som görs? "Fewer grid points used for inversion parameters. Use SBP-preserving interpolation"