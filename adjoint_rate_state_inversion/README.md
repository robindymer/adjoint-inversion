# Testing and verification of numerical methods in Adjoint-based inversion for stress and frictional parameters in earthquake modeling
This archive contains Matlab scripts and functions for testing and verification of the numerical methods presented in the paper 'Adjoint-based inversion for stress and frictional parameters in earthquake modeling' by Stiernström et al.

To run the scripts, first add subdirectories to the Matlab path by running `addSubpaths` in the Matlab command window.

## Reproducing figures from the paper
- Verification of discrete gradient: Run the script `gradientConvergenceFriction.m`. Performs a convergence study of the discrete gradient with displacement and velocity misfits. This is done through successive refinements of a finite difference approximation of the gradient.
- Inverse crime for direct effect: Run the script `direct_effect_inverse_crime.m`. Inverts for the direct effect parameter 'a', using synthetic data recorded from a forward solve on the same computational grid.

## Other useful scripts for testing/verification
   - `gradientTestFriction2D.m` - Compares the adjoint gradient to a single finite difference gradient evaluation
   - `optimize_antiplaneshear.m` - Function used to invert for friction parameters using synthesized data.
   - `fractal_fault_2d.m` - Runs a forward simulation plotting e.g. wavefields.
   - `fractal_fault_generate_data.m` - Script used to generate and save synthetic data from a forward solve.
    
## Parameters used for running tests are found in the subfolder `+pars`.
   - `+pars/rsFriction2DFractalFaultVerification.m` - Parameter set used for gradient verification
   - `+pars/rsFriction2DFractalFaultInversion.m` - Parameter set used for inversions.
   - `+pars/rsFriction2DFractalFaultForward.m` - Parameter set for high-resolution simulation. Used to generate high-resolution synthetic data.

## Further details on the implementation:
   - `AntiplaneShear2DRSFrictionOpt.m` - Optimization class
   - `+elastic/AntiplaneShear2DRSFrictionFwdDiscr.m` - Discretization class of forward problem
   - `+elastic/AntiplaneShear2DRSFrictionAdjDiscr.m` - Discretization class of adjoint problem
   - `+elastic/+friction/+inversion/` - Contains the friction functions (friction law, state evolution). These together with derivative of the functions w.r.t. their parameters are generated through the function `generateFunctions.m` in the same directory.

This code is also available on Sourceforge at https://sourceforge.net/p/elastic-optimization/code/ci/zenodo-archive/tree/, version controlled under mercurical. To retrieve the code from the repository, do a read-only clone and update to the branch `zenodo-archive`.

