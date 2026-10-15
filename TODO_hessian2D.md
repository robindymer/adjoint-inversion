# TODO: 2D Hessian-vector code

Differences between the 2D Hessian path (`AntiplaneShear2DRSFrictionOptHessian` + `hessianVectorTestFriction2D`) and the gradient path (`AntiplaneShear2DRSFrictionOpt` + `gradientTestFriction2D`).

## Blocks running the gradient test's setup

- [x] **Coarse parameter grid (m_p < m) is broken.** The constructor uses `secondOrderForwardDiscr` and `secondOrderAdjointDiscr` (lines 140–141) before they are created (lines 156–157). Fix: delete lines 140–141. They are redundant, since `updateSecondOrderForwardDiscr`/`updateSecondOrderAdjointDiscr` set `rsParams` before every run. Then verify with the Hessian test at m = 51, m_p = 11. The gradient test runs m = 51, m_p = 11; the Hessian test is forced to m = m_p = 26. The `If2c`/`Ic2f` handling of `eps_p` and in the Hessian formula looks consistent, but it is untested.
- [ ] **Adaptive time stepping is unsupported.**
  - [ ] `runSecondOrderForward` uses `tsOpts.k` instead of `forwardTimeIntegrationData.k`. The second-order adjoint already uses the forward steps.
  - [ ] `runSecondOrderForward` overwrites `obj.Ht`.
  - [ ] `runSimulationAdaptive` uses an undefined `R` (line 641; should be `Rec`).
  - [ ] `runSimulationAdaptive` flattens `faultData.V` and `faultData.Psi` when pruning (linear indexing).
  - [ ] The Hessian test lacks the gradient test's "turn off adaptive, reuse k" block before the FD.

## Hard-coded to parameter `a`

- [ ] `updateSecondOrderForwardDiscr` and `hessianVectorFormula` use `F_p{1}`/`G_p{1}` and the explicit `F_a_a`, `G_a_a`, `F_V_a`, `F_Psi_a`, `G_V_a`, `G_Psi_a`. The gradient loops over `inversionPars` through `getPartialDerivativeFun`.
- [ ] `AntiplaneShear2DRSFrictionSecondOrderAdjDiscr` has the `*_a` field names (`F_V_a`, `G_Psi_a`, …) built into its traction and state-evolution terms.
- [ ] `hessianVectorFormula` returns a plain vector, not a struct over `inversionPars` (hence the "Trickery" in the test script).
- [ ] `computeHessianVectorFD` hard-codes `grad.a` and `grad0.a` instead of `grad.(parName)`.
- [ ] Only one parameter, so no mixed p₁ ≠ p₂ blocks: a single `eps_p` vector rather than one ∆p₁ per parameter, and nothing distinguishes p₁ from p₂.

## Correctness issue in an unused branch

- [ ] **Filtering with `filterResidual = false`.** In `updateSecondOrderAdjointDiscr`, the "filter data only" branch computes `approx_i − filtered data` (copied from the first-order adjoint). The second-order source should be just `approx_i` (∆Q, no data), as in the unfiltered branch. Harmless now because `doFilter = false`.

## Parameter-struct simplifications

- [ ] `delta_p` lives in `friction.rsParams`, so it is treated as a friction parameter: evaluated onto the fault, interpolated by `interpParameters`, and stored on every discretization. `obj.delta_p` is never used, and the comment at line 1160 about scaling by it is stale. Remove it, or rename it `fd_step` and keep it outside `rsParams`.
- [ ] `eps_p` is hard-coded as `eye(m_p)(:,k_idx)` in `rsFriction2DFractalFaultVerification`, a unit vector only. Allow any direction ∆p₁ (e.g. smooth or multi-node).

## Debug leftovers and cost

- [ ] Remove debug prints:
  - [ ] `grad0.a` in `computeHessianVectorFD` (prints the `ans = 1.0e+03 * …` blocks).
  - [ ] `M0` in `computeGradientFD` (the gradient class doesn't print it).
  - [ ] `disp` progress messages in `computeHessianVector`.
- [ ] Remove the redundant `obj.runForward()` calls before each `computeGradient()` in `computeHessianVectorFD` (doubles the forward solves).

## Test-script settings (not bugs)

- [ ] Align with `gradientTestFriction2D` once the coarse-grid and adaptive items are fixed: order 6 vs 8, T = 4 vs 5, m = 26/m_p = 26 vs 51/11, and different FD step ranges.
