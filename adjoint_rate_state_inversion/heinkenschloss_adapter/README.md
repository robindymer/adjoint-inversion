Adapters that let `matlab_impl_constr/optimization/newton_cg.m` (M. Heinkenschloss) run on the
rate-and-state inversion. `newton_cg` calls `fval`, `grad`, `Hessvec`, `xnew` and `xprod` by name,
so `adjopt.optimize_general` (option `'newton-cg-ls'`) adds this folder to the path only during
that solve. Do not add it to the path permanently, as `grad` etc. could shadow other functions.
