% Adapter for matlab_impl_constr/optimization/newton_cg.m: Hessian at x applied to v.
function [Hv] = Hessvec(v, x, usr_par)
    Hv = adjopt.hessianVector(x, v, usr_par.optObj) + usr_par.reg.hessVec(x, v);
end
