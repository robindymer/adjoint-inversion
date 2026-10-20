% Adapter for matlab_impl_constr/optimization/newton_cg.m: Euclidean inner product.
function [x1x2] = xprod(x1, x2, usr_par)
    x1x2 = x1'*x2;
end
