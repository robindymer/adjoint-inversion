% Adapter for matlab_impl_constr/optimization/newton_cg.m: misfit at x (forward solve only).
function [fx] = fval(x, usr_par)
    optObj = usr_par.optObj;
    optObj.setParameters(optObj.parametersVectorToNative(x));
    optObj.updateForwardDiscr();
    optObj.runForward(false);
    fx = optObj.computeMisfit();
end
