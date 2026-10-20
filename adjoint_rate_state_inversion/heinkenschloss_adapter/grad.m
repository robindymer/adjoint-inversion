% Adapter for matlab_impl_constr/optimization/newton_cg.m: gradient at x.
% newton_cg always calls fval at x first, so the forward solve is reused and only the
% adjoint is solved. Also logs the iterate, since grad is called once per accepted x.
function [g] = grad(x, usr_par)
    optObj = usr_par.optObj;
    if ~isequal(optObj.parametersNativeToVector(optObj.pars), x)
        fval(x, usr_par);
    end
    optObj.updateAdjointDiscr();
    optObj.runAdjoint(false);
    % Mark first-order solutions as valid at x, so Hessvec does not recompute them
    optObj.firstOrderParsVec = optObj.parametersNativeToVector(optObj.forwardDiscr.friction.rsParams);
    g = optObj.gradientNativeToVector(optObj.gradientFormula());
    usr_par.outputFcn(x, optObj.computeMisfit());
end
