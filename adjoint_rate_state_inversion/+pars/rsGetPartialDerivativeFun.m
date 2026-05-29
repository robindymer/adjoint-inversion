function [F_p, G_p] = rsGetPartialDerivativeFun(funs, param)
    switch param
        case 'a'
            F_p = funs.tau_a;
            G_p = funs.g_a;
        case 'b'
            F_p = @(V,Psi,a) zeros(size(funs.tau_a(V,Psi,a))); % Hacky. F does not depend on b (four our friction law)
            G_p = funs.g_b;
        case 'sigma0'
            error('rsFrictionFunctions does not currently provide sigma0 partial derivatives.');
        case 'f0'
            error('rsFrictionFunctions does not currently provide f0 partial derivatives.');
        case 'V0'
            error('rsFrictionFunctions does not currently provide V0 partial derivatives.');
        case 'D_c'
            error('rsFrictionFunctions does not currently provide D_c partial derivatives.');
        case 'tau0'
            error('rsFrictionFunctions does not currently provide tau0 partial derivatives.');
        otherwise   
            error('Parameter %s not supported.', param);
    end