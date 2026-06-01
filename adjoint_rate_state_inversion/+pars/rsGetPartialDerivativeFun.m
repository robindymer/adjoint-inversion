% NOTE: Currently not in use
function [F_p, G_p, F_V_p, F_Psi_p, G_V_p, G_Psi_p] = rsGetPartialDerivativeFun(funs, param)
    switch param
        case 'a'
            F_p = funs.tau_a;
            G_p = funs.g_a;
            F_V_p = funs.tau_V_a;
            F_Psi_p = funs.tau_Psi_a;
            G_V_p = funs.g_V_a;
            G_Psi_p = funs.g_Psi_a;
        case 'b'
            F_p = funs.tau_b;
            G_p = funs.g_b;
            F_V_p = funs.tau_V_b;
            F_Psi_p = funs.tau_Psi_b;
            G_V_p = funs.g_V_b;
            G_Psi_p = funs.g_Psi_b;
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