function [F_p, G_p] = getPartialDerivativeFun(param)
    switch param
        case 'sigma0'
            F_p = @elastic.friction.inversion.F_sigma0;
            G_p = @elastic.friction.inversion.G_sigma0;
        case 'a'
            F_p = @elastic.friction.inversion.F_a;
            G_p = @elastic.friction.inversion.G_a;
        case 'b'
            F_p = @elastic.friction.inversion.F_b;
            G_p = @elastic.friction.inversion.G_b;
        case 'f0'
            F_p = @elastic.friction.inversion.F_f0;
            G_p = @elastic.friction.inversion.G_f0;
        case 'V0'
            F_p = @elastic.friction.inversion.F_V0;
            G_p = @elastic.friction.inversion.G_V0;
        case 'D_c'
            F_p = @elastic.friction.inversion.F_D_c;
            G_p = @elastic.friction.inversion.G_D_c;
        case 'tau0'
            F_p = @elastic.friction.inversion.F_tau0;
            G_p = @elastic.friction.inversion.G_tau0;
        otherwise   
            error('Parameter %s not supported.', param);
    end