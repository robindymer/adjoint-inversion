function funs = rsFrictionFunctions(friction_law, state_law, params)
    default_arg('friction_law','asinh');
    default_arg('state_law','aging');
    switch friction_law
    case {'erickson-2022'}
        funs.f         = @(V, Psi, a) a.*asinh(abs(V));
        funs.tau       = @(V, Psi, a) a.*asinh(V);
        funs.tau_V     = @(V, Psi, a) a.*1./sqrt(V.^2+1);
        funs.tau_Psi   = @(V, Psi, a) 0*V;
        funs.tau_a     = @(V, Psi, a) asinh(V);
    case {'asinh'}
        syms V Psi a V0 sigma0;
        f = sigma0*a*asinh(abs(V)*exp(Psi/a)/(2*V0));
        tau =  sigma0*a*asinh(V*exp(Psi/a)/(2*V0));
        tau_V = diff(tau,'V');
        tau_Psi = diff(tau,'Psi');
        tau_a = diff(tau,'a');
        % Second order
        tau_V_V = diff(tau_V, 'V');
        tau_V_Psi = diff(tau_V, 'Psi');
        tau_V_a = diff(tau_V, 'a');
        tau_Psi_Psi = diff(tau_Psi, 'Psi');
        tau_Psi_a = diff(tau_Psi, 'a');

        f       = matlabFunction(f,'Vars',{V,Psi,a,V0,sigma0});
        tau     = matlabFunction(tau,'Vars',{V,Psi,a,V0,sigma0});
        tau_V   = matlabFunction(tau_V,'Vars',{V,Psi,a,V0,sigma0});
        tau_Psi = matlabFunction(tau_Psi,'Vars',{V,Psi,a,V0,sigma0});
        tau_a   = matlabFunction(tau_a,'Vars',{V,Psi,a,V0,sigma0});
        tau_V_V = matlabFunction(tau_V_V, 'Vars',{V,Psi,a,V0,sigma0});
        tau_V_Psi = matlabFunction(tau_V_Psi, 'Vars',{V,Psi,a,V0,sigma0});
        tau_V_a = matlabFunction(tau_V_a, 'Vars',{V,Psi,a,V0,sigma0});
        tau_Psi_Psi = matlabFunction(tau_Psi_Psi, 'Vars',{V,Psi,a,V0,sigma0});
        tau_Psi_a = matlabFunction(tau_Psi_a, 'Vars',{V,Psi,a,V0,sigma0});
        clear a V V0 Psi sigma0;

        sigma0 = params.sigma0;
        V0 = params.V0;

        funs.f       =  @(V, Psi, a) f(V, Psi, a, V0, sigma0);
        funs.tau       = @(V, Psi, a) tau(V, Psi, a, V0, sigma0);
        funs.tau_V     = @(V, Psi, a) tau_V(V, Psi, a, V0, sigma0);
        funs.tau_Psi   = @(V, Psi, a) tau_Psi(V, Psi, a, V0, sigma0);
        funs.tau_a     = @(V, Psi, a) tau_a(V, Psi, a, V0, sigma0);
        funs.tau_V_V   = @(V, Psi, a) tau_V_V(V, Psi, a, V0, sigma0);
        funs.tau_V_Psi = @(V, Psi, a) tau_V_Psi(V, Psi, a, V0, sigma0);
        funs.tau_V_a   = @(V, Psi, a) tau_V_a(V, Psi, a, V0, sigma0);
        funs.tau_Psi_Psi = @(V, Psi, a) tau_Psi_Psi(V, Psi, a, V0, sigma0);
        funs.tau_Psi_a = @(V, Psi, a) tau_Psi_a(V, Psi, a, V0, sigma0);
        funs.F = @elastic.friction.asinh;

        % Specific for Erickson2022 interface treatment
        funs.Finv = @elastic.friction.erickson2022.asinh_inv;
        funs.nonlin_solve_fun = @elastic.friction.erickson2022.nonlin_solve_fun_asinh;
        
    case {'linear'}
        funs.f         = @(V, Psi, a) a*(abs(V) + Psi);
        funs.tau       = @(V, Psi, a) a*(V + Psi);
        funs.tau_V     = @(V, Psi, a) a.*ones(size(V));
        funs.tau_Psi   = @(V, Psi, a) a.*ones(size(V));
        funs.tau_a     = @(V, Psi, a) V + Psi;
    end

    switch state_law
    case 'none'
        funs.g       = @(V,Psi,a,b) 0*Psi;
        funs.g_Psi   = @(V,Psi,a,b) 0*Psi;
        funs.g_V     = @(V,Psi,a,b) 0*Psi;
        funs.g_a     = @(V,Psi,a,b) 0*Psi;
        funs.g_b     = @(V,Psi,a,b) 0*Psi;
    case 'aging'
        syms V Psi a b f0 V0 D_c;
        g = (b*V0/D_c)*(exp((f0-Psi)/b) - abs(V)/V0);
        g_V = diff(g,'V');
        g_Psi = diff(g,'Psi');
        g_a = diff(g,'a');
        g_b = diff(g,'b');
        % Second order
        g_V_Psi = diff(g_V, 'Psi');
        g_V_V = diff(g_V, 'V');
        g_V_a = diff(g_V, 'a');
        g_Psi_Psi = diff(g_Psi, 'Psi');
        g_Psi_a = diff(g_Psi, 'a');

        g     = matlabFunction(g,'Vars',{V, Psi, a, b, f0, V0, D_c});
        g_V   = matlabFunction(g_V,'Vars',{V, Psi, a, b, f0, V0, D_c});
        g_Psi = matlabFunction(g_Psi,'Vars',{V, Psi, a, b, f0, V0, D_c});
        g_a   = matlabFunction(g_a,'Vars',{V, Psi, a, b, f0, V0, D_c});
        g_b   = matlabFunction(g_b,'Vars',{V, Psi, a, b, f0, V0, D_c});
        g_V_Psi = matlabFunction(g_V_Psi, 'Vars',{V, Psi, a, b, f0, V0, D_c});
        g_V_V = matlabFunction(g_V_V, 'Vars',{V, Psi, a, b, f0, V0, D_c});
        g_V_a = matlabFunction(g_V_a, 'Vars',{V, Psi, a, b, f0, V0, D_c});
        g_Psi_Psi = matlabFunction(g_Psi_Psi, 'Vars',{V, Psi, a, b, f0, V0, D_c});
        g_Psi_a = matlabFunction(g_Psi_a, 'Vars',{V, Psi, a, b, f0, V0, D_c});
        
        clear V Psi a b f0 V0 D_c;

        f0 = params.f0;
        V0 = params.V0;
        D_c = params.D_c;

        funs.g       = @(V,Psi,a,b) g(V, Psi, a, b, f0, V0, D_c);
        funs.g_V     = @(V,Psi,a,b) g_V(V, Psi, a, b, f0, V0, D_c);
        funs.g_Psi   = @(V,Psi,a,b) g_Psi(V, Psi, a, b, f0, V0, D_c);
        funs.g_a     = @(V,Psi,a,b) g_a(V, Psi, a, b, f0, V0, D_c);
        funs.g_b     = @(V,Psi,a,b) g_b(V, Psi, a, b, f0, V0, D_c);
        funs.g_V_Psi = @(V,Psi,a,b) g_V_Psi(V, Psi, a, b, f0, V0, D_c);
        funs.g_V_V = @(V,Psi,a,b) g_V_V(V, Psi, a, b, f0, V0, D_c);
        funs.g_V_a = @(V,Psi,a,b) g_V_a(V, Psi, a, b, f0, V0, D_c);
        funs.g_Psi_Psi = @(V,Psi,a,b) g_Psi_Psi(V, Psi, a, b, f0, V0, D_c);
        funs.g_Psi_a = @(V,Psi,a,b) g_Psi_a(V, Psi, a, b, f0, V0, D_c);
        funs.G = g;
    case 'linear'
        funs.g       = @(V,Psi,a,b) a*V - b*Psi;
        funs.g_Psi   = @(V,Psi,a,b) -b*ones(size(Psi));
        funs.g_V     = @(V,Psi,a,b) a*ones(size(Psi));
        funs.g_a     = @(V,Psi,a,b) V;
        funs.g_b     = @(V,Psi,a,b) -Psi;
    end
end