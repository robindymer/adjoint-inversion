function generateFunctions(state_evolution)
default_arg('state_evolution','slip_law')

% function for generating functions for the standard asinh rate-state 
% friction coefficient (Rice et al 2011) and aging law state evolution 
% together with partial derivatives to be used in inversion codes.


% Friction coefficient
syms V Psi sigma0 a b f0 V0 D_c tau0;
f = a*asinh(abs(V)*exp(Psi/a)/(2*V0));
f_ss = f0 - (b-a)*log(abs(V)/V0);
F =  sigma0*a*asinh(V*exp(Psi/a)/(2*V0)) - tau0;

F_V = diff(F,'V');
F_Psi = diff(F,'Psi');
F_sigma0 = diff(F,'sigma0');
F_a = diff(F,'a');
F_b = diff(F,'b');
F_f0 = diff(F,'f0');
F_V0 = diff(F,'V0');
F_D_c = diff(F,'D_c');
F_tau0 = diff(F,'tau0');

F_V_V = diff(F_V,'V');
F_V_Psi = diff(F_V,'Psi');
F_Psi_Psi = diff(F_Psi,'Psi');

matlabFunction(f,'File','friction_coeff','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(f_ss,'File','friction_coeff_steady_state','Vars',{V, a, b, f0, V0});
matlabFunction(F,'File','F','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_V,'File','F_V','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_Psi,'File','F_Psi','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_sigma0,'File','F_sigma0','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_a,'File','F_a','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_b,'File','F_b','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_f0,'File','F_f0','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_V0,'File','F_V0','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_D_c,'File','F_D_c','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_tau0,'File','F_tau0','Vars',{V, Psi, a, sigma0, V0, tau0});

matlabFunction(F_V_V,'File','F_V_V','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_V_Psi,'File','F_V_Psi','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_Psi_Psi,'File','F_Psi_Psi','Vars',{V, Psi, a, sigma0, V0, tau0});

% Aging law
switch state_evolution
    case 'aging_law'
        G = (b*V0/D_c)*(exp((f0-Psi)/b) - abs(V)/V0);
    case 'slip_law'
        G = -abs(V)/D_c*(f - f_ss);
    otherwise 
        error('State evolution equation %s not implemented', state_evolution);
end
G_V = diff(G,'V');
G_Psi = diff(G,'Psi');
G_sigma0 = diff(G,'sigma0');
G_a = diff(G,'a');
G_b = diff(G,'b');
G_f0 = diff(G,'f0');
G_V0 = diff(G,'V0');
G_D_c = diff(G,'D_c');
G_tau0 = diff(G,'tau0');

G_V_Psi = diff(G_V,'Psi');
G_V_V = diff(G_V,'V');
G_Psi_Psi = diff(G_Psi,'Psi');

matlabFunction(G,'File','G','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_V,'File','G_V','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_Psi,'File','G_Psi','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_sigma0,'File','G_sigma0','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_a,'File','G_a','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_b,'File','G_b','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_f0,'File','G_f0','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_V0,'File','G_V0','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_D_c,'File','G_D_c','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_tau0,'File','G_tau0','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_V_Psi,'File','G_V_Psi','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_V_V,'File','G_V_V','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_Psi_Psi,'File','G_Psi_Psi','Vars',{V, Psi, a, b, f0, V0, D_c});

% Second derivatives involving the friction parameters, used for Hessian computations.
% Generates F_V_<p>, F_Psi_<p>, G_V_<p>, G_Psi_<p> for each parameter p, and
% F_<p>_<q>, G_<p>_<q> for each pair, named in the order of params below (e.g. F_a_b, not F_b_a).
% See getSecondPartialDerivativeFuns.
F_vars = {V, Psi, a, sigma0, V0, tau0};
G_vars = {V, Psi, a, b, f0, V0, D_c};
params = {sigma0, a, b, f0, V0, D_c, tau0};
for i = 1:numel(params)
    p = params{i};
    p_name = char(p);
    matlabFunction(diff(F_V,p),'File',['F_V_' p_name],'Vars',F_vars);
    matlabFunction(diff(F_Psi,p),'File',['F_Psi_' p_name],'Vars',F_vars);
    matlabFunction(diff(G_V,p),'File',['G_V_' p_name],'Vars',G_vars);
    matlabFunction(diff(G_Psi,p),'File',['G_Psi_' p_name],'Vars',G_vars);
    for j = i:numel(params)
        q = params{j};
        q_name = char(q);
        matlabFunction(diff(diff(F,p),q),'File',sprintf('F_%s_%s',p_name,q_name),'Vars',F_vars);
        matlabFunction(diff(diff(G,p),q),'File',sprintf('G_%s_%s',p_name,q_name),'Vars',G_vars);
    end
end
