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

% NOTE: No \partial_b since f does not depend on b
F_V_V = diff(F_V,'V');
F_V_Psi = diff(F_V,'Psi');
F_V_a = diff(F_V,'a');
F_Psi_Psi = diff(F_Psi,'Psi');
F_Psi_a = diff(F_Psi,'a');
F_a_a = diff(F_a, 'a');

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
matlabFunction(F_V_a,'File','F_V_a','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_Psi_Psi,'File','F_Psi_Psi','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_Psi_a,'File','F_Psi_a','Vars',{V, Psi, a, sigma0, V0, tau0});
matlabFunction(F_a_a,'File','F_a_a','Vars',{V, Psi, a, sigma0, V0, tau0});

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

G_V_Psi = diff(G_V,'V','Psi');
G_V_V = diff(G_V,'V','V');
G_V_a = diff(G_V,'V','a');
G_Psi_Psi = diff(G_Psi,'Psi','Psi');
G_Psi_a = diff(G_Psi,'Psi','a');
G_a_a = diff(G_a, 'a');
G_b_b = diff(G_b, 'b');

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
matlabFunction(G_V_a,'File','G_V_a','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_Psi_Psi,'File','G_Psi_Psi','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_Psi_a,'File','G_Psi_a','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_a_a,'File','G_a_a','Vars',{V, Psi, a, b, f0, V0, D_c});
matlabFunction(G_b_b,'File','G_b_b','Vars',{V, Psi, a, b, f0, V0, D_c});
