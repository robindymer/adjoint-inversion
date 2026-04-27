function [parameters, trueParameters] = rsFrictionTimeDepNonlinear(opts)
    default_arg('opts',struct);
    default_field(opts,'misfitType','velocity');
    
    friction_law = 'asinh';
    state_law = 'aging';

    % Operators
    opset = @sbp.D2Nonequidistant;
    
    % Domain
    L = 2;
    xlims = [-L/2, 0, L/2];
    m = 201;
    T = 1;

    mu_p = 2; % GPa
    mu_m = 1.5; % GPa
    rho = 1.3; %Mg/m3

    
    % TODO: Change to outflow bcs
    bc.type = {'outflow','outflow'};
    bc.data = {[],[]};
    bc.method = 'erickson2022';
    
    % Sources
    sources = [];
    
    % friction    
    a_true = 1.5;
    a = 1;
    b = 0.5;
    V0 = 1.1;
    f0 = 0.6;
    sigma0 = 1.4;
    D_c = 0.1;
    tau0 = f0 * sigma0;

    % External loading parameters
    T_E = 0.1;
    A = 2;
    
    friction.params.a = a;
    friction.params.b = b;
    friction.params.V0 = V0;
    friction.params.f0 = f0;
    friction.params.sigma0 = sigma0;
    friction.params.D_c = D_c;
    friction.params.tau0 = tau0;
    friction.params.A = A;
    friction.params.T_E = T_E;
    
    funs = pars.rsFrictionFunctions(friction_law,state_law,friction.params);
    
    % Set fault interface method.
    friction.method = 'erickson2022';

    % TODO: This should be handled inside the discrs.
    % Add external loading to tau
    tau_L = @(t) A*(1 - exp(-t/T_E)) + tau0;
    switch friction.method
    case 'erickson2022'
        tau = @(t, V, Psi, a) funs.tau(V, Psi, a);
        tauinv = @(t, tau, Psi, a) 2*V0./exp(Psi./a).*sinh(tau./(sigma0*a));
    case 'standard'
        tau = @(t, V, Psi, a) funs.tau(V, Psi, a) - tau_L(t).*ones(size(V));
        tauinv = [];
    end   
    %tau_L = @(t) A./(1 + exp(-(1/T_E)*(t-t_E0))) + tau0;

    % Steady state friction function
    f_ss  = @(V,a,b) f0-(b-a).*log(abs(V)./V0);
    
    % Pack functions
    friction.funs = funs;
    friction.funs.tau = tau;
    friction.funs.tau_L = tau_L;
    friction.funs.f_ss = f_ss;
    friction.funs.tauinv = tauinv;

    
    % initial conditions
    sigma = 1/15;
    x_s = -0.5;
    u0 = @(x) exp(-((x-x_s)/sigma).^2);
    v0 = @(x) 2*(x-x_s)/sigma^2.*u0(x);
    initialconditions.u = {u0, u0};
    initialconditions.v = {v0, v0};
    initialconditions.Psi = 0.5;
    initialconditions.x_s = x_s;
    initialconditions.sigma = sigma;
    
    receivers = struct;
    receivers.x = [-L/8];
    receivers.blockIds = [1];
    receivers.source_fun = [];
    receivers.data = [];


    tsOpts.forwardMethod.order = 3;
    tsOpts.forwardMethod.adaptive = true;
    tsOpts.forwardMethod.rtol = 1e-5; 
    tsOpts.forwardMethod.reportRetry = false;

    tsOpts.adjointMethod.order = 3;
    tsOpts.adjointMethod.adaptive = false;
    tsOpts.adjointMethod.rtol = 1e-5; % Only relevant for adaptive ts
    tsOpts.adjointMethod.reportRetry = false; % Only relevant for adaptive ts
    tsOpts.k = [];

    interpolate_data = tsOpts.adjointMethod.adaptive; % only interpolate if adaptive adjoint timestepping.

    
    % Return struct
    parameters = struct;
    parameters.opset = opset;
    parameters.m = m;
    parameters.xlims = xlims;
    parameters.material.rho = {@(x) 0*x+rho, @(x) 0*x+rho};
    parameters.material.mu = {@(x) 0*x+mu_m, @(x) 0*x+mu_p};
    parameters.T = T;
    parameters.bc = bc;
    parameters.friction = friction;
    parameters.sources = sources;
    parameters.receivers = receivers;
    parameters.initialconditions = initialconditions;
    parameters.misfitType = opts.misfitType;
    parameters.dim = 1;
    parameters.interpolate_data = interpolate_data;
    parameters.tsOpts = tsOpts;
    
    trueParameters = parameters;
    trueParameters.friction.params.a = a_true;
end