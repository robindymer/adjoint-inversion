function [parameters, trueParameters] = rsFrictionOffFaultTransportHessian(opts)
    default_arg('opts',struct);
    default_field(opts,'misfitType','displacement');
    default_field(opts, 'ic_method', 'erickson2022')
    
    friction_law = 'asinh';
    state_law = 'aging';

    % Operators
    opset = @sbp.D2Nonequidistant;
    
    % Domain
    L = 1;
    xlims = [-L/2, 0, L/2];
    m = 201;
    T = 1e-1;

    mu = 30; % GPa
    rho = 10/3; %Mg/m3
    material.rho = {@(x) 0*x+rho, @(x) 0*x+rho};
    material.mu = {@(x) 0*x+mu, @(x) 0*x+mu};

    
    % TODO: Change to outflow bcs
    bc.type = {'outflow','outflow'};
    bc.data = {[],[]};
    bc.method = 'erickson2022';
    
    % Sources
    sources = [];
    secondOrderSources = [];
    
    % Friction parameters
    a_true = 0.01; % Used for generating data
    a_init = 0.008; % Initial parameter

    a = a_init;
    eps_a = 0.01; % pertubation on a
    b = 0.02;
    f0 = 0.6;
    V0 = 1e-6;
    sigma0 = 100;
    D_c = 0.01;

    % External loading parameters
    T_E = 1e-3; 

    % Set initial velocity and estimated final velocity
    % and compute initial state, background traction, 
    Vinit = 1e-12;
    Vmax = 10;
    dtau = 3;
    
    Psi0 = f0 + dtau./sigma0 - b.*log(Vmax./V0);
    tau0 = sigma0.*(a.*log(Vinit./V0) + Psi0);
    A = sigma0.*a.*log(Vmax./Vinit);
    
    % fss = f0-(b-a).*log(abs(Vmax)./V0);
    % Psi0 = a*log(2*V0/Vmax*sinh(1/a*dtau/sigma0+fss/a));
    % tau0 = sigma0*a*asinh(Vinit/(2*V0)*exp(Psi0/a));
    % A = sigma0*a*asinh(Vmax/(2*V0)*exp(Psi0/a)) - tau0;
    % f0 = b*log(Vinit/V0)+Psi0;

    % Psi0 = f0 - b.*log(Vinit./V0); % Psi0 on direct effect curve, i.e g(Psi0,Vinit) = 0
    % tau0 = sigma0*a*asinh(Vinit/(2*V0)*exp(Psi0/a));
    % A = sigma0*a*asinh(Vmax/(2*V0)*exp(Psi0/a)) - tau0;
    
    % Pack parameters
    friction.params.a = a;
    friction.params.b = b;
    friction.params.V0 = V0;
    friction.params.Vinit = Vinit; 
    friction.params.Vmax = Vmax;
    friction.params.f0 = f0;
    friction.params.sigma0 = sigma0;
    friction.params.D_c = D_c;
    friction.params.tau0 = tau0;
    friction.params.dtau = dtau;
    friction.params.A = A;
    friction.params.T_E = T_E;
    friction.params.eps_a = eps_a;
    
    friction.method = opts.ic_method;
    
    % Get friction functions for friction law and state law
    funs = pars.rsFrictionFunctions(friction_law,state_law,friction.params);
     
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
    
    % Initialize data to zero
    friction.data.V = [];
    friction.data.Psi = [];
    
    % Initial conditions
    initialconditions.u = {@(x) 0*x, @(x) 0*x};
    initialconditions.v = {@(x) 0*x - Vinit/2, @(x) 0*x + Vinit/2};
    initialconditions.Psi = Psi0;
    
    receivers = struct;
    receivers.x = [-0.1];
    receivers.blockIds = [1];
    receivers.source_fun = [];
    receivers.data = [];

    % TODO: Necessary? What does this represent?
    secondOrderReceivers = struct;
    secondOrderReceivers.x = [-0.1];
    secondOrderReceivers.blockIds = [1];
    secondOrderReceivers.source_fun = [];
    secondOrderReceivers.data = [];

    tsOpts.forwardMethod.order = 3;
    tsOpts.forwardMethod.adaptive = true;
    tsOpts.forwardMethod.rtol = 1e-6; 
    tsOpts.forwardMethod.reportRetry = true;

    tsOpts.adjointMethod.order = 3;
    tsOpts.adjointMethod.adaptive = false;
    tsOpts.adjointMethod.rtol = 1e-6; 
    tsOpts.adjointMethod.reportRetry = false;

    tsOpts.secondOrderForwardMethod.order = 3;
    tsOpts.secondOrderForwardMethod.adaptive = false;
    tsOpts.secondOrderForwardMethod.rtol = 1e-6; 
    tsOpts.secondOrderForwardMethod.reportRetry = false;

    tsOpts.secondOrderAdjointMethod.order = 3;
    tsOpts.secondOrderAdjointMethod.adaptive = false;
    tsOpts.secondOrderAdjointMethod.rtol = 1e-6; 
    tsOpts.secondOrderAdjointMethod.reportRetry = false;

    %tsOpts.k = 5e-7;
    tsOpts.k = [];

    interp_data = tsOpts.adjointMethod.adaptive;
    
    % Form parameter struct
    parameters = struct;
    parameters.opset = opset;
    parameters.m = m;
    parameters.xlims = xlims;
    parameters.T = T;
    parameters.material = material;
    parameters.bc = bc;
    parameters.friction = friction;
    parameters.sources = sources;
    parameters.secondOrderSources = secondOrderSources;
    parameters.receivers = receivers;
    parameters.secondOrderReceivers = secondOrderReceivers;
    parameters.initialconditions = initialconditions;
    parameters.misfitType = opts.misfitType;
    parameters.interpolate_data = interp_data;
    parameters.tsOpts = tsOpts;
    parameters.dim = 1;
    
    % Parameter struct with true value of a.
    trueParameters = parameters;
    trueParameters.friction.params.a = a_true;
end