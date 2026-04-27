function [parameters, truePars] = rsFrictionLinear(opts)
    default_arg('opts',struct);
    default_field(opts,'friction_law','linear');
    default_field(opts,'state_law','linear');
    default_field(opts,'misfitType','velocity');
    
    friction_law = opts.friction_law;
    state_law = opts.state_law;
    
    % Domain
    L = 2;
    xlims = [-L/2, 0, L/2];
    m = 201;
    T = 1;

    mu = 1; % GPa
    rho = 1; %Mg/m3

    
    % TODO: Change to outflow bcs
    bc.type = {'outflow','outflow'};
    bc.data = {[],[]};
    
    % Sources
    sources = [];
    
    % friction    
    a_true = 2;
    a_perturbed = 1.5;
    b = 1;
        
    [fault_funs, state_funs] = pars.rsFrictionFunctions(friction_law,state_law,[]);
    friction.a = a_perturbed; % Initial perturbed guess for a.
    friction.b = b; 

    fault_funs.tau = @(t, V, Psi, a) fault_funs.tau(V, Psi, a);
    friction.fault_funs = fault_funs;
    friction.state_funs = state_funs;

    % TODO: Distinguis function handles in discr from fault_funs in a good way
    % Fix functions for a_perturbed
    friction.fault_fun = @(t, V, Psi) fault_funs.tau(t, V, Psi, a_perturbed);
    friction.fault_data = []; % Data empty initially.
    friction.state_evo_fun = @(Psi, V) state_funs.g(Psi, V, a_perturbed, b);
    friction.state_evo_data = []; % Data empty initially.
    
    % initial conditions
    sigma = 1/15;
    x_s = -0.5;
    u0 = @(x) exp(-((x-x_s)/sigma).^2);
    initialconditions.u = u0;
    initialconditions.ut = @(x) 2*(x-x_s)/sigma^2.*u0(x);
    initialconditions.Psi = 0.5;
    
    receivers = struct;
    receivers.x = [-L/8];
    receivers.blockIds = [1];
    receivers.source_fun = [];
    receivers.data = [];

    tsOpts.method.order = 3;
    tsOpts.method.adaptive = true;
    %tsOpts.k = 1e-7;
    tsOpts.k = [];
    % Only used for adaptive timestepping
    tsOpts.method.rtol = 1e-1; 
    tsOpts.method.reportRetry = false;
    
    % Return struct
    parameters = struct;
    parameters.m = m;
    parameters.xlims = xlims;
    parameters.material.rho = {@(x) 0*x+rho, @(x) 0*x+rho};
    parameters.material.mu = {@(x) 0*x+mu, @(x) 0*x+mu};
    parameters.T = T;
    parameters.bc = bc;
    parameters.friction = friction;
    parameters.sources = sources;
    parameters.receivers = receivers;
    parameters.initialconditions = initialconditions;
    parameters.misfitType = opts.misfitType;
    parameters.dim = 1;
    parameters.data_opts.interp_data = false;
    parameters.tsOpts = tsOpts;
    
    truePars = parameters;
    truePars.friction.a = a_true;
    truePars.friction.fault_fun = @(t, V, Psi) fault_funs.tau(t, V, Psi, a_true);
    truePars.friction.state_evo_fun = @(Psi,V) state_funs.g(Psi, V, a_true, b);
end