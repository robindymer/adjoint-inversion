function [parameters, trueParameters] = rsFrictionErickson(opts)
    default_arg('opts',struct);
    default_field(opts,'misfitType','displacement');

    friction_law = 'erickson-2022';
    state_law = 'none';

    % Operators
    opset = @sbp.D2Nonequidistant;

    % Domain
    xlims = [-1, 0, 1];
    m = 201;
    T = 1;
    %T = 2.2;
    
    material.rho = {@(x) 0*x+1, @(x) 0*x+1};
    material.mu = {@(x) 0*x+1, @(x) 0*x+1};
    
    bc.type = {'outflow','outflow'};
    bc.data = {[],[]};
    bc.method = 'erickson2022';
    
    % Sources
    sources = [];
    
    % friction
    %a_true = 2;
    a_true = 20;

    a_init = 1.5;
    
    friction.params.a = a_init;
    friction.params.b = 0; % Not used.
    
    friction.funs = pars.rsFrictionFunctions(friction_law,state_law,friction.params);
    % friction.funs.tau = @(t, V, Psi, a) friction.funs.tau(V, Psi, a); % Add time-dependency since this is used by code.
    friction.funs.tau = @(t, V, Psi, a) a.*asinh(V);
    friction.funs.tauinv = @(t, tau, Psi, a) sinh(tau./a);
    
    friction.data = [];
    friction.method = 'erickson2022';
    
    % Initial conditions
    sigma = 1/15;
    x_s = -0.5;
    u0 = @(x) exp(-((x-x_s)/sigma).^2);
    v0 = @(x) 2*(x-x_s)/sigma^2.*u0(x);
    initialconditions.u = {u0, u0};
    initialconditions.v = {v0, v0};
    initialconditions.Psi = 0; % Note: Not used in the simulation.
    initialconditions.x_s = x_s;
    initialconditions.sigma = sigma;
    
    receivers = struct;
    receivers.x = [-0.25];
    receivers.blockIds = [1];
    receivers.source_fun = [];
    receivers.data = [];
    
    tsOpts.forwardMethod.order = 3;
    tsOpts.forwardMethod.adaptive = true;
    tsOpts.forwardMethod.rtol = 1e-6; 
    tsOpts.forwardMethod.reportRetry = false;
    
    tsOpts.adjointMethod.order = 3;
    tsOpts.adjointMethod.adaptive = false;
    % tsOpts.adjointMethod.rtol = 1e-4; 
    % tsOpts.adjointMethod.reportRetry = false;
    tsOpts.k = [];
    interpolate_data = false;
    
    % Return struct
    parameters = struct;
    parameters.opset = opset;
    parameters.m = m;
    parameters.xlims = xlims;
    parameters.material = material;
    parameters.T = T;
    parameters.bc = bc;
    parameters.friction = friction;
    parameters.sources = sources;
    parameters.receivers = receivers;
    parameters.initialconditions = initialconditions;
    parameters.misfitType = opts.misfitType;
    parameters.tsOpts = tsOpts;
    parameters.interpolate_data = interpolate_data;
    parameters.dim = 1;
    
    
    trueParameters = parameters;
    trueParameters.friction.params.a = a_true;
    end