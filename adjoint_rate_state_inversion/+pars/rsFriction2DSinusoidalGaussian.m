function [parameters, trueParameters] = rsFriction2DSinusoidalGaussian(opts)
    default_arg('opts',struct);
    default_field(opts,'misfitType','displacement');
    default_field(opts, 'ic_method', 'erickson2022')
    
    friction_law = 'asinh';
    state_law = 'aging';

    m = 101;

    % Operators
    opset = @sbp.D2Variable;
    
    % Domain
    W = 2;
    H = 3;
    xlims = [-W, 0, W];
    ylims = [-H, 0];
    domain = domains.SinusoidalFault(xlims, ylims);

    T = 1;

    mu = 30; % GPa
    rho = 10/3; %Mg/m3
    material.rho = {@(x,y) 0*x+rho, @(x,y) 0*x+rho};
    material.mu = {@(x,y) 0*x+mu, @(x,y) 0*x+mu};

    
    % TODO: Change to outflow bcs
    b = domain.boundaryGroups.bottom;
    t = domain.boundaryGroups.surface;
    l = domain.boundaryGroups.left;
    r = domain.boundaryGroups.right;
    bc.ids = {l, r, b, t};
    % bc.type = {'traction','dirichlet', 'dirichlet', 'traction'};
    % bc.type = {'traction','traction', 'traction', 'traction'};
    % bc.type = {'dirichlet','dirichlet', 'dirichlet', 'dirichlet'};
    bc.type = {'outflow','outflow', 'outflow', 'outflow'};
    % bc.type = {'outflow', 'outflow', 'outflow', 'dirichlet'};
    bc.data = {[],[],[],[]};
    
    % Sources
    sources = [];
    
    % Friction parameters
    a_true = 1; % Used for generating data
    a_init = 0.5; % Initial parameter

    a = a_init;
    b = 2;
    f0 = 0.6;
    V0 = 1;
    sigma0 = 0.1;
    D_c = 0.01;
    Psi0 = 0.3;

    % Pack parameters
    friction.params.a = a;
    friction.params.b = b;
    friction.params.V0 = V0;
    friction.params.f0 = f0;
    friction.params.sigma0 = sigma0;
    friction.params.D_c = D_c;
    
    % Get friction functions for friction law and state law
    funs = pars.rsFrictionFunctions(friction_law,state_law,friction.params);

    % Pack functions
    friction.funs = funs;
    
    % Initialize data to zero
    friction.data.V = [];
    friction.data.Psi = [];
    
    % Initial conditions
    initialconditions.u = {@(x,y) gaussian2d(x,y,-W/2,-H/2,H/10), @(x,y) 0*gaussian2d(x,y,W/2,-H/2,H/10)};
    initialconditions.v = {@(x,y) 0*x, @(x,y) 0*x};
    initialconditions.Psi = @(x,y) 0*x + Psi0;
    
    receivers = struct;
    receivers.x = {[-W/4, -H/2]};
    receivers.blockIds = [1];
    receivers.source_fun = [];
    receivers.data = [];

    tsOpts.forwardMethod.order = 4;
    tsOpts.forwardMethod.adaptive = false;
    tsOpts.forwardMethod.rtol = 1e-3; 
    tsOpts.forwardMethod.reportRetry = true;

    tsOpts.adjointMethod.order = 4;
    tsOpts.adjointMethod.adaptive = false;
    tsOpts.adjointMethod.rtol = 1e-3; 
    tsOpts.adjointMethod.reportRetry = false;
    tsOpts.k = [];

    interp_data = tsOpts.adjointMethod.adaptive;
    
    % Form parameter struct
    parameters = struct;
    parameters.domain = domain;
    parameters.opset = opset;
    parameters.T = T;
    parameters.material = material;
    parameters.bc = bc;
    parameters.friction = friction;
    parameters.sources = sources;
    parameters.receivers = receivers;
    parameters.initialconditions = initialconditions;
    parameters.misfitType = opts.misfitType;
    parameters.interpolate_data = interp_data;
    parameters.tsOpts = tsOpts;
    parameters.dim = 2;
    parameters.m = m;
    
    % Parameter struct with true value of a.
    trueParameters = parameters;
    trueParameters.friction.params.a = a_true;

end

function f = gaussian2d(x, y, x0, y0, d)
    r = sqrt((x-x0).^2 + (y-y0).^2);
    f = gaussian(r, 0, d);
end