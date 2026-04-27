function [parameters, trueParameters] = rsFriction2D(opts)
    default_arg('opts',struct);
    default_field(opts,'misfitType','displacement');
    
    friction_law = 'asinh';
    state_law = 'aging';

    % Operators
    opset = @sbp.D2Variable;
    m = 51;
    
    % Domain
    W = 2;
    H = 3;
    xlims = [-W, 0, W];
    ylims = [-H, 0];
    domain = domains.VerticalFault(xlims, ylims);

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
    bc.type = {'traction','traction', 'traction', 'traction'};
    % bc.type = {'dirichlet','dirichlet', 'dirichlet', 'dirichlet'};
    % bc.type = {'outflow','outflow', 'outflow', 'outflow'};
    bc.data = {[],[],[],[]};
    
    % Sources
    sources = [];
    
    % Friction parameters
    a_true = 0.01; % Used for generating data
    a_init = 0.009; % Initial parameter

    a = a_init;
    b = 0.02;
    f0 = 0.6;
    V0 = 1e-6;
    sigma0 = 100;
    D_c = 0.01;

    % External loading parameters
    T_E = 1e-4; 

    % Set initial velocity and estimated final velocity
    % and compute initial state, background traction, 
    Vinit = 1e-12;
    Vmax = 10;
    dtau = 3;
    
    % Psi0 = f0 + dtau./sigma0 - b.*log(Vmax./V0);
    % tau0 = sigma0.*(a.*log(Vinit./V0) + Psi0);
    % A = sigma0.*a.*log(Vmax./Vinit);
    
    % fss = f0-(b-a).*log(abs(Vmax)./V0);
    % Psi0 = a*log(2*V0/Vmax*sinh(1/a*dtau/sigma0+fss/a));
    % tau0 = sigma0*a*asinh(Vinit/(2*V0)*exp(Psi0/a));
    % A = sigma0*a*asinh(Vmax/(2*V0)*exp(Psi0/a)) - tau0;
    % f0 = b*log(Vinit/V0)+Psi0;

    Psi0 = f0 - b.*log(Vinit./V0); % Psi0 on direct effect curve, i.e g(Psi0,Vinit) = 0
    tau0 = sigma0*a*asinh(Vinit/(2*V0)*exp(Psi0/a));
    A = sigma0*a*asinh(Vmax/(2*V0)*exp(Psi0/a)) - tau0;
    
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

    
    % Get friction functions for friction law and state law
    funs = pars.rsFrictionFunctions(friction_law,state_law,friction.params);
     
    % Add external loading to tau
    % tau_L = @(t) (A*(1 - exp(-t/T_E)) + tau0);
    tau_L = @(t, x, y) (1 + gaussian2d(x,y,-W/20,-H/3,H/5)) * (A*(1 - exp(-t/T_E))) + tau0;

    % Steady state friction function
    f_ss  = @(V,a,b) f0-(b-a).*log(abs(V)./V0);

    % Pack functions
    friction.funs = funs;
    friction.funs.tau_L = tau_L;
    friction.funs.f_ss = f_ss;
    
    % Initialize data to zero
    friction.data.V = [];
    friction.data.Psi = [];
    
    % Initial conditions
    initialconditions.u = {@(x,y) 0*gaussian2d(x,y,-W/2,-H/2,H/10), @(x,y) 0*gaussian2d(x,y,W/2,-H/2,H/10)};
    % initialconditions.v = {@(x,y) 0*x - Vinit/2, @(x,y) 0*x + Vinit/2};
    initialconditions.v = {@(x,y) 0*x, @(x,y) 0*x};
    initialconditions.Psi = @(x,y) 0*x + Psi0;
    
    receivers = struct;
    receivers.x = {[-W/4, -H/2]};
    receivers.blockIds = [1];
    receivers.source_fun = [];
    receivers.data = [];

    tsOpts.forwardMethod.order = 5;
    tsOpts.forwardMethod.adaptive = true;
    tsOpts.forwardMethod.rtol = 1e-3; 
    tsOpts.forwardMethod.reportRetry = true;

    tsOpts.adjointMethod.order = 4;
    tsOpts.k = [];

    
    % Form parameter struct
    parameters = struct;
    parameters.m = m;
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
    parameters.tsOpts = tsOpts;
    parameters.dim = 2;
    
    % Parameter struct with true value of a.
    trueParameters = parameters;
    trueParameters.friction.params.a = a_true;

end

function f = gaussian2d(x, y, x0, y0, d)
    r = sqrt((x-x0).^2 + (y-y0).^2);
    f = gaussian(r, 0, d);
end