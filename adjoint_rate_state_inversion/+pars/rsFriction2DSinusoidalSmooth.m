function [parameters, trueParameters] = rsFriction2DSinusoidalSmooth(opts)
    default_arg('opts',struct);
    default_field(opts,'misfitType','displacement');

    %% SBP operators and grid points
    opset = @sbp.D2Nonequidistant;
    % Grid points along the fault. The total number of gridpoints is approximately
    % [2m,m]^2
    m = 51;
    m_p = 17;
    
    %% Domain
    W = 2; % Width of single block
    H = 3; % Height
    T = 2;  % Final time
    xlims = [-W, 0, W];
    ylims = [-H, 0];
    domain = domains.SinusoidalFault(xlims, ylims);

    %% Domain material parameters
    % material.rho - cell array with function handle for density for each block
    % material.mu - cell array with function handle for bulk modulus for each block
    mu = 2; % GPa
    rho = 1.3; %Mg/m3
    material.rho = {@(x,y) 0*x+rho, @(x,y) 0*x+rho};
    material.mu = {@(x,y) 0*x+mu, @(x,y) 0*x+mu};
    
    %% Boundary conditions
    % Store boundary ids
    b = domain.boundaryGroups.bottom;
    t = domain.boundaryGroups.surface;
    l = domain.boundaryGroups.left;
    r = domain.boundaryGroups.right;
    bc.ids = {l, r, b, t};

    % Set types of conditions and potential data
    % bc.type - cell array with entries 'dirichlet', 'traction', or 'outflow'
    %           for each boundary. Ex: See below.
    % bc.data - cell array with function handles to time-dependent data. 
    %           May be specified for dirichlet or traction conditions.
    %           ex: bc.type = {'outflow', 'outflow', 'outflow', 'traction'};
    %               bc.data = {[],[],[],@(t) g(t)};
    bc.type = {'outflow', 'outflow', 'outflow', 'traction'};
    bc.data = {[],[],[],[]};
    
    %% Point sources
    % Set to empty if no point sources are included. Otherwise specify
    % sources.x - cell array with coordinate vectors for receiver positions
    %             ex: sources.x = {[x1,y1], [x2,y2], ...}
    % sources.blockIds - vector with grid block ids for each source. 
    %                    See domain for detail on block ids.
    %                    ex: sources.blockIds = [blockId1, blockId2, ...] 
    % sources.funs - cell array with function handles for time dependent part of sources
    %                ex: sources.funs = {@(t) f1(t), @(t) f2(t), ...}
    sources = [];

    %% Receivers
    % receivers.x - cell array with coordinate vectors for receiver positions
    %               or character string 'surface' to use all grid points on the surface.
    %               ex: receivers.x = {[x1,y1], [x2,y2], ...}
    % receivers.blockIds - vector with grid block ids for each source. Set to []
    %                      if receivers.x = 'surface'
    %                      See domain for detail on block ids.
    %                      ex: sources.blockIds = [blockId1, blockId2, ...] 
    receivers = struct;
    receivers.x = {[-W/4, 0], [-W/4, -H/4], [-W/4, -H/2], [-W/4, -3*H/4], [-W/4, -H], ...
                   [W/4, 0], [W/4, -H/4], [W/4, -H/2], [W/4, -3*H/4], [W/4, -H]};
    receivers.blockIds = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2];
    
    %% Friction
    % Parameters
    % Must set a, b, V0, f0, sigma0, D_c, Psi0
    % These are also the variables which may be inverted for
    % The parameters can be either scalars or function handles of
    % (x,y)
    a = @(x,y) 1*(y > -H/2) + 1.2*(y <= -H/2);
    b = 0.5;
    V0 = 1.1;
    f0 = 0.6;
    sigma0 = 1.4;
    D_c = 0.1;
    Psi0 = 0.5;
    tau0 = f0*sigma0;

    % 'True' parameter values, used to generate data for misfit
    % Note: Remember to set fields in true_params below.
    a_true = @(x,y) a(x,y) + 0.5;
    % b_true = 0.75;
    %tau0_true = 1.1*f0*sigma0;

    % External loading parameters
    A = 2; % Amplitud
    T_E = 1e-1; % Rise time
    
    
    % Set functions for friction law and state law
    % These functions are defined in elastic.friction.inversion.
    % See elastic.friction.inversion.generateFunctions for details
    % on how to generate the friction law and state law of interest.
    F = @elastic.friction.inversion.F;
    F_V = @elastic.friction.inversion.F_V;
    F_Psi = @elastic.friction.inversion.F_Psi;
    G = @elastic.friction.inversion.G;
    G_V = @elastic.friction.inversion.G_V;
    G_Psi = @elastic.friction.inversion.G_Psi;
    % Functions used to solve for characteristic variables
    % Note: Must match with F above! 
    % TODO: Should probably be generated together with above functions
    Finv = @elastic.friction.erickson2022.asinh_inv;
    nonlin_solve_fun = @elastic.friction.erickson2022.nonlin_solve_fun_asinh;
    % External loading
    % TODO: Include tau0 in F such that it can be inverted for
    tau_L = @(t, x, y) (1 + gaussian2d(x,y,-W/20,-H/3,H/5)) * (A*(1 - exp(-t/T_E)));

    % Steady state friction function
    f_ss  = @(V,a,b) f0-(b-a).*log(abs(V)./V0);

    % Pack parameters
    friction.rsParams.a = a;
    friction.rsParams.b = b;
    friction.rsParams.V0 = V0;
    friction.rsParams.f0 = f0;
    friction.rsParams.sigma0 = sigma0;
    friction.rsParams.D_c = D_c;
    friction.rsParams.tau0 = tau0;

    friction.loadingParams.A = A;
    friction.loadingParams.T_E = T_E;

    % Pack functions
    friction.funs.F = F;
    friction.funs.F_V = F_V;
    friction.funs.F_Psi = F_Psi;
    friction.funs.G = G;
    friction.funs.G_V = G_V;
    friction.funs.G_Psi = G_Psi;
    friction.funs.Finv = Finv;
    friction.funs.nonlin_solve_fun = nonlin_solve_fun;
    friction.funs.tau_L = tau_L;
    friction.funs.f_ss = f_ss;
    
    %% Initial conditions
    % Initial values of fields u, v and Psi as a cell array
    % of grid functions for each block.
    % Set initialconditions = [] for zero initial conditions
    initialconditions.u = {@(x,y) 0*gaussian2d(x,y,-W/2,-H/2,H/10), @(x,y) 0*gaussian2d(x,y,W/2,-H/2,H/10)};
    initialconditions.v = {@(x,y) 0*x, @(x,y) 0*x};
    initialconditions.Psi = @(x,y) 0*x + Psi0;
    
    %% ERK time-stepper options.
    % For gradient computations 
    % Either use order = 4, adaptive_fwd = false;
    % or         order = 3, adaptive_fwd = true;
    order = 4;
    adaptive_fwd = false;
    tol = 1e-5; % % Only used if adaptive_fwd = true;
    reportRetry = false; % Only used if adaptive_fwd = true;
    k = []; % Time step. Will be aligned to final time. Set to empty if default CFL conditions are used
            % (specified in the AntiplaneShear2DRS<Fwd/Adj>Discr.m)

    % Pack tsOpts struct
    tsOpts.forwardMethod.order = order;
    tsOpts.forwardMethod.adaptive = adaptive_fwd;
    tsOpts.forwardMethod.rtol = tol;
    tsOpts.forwardMethod.reportRetry = reportRetry;
    tsOpts.adjointMethod.order = order; 
    tsOpts.k = k;
    
    %% Pack parameters struct
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
    parameters.tsOpts = tsOpts;
    parameters.dim = 2;
    parameters.m = m;
    parameters.m_p = m_p;
    
    % Parameter struct with true value of a.
    trueParameters = parameters;
    trueParameters.friction.rsParams.a = a_true;
    % trueParameters.friction.rsParams.b = b_true;
    %trueParameters.friction.rsParams.tau0 = tau0_true;

end

function f = gaussian2d(x, y, x0, y0, d)
    r = sqrt((x-x0).^2 + (y-y0).^2);
    f = gaussian(r, 0, d);
end