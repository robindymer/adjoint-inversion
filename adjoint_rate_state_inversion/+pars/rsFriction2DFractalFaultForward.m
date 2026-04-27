function [parameters, trueParameters] = rsFriction2DFractalFaultForward(opts)
    default_arg('opts',struct);
    default_field(opts,'misfitType','displacement');

    %% SBP operators and grid points
    opset = @sbp.D2Nonequidistant;
    % Grid points along the fault. The total number of gridpoints is
    % approximately [m,m]^2
    m = 1001;
    m_p = 1001;

    %% Domain
    W = 20; % Width km
    H = 20; % Height km
    T = 6;  % Final time
    xlims = [-W/2, W/2];
    ylims = [-H/2, H/2];
    domain = domains.FractalFault(xlims, ylims);

    %% Domain material parameters
    % material.rho - cell array with function handles of (x,y) for density for each block
    % material.mu - cell array with function handles (x,y) for bulk modulus for each block
    cs = 3.464; % km/s
    rho = 2.67; % g/cm3
    mu = cs^2*rho; % GPa

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
    %           for each boundary.
    % bc.data - cell array with function handles to time-dependent data. 
    %           May be specified for dirichlet or traction conditions.
    %           For zero data set entries to [].
    %           ex: bc.type = {'outflow', 'outflow', 'outflow', 'traction'};
    %               bc.data = {[],[],[],@(t) g(t)};
    bc.type = {'outflow', 'outflow', 'outflow', 'outflow'};
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
    
    receiverSpacing = 1;
    boundaries.outer.x = [-9, 9];
    boundaries.outer.y = [-9, 9];
    boundaries.inner.y = [-2, 2];
    if receiverSpacing == 1
        boundaries.inner.x = [-6, 7];
    elseif receiverSpacing == 2
        boundaries.inner.x = [-7, 7];
    end
    receivers = pars.placeReceivers(boundaries,receiverSpacing);
    
    
    %% Friction
    % Parameters
    % Must set a, b, V0, f0, sigma0, D_c, Psi0
    % These are also the variables which may be inverted for
    % The parameters can be either scalars or function handles of
    % (x,y).

    if_velocity_weakening = @(x,y) (x >= - 5) & (x <= 6);
    if_velocity_strengthening = @(x,y) 1-if_velocity_weakening(x,y);

    % 'True' parameter values, used to generate data for misfit
    % Note: Remember to set fields in true_params below.
    a_true =  @(x,y) 0.009*if_velocity_weakening(x,y) + 0.013*if_velocity_strengthening(x,y);
    %a = @(x,y) 1.1*a_true(x,y); % Initial guess. Perturbed a.
    a = 0.010;
    b = 0.011;
    sigma0 = 120;
    f0 = 0.6;
    V0 = 1e-6;
    D_c = @(x,y) (0.2*if_velocity_weakening(x,y) + 1*if_velocity_strengthening(x,y));
    Psi0 = 0.724339595021678;
    tau0 = 72; % Note: Really sigma_yz. Will be dotted with the normal component on the fault!
    

    % Set initial velocity and estimated final velocity
    % and compute initial state, and background traction
    A = 25;
    Vinit = 1e-12;
    

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
    x0 = 3;
    %y0 = -8.955116271972656e-01;
    std_dev = 2;
    tau_L = @(t, x, y) A*gaussian1d(x,x0,std_dev);
    
    % Pack parameters
    friction.rsParams.a = a;
    friction.rsParams.b = b;
    friction.rsParams.sigma0 = sigma0;
    friction.rsParams.f0 = f0;
    friction.rsParams.V0 = V0;
    friction.rsParams.D_c = D_c;
    friction.rsParams.tau0 = tau0;

    friction.loadingParams.Vinit = Vinit; 
    friction.loadingParams.A = A;

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
        
   %% Initial conditions
    % Initial values of fields u, v and Psi as a cell array
    % of grid functions for each block.
    % Set initialconditions = [] for zero initial conditions
    initialconditions.u = {@(x,y) 0*x, @(x,y) 0*x};
    initialconditions.v = {@(x,y) 0*x - Vinit/2, @(x,y) 0*x + Vinit/2};
    initialconditions.Psi = @(x,y) 0*x + Psi0;
    
    %% ERK time-stepper options.
    % For gradient computations 
    % Either use order = 4, adaptive_fwd = false;
    % or         order = 3, adaptive_fwd = true/false;
    order = 4;
    adaptive_fwd = false;
    tol = 1e-5; % % Only used if adaptive_fwd = true;
    reportRetry = false; % Only used if adaptive_fwd = true;
    k = 0.0005; % Time step. Will be aligned to final time. Set to empty if default CFL conditions are used
            % (specified in the AntiplaneShear2DRS<Fwd/Adj>Discr.m)

    % Pack tsOpts struct
    tsOpts.forwardMethod.order = order;
    tsOpts.forwardMethod.adaptive = adaptive_fwd;
    tsOpts.forwardMethod.rtol = tol;
    tsOpts.forwardMethod.reportRetry = reportRetry;
    tsOpts.adjointMethod.order = order; 
    tsOpts.k = k;
    
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
    parameters.tsOpts = tsOpts;
    parameters.dim = 2;
    parameters.m = m;
    parameters.m_p = m_p;
    parameters.filterOpts = [];
    
    % Parameter struct with true value of a.
    trueParameters = parameters;
    trueParameters.friction.rsParams.a = a_true;

end

% TODO: Check that gaussian matches Erics impl.
function f = gaussian1d(x, x0, d)
    f = exp(-0.5*(x-x0).^2/d^2);
end

function f = gaussian2d(x, y, x0, y0, d)
    f = exp(-0.5*((x-x0).^2 + (y-y0).^2)/d^2);
end