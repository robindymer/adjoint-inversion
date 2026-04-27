classdef AntiplaneShearRSFrictionInterfaceDiscr < noname.Discretization
% Discretizes the 1D anti-plane shear problem
% 
% u_tt = 1/rho ( (mu*u_x)_x + S)
% Psi_t = g(V, Psi)
% 
% for velocity potential u, and state variable Psi, where
% S is a dirac delta source term, and g a state evolution law.
properties
    name         = "AntiplaneShearRSFrictionInterfaceDiscr"
    description  = 'Discretizes the 1D wave equation for antiplane shear with a rate- and state friction interface. Solves for displacement u, velocity v and state variable Psi'
    order                % Order of accuracy
    grid                 % Grid
    material             % Material parameters
    interpolate_data     % Flag specifying whether data should be interpolated
    ic_method            % String: 'erickson2022' or 'standard'
    
    % Discretization matrix. Formulate 1st-order in time system
    % v = u_t
    % w = [u; v; Psi]
    % such that
    % D = [0      I     0;
    %      1/rho D2(mu) 0;
    %      0      0     0];  
    % F = [0; 1/rho(S + SAT); g]^T
    % w_t = D*w + F(w)
    D 

    % Quadratures & norms
    H  % Spatial quadature of the SBP operator

    % Component operators for first-order in time system
    % Eu*w = u, Eu'*u = [u; 0; 0]; 
    % Ev*w = v, Ev'*v = [0; v; 0];
    % EPsi*w = Psi, EPsi'*Psi = [0; 0; Psi];
    Eu % Eu*w = u
    Ev % Ev*w = v
    EPsi % Epsi*w = Psi
    Eus
    Eusim
    Eusip
    V_star

    % Initial conditions
    w0 % w0 = [u0; v0; Psi0]
        
    % Boundary functions and operators
    boundary_data_fun % Boundary data forcing

    % Friction and operators
    friction             % Struct holding parameters and friction functions/data
    e_fault              % Fault restriciton operator
    penalty_fault        % Fault data penalty operator
    fault_traction_fwd   % Forward fault traction forcing
    fault_traction_adj   % Adjoint fault traction forcing
    state_evo_fwd        % Forward state evolution forcing
    state_evo_adj        % Adjoint state evolution forcing

    % Sources
    sources         % Struct holding parameters source functions/data
    S_cont          % Continuous source term
    S_discr         % Discrete source term
    S_friction
    dirac_deltas    % dirac delta functions

end

methods
    
    %   m          --      Number of grid points for each domain
    %   lims       --      Limits [x_l, x_i, x_r], specifying left/right boundaries and interface positions.
    %   order      --      order of accuracy
    %   material   --      struct for material properties
    %       material.rho   cell array with function rho(x) for density in each domain                 
    %       material.mu    cell array with function mu(x) for shear modulus in each domain                 
    %   bc         --      struct for boundary conditions
    %       bc.type             --  cell array with boundary conditions for left/right boundary
    %                               'dirichlet': Dirichlet conditions
    %                               'neumann': neumann conditions
    %                               'outflow bc': u_t +/- c*u_x = 0
    %       bc.data             --  cell array with boundary data for left/right boundary
    %                               empty for homogeneous data. Otherwise function data(t,u,ut).
    %   friction   --       struct for friction
    %       friction.fault_fun      -- function @(t,V,Psi)
    %       friction.fault_data     -- data.ut data for all time levels to be multiplied by velocity ut
    %                               -- data.Psi data for all time levels to be multiplied by velocity state variable Psi
    %       friction.state_evo_fun  -- function @(Psi,ut)
    %       friction.state_evo_data -- data.Psi data for all time levels to be multiplied by state variable Psi
    %                               -- data.ut data for all time levels to be multiplied by velocity ut
    %   sources    --      struct for point sources
    %       sources.x           -- array with source positions
    %       sources.source_fun  -- cell array with source functions @(t)
    %       sources.data        -- data for all time levels.
    %
    %   initialcondition    --      struct for initial conditions.
    %       initialcondition.u    -- function @(x,y) for potential condition
    %       initialcondition.ut   -- function @(x,y) for velocity condition
    %
    %   bc_method           --      string for numerical method used to impose boundary conditions
    %                                'standard' or 'erickson2022'
    %
    %   ic_method           --      string for numerical method used to impose nonlinear interface condition
    %                                'standard' or 'erickson2022'

    function obj = AntiplaneShearRSFrictionInterfaceDiscr(opset, m, lims, order, material, bc, friction,...
                                                         sources, initialcondition, interpolate_data,...
                                                         bc_method, ic_method)
        
        % ---- Default values --------------
        default_arg('opset',@sbp.D2Variable);
        default_arg('m',101);
        default_arg('lims',[-1, 0, 1]);
        default_arg('order',4);
        
        default_arg('material', struct);
        default_field(material, 'rho', {@(x) 0*x + 1, @(x) 0*x + 1});
        default_field(material, 'mu',  {@(x) 0*x + 1, @(x) 0*x + 1});
        
        
        default_arg('bc', struct);
        default_field(bc, 'type', {'outflow','outflow'});
        default_field(bc, 'data', {[],[]});
        
        default_arg('friction', struct);
        default_field(friction, 'params', []);
        default_field(friction, 'funs', struct);
        default_field(friction.funs, 'tau', []);
        default_field(friction.funs, 'g', []);
         
        default_arg('sources',[]);
        default_arg('initialcondition',[]);
        
        default_arg('interpolate_data', false);
        default_arg('bc_method', 'standard');
        default_arg('ic_method', 'standard');
        
        domain = elastic.AntiplaneShearRSFrictionDomain(lims,{'-','+'});
        
        % Get grid and multiblock diff op
        if isequal(opset,@sbp.D2Nonequidistant) % Boundary-optimized operators 
            mbGrid = domain.getGrid(m,'boundaryopt',order,'acc');
        else
            mbGrid = domain.getGrid(m,'equidist');
        end
        scheme_params_m = {@(x) 1./material.rho{1}(x), material.mu{1}, opset};
        scheme_params_p = {@(x) 1./material.rho{2}(x), material.mu{2}, opset};
        scheme_params = {scheme_params_m, scheme_params_p};

        mbDiffOp = multiblock.DiffOp(@scheme.Laplace1dVariable, mbGrid, order, scheme_params);

        % ---- Boundary conditions  -----

        % Initialize struct for stuff related to erickson2022 type BC
        erickson2022_struct = struct;
        erickson2022_struct.R = {};
        erickson2022_struct.Z = {};
        erickson2022_struct.gamma = {};
        erickson2022_struct.tau_star_penalty = {};
        erickson2022_struct.u_star_penalty = {};
        erickson2022_struct.tau_op = {};
        erickson2022_struct.e = {};

        % Initialize before looping over bc
        bc.ids = {domain.boundaryGroups.left, domain.boundaryGroups.right};
        closure_bc = 0*mbDiffOp.D;
        boundary_data_fun = @(t,u,ut) sparse(length(mbDiffOp.D), 1);

        % Loop over boundaries to impose bc on
        for i = 1:numel(bc.ids)
            bid = bc.ids{i};
            type = bc.type{i};
            e = mbDiffOp.getBoundaryOperator('e',bid);
            x = mbGrid.getBoundary(bid);

            switch bc_method
            case 'standard'
                switch type
                case 'outflow'
                    mu = e'*multiblock.evalOn(mbGrid, material.mu);
                    rho = e'*multiblock.evalOn(mbGrid, material.rho);
                    c = sqrt(mu/rho);
                    [closure, penalty] = mbDiffOp.boundary_condition(bid,'neumann');
                    boundary_data_fun = @(t,u,ut) boundary_data_fun(t,u,ut) -1/c*penalty*e'*ut;

                otherwise % Dirichlet or Neumann conditions
                    [closure, penalty] = mbDiffOp.boundary_condition(bid,type);
                    data = bc.data{1};
                    if ~isempty(data)
                        boundary_data_fun = @(t,u,ut) boundary_data_fun(t,u,ut) + penalty*data(t,u,ut);
                    end
                end

            case 'erickson2022'

                % Get closure only 
                [closure, ~] = mbDiffOp.boundary_condition(bid,'erickson2022');

                % Get operators
                rho_vec = multiblock.evalOn(mbGrid, material.rho);
                rho_mat = spdiag(rho_vec);
                mu_vec = multiblock.evalOn(mbGrid, material.mu);
                mu_boundary = e'*mu_vec;
                rho_boundary = e'*rho_vec;
                Z = spdiag(sqrt(mu_boundary.*rho_boundary));
                mu_boundary = diag(mu_boundary);
                rho_boundary = diag(rho_boundary);
                H = mbDiffOp.H;

                tau_op = mbDiffOp.getBoundaryOperator('d', bid)*mu_boundary;
                % Sign for normal derivate. Hacky solution :(
                if i == 1
                    tau_op = -tau_op;
                end

                % Set up penalties in front of fluxes
                tau_star_penalty = (rho_mat*H)\e;
                u_star_penalty = -(rho_mat*H)\tau_op;

                % Compute alpha in BC: tau = - alpha * u_t
                % Then compute reflection coefficient R from alpha
                % (or set R directly).
                I = eye(size(rho_boundary));
                switch type
                case 'outflow'
                    alpha = sqrt(rho_boundary\mu_boundary);
                    R = (I+alpha)\(I - alpha);
                case 'dirichlet'
                    R = -1*I;
                case 'neumann'
                    R = I;
                end
                R = sparse(R);

                % Penalty strength gamma
                gamma = 1.1/mbDiffOp.diffOps{1}.gamm*mu_boundary;

                erickson2022_struct.R{end+1} = R;
                erickson2022_struct.Z{end+1} = Z;
                erickson2022_struct.gamma{end+1} = gamma;
                erickson2022_struct.tau_star_penalty{end+1} = tau_star_penalty;
                erickson2022_struct.u_star_penalty{end+1} = u_star_penalty;
                erickson2022_struct.tau_op{end+1} = tau_op;
                erickson2022_struct.e{end+1} = e;
            end
            closure_bc = closure_bc + closure;
        end
        
        % ---- Sources -----
        dirac_deltas = [];
        if ~isempty(sources)
            ns = length(sources.x);
            dirac_deltas = cell(ns,1);
            % Create dirac delta functions
            for i = 1:ns
                x_i  = sources.x(i);
                blockId = sources.blockIds(i);
                grd = mbGrid.grids{i};
                H = mbDiffOp.diffOps{i}.H;
                delta_fun_local = diracDiscr(grd, x_i, order, 0, H);
                delta_fun = mbGrid.expandFunc(delta_fun_local, blockId);
                dirac_deltas{i} = delta_fun;
            end
        end

        % ---- Friction  -----
        bid = domain.boundaryGroups.interface;
        e_fault = mbDiffOp.getBoundaryOperator('e',bid);
        switch ic_method
        case 'standard'
            [closure_fault, penalty_fault] = mbDiffOp.boundary_condition(bid, 'traction');
        case 'erickson2022'
            % Penalty is not used, is a zero matrix.
            [closure_fault, penalty_fault] =  mbDiffOp.boundary_condition(bid, 'erickson2022');
        end

        % ---- Determine total number of unknowns -------%
        n = mbGrid.nPoints;
        switch bc_method
        case 'standard'
            ns = 0;
        case 'erickson2022'
            ns = length(erickson2022_struct.R);
        end

        switch ic_method
        case 'standard'
            nsip = 0;
            nsim = 0;
        case 'erickson2022'
            bidm = domain.boundaryGroups.interface{1};
            bidp = domain.boundaryGroups.interface{2};
            em = mbDiffOp.getBoundaryOperator('e', bidm);
            ep = mbDiffOp.getBoundaryOperator('e', bidp);
            [~, nsim] = size(em);
            [~, nsip] = size(ep);
        end

        % ---- 1st order system -----
        % u = [u_m; u_p];
        % v = u_t

        % -- Standard:
        % u = [u_m; u_p];
        % v = u_t
        % w = [u; v; Psi]
        % D = [0         I 0;
        %      D_laplace 0 0;
        %      0         0 0];
        % F = [0; S + penalty_bc*bc_data + penalty_fault*fault_fun; g]
        % w_t = D*w + F(w)

        % -- Erickson2022:
        % us: boundary fluxes. The ODE for us is linear and can be built into D.
        % usim, usip: interface fluxes (minus and plus)
        % The ODEs for the interface fluxes have both linear and nonlinear parts.
        %
        % w = [u; v; Psi; us; usim; usip]
        % F = [0; S + penalty_bc*bc_data + penalty_fault*fault_fun; g; 0; nonlinear part of evolution for interface u*]
        % w_t = D*w + F(w)
            

        % ---- Component operators -----   
        I = speye(n);
        Eu = [I, sparse(n,n), sparse(n,1), sparse(n,ns), sparse(n,nsim), sparse(n,nsip)];
        Ev = [sparse(n,n), I, sparse(n,1), sparse(n,ns), sparse(n,nsim), sparse(n,nsip)];
        EPsi = [sparse(1,n), sparse(1,n), 1, sparse(1,ns), sparse(1,nsim), sparse(1,nsip)];
        Eus = [sparse(ns,n), sparse(ns,n), sparse(ns,1), speye(ns), sparse(ns,nsim), sparse(ns, nsip)];
        Eusim = [sparse(nsim,n), sparse(nsim,n), sparse(nsim,1), sparse(nsim,ns), speye(nsim), sparse(nsim,nsip)];
        Eusip = [sparse(nsip,n), sparse(nsip,n), sparse(nsip,1), sparse(nsip,ns), sparse(nsip,nsim), speye(nsip)];

        % ---- Construct scheme for multiblock Laplace including the closures -----
        D_laplace = mbDiffOp.D + closure_bc + closure_fault;

        %---- Build system matrix-----
        D = Ev'*D_laplace*Eu;
        D = D + Eu'*I*Ev;   

        % ---- Initial conditions -----
        if isempty(initialcondition)
            u0 = zeros(n,1);
            v0 = zeros(n,1);
            Psi0 = 0;
        else
            u0 = multiblock.evalOn(mbGrid,initialcondition.u);
            v0 = multiblock.evalOn(mbGrid,initialcondition.v);
            Psi0 = initialcondition.Psi;
        end
        us0 = zeros(ns,1);
        usim0 = zeros(nsim,1);
        usip0 = zeros(nsip,1);

        w0 = [u0; v0; Psi0; us0; usim0; usip0];
        % ----------------------------------------------

        %--- Add special erickson2022 contributions to BC
        switch bc_method
        case 'erickson2022'
            %---- Unpack erickson struct ------
            % Block-diagonal matrices
            R = elastic.helpers.cell_row_to_diag_blockmatrix(erickson2022_struct.R);
            Z = elastic.helpers.cell_row_to_diag_blockmatrix(erickson2022_struct.Z);
            gamma = elastic.helpers.cell_row_to_diag_blockmatrix(erickson2022_struct.gamma);

            % Tall and skinny matrices, horizontally stacked
            tau_star_penalty = blockmatrix.toMatrix(erickson2022_struct.tau_star_penalty);
            u_star_penalty = blockmatrix.toMatrix(erickson2022_struct.u_star_penalty);
            tau_op = blockmatrix.toMatrix(erickson2022_struct.tau_op);
            e = blockmatrix.toMatrix(erickson2022_struct.e);

            ns = length(R);
            IR = eye(size(R));
            %------------------------------------

            % ----Build operators for fluxes ----

            % Contributions to tau* from u, u_t and u*.
            taus_u = (IR - R)/2*tau_op' + gamma/2*(R - IR)*e';
            taus_v = Z*(R - IR)/2*e';
            taus_us = -gamma/2*(R - IR);

            % Contributions to u*_t from u, u_t and u_*
            ust_u = -Z\((IR + R)/2*tau_op' - gamma/2*(IR + R)*e');
            ust_v = (R + IR)/2*e';
            ust_us = -Z\(gamma/2*(IR + R));
            
            %- Add contributions to system matrix-         
            % SATs with tau*
            D = D + Ev'*tau_star_penalty*(taus_u*Eu + taus_v*Ev + taus_us*Eus);

            % SATs with u*
            D = D + Ev'*u_star_penalty*Eus;

            % Evolution of u*
            D = D + Eus'*(ust_u*Eu + ust_v*Ev + ust_us*Eus);
        end

        % ----  Add special erickson2022 contributions to fricitonal interface -----
        switch ic_method
        case 'erickson2022'

            % Compute impedances ( Z = sqrt(rho*mu) ) at interface
            rho_vec = multiblock.evalOn(mbGrid, material.rho);
            mu_vec = multiblock.evalOn(mbGrid, material.mu);
            rhom = em'*rho_vec;
            mum = em'*mu_vec;
            Zm = sqrt(rhom.*mum);
            rhop = ep'*rho_vec;
            mup = ep'*mu_vec;
            Zp = sqrt(rhop.*mup);
            eta = Zp.*Zm./(Zp + Zm);

            %--- Setup nonlinear function that needs to be solved for V* ----
            F = @(t, V, Psi) friction.funs.tau(t, V, Psi, friction.params.a);
            Finv = @(t, tau, Psi) friction.funs.tauinv(t, tau, Psi, friction.params.a);
            % Last term is F(V), but extra layers of function calls slow down execution.
            nonlin_solve_fun = @(t, V, tau_l, Psi) -eta*V + tau_l - friction.funs.tau(t, V, Psi, friction.params.a);
            %----------------------------------------------------------------

            % Bracket for rootfinding
            bracket = @(t, tau_l, Psi) (tau_l >=0 )*[0, Finv(t, tau_l, Psi)] + (tau_l < 0)*[Finv(t, tau_l, Psi), 0]; 

            % Function that computes V* given tau_l and Psi (via nonlinear solve). 
            V_star_from_tau_l = @(t, tau_l, Psi) elastic.helpers.vectorBisection(@(V) nonlin_solve_fun(t, V, tau_l, Psi), bracket(t, tau_l, Psi), 1e-6);

            %---- Operators that compute characteristic variables ----
            % Traction operators, with sign correction
            Tm = (-1 + 2*(bidm{2}=='r')) * mbDiffOp.getBoundaryOperator('d', bidm)*mum;
            Tp = (-1 + 2*(bidp{2}=='r')) * mbDiffOp.getBoundaryOperator('d', bidp)*mup;

            % Penalty strengths
            gammam = 1.1/mbDiffOp.diffOps{1}.gamm*mum;
            gammap = 1.1/mbDiffOp.diffOps{1}.gamm*mup;
            
            Wm = Zm*em'*Ev - Tm'*Eu - gammam*(Eusim - em'*Eu);
            Wp = Zp*ep'*Ev - Tp'*Eu - gammap*(Eusip - ep'*Eu);
            %---------------------------------------------------------

            % Function that computes tau_l given the entire solution vector
            tau_l_fun = @(t, U) (Zp.*(Wm*U) - Zm.*(Wp*U)) ./ (Zp + Zm);
            if isfield(friction.funs, 'tau_L')
                tau_l_fun = @(t, U) (Zp.*(Wm*U) - Zm.*(Wp*U)) ./ (Zp + Zm) - friction.funs.tau_L(t);
            end

            V_star = @(t, U) V_star_from_tau_l(t, tau_l_fun(t, U), EPsi*U);
            obj.V_star = V_star;
            taus_p = @(t, U) F(t, V_star(t, U), EPsi*U); % Use that tau_m = -tau_p to avoid unnecessary solves.
            if isfield(friction.funs, 'tau_L')
                taus_p = @(t, U) F(t, V_star(t, U), EPsi*U) + friction.funs.tau_L(t);
            end


            %---- ODEs for us_p and us_m --------------------------------
            % Linear parts
            D = D + Eusim'*(em'*Ev - Zm\(Tm'*Eu + gammam*(Eusim - em'*Eu)));
            D = D + Eusip'*(ep'*Ev - Zp\(Tp'*Eu + gammap*(Eusip - ep'*Eu)));

            % Nonlinear parts
            interface_ode = (Eusip'*inv(Zp) - Eusim'*inv(Zm)); % To be multiplied by tau*_+(t, U)
            %--------------------------------------------------------------

            %---- Penalties ------
            H = mbDiffOp.H;
            rho_mat = spdiag(rho_vec);

            % Terms with u* are linear
            D = D - Ev' * inv(rho_mat*H) * Tm*Eusim;
            D = D - Ev' * inv(rho_mat*H) * Tp*Eusip;

            % Terms with tau* are nonlinear
            interface_penalties = Ev'* inv(rho_mat*H)*(ep-em); % To be multiplied by tau*_+(t, U)
            %------------------------------------

            obj.S_friction = @(t, U) (interface_ode + interface_penalties)*taus_p(t, U);
        end

        % ---- Set properties -----
        obj.grid = mbGrid;
        obj.order = order;
        obj.material = material;
        obj.interpolate_data = interpolate_data;
        obj.ic_method = ic_method;
        
        % Operators for full domain
        obj.D = D;
        obj.H = mbDiffOp.H;
        obj.Eu = Eu;
        obj.Ev = Ev;
        obj.EPsi = EPsi;
        obj.Eus = Eus;
        obj.Eusim = Eusim;
        obj.Eusip = Eusip;


        % Initial condition
        obj.w0 = w0;
        
        % Outer boundary
        obj.boundary_data_fun = boundary_data_fun;
        
        % Friction properties
        obj.friction = friction;
        obj.e_fault = e_fault;
        obj.penalty_fault = penalty_fault;
        
        % Point sources
        obj.sources = sources;
        obj.dirac_deltas = dirac_deltas;
        
    end
    
    % Prints some info about the discretisation
    function printInfo(obj)
        fprintf('Name: %s\n',obj.name);
        fprintf('Size: %d\n',obj.size());
    end
    
    % Return the number of DOF
    function n = size(obj)
        n = obj.grid.nPoints;
    end
    
    function V = slipRate(obj, v)
        V = diff(obj.e_fault'*v);
    end

    function obj = setFwdFaultTraction(obj)
        penalty = obj.penalty_fault;
        tau = obj.friction.funs.tau;
        a = obj.friction.params.a;
        
        obj.fault_traction_fwd = @(t, V, Psi) penalty*[ tau(t,  V, Psi, a); ...
                                                        -tau(t, V, Psi, a)];
    end
    
    function obj = setAdjFaultTraction(obj)
        penalty = obj.penalty_fault;
        tau_V = obj.friction.data.tau_V;
        g_V = obj.friction.data.g_V;

        obj.fault_traction_adj = @(i_t, V_adj, Psi_adj) penalty*[ tau_V(i_t)*V_adj + g_V(i_t)*Psi_adj; ...
                                                                 -tau_V(i_t)*V_adj - g_V(i_t)*Psi_adj];
    end
    
    function obj = setInterpAdjFaultTraction(obj)
        penalty = obj.penalty_fault;
        tau_V = obj.friction.data.tau_V;
        g_V = obj.friction.data.g_V;
        t = obj.friction.data.t;
        
        tau_V_pp = spline(t, tau_V);
        g_V_pp = spline(t, g_V);
        
        obj.fault_traction_adj = @(t, V_adj, Psi_adj) penalized_fault_traction(t, V_adj, Psi_adj, tau_V_pp, g_V_pp, penalty);
        
        % Helper function for constructing the penalized data
        function r = penalized_fault_traction(t, V_adj, Psi_adj, tau_V_pp, g_V_pp, penalty)
            val = ppval(tau_V_pp, t)*V_adj + ppval(g_V_pp, t)*Psi_adj;
            r = penalty*[val; -val];
        end
    end
    
    function obj = setFwdStateEvolution(obj)
        g = obj.friction.funs.g;
        a = obj.friction.params.a;
        b = obj.friction.params.b;
        obj.state_evo_fwd = @(V, Psi) g(V, Psi, a, b);
    end
    
    function obj = setAdjStateEvolution(obj)
        tau_Psi = obj.friction.data.tau_Psi;
        g_Psi = obj.friction.data.g_Psi;
        obj.state_evo_adj = @(i_t, V_adj, Psi_adj) g_Psi(i_t)*Psi_adj + tau_Psi(i_t)*V_adj;
    end
    
    function obj = setInterpAdjStateEvolution(obj)
        tau_Psi = obj.friction.data.tau_Psi;
        g_Psi = obj.friction.data.g_Psi;
        t = obj.friction.data.t;
        
        tau_Psi_pp = spline(t, tau_Psi);
        g_Psi_pp = spline(t, g_Psi);

        obj.state_evo_adj = @(t, V_adj, Psi_adj) ppval(g_Psi_pp, t)*Psi_adj + ppval(tau_Psi_pp, t)*V_adj;
    end
    
    function obj = setContPointSources(obj)
        if isempty(obj.sources)
            F = [];
        else
            source_funs = obj.sources.funs;
            deltas = obj.dirac_deltas;
            ns = numel(deltas);
            rho = multiblock.evalOn(obj.grid,obj.material.rho);
            Rho_inv = spdiag(1./rho);
            F = @(t) 0;
            for i = 1:ns
                F_i = source_funs{i};
                F = @(t) F(t) + Rho_inv*F_i(t)*deltas{i};
            end
        end
        obj.S_cont = F;
    end
    
    function obj = setDiscrPointSources(obj)
        if isempty(obj.sources)
            F = [];
        else
            source_data = obj.sources.data;
            deltas = obj.dirac_deltas;
            ns = numel(deltas);
            rho = multiblock.evalOn(obj.grid,obj.material.rho);
            Rho_inv = spdiag(1./rho);
            F = 0;
            for i = 1:ns
                data_i = source_data{i};
                if iscolumn(data_i)
                    data_i = transpose(data_i);
                end
                F = F + Rho_inv*kron(data_i, deltas{i});
            end
            F = @(i_t) F(:,i_t);
        end
        obj.S_discr = F;
    end
        
    % Assembles the time-continuous forcing terms to a single function F_cont(t,w)
    function F_cont = assembleContForcing(obj)
        F_cont = [];
        V = @obj.slipRate;
        % Sources
        if ~isempty(obj.S_cont)
            F_cont = @(t, w) obj.Ev'*obj.S_cont(t);
        end
        % Fwd fault traction
        if ~isempty(obj.fault_traction_fwd)
            if ~isempty(F_cont)
                F_cont = @(t, w) F_cont(t, w) + obj.Ev'*obj.fault_traction_fwd(t, V(obj.Ev*w), obj.EPsi*w);
            else
                F_cont = @(t, w) obj.Ev'*obj.fault_traction_fwd(t, V(obj.Ev*w), obj.EPsi*w);
            end
        end
        % Fwd erickson2022 friction
        if ~isempty(obj.S_friction)
            if ~isempty(F_cont)
                F_cont = @(t, w) F_cont(t, w) + obj.S_friction(t, w);
            else
                F_cont = @(t, w) obj.S_friction(t, w);
            end
        end
        % BC
        if ~isempty(obj.boundary_data_fun)
            if ~isempty(F_cont)
                F_cont = @(t, w) F_cont(t, w) + obj.Ev'*obj.boundary_data_fun(t, obj.Eu*w, obj.Ev*w);
            else
                F_cont = @(t, w)  obj.Ev'*obj.boundary_data_fun(t, obj.Eu*w, obj.Ev*w);
            end
        end
        % Fwd state evolution
        if ~isempty(obj.state_evo_fwd)
            % Erickson2022 interface
            if ~isempty(obj.V_star)
                Vs = obj.V_star;
                if ~isempty(F_cont)
                    F_cont = @(t,w) F_cont(t,w) + obj.EPsi'*obj.state_evo_fwd(Vs(t, w), obj.EPsi*w);
                else
                    F_cont = @(t,w) obj.EPsi'*obj.state_evo_fwd(Vs(t, w), obj.EPsi*w);
                end

            % standard interface
            else
                if ~isempty(F_cont)
                    F_cont = @(t,w) F_cont(t,w) + obj.EPsi'*obj.state_evo_fwd(V(obj.Ev*w), obj.EPsi*w);
                else
                    F_cont = @(t,w) obj.EPsi'*obj.state_evo_fwd(V(obj.Ev*w), obj.EPsi*w);
                end
            end
        end

        % Interpolated forcing is treated as time-continuous
        if obj.interpolate_data 
            F_discr_interp = obj.assembleDiscrForcing();
            if ~isempty(F_discr_interp)
                if ~isempty(F_cont)
                    F_cont = @(t,w) F_cont(t,w) + F_discr_interp(t,w);
                else
                    F_cont = F_discr_interp;
                end
            end
        end
    end

    % Assembles the time-discrete forcing terms to a single function F_discr(i_t,w)
    function F_discr = assembleDiscrForcing(obj)
        F_discr = [];
        V = @obj.slipRate;
        % Sources
        if ~isempty(obj.S_discr)
            F_discr =  @(i_t,w) obj.Ev'*obj.S_discr(i_t);
        end
       % Adjoint fault traction
        if ~isempty(obj.fault_traction_adj)
            if ~isempty(F_discr)
                F_discr = @(i_t,w) F_discr(i_t,w) + obj.Ev'*obj.fault_traction_adj(i_t, V(obj.Ev*w), obj.EPsi*w);
            else
                F_discr = @(i_t,w) obj.Ev'*obj.fault_traction_adj(i_t, V(obj.Ev*w), obj.EPsi*w);
            end
        end
        % Adjoint state evolution
        if ~isempty(obj.state_evo_adj)
            if ~isempty(F_discr)
                F_discr = @(i_t,w) F_discr(i_t,w) + obj.EPsi'*obj.state_evo_adj(i_t, V(obj.Ev*w), obj.EPsi*w);
            else
                F_discr = @(i_t,w) obj.EPsi'*obj.state_evo_adj(i_t, V(obj.Ev*w), obj.EPsi*w);
            end
        end
    end
    
    % Returns a timestepper for integrating the discretisation in time
    %     method is a string that states which timestepping method should be used.
    %          The implementation should switch on the string and deliver
    %          the appropriate timestepper. It should also provide a default value.
    %     time_align is a time that the timesteps should align with so that for some
    %                integer number of timesteps we end up exactly on time_align
    function [ts, N] = getTimestepper(obj, method, time_align, k)
        default_arg('method','struct');
        default_field(method,'adaptive',false);
        default_field(method,'order',4);
        default_arg('time_align',[]);
        default_arg('k',[]);
        
        t = 0;
        if method.adaptive
            switch method.order
            case 3
                cfl = 0.2;
            case 5
                cfl = 0.2;
            otherwise
                error('Adaptive RK of order %d not implemented',method.order);
            end
        else
            switch method.order
            case 3
                cfl = 0.2;
            case 4
                cfl = 0.2;
            case 5
                cfl = 0.2;
            case 6
                cfl = 0.1;
            otherwise
                error('RK of order %d not implemented',method.order);
            end
        end
        
        if isempty(k)
            k = obj.getTimestep(method,cfl);
        end
        
        if ~isempty(time_align)
            [k, N] = alignedTimestep(k, time_align);
        end
        F_cont = obj.assembleContForcing();
        F_discr = obj.assembleDiscrForcing();
        
        if method.adaptive
            default_error_check = @(vNew, vStar) norm(vNew - vStar,inf);
            default_field(method,'errorCheckCallback',default_error_check);
            default_field(method,'reportRetry',false);
            if ~isempty(F_cont)
                F = @(t,w) obj.D*w + F_cont(t,w);
            else
                F = @(t,w) obj.D*w;
            end
            ts = time.EmbeddedRungeKutta(F, k, t, obj.w0, method.order, method.rtol, method.errorCheckCallback, [], method.reportRetry);
        else
            ts = time.ExplicitRungeKuttaDiscreteData(...
            obj.D, F_cont, F_discr, k, t, obj.w0, method.order);
        end
    end
    
    
    function k = getTimestep(obj, method, cfl)
        default_arg('cfl',[]);
        if isempty(cfl)
            error('Specify CFL')
        end
        
        h  = min(obj.grid.grids{1}.scaling(),obj.grid.grids{2}.scaling());
        k =  cfl*h;
    end
    
    function r = getTimeSnapshot(obj, ts)
        if ts == 0
            r.t = 0;
            r.u = obj.Eu*obj.w0;
            r.v = obj.Ev*obj.w0;
            r.Psi = obj.EPsi*obj.w0;
            r.us = obj.Eus*obj.w0;
            r.usim = obj.Eusim*obj.w0;
            r.usip = obj.Eusip*obj.w0;
            % r.Vs = obj.V_star(obj.w0);
            return
        end
        r.t = ts.t;
        w = ts.getV();
        r.u = obj.Eu*w;
        r.v = obj.Ev*w;
        r.Psi = obj.EPsi*w;
        r.usim = obj.Eusim*w;
        r.usip = obj.Eusip*w;
        % r.Vs = obj.V_star(w);
    end
    
    % Sets up a plot of the discretisation
    %     update is a function_handle accepting a timestepper that updates the plot to the
    %            state of the timestepper
    function [update,figure_handle] = setupPlot(obj, type)
        figure_handle = figure();
        x = obj.grid.points();
        default_arg('type',struct);
        default_field(type,'plot_variables','trajectory')
        default_field(type,'axis_u',[x(1) x(end), -1,1]);
        default_field(type,'axis_v',[x(1) x(end), -1,1]);
        
        function plot_comp(comp,axis_lim)
            plot(x,comp,'lineWidth',2);
            %axis(axis_lim);
        end
        
        function plot_uv(r)
            subplot(2,1,1);
            plot_comp(r.u,type.axis_u);
            xlabel('$x$','interpreter','latex');
            ylabel('$u$','interpreter','latex');
            subplot(2,1,2);
            plot_comp(r.v,type.axis_v);
            xlabel('$x$','interpreter','latex');
            ylabel('$v$','interpreter','latex');
        end
        
        t_plot  = [];
        Psi_plot  = [];
        V_plot  = [];
        usim_plot = [];
        Vs_plot = [];
        function plot_all(r)
            t_plot = [t_plot, r.t];
            Psi_plot = [Psi_plot, r.Psi];
            V_plot = [V_plot, obj.slipRate(r.v)];
            subplot(4,1,1);
            plot_comp(r.u,type.axis_u);
            xlabel('$x$','interpreter','latex');
            ylabel('$u$','interpreter','latex');
            subplot(4,1,2);
            plot_comp(r.v,type.axis_v);
            xlabel('$x$','interpreter','latex');
            ylabel('$v$','interpreter','latex');
            subplot(4,1,3);
            plot(t_plot, Psi_plot, 'linewidth',  2);
            xlabel('$t$','interpreter','latex');
            ylabel('$\Psi$','interpreter','latex');
            subplot(4,1,4);
            plot(t_plot, V_plot, 'linewidth',  2);
            xlabel('$t$','interpreter','latex');
            ylabel('$V$','interpreter','latex');
        end

        function plot_with_traction_trajectory(r)
            t_plot = [t_plot, r.t];
            Psi_plot = [Psi_plot, r.Psi];
            V_plot = [V_plot, obj.slipRate(r.v)];
            log_abs_V_plot = log(abs(V_plot));
            a = obj.friction.params.a;
            b = obj.friction.params.b;
            sigma0 = obj.friction.params.sigma0;
            switch obj.ic_method
            case 'standard'
                tau_plot = 1/sigma0.*(obj.friction.funs.tau(t_plot,V_plot,Psi_plot,a) + obj.friction.funs.tau_L(t_plot));
            case 'erickson2022'
                tau_plot = 1/sigma0.*obj.friction.funs.tau(t_plot,V_plot,Psi_plot,a);
            end
            f_ss_plot = obj.friction.funs.f_ss(V_plot,a,b);
            subplot(5,1,1);
            plot_comp(r.u,type.axis_u);
            xlabel('$x$','interpreter','latex');
            ylabel('$u$','interpreter','latex');
            subplot(5,1,2);
            plot_comp(r.v,type.axis_v);
            xlabel('$x$','interpreter','latex');
            ylabel('$v$','interpreter','latex');
            subplot(5,1,3);
            plot(t_plot, Psi_plot, 'linewidth',  2);
            xlabel('$t$','interpreter','latex');
            ylabel('$\Psi$','interpreter','latex');
            subplot(5,1,4);
            plot(t_plot, V_plot, 'linewidth',  2);
            xlabel('$t$','interpreter','latex');
            ylabel('$V$','interpreter','latex');
            subplot(5,1,5);
            plot(log_abs_V_plot,tau_plot,log_abs_V_plot, f_ss_plot , 'linewidth',  2);
            xlabel('$ln(|V|)$','interpreter','latex');
            ylabel('$\tau/\sigma_0$','interpreter','latex');
        end

        
        switch type.plot_variables
        case 'u'
            update = @(r)plot_comp(r.u,type.axis_u);
        case 'v'
            update = @(r)plot_comp(r.v,type.axis_v);
        case 'uv'
            update = @(r)plot_uv(r);
        case 'all'
            update = @(r)plot_all(r);
        case 'trajectory'
            update = @(r)plot_with_traction_trajectory(r);
        end
    end
    
    % Compare two functions u and v in the discrete l2 norm.
    function e = compareSolutions(obj, u, v)
        e = [];
        error('Not implemented')
    end
    
end
end