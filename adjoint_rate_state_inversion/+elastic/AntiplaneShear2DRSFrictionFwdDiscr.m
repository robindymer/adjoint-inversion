classdef AntiplaneShear2DRSFrictionFwdDiscr < noname.Discretization
% Discretizes the 1D anti-plane shear problem
% 
% u_tt = 1/rho ( (mu*u_x)_x + S)
% Psi_t = g(V, Psi)
% 
% for velocity potential u, and state variable Psi, where
% S is a dirac delta source term, and g a state evolution law.
properties
    name   = "AntiplaneShear2DRSFrictionFwdDiscr"
    description  = 'Discretizes the 1D wave equation for antiplane shear with a rate- and state friction interface. Solves for displacement u, velocity v and state variable Psi'
    order                % Order of accuracy
    grid                 % Grid
    material             % Material parameters
    
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

    % DiffOp
    mbDiffOp

    % Component operators for first-order in time system
    % E.u*w = u, Eu'*u = [u; 0; 0]; 
    % E.v*w = v, Ev'*v = [0; v; 0];
    % E.Psi*w = Psi, EPsi'*Psi = [0; 0; Psi];
    E

    % Initial conditions
    w0 % w0 = [u0; v0; Psi0]
        
    % Boundary functions and operators
    penalized_data % Boundary data forcing

    % Friction and operators
    friction             % Struct holding parameters and friction functions/data
    e_fault_m            % Fault restriciton operator on minus side
    e_fault_p            % Fault restriciton operator on plus side
    penalty_fault        % Fault data penalty operator
    fault_traction       % Fault traction forcing
    state_evo            % State evolution forcing
    V_star
    tau_l_fun
    eta
    tau_L

    % Sources
    sources         % Struct holding parameters source functions/data
    S               % Continuous source term
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
    %       friction.state_evo_fun  -- function @(Psi,ut)
    %   sources    --      struct for point sources
    %       sources.x           -- array with source positions
    %       sources.source_fun  -- cell array with source functions @(t)
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

    function obj = AntiplaneShear2DRSFrictionFwdDiscr(opset, domain, m, order, material, bc, friction,...
                                                      sources, initialcondition)
        
        % ---- Default values --------------
        default_arg('opset',@sbp.D2Variable);
        default_arg('m',101);
        default_arg('order',4);
        
        default_arg('material', struct);
        default_field(material, 'rho', {@(x,y) 0*x + 1, @(x,y) 0*x + 1});
        default_field(material, 'mu',  {@(x,y) 0*x + 1, @(x,y) 0*x + 1});
        
        default_arg('sources',[]);
        default_arg('initialcondition',[]);
    
        
        % Get grid and multiblock diff op
        if isequal(opset,@sbp.D2Nonequidistant) % Boundary-optimized operators 
            mbGrid = domain.getGrid(m, 'boundaryopt', order, 'acc');
        else
            mbGrid = domain.getGrid(m, 'equidist');
        end

        nBlocks = mbGrid.nBlocks();
        doParam = cell(nBlocks, 1);
        for i = 1:nBlocks
            if iscell(material.rho)
                doParam{i} = {@(x,y) 1./material.rho{i}(x,y), material.mu{i}, opset};
            elseif isa(rho, 'function_handle')
                doParam{i} = {@(x,y) 1./material.rho(x,y), material.mu, opset};
            else
                error('Inconsistent format for mu and rho.');
            end
        end

        mbDiffOp = multiblock.DiffOp(@scheme.LaplaceCurvilinearNew, mbGrid, order, doParam);
        
        % Spatial operators
        [D, E, bc_struct, ic_struct] = elastic.discrs.antiplaneshear.discrOps2D(domain, mbGrid, mbDiffOp, material, bc);
        
        % Sources
        dirac_deltas = elastic.discrs.antiplaneshear.diracDelta2D(sources, mbGrid, order, {opset, opset});
        
        % Problem sizes
        [n, ns, nsip, nsim] = elastic.discrs.antiplaneshear.nunknowns2D(domain, mbGrid, mbDiffOp, bc);
        
        % ---- Fault interface -----
        bidm = domain.boundaryGroups.fault_minus;
        X_fault = mbGrid.getBoundary(bidm);

         % ---- Functions for non-linear solve for target slip velocity V* -----
        Zm = ic_struct.Zm;
        Zp = ic_struct.Zp;
        Wm = ic_struct.Wm;
        Wp = ic_struct.Wp;
        eta = Zp.*Zm./(Zp + Zm);
        % Function that computes tau_l given the entire solution vector
        tau_l_fun = @(t, U) (Zp.*(Wm*U) - Zm.*(Wp*U)) ./ (Zp + Zm);
        if isfield(friction.funs, 'tau_L')
            tau_l_fun = @(t, U) (Zp.*(Wm*U) - Zm.*(Wp*U)) ./ (Zp + Zm) - friction.funs.tau_L(t, X_fault(:,1), X_fault(:,2));
            % Store tau_L for setting up RHS
            tau_L = @(t) friction.funs.tau_L(t, X_fault(:,1), X_fault(:,2));
        else
            tau_L = [];
        end

        % ---- Rate-state friction parameters ----- %
        names = fieldnames(friction.rsParams);
        for i = 1:numel(names)
            name = names{i};
            p = friction.rsParams.(name);
            if isa(p,'function_handle')
                p_fun = p;
            elseif isscalar(p)
                p_fun = @(x,y) x*0 + p;
            else
                error('Friction parameters must be either scalars of function handles');
            end
            friction.rsParams.(name) = elastic.helpers.evalOnLine(X_fault, p_fun);
        end
        
        % TODO: Compute from stress tensor.
        gidm = bidm{1}{1};
        boundary = bidm{1}{2};
        switch boundary
        case 'e'
            n = mbDiffOp.diffOps{gidm}.n_e{1};
        case 'w'
            n = mbDiffOp.diffOps{gidm}.n_w{1};
        case 'n'
            n = mbDiffOp.diffOps{gidm}.n_n{2};
        case 's'
            n = mbDiffOp.diffOps{gidm}.n_s{2};
        end
        friction.rsParams.tau0 = n*friction.rsParams.tau0;


         % ---- Initial conditions -----
         if isempty(initialcondition)
            u0 = zeros(n,1);
            v0 = zeros(n,1);
            Psi0 = zeros(nsim, 1);
        else
            u0 = multiblock.evalOn(mbGrid,initialcondition.u);
            v0 = multiblock.evalOn(mbGrid,initialcondition.v);
            Psi0 = elastic.helpers.evalOnLine(X_fault, initialcondition.Psi);
        end
        us0 = zeros(ns,1);
        usim0 = zeros(nsim,1);
        usip0 = zeros(nsip,1);
        %Vs0 = zeros(nsim,1);

        w0 = [u0; v0; Psi0; us0; usim0; usip0];%; Vs0];

        % ---- Set properties -----
        obj.grid = mbGrid;
        obj.order = order;
        obj.material = material;
        
        % Operators for full domain
        obj.D = D;
        obj.H = mbDiffOp.H;
        obj.E = E;
        obj.mbDiffOp = mbDiffOp;

        % Initial condition
        obj.w0 = w0;
        
        % Outer boundary
        obj.penalized_data = bc_struct.penalized_data;
        
        % Friction properties
        obj.friction = friction;
        obj.e_fault_m = ic_struct.em;
        obj.e_fault_p = ic_struct.ep;
        obj.penalty_fault = ic_struct.penalty;
        obj.tau_l_fun = tau_l_fun;
        obj.eta = eta;
        obj.tau_L = tau_L;
        
        % Point sources
        obj.sources = sources;
        obj.dirac_deltas = dirac_deltas;

        % Set fault interface, state evolution and point source
        % forcing properties here to avoid partial initialization.
        obj = obj.setFaultTraction();
        obj = obj.setStateEvolution();
        obj = obj.setPointSources();
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
    
    function V = fault_jump(obj, v)
        V = (obj.e_fault_p-obj.e_fault_m)'*v;
    end

    % TODO: 
    % One could consider having cell arrays for friction law parameter and
    % state eq parameters.
    % The friction law would have the interface F(V,Psi,param1,param2,...)
    % which in the code would then be evaluated as F(V,Psi,param_carray{:}).
    % Similarly for the state eq. In this way the functions setting fault 
    % traction and state eq would have general interfaces
    function obj = setFaultTraction(obj)  
        F = obj.friction.funs.F;
        rs = obj.friction.rsParams;
        Finv = obj.friction.funs.Finv; % NOTE! Does not include tau0. TODO: Fix dep on tau0!
        nonlin_solve_fun = obj.friction.funs.nonlin_solve_fun; % NOTE! Does not include tau0. TODO: Fix dep on tau0!
        tol = 1e-13;
        
        % Operators 
        penalty = obj.penalty_fault;
        E = obj.E;
        eta = obj.eta;

        %--- Setup nonlinear function that needs to be solved for V* ----
        % Function that computes V* given tau_l and Psi (via nonlinear solve). 
        V_star_from_tau_l = @(tau_l, Psi) elastic.helpers.vectorBisection( @(V) nonlin_solve_fun(V, Psi, eta, tau_l, rs.a, rs.sigma0, rs.V0), ...
                                                                                -Finv(tau_l, Psi, rs.a, rs.sigma0, rs.V0) .* [tau_l >=0, tau_l < 0], ... % Bracket
                                                                                tol); % Tolerance
        obj.V_star = @(t, U) V_star_from_tau_l(obj.tau_l_fun(t, U) - rs.tau0, E.Psi*U); % NOTE: here we subtract by tau0 to get the correct behaviour. TODO: Fix dep on tau0
        % Fault traction function
        if ~isempty(obj.tau_L)
            obj.fault_traction = @(t, U, Vs) penalty*(F(Vs, E.Psi*U, rs.a, rs.sigma0, rs.V0, rs.tau0) - obj.tau_L(t));
        else 
            obj.fault_traction = @(t, U, Vs) penalty*(F(Vs, E.Psi*U, rs.a, rs.sigma0, rs.V0, rs.tau0));
        end
    end
        
    function obj = setStateEvolution(obj)
        G = obj.friction.funs.G;
        rs = obj.friction.rsParams;
        E = obj.E;
        obj.state_evo = @(U, Vs) E.Psi'*G(Vs, E.Psi*U, rs.a, rs.b, rs.f0, rs.V0, rs.D_c);
    end
    
    function obj = setPointSources(obj)
        if isempty(obj.sources)
            F = [];
        else
            source_funs = obj.sources.funs;
            deltas = obj.dirac_deltas;
            ns = numel(deltas);
            rho = multiblock.evalOn(obj.grid,obj.material.rho);
            Rho_inv = spdiag(1./rho);
            E = obj.E;
            F = @(t) 0;
            for i = 1:ns
                F_i = source_funs{i};
                F = @(t) F(t) + Rho_inv*F_i(t)*deltas{i};
            end
            F = @(t) E.v'*F(t);
        end
        obj.S = F;
    end
            
    % Assembles the forcing terms to a single function F(t, U)
    function F_tot = assembleForcing(obj)
        has_source = ~isempty(obj.S);
        has_data = ~isempty(obj.penalized_data);
        
        % Assemble rhs F(t, U, Vs)
        if has_source && has_data
            F = @(t, U, Vs) obj.fault_traction(t, U, Vs) + obj.state_evo(U, Vs) + obj.S(t) + obj.penalized_data(t);
        elseif has_source
            F = @(t, U, Vs) obj.fault_traction(t, U, Vs) + obj.state_evo(U, Vs) + obj.S(t);
        elseif has_data
            F = @(t, U, Vs) obj.fault_traction(t, U, Vs) + obj.state_evo(U, Vs) + obj.penalized_data(t);
        else
            F = @(t, U, Vs) obj.fault_traction(t, U, Vs) + obj.state_evo(U, Vs);
        end
        
        % % Compute V_star once, pass along to other forcing terms
        % % and store in returned vector r
        % function r = F_precomp_Vs(t,U)
        %     Vs = obj.V_star(t, U);
        %     r = obj.E.Vs'*Vs + F(t, U, Vs);
        % end
        % F_tot = @F_precomp_Vs;
        % Compute V_star once, pass along to other forcing terms
        F_tot = @(t, U) F(t, U, obj.V_star(t, U));
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
                cfl = 0.1;
            case 5
                cfl = 0.1;
            otherwise
                error('Adaptive RK of order %d not implemented',method.order);
            end
        else
            switch method.order
            case 3
                cfl = 0.1;
            case 4
                cfl = 0.1;
            case 5
                cfl = 0.1;
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
        F = obj.assembleForcing();
        
        if method.adaptive
            default_error_check = @(vNew, vStar) norm(vNew - vStar,inf);
            default_field(method,'errorCheckCallback',default_error_check);
            default_field(method,'reportRetry',false);
            if ~isempty(F)
                F = @(t,w) obj.D*w + F(t,w);
            else
                F = @(t,w) obj.D*w;
            end
            ts = time.EmbeddedRungeKutta(F, k, t, obj.w0, method.order, method.rtol, method.errorCheckCallback, [], method.reportRetry);
        else
            ts = time.ExplicitRungeKuttaDiscreteData(obj.D, F, [], k, t, obj.w0, method.order);
        end
    end
    
    
    function k = getTimestep(obj, method, cfl)
        default_arg('cfl',[]);
        if isempty(cfl)
            error('Specify CFL')
        end
        
        nBlocks = obj.grid.nBlocks;
        k = inf;
        for i = 1:nBlocks
            a = diag(obj.mbDiffOp.diffOps{i}.a);
            b = diag(obj.mbDiffOp.diffOps{i}.b);

            % Physical wave speed
            v = sqrt(a.*b);

            % Reference grid spacing
            h = obj.mbDiffOp.diffOps{i}.h;
            h_xi = h(1);
            h_eta = h(2);

            x_xi = obj.mbDiffOp.diffOps{i}.x_u;
            x_eta = obj.mbDiffOp.diffOps{i}.x_v;
            y_xi = obj.mbDiffOp.diffOps{i}.y_u;
            y_eta = obj.mbDiffOp.diffOps{i}.y_v;

            % Approximate physical grid spacings
            h1 = sqrt( (x_xi*h_xi).^2 + (y_xi*h_xi).^2 );
            h2 = sqrt( (x_eta*h_eta).^2 + (y_eta*h_eta).^2 );

            k1 = cfl*min(h1./v);
            k2 = cfl*min(h2./v);

            % Time-step suggested by current block
            k_temp = min(k1,k2);

            % Use smallest time steps of all blocks
            k = min(k, k_temp);
        end
    end
    
    function r = getTimeSnapshot(obj, ts)
        if ts == 0
            r.t = 0;
            r.u = obj.E.u*obj.w0;
            r.v = obj.E.v*obj.w0;
            r.Psi = obj.E.Psi*obj.w0;
            r.us = obj.E.us*obj.w0;
            r.usim = obj.E.usim*obj.w0;
            r.usip = obj.E.usip*obj.w0;
            %r.Vs = obj.E.Vs*obj.w0;
            return
        end
        r.t = ts.t;
        w = ts.getV();
        r.u = obj.E.u*w;
        r.v = obj.E.v*w;
        r.Psi = obj.E.Psi*w;
        r.usim = obj.E.usim*w;
        r.usip = obj.E.usip*w;
        r.us = obj.E.us*w;
        %r.Vs = obj.E.Vs*w;
    end
    
    % Sets up a plot of the discretisation
    %     update is a function_handle accepting a timestepper that updates the plot to the
    %            state of the timestepper
    function [update,figure_handle] = setupPlot(obj, type)
        figure_handle = figure();
        x = obj.grid.points();
        default_arg('type',struct);
        default_field(type,'plot_variables','uv')
        default_field(type,'axlims', []);

        
        function S = plot_comp(comp)
            S = multiblock.Surface(obj.grid, comp);
            xlabel('$x$','interpreter','latex');
            ylabel('$y$','interpreter','latex');
            elastic.cmap;
            colorbar;
            shading interp;
        end

        function Surs = setup_v(r)
            Surs = plot_comp(r.v);
            title(sprintf('v, t = %.2e', r.t));
            if ~isempty(type.axlims)
                axis(type.axlims)
            end
        end

        function update_v(r, Surs)
            Surs.ZData = r.v;
            Surs.CData = r.v;
            title(sprintf('v, t = %.2e', r.t));
        end
        
        function Surs = setup_uv(r)
            subplot(2,1,1);
            Sur1 = plot_comp(r.u);
            title(sprintf('u, t = %.2e', r.t));
            axis 
            
            subplot(2,1,2);
            Sur2 = plot_comp(r.v);
            title(sprintf('v'));

            Surs = {Sur1, Sur2};
        end

        function update_uv(r, Surs)
            Surs{1}.ZData = r.u;
            Surs{1}.CData = r.u;
            Surs{2}.ZData = r.v;
            Surs{2}.CData = r.v;
            subplot(2,1,1);
            title(sprintf('u, t = %.2e', r.t));
            if ~isempty(type.axlims)
                axis(type.axlims)
                subplot(2,1,2);
                axis(type.axlims)
            end
        end
        
        function [Surs, Lines] = setup_all(r)
            % V_plot = [V_plot, obj.slipRate(r.v)];
            subplot(4,1,1);
            Sur1 = plot_comp(r.u);
            title(sprintf('u, t = %.2e', r.t));
          
            
            subplot(4,1,2);
            Sur2 = plot_comp(r.v);
            title(sprintf('v'));
            

            subplot(4,1,3);
            line1 = plot(r.Psi, 'linewidth',  2);
            xlabel('number along fault');
            ylabel('$\Psi$','interpreter','latex');

            subplot(4,1,4);
            ylab = '$V^*$';
            %ylab = '$\log(|V^*|)$';
            w = [r.u; r.v; r.Psi; r.us; r.usim; r.usip];% r.Vs];
            Vs = obj.V_star(r.t, w);
            %line2 = plot(log10(abs(Vs)), 'linewidth',  2);
            line2 = plot(Vs, 'linewidth',  2);
            %line2 = plot(log10(abs(r.Vs)), 'linewidth',  2);
            xlabel('number along fault');
            ylabel(ylab,'interpreter','latex');

            Surs = {Sur1, Sur2};
            Lines = {line1, line2};
        end

        function update_all(r, Surs, Lines)
            subplot(4,1,1);
            Surs{1}.ZData = r.u;
            Surs{1}.CData = r.u;
            cmax = max(abs(r.u)) + 1e-10;
            caxis([-cmax, cmax]);
            title(sprintf('u, t = %.2e', r.t));
            if ~isempty(type.axlims)
                axis(type.axlims)
            end

            subplot(4,1,2);
            Surs{2}.ZData = r.v;
            Surs{2}.CData = r.v;
            cmax = max(abs(r.v)) + 1e-10;
            caxis([-cmax, cmax]);
            if ~isempty(type.axlims)
                axis(type.axlims)
            end

            subplot(4,1,3);
            Lines{1}.YData = r.Psi; 

            subplot(4,1,4);
            w = [r.u; r.v; r.Psi; r.us; r.usim; r.usip];% r.Vs];
            Vs = obj.V_star(r.t, w);
            %Lines{2}.YData = log10(abs(Vs));
            Lines{2}.YData = Vs;
            % Lines{2}.YData = log10(abs(r.Vs));
        end

        function plot_with_traction_trajectory(r)
            error('Not implemented');
        end

        
        switch type.plot_variables
        case 'v'
            r = obj.getTimeSnapshot(0);
            Surs = setup_v(r);
            update = @(r)update_v(r, Surs);
        case 'uv'
            r = obj.getTimeSnapshot(0);
            Surs = setup_uv(r);
            update = @(r)update_uv(r, Surs);
        case 'all'
            r = obj.getTimeSnapshot(0);
            [Surs, Lines] = setup_all(r);
            update = @(r)update_all(r, Surs, Lines);
        end
    end
    
    % Compare two functions u and v in the discrete l2 norm.
    function e = compareSolutions(obj, u, v)
        e = [];
        error('Not implemented')
    end
    
end
end