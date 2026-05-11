% TODO: Update
classdef AntiplaneShearRSFrictionSecondOrderFwdDiscr < noname.Discretization
% Discretizes the 1D anti-plane shear problem
% 
% u_tt = 1/rho ( (mu*u_x)_x + S)
% Psi_t = g(V, Psi)
% 
% for velocity potential u, and state variable Psi, where
% S is a dirac delta source term, and g a state evolution law.
properties
    name   = "AntiplaneShearRSFrictionSecondOrderFwdDiscr"
    description  = 'Discretizes the 1D wave equation for antiplane shear with a rate- and state friction interface. Solves for displacement u, velocity v and state variable Psi'
    order                % Order of accuracy
    grid                 % Grid
    material             % Material parameters
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
    % E.u*w = u, Eu'*u = [u; 0; 0]; 
    % E.v*w = v, Ev'*v = [0; v; 0];
    % E.Psi*w = Psi, EPsi'*Psi = [0; 0; Psi];
    E

    % Initial conditions
    w0 % w0 = [u0; v0; Psi0]
        
    % Boundary functions and operators
    boundary_data_fun % Boundary data forcing

    % Friction and operators
    friction             % Struct holding parameters and friction functions/data
    e_fault              % Fault restriciton operator
    penalty_fault        % Fault data penalty operator
    fault_traction       % Fault traction forcing
    state_evo            % State evolution forcing
    V_star
    tau_l_fun
    eta

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

    function obj = AntiplaneShearRSFrictionSecondOrderFwdDiscr(opset, m, lims, order, material, bc, friction,...
                                                    sources, initialcondition)
        
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
        default_field(bc, 'method', 'standard');
        
        default_arg('friction', struct);
        default_field(friction, 'params', []);
        default_field(friction, 'funs', struct);
        default_field(friction.funs, 'tau', []);
        default_field(friction.funs, 'g', []);
        default_field(friction, 'method', 'standard');
         
        default_arg('sources',[]);
        default_arg('initialcondition',[]);
    
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
        
        % Spatial operators
        [D, E, bc_struct, ic_struct] = elastic.discrs.antiplaneshear.discrOps1D(domain, mbGrid, mbDiffOp, material, bc, friction.method);
        
        % Sources
        dirac_deltas = elastic.discrs.antiplaneshear.diracDelta1D(sources, mbGrid, mbDiffOp, order);
        
        % ---- Boundary data -----
        if ~isempty(bc_struct.boundary_data_fun)
            boundary_data_fun = @(t, U) F(t, U) + E.v'*obj.boundary_data_fun(t, E.u*U, E.v*U);
        else
            boundary_data_fun = [];
        end
        
        % ---- Fault interface -----
        switch friction.method
        case 'standard'
            eta = [];
            tau_l_fun = [];
        case 'erickson2022'
            Zm = ic_struct.Zm;
            Zp = ic_struct.Zp;
            Wm = ic_struct.Wm;
            Wp = ic_struct.Wp;
            eta = Zp.*Zm./(Zp + Zm);
            if isfield(friction.funs, 'tau_L')
                tau_l_fun = @(t,U) (Zp.*(Wm*U) - Zm.*(Wp*U)) ./ (Zp + Zm) - friction.funs.tau_L(t);
            else
                tau_l_fun = @(t,U) (Zp.*(Wm*U) - Zm.*(Wp*U)) ./ (Zp + Zm);
            end
        end

         % ---- Initial conditions -----
         [n, ns, nsip, nsim] = elastic.discrs.antiplaneshear.nunknowns(domain, mbGrid, mbDiffOp, bc.method, friction.method);
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

        % ---- Set properties -----
        obj.grid = mbGrid;
        obj.order = order;
        obj.material = material;
        
        % Operators for full domain
        obj.D = D;
        obj.H = mbDiffOp.H;
        obj.E = E;

        % Initial condition
        obj.w0 = w0;
        
        % Outer boundary
        obj.boundary_data_fun = boundary_data_fun;
        
        % Friction properties
        obj.friction = friction;
        obj.e_fault = ic_struct.e;
        obj.penalty_fault = ic_struct.penalty;
        obj.tau_l_fun = tau_l_fun;
        obj.eta = eta;
        
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
        V = diff(obj.e_fault'*v);
    end

    function obj = setFaultTraction(obj)
        switch obj.friction.method
        case 'standard'
            penalty = obj.E.v'*obj.penalty_fault;
            F = obj.friction.funs.tau;
            a = obj.friction.params.a;
            E = obj.E;
            obj.fault_traction = @(t, U) penalty*[ F(t, obj.fault_jump(E.v*U), E.Psi*U, a); ...
                                                  -F(t, obj.fault_jump(E.v*U), E.Psi*U, a)];
        case 'erickson2022'
            penalty = obj.penalty_fault;
            F = obj.friction.funs.tau;
            Finv = obj.friction.funs.tauinv;
            a = obj.friction.params.a;
            E = obj.E;
            eta = obj.eta;

            %--- Setup nonlinear function that needs to be solved for V* ----
            nonlin_solve_fun = @(t, V, tau_l, Psi) eta*V + F(t, V, Psi, a) + tau_l;
            % Bracket for rootfinding
            bracket = @(t, tau_l, Psi) (tau_l >=0 )*[-Finv(t, tau_l, Psi, a), 0] + (tau_l < 0)* [0, -Finv(t, tau_l, Psi, a)];
            % Function that computes V* given tau_l and Psi (via nonlinear solve). 
            V_star_from_tau_l = @(t, tau_l, Psi) elastic.helpers.vectorBisection(@(V) nonlin_solve_fun(t, V, tau_l, Psi), bracket(t, tau_l, Psi), 1e-13);
            obj.V_star = @(t, U) V_star_from_tau_l(t, obj.tau_l_fun(t, U), E.Psi*U);
            % Fault traction function
            if isfield(obj.friction.funs, 'tau_L')
                obj.fault_traction = @(t, U) -penalty*(F(t, obj.V_star(t, U), E.Psi*U, a) - obj.friction.funs.tau_L(t));
            else 
                obj.fault_traction = @(t, U) -penalty*F(t, obj.V_star(t, U), E.Psi*U, a);
            end
        end
    end
        
    function obj = setStateEvolution(obj)
        G = obj.friction.funs.g;
        a = obj.friction.params.a;
        b = obj.friction.params.b;
        E = obj.E;
        switch obj.friction.method
        case 'standard'
            obj.state_evo = @(t, U) E.Psi'*G(obj.fault_jump(E.v*U), E.Psi*U, a, b);
        case 'erickson2022'
            obj.state_evo = @(t, U) E.Psi'*G(obj.V_star(t, U), E.Psi*U, a, b);
        end
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
    function F = assembleForcing(obj)
        F = [];
        % Sources
        if ~isempty(obj.S)
            F = @(t, U) obj.S(t);
        end
        % Fault traction
        if ~isempty(obj.fault_traction)
            if ~isempty(F)
                F = @(t, U) F(t, U) + obj.fault_traction(t, U);
            else
                F = obj.fault_traction;
            end
        end
        % BC
        if ~isempty(obj.boundary_data_fun)
            if ~isempty(F)
                F = @(t, U) F(t, U) + obj.boundary_data_fun(t, U);
            else
                F = obj.boundary_data_fun(t, U);
            end
        end
        % State evolution
        if ~isempty(obj.state_evo)
            if ~isempty(F)
                F = @(t, U) F(t, U) + obj.state_evo(t, U);
            else
                F = obj.state_evo;
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
            ts = time.ExplicitRungeKuttaDiscreteData(...
            obj.D, F, [], k, t, obj.w0, method.order);
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
            r.u = obj.E.u*obj.w0;
            r.v = obj.E.v*obj.w0;
            r.Psi = obj.E.Psi*obj.w0;
            r.us = obj.E.us*obj.w0;
            r.usim = obj.E.usim*obj.w0;
            r.usip = obj.E.usip*obj.w0;
            % r.Vs = obj.V_star(obj.w0);
            return
        end
        r.t = ts.t;
        w = ts.getV();
        r.u = obj.E.u*w;
        r.v = obj.E.v*w;
        r.Psi = obj.E.Psi*w;
        r.usim = obj.E.usim*w;
        r.usip = obj.E.usip*w;
        % r.Vs = obj.V_star(w);
    end
    
    % Sets up a plot of the discretisation
    %     update is a function_handle accepting a timestepper that updates the plot to the
    %            state of the timestepper
    function [update,figure_handle] = setupPlot(obj, type)
        figure_handle = figure();
        x = obj.grid.points();
        default_arg('type',struct);
        default_field(type,'plot_variables','all')
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
            V_plot = [V_plot, obj.fault_jump(r.v)];
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
            V_plot = [V_plot, obj.fault_jump(r.v)];
            log_abs_V_plot = log(abs(V_plot));
            a = obj.friction.params.a;
            b = obj.friction.params.b;
            sigma0 = obj.friction.params.sigma0;
            switch obj.friction.method
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