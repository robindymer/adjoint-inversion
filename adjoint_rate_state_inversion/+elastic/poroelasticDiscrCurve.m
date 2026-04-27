classdef poroelasticDiscrCurve < noname.Discretization
    properties
        name         = 'poroelastic 2D'
        description  = 'poroelastic equations'
        order        %Order of accuracy
        
                     % Solution: v = [u1, u2, p];
                     % Displacements and pressure.
        A            %Matrix for v_t
        B            %Matrix for v
        E_u, E_u1, E_u2  %Pick out displacements 
        E_p          %Picks out pressure
        Extract_p    % Matrix to compute p from u and p'
        Extract_pprime    % Matrix to compute p' from u and p
        H            %Total quadrature, for all components           
        S            %Function handle, A*v_t = B*v + S(t)
        v0           %Initial data
        diffOpElastic    %Various operators 
        diffOpHeat       %Various operators 
        Div          %Divergence operator
        grid         %Mulitblock grid in general

        % Metric coefficients
        J, Ji
        b

        M_u, M_p, F_u, F_p % Matrices for the different components.`
        M_bu, M_bp, F_b % Boundary data matrices
        
    end

    methods
        
        % boundary_data should be a struct containing
        %       - types: e.g. {'d','d'; 'free','dirichlet'; 'd','n'}
        %       - funcs: cell array of function handles for boundary data, e.g
        %                {g_u1S, g_u1N; g_u2S, g_u2N; g_pS, g_pN}
        %
        % The ordering is {u1S, u1N; ...
        %                  u2S, u2N; ...
        %                   pS,  pN}; 

        function obj = poroelasticDiscrCurve(g, order, coeff_fun, F, boundary_data, opSet, BCstrength)
            default_arg('BCstrength',[]);
            default_arg('opSet',[]);
            default_arg('F',[]);
            default_arg('boundary_data',[]);
            default_arg('order',4);
            dim = 2;

            K_fun = coeff_fun.K;
            M_fun = coeff_fun.M;
            alpha_fun = coeff_fun.alpha;
            kappa_fun = coeff_fun.kappa;
            mu_fun = coeff_fun.mu;
            lambda_fun = coeff_fun.lambda;

            default_arg('K_fun', @(x,y) 0*x + 1);
            default_arg('M_fun', @(x,y) 0*x + 1);
            default_arg('alpha_fun', @(x,y) 0*x + 1);
            default_arg('kappa_fun',[]);
            default_arg('mu_fun',[]);
            default_arg('lambda_fun',[]);

            %--- Grid ---
            obj.grid = g;
            N = obj.grid.N();
            m = obj.grid.size();
            % --------------------------------

            % Boundary conditions
            BC = boundary_data.types;

            % DiffOp building blocks
            rho_fun = [];
            diffOpElastic = scheme.Elastic2dCurvilinear(g, order, lambda_fun, mu_fun, rho_fun, opSet);
            diffOpHeat = scheme.Heat2dCurvilinear(g ,order, kappa_fun, opSet);
            H = diffOpHeat.H;
            J = diffOpHeat.J;
            Ji = diffOpHeat.Ji;
            b = diffOpHeat.b;
            beta = diffOpHeat.beta;
            obj.H = kron(H,speye(dim+1));
            obj.J = kron(J,speye(dim+1));
            obj.Ji = kron(Ji,speye(dim+1));
            obj.b = b;
            Laplace_kappa = diffOpHeat.D;
            D_elastic = diffOpElastic.D;

            % Coeffcients
            alpha = spdiag(grid.evalOn(g, alpha_fun));
            M = spdiag(grid.evalOn(g, M_fun));
            K = spdiag(grid.evalOn(g, K_fun));
            G = diffOpElastic.MU;
            Mi = inv(M);

            % Forcing function F
            F_comb = @(x,y,t) [F{1}(x,y,t); F{2}(x,y,t); F{3}(x,y,t)];
            function F = vectorTimeForcing(t, F_comb, g)
                F = [];
                for i = 1:length(t)
                    F = [F; grid.evalOn(g, @(x,y)F_comb(x,y,t(i)))];
                end
            end

            if ~isempty(F)
                Ft = @(t) vectorTimeForcing(t, F_comb, g);
            else
                Ft = [];
            end

            % Divergence operator
            E = diffOpElastic.E;
            D1 = diffOpElastic.D1;
            Div = sparse(N,N*dim);

            for j = 1:dim
                for k = 1:dim
                    Div = Div + b{j,k}*D1{k}*E{j}';
                end
            end

            % Operators for special BC
            Hi = diffOpElastic.Hi;
            H_boundary_l = diffOpHeat.H_boundary_l;
            H_boundary_r = diffOpHeat.H_boundary_r;
            eS = diffOpElastic.get_boundary_operator('e', 's');
            eN = diffOpElastic.get_boundary_operator('e', 'n');
            kappa = diffOpHeat.KAPPA;
            kappaS = eS'*kappa*eS;
            kappaN = eN'*kappa*eN;

            %--- Fluid flow equation ----
            BC_p_S = BC{3,1};
            BC_p_N = BC{3,2};
            [closure_S_heat, penalty_S_heat] = diffOpHeat.boundary_condition('s',BC_p_S,true,BCstrength);
            [closure_N_heat, penalty_N_heat] = diffOpHeat.boundary_condition('n',BC_p_N,true,BCstrength);
            Laplace_kappa = Laplace_kappa + closure_N_heat + closure_S_heat;

            % Matrix for p
            F_p = M*Laplace_kappa;

            % Matrix for displacements
            F_u = sparse(N,2*N); % Zero before change of variables

            % Boundary data object for pressure
            penalty_S_heat = M*penalty_S_heat;
            penalty_N_heat = M*penalty_N_heat;

            % If flux BC, multiply penalty by kappa because pressure data will be used.
            switch BC_p_S
            case {'N','n','Neumann','neumann'}
                penalty_S_heat = penalty_S_heat*kappaS;
            end
            switch BC_p_N
            case {'N','n','Neumann','neumann'}
                penalty_N_heat = penalty_N_heat*kappaN;
            end
            F_b.S = penalty_S_heat;
            F_b.N = penalty_N_heat;
            %-------------------------------------

            %--- Mechanical equilibrium equation ----
            BC_displ_S = {BC{1,1}, BC{2,1}};
            BC_displ_N = {BC{1,2}, BC{2,2}};

            [closure_S_elastic, penalty_S_elastic] = ...
                                diffOpElastic.boundary_condition('s', BC_displ_S, BCstrength);

            [closure_N_elastic, penalty_N_elastic] = ...
                                diffOpElastic.boundary_condition('n',BC_displ_N, BCstrength);

            D_elastic = D_elastic + closure_S_elastic + closure_N_elastic;

            % Matrix for displacements
            M_u = D_elastic;

            M_bu.S = [penalty_S_elastic{1}, penalty_S_elastic{2}];
            M_bu.N = [penalty_N_elastic{1}, penalty_N_elastic{2}];

            % Matrix for p
            M_p = sparse(N*dim, N);
            for i = 1:dim
                for k = 1:dim
                    M_p = M_p - E{i}*b{i,k}*D1{k}*alpha; 
                end
            end

            % If traction BC for displacements, modify pressure BC too
            % --- South ---- %
            n = -1;
            for k = 2
                switch BC_displ_S{k}
                case {'T','t','Traction','traction'};
                    M_p = M_p + penalty_S_elastic{k}*n*eS'*alpha; 
                end
            end
            % ---------------- %

            % --- North --- -%
            n = 1;
            for k = 2
                switch BC_displ_N{k}
                case {'T','t','Traction','traction'};
                    M_p = M_p + penalty_N_elastic{k}*n*eN'*alpha;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Penalties in mechanical equilibrium equation depending on
            % Pressure BC
            M_bp = struct;
            M_bp.W = sparse(dim*N, m(2));
            M_bp.E = sparse(dim*N, m(2));
            M_bp.S = sparse(dim*N, m(1));
            M_bp.N = sparse(dim*N, m(1));

            j = 2; % Boundary number
            % ----- South ---- %
            n = -1;
            H_gamma = H_boundary_l{j};
            switch BC_p_S
            case {'D','d','Dirichlet','dirichlet'}
                for i = 1:dim
                    M_p = M_p + n*E{i}*Ji*Hi*inv(beta{j})*J*b{i,j}*alpha*eS*H_gamma*eS';
                    M_bp.S = M_bp.S - n*E{i}*Ji*Hi*inv(beta{j})*J*b{i,j}*alpha*eS*H_gamma;
                end
            otherwise
                M_bp.S = 0*E{j}*Hi*alpha*eS*H_gamma;
            end

            % ---- North ----- %
            n = 1;
            H_gamma = H_boundary_r{j};
            switch BC_p_N
            case {'D','d','Dirichlet','dirichlet'}
                for i = 1:dim
                    M_p = M_p + n*E{i}*Ji*Hi*inv(beta{j})*J*b{i,j}*alpha*eN*H_gamma*eN';
                    M_bp.N = M_bp.N - n*E{i}*Ji*Hi*inv(beta{j})*J*b{i,j}*alpha*eN*H_gamma;
                end
            otherwise
                M_bp.N =  0*E{j}*Hi*alpha*eN*H_gamma;    
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Ordering: [u1; u2; p].
            E = cell(dim,1);
            I = speye(N,N);
            for i = 1:dim+1
                e = sparse(dim+1,1);
                e(i) = 1;
                E{i} = kron(I,e);
            end
            E_p = E{3};
            e_u = sparse(dim+1,2);
            e_u(1,1) = 1;
            e_u(2,2) = 1;
            E_u = kron(I,e_u);
            % E_u = [E{1}, E{2}];
            obj.E_u = E_u;
            obj.E_u1 = E{1};
            obj.E_u2 = E{2};
            obj.E_p = E_p;

            %%%%%%%% Change of variables %%%%%%%%%%%%%%%%%%%
            % p = p' - M * alpha * d_j u_j
            F_u = F_u - F_p*M*alpha*Div;
            M_u = M_u - M_p*M*alpha*Div;

            % Matrices for converting between p' and p.
            obj.Extract_p = E_p' - M*alpha*Div*E_u'; 
            obj.Extract_pprime = E_p' + M*alpha*Div*E_u';   

            % ODE on form A*v_t = B*v + f;
            % v = [u; p];
            % A = {0*M_u, 0*M_p; 0*F_u, I};
            % B = {M_u, M_p; F_u, F_p};
            A = E_p*E_p';
            B = E_u*M_u*E_u' + E_u*M_p*E_p' + ...
                E_p*F_u*E_u' + E_p*F_p*E_p';
            obj.A = A;
            obj.B = B;

            % Forcing and data
            bc_data = boundary_data.funcs;
            penalty_S = {E_u*penalty_S_elastic{1}; E_u*penalty_S_elastic{2}; ...
                        E_p*penalty_S_heat + E_u*M_bp.S};
            penalty_N = {E_u*penalty_N_elastic{1}; E_u*penalty_N_elastic{2}; ...
                        E_p*penalty_N_heat + E_u*M_bp.N};
            data_funcs = cell(4,2);
            data_funcs{1,1} = Ft;

            X = g.points();
            x = X(:,1);
            y = X(:,2);

            % South boundary
            for k = 1:3
                data_funcs{1+k,1} = @(t) penalty_S{k}*bc_data{k,1}(eS'*x, eS'*y, t);
            end
            % North boundary
            for k = 1:3
                data_funcs{1+k,2} = @(t) penalty_N{k}*bc_data{k,2}(eN'*x, eN'*y, t);
            end

            % Create data function S(t)
            S = [];
            for i = 1:numel(data_funcs)
                data = data_funcs{i};
                if ~isempty(data)
                    if(isempty(S))
                        S = data;
                    else
                        S = @(t) S(t) + data(t);
                    end
                end
            end
            obj.S = S;

            % Initial data, two components of displacement and pressure
            x0 = 0.5;
            y0 = 0.5;
            sigma = 0.5;
            v0_fun = @(x,y) [(exp(-((x-x0).^2+(y-y0).^2)/sigma^2)); ...
                            0*x;...
                            0*x];

            obj.v0 = grid.evalOn(g, v0_fun);

            % Misc.
            obj.M_u = M_u;
            obj.M_p = M_p;
            obj.F_u = F_u;
            obj.F_p = F_p;
            obj.M_bu = M_bu;
            obj.M_bp = M_bp;
            obj.F_b = F_b;
            obj.diffOpElastic = diffOpElastic;
            obj.diffOpHeat = diffOpHeat;
            obj.order = order;

            obj.Div = Div;


        end
        % Prints some info about the discretisation
        function printInfo(obj)
            fprintf('Name: %s\n',obj.name);
            fprintf('Size: %d\n',obj.size());
        end

        % Return the number of DOF
        function n = size(obj)
            n = length(obj.v0);
        end

        % Returns a timestepper for integrating the discretisation in time
        %     method is a string that states which timestepping method should be used.
        %          The implementation should switch on the string and deliver
        %          the appropriate timestepper. It should also provide a default value.
        %     time_align is a time that the timesteps should align with so that for some
        %                integer number of timesteps we end up exactly on time_align
        function [ts, N] = getTimestepper(obj,method,time_align) %% ???
            default_arg('method','sbp');
            default_arg('time_align',[]);
            switch method
                case 'sbp'
                    cfl = 2;
                    k = obj.getTimestep(method,cfl);
                    
                    if ~isempty(time_align)
                        [k, N] = alignedTimestep(k, time_align);
                    end

                    t = 0;
                    ts = time.SBPInTimeImplicitFormulation(obj.A, obj.B, obj.S, k, t, obj.v0);
                otherwise
                    error('Timestepping method ''%s'' not supported',method);
            end
        end

        function k = getTimestep(obj, method, cfl)
            default_arg('cfl',[]);
            if isempty(cfl)   
                error('Specify cfl')
            else
                kappa = max(max(obj.diffOpHeat.KAPPA));
                h = min(obj.grid.h);
                if(kappa > 0)
                    k = cfl/kappa*h;
                else
                    k = cfl*h;
                end
            end
        end

        function r = getTimeSnapshot(obj, ts)
            if ts == 0
                r.t = 0;
                r.v = obj.v0;
                return
            end
            r.t = ts.t;
            r.v = ts.getV();
        end

        % Sets up movie recording to a given file.
        %     saveFrame is a function_handle with no inputs that records the current state
        %               as a frame in the moive.
        function saveFrame = setupMov(obj, file)
            error('not implemented');
        end

        % Sets up a plot of the discretisation
        %     update is a function_handle accepting a timestepper that updates the plot to the
        %            state of the timestepper
        function [update,figure_handle] = setupPlot(obj, type)
            
            g = obj.grid;
            figure_handle = figure();
            v0 = obj.v0;
            v0_1 = obj.E_u1'*v0;
            v0_2 = obj.E_u2'*v0;
            v0_3 = obj.E_p'*v0;

            v0_1 = grid.funcToPlotMatrix(g, v0_1);
            v0_2 = grid.funcToPlotMatrix(g, v0_2);
            v0_3 = grid.funcToPlotMatrix(g, v0_3);

            X = g.matrices();
            h1 = subplot(3,1,1);
            Sur1 = surf(X{1}',X{2}',v0_1);
            view(0,90)
            shading interp
            % axis equal
            xlabel('x')
            ylabel('y')
            colorbar
            xlims1 = xlim;
            ylims1 = ylim;

            h2 = subplot(3,1,2);
            Sur2 = surf(X{1}',X{2}',v0_2);
            view(0,90)
            % axis equal
            shading interp
            xlabel('x')
            ylabel('y')
            colorbar
            xlims2 = xlim;
            ylims2 = ylim;

            h3 = subplot(3,1,3);
            Sur3 = surf(X{1}',X{2}',v0_3);
            view(0,90)
            % axis equal
            shading interp
            xlabel('x')
            ylabel('y')
            colorbar
            xlims3 = xlim;
            ylims3 = ylim;

            a = gca;

            function update_fun(r,E_u1,E_u2,E_p,...
                                xlims1,ylims1,xlims2,ylims2,xlims3,ylims3)
                t = r.t;
                v1 = E_u1'*r.v;
                v1 = grid.funcToPlotMatrix(g, v1);
                v2 = E_u2'*r.v;
                v2 = grid.funcToPlotMatrix(g, v2);
                v3 = E_p'*r.v;
                v3 = grid.funcToPlotMatrix(g, v3);
                if ishandle(a)
                    title(a,sprintf('T = %.3f',t))
                end
                Sur1.ZData = v1;
                Sur1.CData = v1;
                Sur2.ZData = v2;
                Sur2.CData = v2;
                Sur3.ZData = v3;
                Sur3.CData = v3;
                xlim(h1, xlims1);
                ylim(h1, ylims1);
                xlim(h2, xlims2);
                ylim(h2, ylims2);
                xlim(h3, xlims3);
                ylim(h3, ylims3);
                % caxis(h1, [-1,1]);
                % caxis(h2, [-1,1]);
                % caxis(h3, [1,3]);
            end
            update = @(r)update_fun(r,obj.E_u1,obj.E_u2,obj.E_p,...
                    xlims1,ylims1,xlims2,ylims2,xlims3,ylims3);
        end
        
        
        % Compare two functions u and v in the discrete l2 norm.
        function e = compareSolutions(obj, u, v)
            g = obj.grid;
            H = obj.H;
            J = obj.J;
            
            evec = u - v;
            e = sqrt(evec'*H*J*evec);
        end

        % Compare the grid function u to the analytical function g in the discrete l2 norm.
        function e = compareSolutionsAnalytical(obj, u, g)
            gr = obj.grid;
            v = grid.evalOn(gr, g);
            e = obj.compareSolutions(u, v);
        end

    end

    methods(Static)
        
    end
end