classdef elasticDiscrCurve < noname.Discretization
    properties
        name         = 'elastic 2D curvlilinear'
        description  = 'rho*u_{i,tt} = d_j a d_i u_j + d_j a d_j u_i + F(x,y,t)'
        order        %Order of accuracy

        D            %Discretization matrix including BC
        E            %Matrix for u_t
        H            %Total quadrature, for both components
        S            %Function handle, v_tt = D*v + E*v_t + S(t)
        v0           %Initial data
        v0t          %Initial data for time derivative
        diffOp       %Various operators
        grid         %Mulitblock grid in general

    end

    methods

        % boundary_data should be a struct containing
        %       - names: e.g. {'w','e','n','s'}
        %       - types: e.g. {{'d','d'},{'free','dirichlet'},{'free','free'},{'d','free'}}
        %       - funcs: cell array of function handles for boundary data, e.g
        %                {{g1W, g2W}, {g1E, g2E}, {g1S, g2S}, {g1N, g2N}}
        function obj = elasticDiscrCurve(g, order, lambda_fun, mu_fun, rho_fun, F, boundary_data, opSet, delta)
            default_arg('delta', []);
            default_arg('opSet',{@sbp.D2Variable, @sbp.D2Variable});
            default_arg('F',[]);
            default_arg('boundary_data',[]);
            default_arg('rho_fun',[]);
            default_arg('mu_fun',[]);
            default_arg('lambda_fun',[]);
            default_arg('order',4);
            dim = 2;

            obj.grid = g;
            diffOp = scheme.Elastic2dCurvilinear(g, order, lambda_fun, mu_fun, rho_fun, opSet);

            H = diffOp.H;
            obj.H = kron(H,speye(2));
            RHO = diffOp.RHO;
            RHO_kron = kron(RHO,speye(2));

            % Forcing function F
            F_comb = @(x,y,t) [F{1}(x,y,t); F{2}(x,y,t)];
            if ~isempty(F)
                Ft = @(t) RHO_kron\grid.evalOn(g, @(x,y)F_comb(x,y,t));
            else
                Ft = [];
            end

            % Delta function
            if ~isempty(delta)
                E = diffOp.E;
                delta = @(t) RHO_kron\(E{1}*delta{1}(t) + E{2}*delta{2}(t));
            else
                delta = [];
            end

            % Set BC
            X = g.points();
            x = X(:,1);
            y = X(:,2);
            D = diffOp.D;
            if isempty(boundary_data)
                boundaries = {'w','e','n','s','w','e','n','s'};
                bc_types = { {1, 't'}, {1, 't'}, {1, 't'}, {1, 't'},...
                             {2, 't'}, {2, 't'}, {2, 't'}, {2, 't'} };
                boundary_data_funcs = cell(numel(boundaries), 1);
            else
                boundaries = boundary_data.names;
                boundary_data_funcs = boundary_data.funcs;
                bc_types = boundary_data.types;
            end

            nb = length(boundaries);
            boundary_data_tot = cell(nb,1);

            % Loop over boundary conditions
            for j = 1:nb
                b = boundaries{j};
                data = boundary_data_funcs{j};

                [closure, penalty] = diffOp.boundary_condition(b,bc_types{j});
                D = D + closure;

                % Create boundary data functions
                eB = diffOp.get_boundary_operator('e', b);
                boundary_data_tot{j} = @(t) penalty*data(eB'*x,eB'*y,t) ;
            end

            % Create data function S(t)
            data_funcs = cell(nb+2,1);
            data_funcs{1} = Ft;
            data_funcs{2} = delta;
            for j = 1:nb
                data_funcs{2+j} = boundary_data_tot{j};
            end

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

            % Initial data, two components of displacement
            x0 = 0.5;
            y0 = 0.5;
            sigma = 0.5;
            v0_fun = @(x,y) [(exp(-((x-x0).^2+(y-y0).^2)/sigma^2)); ...
                            0*x];
            v0t_fun = @(x,y) [0*x; 0*y];

            obj.v0 = grid.evalOn(g, v0_fun);
            obj.v0t = grid.evalOn(g, v0t_fun);

            obj.diffOp = diffOp;
            obj.D = D;
            obj.E = [];
            obj.order = order;

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
            default_arg('method','rk4');
            default_arg('time_align',[]);
            switch method
                case 'rk4'
                    cfl = 0.25;
                    k = obj.getTimestep(method,cfl);

                    if ~isempty(time_align)
                        [k, N] = alignedTimestep(k, time_align);
                    end

                    t = 0;
                    ts = time.Rungekutta4SecondOrder(obj.D, obj.E, obj.S, k, t, obj.v0, obj.v0t);
                otherwise
                    error('Timestepping method ''%s'' not supported',method);
            end
        end

        function k = getTimestep(obj, method, cfl)
            default_arg('cfl',[]);
            if isempty(cfl)
%             try
%                 k = 0.5*time.rk4.get_rk4_time_step(eigs(obj.D,1));
%             catch
                k = 0.5*time.rk4.get_rk4_time_step(max(abs(eig(full(obj.D)))));
%             end
            else
                RHO = diag(obj.diffOp.RHO);
                mu = diag(obj.diffOp.MU);
                lambda = diag(obj.diffOp.LAMBDA);

                v_s = sqrt(max( mu./RHO ));
                v_p = sqrt(max( (lambda + 2*mu)./RHO ));
                v = max(v_p,v_s);

                try
                    h = min(obj.grid.h);
                catch
                    h = min(obj.grid.scaling());
                end

                k = cfl/v*h;
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
            r.vt = ts.getVt();
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
            E = obj.diffOp.E;
            v0 = obj.v0;
            v0_1 = E{1}'*v0;
            v0_2 = E{2}'*v0;

            v0_1 = grid.funcToPlotMatrix(g, v0_1);
            v0_2 = grid.funcToPlotMatrix(g, v0_2);

            p = g.points();
            X = grid.funcToPlotMatrix(g, p(:,1));
            Y = grid.funcToPlotMatrix(g, p(:,2));
            h1 = subplot(2,1,1);
            Sur1 = surf(X,Y,v0_1);
            view(0,90)
            % shading interp
            axis equal
            xlabel('x')
            ylabel('y')
            colorbar
            xlims1 = xlim;
            ylims1 = ylim;

            h2 = subplot(2,1,2);
            Sur2 = surf(X,Y,v0_2);
            view(0,90)
            axis equal
            % shading interp
            xlabel('x')
            ylabel('y')
            colorbar
            xlims2 = xlim;
            ylims2 = ylim;
            a = gca;

            function update_fun(r,E,xlims1,ylims1,xlims2,ylims2)
                t = r.t;
                v1 = E{1}'*r.v;
                v1 = grid.funcToPlotMatrix(g, v1);
                v2 = E{2}'*r.v;
                v2 = grid.funcToPlotMatrix(g, v2);
                if ishandle(a)
                    title(a,sprintf('T = %.3f',t))
                end
                Sur1.ZData = v1;
                Sur1.CData = v1;
                Sur2.ZData = v2;
                Sur2.CData = v2;
                xlim(h1, xlims1);
                ylim(h1, ylims1);
                xlim(h2, xlims2);
                ylim(h2, ylims2);
                caxis(h1, [-1,1]);
                caxis(h2, [-1,1]);
            end
            update = @(r)update_fun(r,E,xlims1,ylims1,xlims2,ylims2);
        end


        % Compare two functions u and v in the discrete l2 norm.
        function e = compareSolutions(obj, u, v)
            g = obj.grid;
            H = obj.diffOp.H_kron;
            J = obj.diffOp.J_kron;
            H = H*J;

            evec = u - v;
            e = sqrt( evec'*H*evec );
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