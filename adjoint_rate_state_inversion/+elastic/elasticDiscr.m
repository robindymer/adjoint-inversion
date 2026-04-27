classdef elasticDiscr < noname.Discretization
    properties
        name         = 'elastic 2D with point source'
        description  = 'rho*u_{i,tt} = d_j a d_i u_j + d_j a d_j u_i + F(x,y,t)'
        order        %Order of accuracy

        D            %Discretization matrix including BC
        E            %Matrix for u_t
        H            %Total quadrature, for both components

        % v_tt = D*v + E*v_t + S_cont(t) + S_discr
        S_cont       %Function handle, v_tt = D*v + E*v_t + S_cont(t) + S_discr
        S_discr      %Matrix of time-dependent data, one column per time-stepping stage.

        v0           %Initial data
        v0t          %Initial data for time derivative
        diffOp       %Various operators
        grid         %Mulitblock grid in general

    end

    methods

        % boundary_data should be a struct containing
        %       - names: e.g. {'w','e','n','s','w','e','n','s'}
        %       - types: Coordinate direction and kind of BC (displacement or traction)
        %                e.g. { {1, 't'}, {1, 't'}, {1, 'd'}, {1, 't'},...
        %                       {2, 'd'}, {2, 't'}, {2, 't'}, {2, 't'} };
        %       - funcs: cell array of function handles for boundary data, e.g
        %                {g1W, g1E, g1N, g1S, g2W, g2E, g2N, g2S}

        % point_sources should be a struct containing cell arrays with one element per source
        %       - x:    cell array of coordinate vectors specifying source locations,
        %               e.g. x = {[x1, y1], [x2, y2]}
        %       - g:    cell array of cell arrays of scalar source time function handles. An array {gx, gy}
        %               specifies source strength gx in x-dir and gy in y-dir. Example for two sources:
        %               g = {{g1x, g1y}, {g2x, g2y}}
        function obj = elasticDiscr(m, order, xlim, ylim, lambda, mu, rho,...
                                               F, boundary_data, point_sources, opSet, optFlag)
            default_arg('optFlag',[]);
            default_arg('opSet',[]);
            default_arg('point_sources',[]);
            default_arg('F',[]);
            default_arg('boundary_data',[]);
            default_arg('rho',[]);
            default_arg('mu',[]);
            default_arg('lambda',[]);
            default_arg('xlim',{0,1});
            default_arg('ylim',{0,1});
            default_arg('order',4);
            dim = 2;

            %--- Domain definition, grid, diffOp ---
            % Correct grid end points if periodic.
            if ~isempty(opSet)
                ops = opSet{1}(m(1), xlim, order);
                x = ops.x;
                hx = ops.h;

                ops = opSet{2}(m(2), ylim, order);
                y = ops.x;
                hy = ops.h;

                g = grid.Cartesian(x,y);
            else
                g = grid.equidistant(m, xlim, ylim);
                hx = (xlim{2}-xlim{1})/(m(1) - 1);
                hy = (ylim{2}-ylim{1})/(m(2) - 1);
            end
            g.h = [hx, hy];
            g.lim = {xlim, ylim};
            obj.grid = g;
            % --------------------------------

            diffOp = scheme.Elastic2dVariable(g, order, lambda, mu, rho, opSet, optFlag);

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

            % Point sources
            if ~isempty(point_sources)

                x_s = point_sources.x;
                n_s = length(x_s);
                E = diffOp.E;       % E{i}' picks out component i
                H_1D = diffOp.H_1D; % 1D quadratures

                % Initialize
                source_func_cont = @(t) 0*t;
                source_func_discr = 0;

                % Loop over sources
                for s = 1:n_s

                    % Delta function corresponding to one component
                    delta_fun_local = diracDiscr(x_s{s}, obj.grid.x, order, 0, H_1D);

                    for i = 1:dim
                        % Project source for component i to full grid size and accumulate different sources.
                        if isa(point_sources.g{s}{i}, 'function_handle')
                            source_func_cont = @(t) source_func_cont(t) + ...
                                               RHO_kron\(point_sources.g{s}{i}(t)*E{i}*delta_fun_local);
                        else
                            source_vec = point_sources.g{s}{i};
                            if iscolumn(source_vec)
                                source_vec = transpose(source_vec);
                            end
                            source_func_discr = source_func_discr + ...
                                               RHO_kron\kron(source_vec, E{i}*delta_fun_local);
                        end
                    end
                end
                % If there were no discrete source functions, make it empty
                if isscalar(source_func_discr)
                    source_func_discr = [];
                end
                % If there were no continuous source functions, make it empty
                if isscalar(source_func_cont(0))
                    source_func_cont = [];
                end
            else
                source_func_cont = [];
                source_func_discr = [];
            end
            obj.S_discr = source_func_discr;

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
                [closure, penalty] = diffOp.boundary_condition(b,bc_types{j});
                D = D + closure;

                data = boundary_data_funcs{j};
                if isempty(data)
                    boundary_data_tot{j} = [];
                else
                    % Create boundary data functions
                    eB = diffOp.get_boundary_operator('e', b);
                    boundary_data_tot{j} = @(t) penalty*data(eB'*x,eB'*y,t) ;
                end
            end

            % Create data function S(t)
            data_funcs = cell(nb+2,1);
            data_funcs{1} = Ft;
            data_funcs{2} = source_func_cont;
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
            obj.S_cont = S;

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
        function [ts, N] = getTimestepper(obj, method, time_align, k) %% ???
            default_arg('method','rk4');
            default_arg('time_align',[]);
            default_arg('k',[]);
            switch method
                case 'rk4'
                    cfl = 0.25;
                    if isempty(k)
                        k = obj.getTimestep(method,cfl);
                    end

                    if ~isempty(time_align)
                        [k, N] = alignedTimestep(k, time_align);
                    end

                    t = 0;
                    ts = time.ExplicitRungeKuttaSecondOrderDiscreteData(...
                                            obj.D, obj.E, obj.S_cont, obj.S_discr, k, t, obj.v0, obj.v0t);
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

                h = min(obj.grid.h);
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

            X = g.matrices();
            h1 = subplot(2,1,1);
            Sur1 = surf(X{1}',X{2}',v0_1);
            view(0,90)
            shading interp
            axis equal
            xlabel('x')
            ylabel('y')
            colorbar
            xlims1 = xlim;
            ylims1 = ylim;

            h2 = subplot(2,1,2);
            Sur2 = surf(X{1}',X{2}',v0_2);
            view(0,90)
            axis equal
            shading interp
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

        % Generates delta functions.
        % Source coordinates should be a cell array of coordinate vectors.
        % Returns a (nSources x dim) cell array of delta functions.
        function deltaFunctions = generateDeltaFunctions(obj, sourceCoord, nMoment, nSmooth)
            default_arg('nMoment', obj.order);
            default_arg('nSmooth', 0);

            nSources = numel(sourceCoord);
            deltaFunctions = cell(nSources, obj.diffOp.dim);

            E = obj.diffOp.E;       % E{i}' picks out component i
            H_1D = obj.diffOp.H_1D;

            for i = 1:nSources
                deltaLocal = diracDiscr(sourceCoord{i}, obj.grid.x, nMoment, nSmooth, H_1D);
                for d = 1:obj.diffOp.dim
                    deltaFunctions{i,d} = E{d}*deltaLocal;
                end
            end
        end

        % Compare two functions u and v in the discrete l2 norm.
        function e = compareSolutions(obj, u, v)
            g = obj.grid;
            h = obj.grid.h;

            evec = u - v;
            e = sqrt(prod(h))*norm(evec);
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