classdef acousticDiscrCurve < noname.Discretization
    properties
        name         = 'wave 2D curvlilinear multiblock'
        description  = 'u_{tt} = a div( b grad u ) + f'
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
        grid         %Multiblock grid in general

        A            %Coefficient matrix for a, u_tt = a*grad*b*div(u)

    end

    methods

        function obj = acousticDiscrCurve(g, order, a, b, F, bc, pointSources, opSet)
            default_arg('opSet', []);
            default_arg('pointSources', []);
            default_arg('F', []);
            default_arg('a', []);
            default_arg('b', []);
            default_arg('order', 4);
            dim = 2;

            obj.grid = g;
            nBlocks = g.nBlocks;

            % Create diffOpPar cell array
            if iscell(a) && numel(a) > 2
                doPar = cell(nBlocks,1);
                for i = 1:nBlocks
                    doPar{i} = {a{i}, b{i}, opSet};
                end
            else
                doPar = {a,b,opSet};
            end

            diffOp = multiblock.DiffOp(@scheme.LaplaceCurvilinearNew, g, order, doPar);

            H = diffOp.H;
            obj.H = H;
            D = diffOp.D;

            A = cell(nBlocks, nBlocks);
            for i = 1:nBlocks
                A{i,i} = diffOp.diffOps{i}.a;
            end
            A = blockmatrix.toMatrix(A);
            obj.A = A;

            % Forcing function F
            if ~isempty(F)
                Ft = @(t) grid.evalOn(g, @(x,y) F(t,x,y));
            else
                Ft = [];
            end

            % Set BC
            [closure, penalty] = scheme.bcSetup(diffOp, bc);
            D = D + closure;

            %----- Point sources ----------------
            discrSource = false;
            contSource = false;
            if ~isempty(pointSources)

                x_s = pointSources.x;
                blockIds = pointSources.blockIds;
                n_s = length(x_s);

                % Initialize
                source_func_cont = @(t) 0*t;
                source_func_discr = 0;

                for s = 1:n_s

                    % Check which grid block source is in
                    blockId = blockIds(s);

                    % Delta function corresponding to local grid
                    delta_fun_local = elastic.diracDiscrCurve(x_s{s}, g.grids{blockId}, order, 0);

                    % Extend delta fun to global grid.
                    delta_fun = g.expandFunc(delta_fun_local, blockId);

                    % Accumulate different sources.
                    if isa(pointSources.g{s}, 'function_handle')
                        source_func_cont = @(t) source_func_cont(t) + ...
                                           obj.A * (pointSources.g{s}(t)*delta_fun);
                        contSource = true;
                    else
                        source_vec = pointSources.g{s};
                        if iscolumn(source_vec)
                            source_vec = transpose(source_vec);
                        end
                        source_func_discr = source_func_discr + ...
                                           obj.A * kron(source_vec, delta_fun);

                        discrSource = true;
                    end
                end
            end

            if discrSource
                ps_discr = source_func_discr;
            else
                ps_discr = [];
            end

            if contSource
                ps_cont = source_func_cont;
            else
                ps_cont = [];
            end
            % -----------------------------------

            % Create data function S(t)
            data_funcs = {Ft, penalty, ps_cont};
            S = elastic.addFunctionHandles(data_funcs);
            obj.S_cont = S;
            obj.S_discr = ps_discr;

            % Zero initial data
            v0_fun = @(x,y) 0*x;
            v0t_fun = @(x,y) 0*x;

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
        function [ts, N] = getTimestepper(obj,method,time_align,cfl)
            default_arg('method','rk4');
            default_arg('time_align',[]);
            switch method
                case 'rk4'
                    switch obj.order
                    case 2
                        default_arg('cfl', 0.5);
                    case 4
                        default_arg('cfl', 0.35);
                    case 6
                        default_arg('cfl', 0.27);
                    end
                    k = obj.getTimestep(method,cfl);

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
                switch obj.order
                case 2
                    cfl = 0.5;
                case 4
                    cfl = 0.35;
                case 6
                    cfl = 0.27;
                end
            end

            nBlocks = obj.grid.nBlocks;
            k = inf;
            for i = 1:nBlocks
                a = diag(obj.diffOp.diffOps{i}.a);
                b = diag(obj.diffOp.diffOps{i}.b);

                % Physical wave speed
                v = sqrt(a.*b);

                % Reference grid spacing
                h = obj.diffOp.diffOps{i}.h;
                h_xi = h(1);
                h_eta = h(2);

                x_xi = obj.diffOp.diffOps{i}.x_u;
                x_eta = obj.diffOp.diffOps{i}.x_v;
                y_xi = obj.diffOp.diffOps{i}.y_u;
                y_eta = obj.diffOp.diffOps{i}.y_v;

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

            h = multiblock.Surface(g, v0);
            shading interp

            axis equal
            view(0,90)
            % shading interp
            axis equal
            xlabel('x')
            ylabel('y')
            colorbar
            xlims1 = xlim;
            ylims1 = ylim;

            a = gca;

            function update_fun(r,xlims1,ylims1,h,g)
                t = r.t;
                v = r.v;
                if ishandle(a)
                    title(a,sprintf('T = %.3f',t))
                end
                h.ZData = v;
                h.CData = v;
                xlim(a, xlims1);
                ylim(a, ylims1);
                caxis(a, [-1,1]);
                drawnow;
            end
            update = @(r)update_fun(r,xlims1,ylims1,h,g);
        end


        % Compare two functions u and v in the discrete l2 norm.
        function e = compareSolutions(obj, u, v)
            H = obj.diffOp.H;

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