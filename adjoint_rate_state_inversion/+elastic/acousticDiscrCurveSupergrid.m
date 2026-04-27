classdef acousticDiscrCurveSupergrid < noname.Discretization
    properties
        name         = 'wave 2D curvilinear multiblock with absorbing boundary layers'
        description  = 'u_{tt} = a div( b grad u ) + f + damping'
        order        %Order of accuracy

        D            %Discretization matrix including BC
        E            %Matrix for u_t
        H            %Quadrature
        % v_tt = D*v + E*v_t + S_cont(t) + S_discr
        S_cont       %Function handle, v_tt = D*v + E*v_t + S_cont(t) + S_discr
        S_discr      %Matrix of time-dependent data, one column per time-stepping stage.

        v0           %Initial data
        v0t          %Initial data for time derivative
        diffOp       %Various operators
        grid         %Mulitblock grid in general

        A            %Coefficient matrix, u_tt = a*grad*b*div(u)
        superGrid    %Super grid damping details

        Dx, Dy       %Physical derivatives

    end

    methods

        function obj = acousticDiscrCurveSupergrid(g, order, a, b, F, bc, pointSources, opSet, superGrid)
            default_arg('superGrid', [])
            default_arg('opSet', []);
            default_arg('pointSources', []);
            default_arg('F', []);
            default_arg('bc', []);
            default_arg('a', 1);
            default_arg('b', 1);
            default_arg('order', 4);
            dim = 2;

            obj.grid = g;
            nBlocks = g.nBlocks;

            % Create diffOpPar cell array
            if iscell(a)
                doPar = cell(nBlocks,1);
                for i = 1:nBlocks
                    doPar{i} = {a{i}, b{i}, opSet};
                end
            else
                doPar = {a,b,opSet};
            end

            diffOp = multiblock.DiffOp(@scheme.LaplaceCurvilinearNew, g, order, doPar);

            obj.grid = g;
            H = diffOp.H;
            obj.H = H;
            D = diffOp.D;

            A = cell(nBlocks, nBlocks);
            B = cell(nBlocks, nBlocks);
            Dx = cell(nBlocks, nBlocks);
            Dy = cell(nBlocks, nBlocks);
            for i = 1:nBlocks
                A{i,i} = diffOp.diffOps{i}.a;
                B{i,i} = diffOp.diffOps{i}.b;
                Dx{i,i} = diffOp.diffOps{i}.Dx;
                Dy{i,i} = diffOp.diffOps{i}.Dy;
            end
            A = blockmatrix.toMatrix(A);
            B = blockmatrix.toMatrix(B);
            Dx = blockmatrix.toMatrix(Dx);
            Dy = blockmatrix.toMatrix(Dy);
            obj.A = A;
            obj.Dx = Dx;
            obj.Dy = Dy;

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

            % ------- Absorbing boundary layers ------- %
            defaultSuperGrid.gamma = {0, 1e-1, 0};
            defaultSuperGrid.DOI_IDs = 1:g.nBlocks();
            default_struct('superGrid', defaultSuperGrid);
            [gamma, DOI_IDs] = dealStruct(superGrid, {'gamma', 'DOI_IDs'});

            av = spdiags(A,0);
            bv = spdiags(B,0);
            cv = sqrt(av.*bv);
            c = spdiag(cv);

            % Build operators
            nBlocks = g.nBlocks;
            D1xi = cell(nBlocks, nBlocks);
            D1eta = cell(nBlocks, nBlocks);
            D2xi = cell(nBlocks, nBlocks);
            D2eta = cell(nBlocks, nBlocks);
            D4xi = cell(nBlocks, nBlocks);
            D4eta = cell(nBlocks, nBlocks);
            Hxi = cell(nBlocks, nBlocks);
            Heta = cell(nBlocks, nBlocks);

            for i = 1:nBlocks

                gi = g.grids{i};
                mi = gi.size();

                hxi = 1/(mi(1) - 1);
                heta = 1/(mi(2) - 1);

                ops_xi = sbp.D2Variable(mi(1), {0,1}, 2);
                ops_eta = sbp.D2Variable(mi(2), {0,1}, 2);

                d1xi = ops_xi.D1;
                d1eta = ops_eta.D1;
                d2xi = ops_xi.D2(ones(mi(1),1));
                d2eta = ops_eta.D2(ones(mi(2),1));
                Ixi = speye(mi(1));
                Ieta = speye(mi(2));

                Hxi{i,i} = hxi*kron(Ixi, Ieta);
                Heta{i,i} = heta*kron(Ixi, Ieta);

                flag = ~ismember(i, DOI_IDs);
                D1xi{i,i} = flag*kron(d1xi, Ieta);
                D1eta{i,i} = flag*kron(Ixi, d1eta);
                D2xi{i,i} = flag*kron(d2xi, Ieta);
                D2eta{i,i} = flag*kron(Ixi, d2eta);

                D4xi{i,i} = flag*D2xi{i,i}*D2xi{i,i};
                D4eta{i,i} = flag*D2eta{i,i}*D2eta{i,i};
            end

            % Compute "reference Jacobian", i.e value of Jacobian in DOI
            Jref = obj.diffOp.diffOps{DOI_IDs(1)}.J;
            Jref = mean(diag(Jref));

            D1xi = blockmatrix.toMatrix(D1xi);
            D1eta = blockmatrix.toMatrix(D1eta);
            D2xi = blockmatrix.toMatrix(D2xi);
            D2eta = blockmatrix.toMatrix(D2eta);
            D4xi = blockmatrix.toMatrix(D4xi);
            D4eta = blockmatrix.toMatrix(D4eta);
            Hxi = blockmatrix.toMatrix(Hxi);
            Heta = blockmatrix.toMatrix(Heta);

            E = - gamma{1} * Jref^(-1/dim) * ( Hxi \ (Hxi.^2 * (obj.H \ (D1xi' * obj.H * c * D1xi)))) ...
                - gamma{2} * Jref^(-1/dim) * ( Hxi \ (Hxi.^4 * (obj.H \ (D2xi' * obj.H * c * D2xi)))) ...
                - gamma{3} * Jref^(-1/dim) * ( Hxi \ (Hxi.^8 * (obj.H \ (D4xi' * obj.H * c * D4xi)))) ...
                - gamma{1} * Jref^(-1/dim) * ( Heta \ (Heta.^2 * (obj.H \ (D1eta' * obj.H * c * D1eta)))) ...
                - gamma{2} * Jref^(-1/dim) * ( Heta \ (Heta.^4 * (obj.H \ (D2eta' * obj.H * c * D2eta)))) ...
                - gamma{3} * Jref^(-1/dim) * ( Heta \ (Heta.^8 * (obj.H \ (D4eta' * obj.H * c * D4eta)))) ;
            %------------------------------------------ %

            obj.E = E;
            obj.superGrid = superGrid;


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

            DOI_IDs = obj.superGrid.DOI_IDs;
            k = inf;
            for i = DOI_IDs
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

            % Get total grid and grid for DOI only
            g = obj.grid;
            DOI_IDs = obj.superGrid.DOI_IDs;
            nDOIGrids = numel(DOI_IDs);
            if nDOIGrids > 1
                g_DOI = g.grids(DOI_IDs);
                conn = g.connections(DOI_IDs, DOI_IDs);
                g_DOI = multiblock.Grid(g_DOI, conn);
            else
                g_DOI = g.grids{DOI_IDs};
            end

            figure_handle = figure();
            v0 = obj.v0;

            % Loop through grids and compress them
            % Increase magnitude of solution inside SG layer
            factor = 1000;
            v_split = g.splitFunc(v0);
            for i = 1:g.nBlocks

                if ismember(i, obj.superGrid.W_IDs)
                    X = g.grids{i}.coords(:,1);
                    xr = max(X);
                    X = xr - (xr-X)/factor;
                    g.grids{i}.coords(:,1) = X;
                end

                if ismember(i, obj.superGrid.N_IDs)
                    Y = g.grids{i}.coords(:,2);
                    yl = min(Y);
                    Y = yl + (Y-yl)/factor;
                    g.grids{i}.coords(:,2) = Y;
                end

                if ismember(i, obj.superGrid.S_IDs)
                    Y = g.grids{i}.coords(:,2);
                    yr = max(Y);
                    Y = yr - (yr-Y)/factor;
                    g.grids{i}.coords(:,2) = Y;
                end

                if ismember(i, obj.superGrid.E_IDs)
                    X = g.grids{i}.coords(:,1);
                    xl = min(X);
                    X = xl + (X-xl)/factor;
                    g.grids{i}.coords(:,1) = X;
                end

                if ~ismember(i, DOI_IDs)
                    v_split{i} = 1e3*v_split{i};
                end

            end
            v = cell2mat(v_split);

            % Full domain
            figure
            Sur_full = multiblock.Surface(obj.grid, v);
            h_full = gca;
            xlabel('x')
            ylabel('y')
            shading interp
            colorbar
            xlims_full = xlim;
            ylims_full = ylim;
            axis equal
            a = gca;

            % Figure 2, DOI only
            figure
            v_split = g.splitFunc(v0);
            v_DOI = cell2mat(v_split(DOI_IDs));
            if nDOIGrids > 1
                Sur_DOI = multiblock.Surface(g_DOI, v_DOI);
                xlabel('x')
                ylabel('y')
                view(0,90)
                shading interp
                colorbar
                xlims_DOI = xlim;
                ylims_DOI = ylim;
                axis equal
            else
                p = g_DOI.points();
                X = grid.funcToPlotMatrix(g_DOI, p(:,1));
                Y = grid.funcToPlotMatrix(g_DOI, p(:,2));
                Z = grid.funcToPlotMatrix(g_DOI, v_DOI);

                Sur_DOI = surf(X,Y,Z);
                xlabel('x')
                ylabel('y')
                view(0,90)
                shading interp
                colorbar
                xlims_DOI = xlim;
                ylims_DOI = ylim;
                axis equal
            end
            h_DOI = gca;

            function update_fun(r)
                t = r.t;
                v = r.v;
                if ishandle(a)
                    title(a,sprintf('T = %.3f',t))
                end

                v_split = g.splitFunc(v);
                for i = 1:g.nBlocks
                    if ~ismember(i, DOI_IDs)
                        v_split{i} = 1e3*v_split{i};
                    end
                end
                v = cell2mat(v_split);

                Sur_full.ZData = v;
                Sur_full.CData = v;
                xlim(h_full, xlims_full);
                ylim(h_full, ylims_full);

                cmax = 1;
                caxis(h_full, [-cmax, cmax]);

                % Figure 2, DOI only
                v_split = g.splitFunc(v);
                v_DOI = cell2mat(v_split(DOI_IDs));
                if nDOIGrids > 1
                    Sur_DOI.ZData = v_DOI;
                    Sur_DOI.CData = v_DOI;
                else
                    Sur_DOI.ZData = grid.funcToPlotMatrix(g_DOI, v_DOI);
                    Sur_DOI.CData = grid.funcToPlotMatrix(g_DOI, v_DOI);
                end
                xlim(h_DOI, xlims_DOI);
                ylim(h_DOI, ylims_DOI);

                cmax = 1e-1;
                caxis(h_DOI, [-cmax, cmax]);
            end
            update = @(r)update_fun(r);
        end

        % Compare two functions u and v in the discrete l2 norm.
        function e = compareSolutions(obj, u, v)
            h = obj.grid.grids{1}.scaling;
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