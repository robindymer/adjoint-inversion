classdef elasticAnisotropicSupergridMixedDiscr < noname.Discretization
    properties
        name         = 'Anisotropic elastic 2D multiblock curvilinear with supergrid'
        description  = 'rho*u_{j,tt} = d_i C_ijkl d_k u_l + F(x,y,t)'
        order        %Order of accuracy
        dim          %Number of spatial dimensions

        D            %Discretization matrix including BC
        E            %Matrix for u_t
        H            %Total quadrature, for both components

        % v_tt = D*v + E*v_t + S_cont(t) + S_discr
        S_cont       %Function handle, v_tt = D*v + E*v_t + S_cont(t) + S_discr
        S_discr      %Matrix of time-dependent data, one column per time-stepping stage.

        v0           %Initial data
        v0t          %Initial data for time derivative
        diffOp, diffOps       %Various operators
        grid         %Multiblock grid in general

        RHO          %Density, diagonal matrix
        RHO_kron     %Kroneckered density
        C            %Elastic stiffness tensor
        sigma        %Elastic stress operator

        Ecomp        % Ecomp{1}' picks out first component, Ecomp{2}' picks out 2nd comp.
        superGrid
    end

    methods

        % bc should be a cell array of structs, one struct per BC.
        % Each struct contains the fields
        %       - boundary: e.g. {1, 'w'} for west boundary, 1st block
        %       - type:     e.g. {2, 'd'} for displacement condition on the 2nd component
        %       - data:     boundary data function handle
        function obj = elasticAnisotropicSupergridMixedDiscr(g, order, C_upwind, C_narrow, C_wide, rho, F,...
                                             bc, pointSources, opSet, superGrid, optFlag)
            default_arg('optFlag',[]);
            default_arg('superGrid',[]);
            default_arg('pointSources',[]);
            default_arg('opSet',[]);
            default_arg('F',[]);
            default_arg('bc',[]);
            default_arg('rho',[]);
            default_arg('C_upwind',[]);
            default_arg('C_narrow',[]);
            default_arg('C_wide',[]);
            default_arg('xlim',[-1,0,1]);
            default_arg('ylim',[1,0]);
            default_arg('order',4);
            dim = 2;
            obj.dim = dim;

            % Domain definition, grid, diffOp
            obj.grid = g;

            nBlocks = g.nBlocks();
            doParam = cell(nBlocks, 1);

            if isa(rho, 'double') && ~isempty(rho);
                rho = g.splitFunc(rho);
            end

            if isempty(C_upwind)
                opSetUpwind = {@sbp.D1Upwind, @sbp.D1Upwind};
            else
                opSetUpwind = {@sbp.D1UpwindCompatible, @sbp.D1UpwindCompatible};
            end
            opSetWide = {@sbp.D2Standard, @sbp.D2Standard};

            for i = 1:nBlocks
                if iscell(rho)
                    doParamUpwind{i} = {rho{i}, C_upwind{i}, opSetUpwind, optFlag};
                    doParamNarrow{i} = {rho{i}, C_narrow{i}, opSet, optFlag};
                    doParamWide{i} = {rho{i}, C_wide{i}, opSetWide, optFlag};
                elseif isa(rho, 'function_handle')
                    doParamUpwind{i} = {rho, C_upwind, opSetUpwind, optFlag};
                    doParamNarrow{i} = {rho, C_narrow, opSet, optFlag};
                    doParamWide{i} = {rho, C_wide, opSetWide, optFlag};
                else
                    error('Inconsistent format for C and rho.');
                end
            end

            diffOps = struct;
            diffOps.narrow = multiblock.DiffOp(@scheme.Elastic2dCurvilinearAnisotropic, g, order, doParamNarrow);
            diffOps.upwind = multiblock.DiffOp(@scheme.Elastic2dCurvilinearAnisotropicUpwind, g, order, doParamUpwind);
            diffOps.wide = multiblock.DiffOp(@scheme.Elastic2dCurvilinearAnisotropicUpwind, g, order, doParamWide);

            obj.diffOp = diffOps.narrow;

            DNarrow = diffOps.narrow.D;
            DUpwind = diffOps.upwind.D;
            DWide = diffOps.wide.D;

            H = diffOps.narrow.H;
            obj.H = kron(H,speye(2));

            RHO = cell(nBlocks, nBlocks);
            RHO_kron = cell(nBlocks, nBlocks);
            E1 = cell(nBlocks, nBlocks);
            E2 = cell(nBlocks, nBlocks);

            C = cell(dim, dim, dim, dim);
            sigma = cell(dim, dim);
            for i = 1:dim
                for j = 1:dim
                    sigma{i,j} = cell(nBlocks, nBlocks);
                    for k = 1:dim
                        for l = 1:dim
                            C{i,j,k,l} = cell(nBlocks, nBlocks);
                        end
                    end
                end
            end

            I_dim = speye(2,2);
            for b = 1:nBlocks
                RHO_kron{b,b} = kron(diffOps.narrow.diffOps{b}.RHO, I_dim);
                RHO{b,b} = diffOps.narrow.diffOps{b}.RHO;
                E1{b,b} = diffOps.narrow.diffOps{b}.E{1};
                E2{b,b} = diffOps.narrow.diffOps{b}.E{2};
                for i = 1:dim
                    for j = 1:dim
                        sigma{i,j}{b,b} = diffOps.narrow.diffOps{b}.sigma{i,j} + diffOps.upwind.diffOps{b}.sigma{i,j} + diffOps.wide.diffOps{b}.sigma{i,j};
                        for k = 1:dim
                            for l = 1:dim
                                C{i,j,k,l}{b,b} = diffOps.narrow.diffOps{b}.C{i,j,k,l} + diffOps.upwind.diffOps{b}.C{i,j,k,l} + diffOps.wide.diffOps{b}.C{i,j,k,l};
                            end
                        end
                    end
                end
            end
            obj.RHO_kron = blockmatrix.toMatrix(RHO_kron);
            obj.RHO = blockmatrix.toMatrix(RHO);

            obj.C = cell(dim, dim, dim, dim);
            obj.sigma = cell(dim, dim);
            for i = 1:dim
                for j = 1:dim
                    obj.sigma{i,j} = blockmatrix.toMatrix(sigma{i,j});
                    for k = 1:dim
                        for l = 1:dim
                            obj.C{i,j,k,l} = blockmatrix.toMatrix(C{i,j,k,l});
                        end
                    end
                end
            end

            obj.Ecomp{1} = blockmatrix.toMatrix(E1);
            obj.Ecomp{2} = blockmatrix.toMatrix(E2);

            RHO_kron = obj.RHO_kron;
            C = obj.C;
            RHO = obj.RHO;
            Ec = obj.Ecomp;

            %---- Forcing function F ----
            % Put the components of F in a vector, for each block
            F_comb = cell(nBlocks,1);
            for i = 1:nBlocks
                F_comb{i} = @(t,x,y) [F{i}{1}(t,x,y); F{i}{2}(t,x,y)];
            end

            % Evaluate F
            if ~isempty(F)
                Ft = @(t) RHO_kron\multiblock.evalOn(g, F_comb, t);
            else
                Ft = [];
            end
            %---------------------------

            % Set BC
            penalty = {};
            for i = 1:numel(bc)

                [closureNarrow, penaltyNarrow] = scheme.bcSetup(diffOps.narrow, {bc{i}});
                [closureUpwind, penaltyUpwind] = scheme.bcSetup(diffOps.upwind, {bc{i}});
                [closureWide, penaltyWide] = scheme.bcSetup(diffOps.wide, {bc{i}});

                type = bc{i}.type{2};
                switch type
                    case {'t', 'traction'}
                        penalty{end+1} = penaltyNarrow;
                    case {'d', 'displacement'}
                        penalty{end+1} = penaltyNarrow;
                        penalty{end+1} = penaltyUpwind;
                        penalty{end+1} = penaltyWide;
                end
                DNarrow = DNarrow + closureNarrow;
                DUpwind = DUpwind + closureUpwind;
                DWide = DWide + closureWide;
            end

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

                    % Loop over components
                    for i = 1:dim

                        % Project source for component i to full grid size and accumulate different sources.
                        if isa(pointSources.g{s}{i}, 'function_handle')
                            source_func_cont = @(t) source_func_cont(t) + ...
                                               RHO_kron \ (pointSources.g{s}{i}(t)*Ec{i}*delta_fun);
                            contSource = true;
                        else
                            source_vec = pointSources.g{s}{i};
                            if iscolumn(source_vec)
                                source_vec = transpose(source_vec);
                            end
                            source_func_discr = source_func_discr + ...
                                               RHO_kron \ kron(source_vec, Ec{i}*delta_fun);

                            discrSource = true;
                        end
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
            data_funcs = {Ft, penalty{:}, ps_cont};
            S = elastic.addFunctionHandles(data_funcs);
            obj.S_cont = S;
            obj.S_discr = ps_discr;

            % Zero initial data, two components of displacement
            v0_fun = @(x,y) [0*x; 0*y];
            v0t_fun = @(x,y) [0*x; 0*y];

            obj.v0 = grid.evalOn(g, v0_fun);
            obj.v0t = grid.evalOn(g, v0t_fun);

            % Misc
            diffOps.D = DNarrow + DUpwind + DWide;
            obj.diffOps = diffOps;
            obj.D = diffOps.D;
            obj.order = order;

            % ------- Supergrid ------- %
            defaultSuperGrid.gamma = {0, 1e-1, 0};
            defaultSuperGrid.DOI_IDs = 1:g.nBlocks();
            defaultSuperGrid.W_IDs = [];
            defaultSuperGrid.E_IDs = [];
            defaultSuperGrid.S_IDs = [];
            defaultSuperGrid.N_IDs = [];
            default_struct('superGrid', defaultSuperGrid);
            [gamma, DOI_IDs] = dealStruct(superGrid, {'gamma', 'DOI_IDs'});

            %--- Estimate max wave speed ----
            % Compute speed in a few different directions
            theta_range = linspace(0,pi,5);
            v_max = zeros(length(RHO), length(theta_range));
            for th = 1:length(theta_range)
                theta = theta_range(th);
                y = [cos(theta), sin(theta)];

                R = cell(dim,dim);
                for j = 1:dim
                    for l = 1:dim
                        R{j,l} = 0*C{1,1,1,1};
                    end
                end

                % R_jl = y_i C_ijkl y_k
                for i = 1:dim
                    for j = 1:dim
                        for k = 1:dim
                            for l = 1:dim
                                R{j,l} = R{j,l} + y(i)*C{i,j,k,l}*y(k);
                            end
                        end
                    end
                end

                % Solve quadratic for v^2.
                p = R{1,1} + R{2,2};
                q = R{1,1}*R{2,2} - R{1,2}*R{2,1};
                p = diag(p);
                q = diag(q);
                rho = diag(RHO);
                v_max_2 = 1./rho .* (p/2 + sqrt(1/4*p.^2 - q) );
                v_max(:,th) = sqrt(max(v_max_2));
            end
            v = max(v_max, [], 2);
            v = spdiag(v);

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
            Jref = obj.diffOps.narrow.diffOps{DOI_IDs(1)}.J;
            Jref = mean(diag(Jref));

            D1xi = blockmatrix.toMatrix(D1xi);
            D1eta = blockmatrix.toMatrix(D1eta);
            D2xi = blockmatrix.toMatrix(D2xi);
            D2eta = blockmatrix.toMatrix(D2eta);
            D4xi = blockmatrix.toMatrix(D4xi);
            D4eta = blockmatrix.toMatrix(D4eta);
            Hxi = blockmatrix.toMatrix(Hxi);
            Heta = blockmatrix.toMatrix(Heta);

            H = diffOps.narrow.H;

            E = - gamma{1} * Jref^(-1/dim) * ( Hxi \ (Hxi.^2 * (H\ (D1xi' * H* v * D1xi)))) ...
                - gamma{2} * Jref^(-1/dim) * ( Hxi \ (Hxi.^4 * (H\ (D2xi' * H* v * D2xi)))) ...
                - gamma{3} * Jref^(-1/dim) * ( Hxi \ (Hxi.^8 * (H\ (D4xi' * H* v * D4xi)))) ...
                - gamma{1} * Jref^(-1/dim) * ( Heta \ (Heta.^2 * (H\ (D1eta' * H* v * D1eta)))) ...
                - gamma{2} * Jref^(-1/dim) * ( Heta \ (Heta.^4 * (H\ (D2eta' * H* v * D2eta)))) ...
                - gamma{3} * Jref^(-1/dim) * ( Heta \ (Heta.^8 * (H\ (D4eta' * H* v * D4eta)))) ;
            %------------------------------------------ %

            obj.E = Ec{1}*E*Ec{1}' + Ec{2}*E*Ec{2}';
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
        function [ts, N] = getTimestepper(obj, method, time_align, k)
            default_arg('method','rk4');
            default_arg('time_align',[]);
            default_arg('k',[]);
            switch method
                case 'rk4'
                    if isempty(k)
                        k = obj.getTimestep(method);
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
                switch obj.order
                case 2
                    cfl = 0.8;
                case 4
                    cfl = 0.6;
                case 6
                    cfl = 0.4;
                end
            end
            DOI_IDs = obj.superGrid.DOI_IDs;
            nBlocks = length(DOI_IDs);
            dt = zeros(nBlocks, 1);
            for ii = 1:length(DOI_IDs)
                b = DOI_IDs(ii);

                % Get transformed density and stiffness tensor
                rho = diag(obj.diffOps.narrow.diffOps{b}.refObj.RHO);
                CNarrow = obj.diffOps.narrow.diffOps{b}.refObj.C;
                CUpwind = obj.diffOps.upwind.diffOps{b}.refObj.C;
                CWide = obj.diffOps.wide.diffOps{b}.refObj.C;

                dim = 2;
                for i = 1:dim
                    for j = 1:dim
                        for k = 1:dim
                            for l = 1:dim
                                C{i,j,k,l} = diag(CNarrow{i,j,k,l} + CUpwind{i,j,k,l} + CWide{i,j,k,l});
                            end
                        end
                    end
                end

                % Compute speed in a few different directions
                theta_range = linspace(0,pi/2,5);
                dt_local = zeros(size(theta_range));
                for th = 1:length(theta_range)
                    theta = theta_range(th);
                    y = [cos(theta), sin(theta)];

                    R = cell(dim,dim);
                    for j = 1:dim
                        for l = 1:dim
                            R{j,l} = 0*C{1,1,1,1};
                        end
                    end

                    % R_jl = y_i C_ijkl y_k
                    for i = 1:dim
                        for j = 1:dim
                            for k = 1:dim
                                for l = 1:dim
                                    R{j,l} = R{j,l} + y(i)*C{i,j,k,l}*y(k);
                                end
                            end
                        end
                    end

                    % Solve quadratic for v^2.
                    p = R{1,1} + R{2,2};
                    q = R{1,1}.*R{2,2} - R{1,2}.*R{2,1};
                    v_max_2 = 1./rho .* (p/2 + sqrt( 1/4*p.^2 - q ) );
                    v_max = sqrt(max(v_max_2));

                    m = obj.grid.grids{b}.logic.m;
                    hxi = 1/(m(1)-1);
                    heta = 1/(m(2)-1);
                    h = min(hxi,heta);

                    dt_local(th) = cfl/v_max*h;
                end
                dt(ii) = min(dt_local);
            end
            dt = min(dt);
            k = dt;
        end

        function r = getTimeSnapshot(obj, ts)
            if ts == 0
                r.t = 0;
                r.v = obj.v0;
                r.vt = obj.v0t;
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
            elastic.cmap;
            v0 = obj.v0t;
            v0_1 = obj.Ecomp{1}'*v0;
            v0_2 = obj.Ecomp{2}'*v0;

            % Loop through grids and compress them
            % Increase magnitude of solution inside SG layer
            factor = 100;
            v1_split = g.splitFunc(v0_1);
            v2_split = g.splitFunc(v0_2);
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
                    v1_split{i} = 1e3*v1_split{i};
                    v2_split{i} = 1e3*v2_split{i};
                end

            end
            v1 = cell2mat(v1_split);
            v2 = cell2mat(v2_split);

            % u1, full domain
            h1 = subplot(2,1,1);
            Sur1 = multiblock.Surface(g, v1);
            xlabel('x')
            ylabel('y')
            shading interp
            colorbar
            xlims1 = xlim;
            ylims1 = ylim;
            axis equal
            a = gca;

            % u2, full domain
            h2 = subplot(2,1,2);
            Sur2 = multiblock.Surface(g, v2);
            xlabel('x')
            ylabel('y')
            shading interp
            colorbar
            xlims2 = xlim;
            ylims2 = ylim;
            axis equal;

            % Figure 2, DOI only
            figure;
            elastic.cmap;
            v1_split = g.splitFunc(v0_1);
            v1_DOI = cell2mat(v1_split(DOI_IDs));
            v2_split = g.splitFunc(v0_2);
            v2_DOI = cell2mat(v2_split(DOI_IDs));
            if nDOIGrids > 1
                h1_DOI = subplot(2,1,1);
                Sur1_DOI = multiblock.Surface(g_DOI, v1_DOI);
                xlabel('x')
                ylabel('y')
                view(0,90)
                shading interp
                colorbar
                xlims1_DOI = xlim;
                ylims1_DOI = ylim;
                axis equal

                h2_DOI = subplot(2,1,2);
                Sur2_DOI = multiblock.Surface(g_DOI, v2_DOI);
                xlabel('x')
                ylabel('y')
                view(0,90)
                shading interp
                colorbar
                xlims2_DOI = xlim;
                ylims2_DOI = ylim;
                axis equal
            else
                p = g_DOI.points();
                X = grid.funcToPlotMatrix(g_DOI, p(:,1));
                Y = grid.funcToPlotMatrix(g_DOI, p(:,2));
                Z1 = grid.funcToPlotMatrix(g_DOI, v1_DOI);
                Z2 = grid.funcToPlotMatrix(g_DOI, v2_DOI);

                h1_DOI = subplot(2,1,1);
                Sur1_DOI = surf(X,Y,Z1);
                xlabel('x')
                ylabel('y')
                view(0,90)
                shading interp
                colorbar
                xlims1_DOI = xlim;
                ylims1_DOI = ylim;
                axis equal

                h2_DOI = subplot(2,1,2);
                Sur2_DOI = surf(X,Y,Z2);
                xlabel('x')
                ylabel('y')
                view(0,90)
                shading interp
                colorbar
                xlims2_DOI = xlim;
                ylims2_DOI = ylim;
                axis equal
            end

            function update_fun(r,E,xlims1,ylims1,xlims2,ylims2)
                t = r.t;
                v1 = E{1}'*r.vt;
                v2 = E{2}'*r.vt;
                if ishandle(a)
                    title(a,sprintf('T = %.3f',t))
                end

                v1_split = g.splitFunc(v1);
                v2_split = g.splitFunc(v2);
                for i = 1:g.nBlocks
                    if ~ismember(i, DOI_IDs)
                        v1_split{i} = 1e3*v1_split{i};
                        v2_split{i} = 1e3*v2_split{i};
                    end
                end
                v1 = cell2mat(v1_split);
                v2 = cell2mat(v2_split);

                Sur1.ZData = v1;
                Sur1.CData = v1;
                Sur2.ZData = v2;
                Sur2.CData = v2;
                xlim(h1, xlims1);
                ylim(h1, ylims1);
                xlim(h2, xlims2);
                ylim(h2, ylims2);

                cmax = 1;
                caxis(h1, [-cmax, cmax]);
                caxis(h2, [-cmax, cmax]);

                % Figure 2, DOI only
                v1_split = g.splitFunc(v1);
                v1_DOI = cell2mat(v1_split(DOI_IDs));
                v2_split = g.splitFunc(v2);
                v2_DOI = cell2mat(v2_split(DOI_IDs));
                if nDOIGrids > 1
                    Sur1_DOI.ZData = v1_DOI;
                    Sur2_DOI.ZData = v2_DOI;
                    Sur1_DOI.CData = v1_DOI;
                    Sur2_DOI.CData = v2_DOI;
                else
                    Sur1_DOI.ZData = grid.funcToPlotMatrix(g_DOI, v1_DOI);
                    Sur2_DOI.ZData = grid.funcToPlotMatrix(g_DOI, v2_DOI);
                    Sur1_DOI.CData = grid.funcToPlotMatrix(g_DOI, v1_DOI);
                    Sur2_DOI.CData = grid.funcToPlotMatrix(g_DOI, v2_DOI);
                end
                xlim(h1_DOI, xlims1_DOI);
                ylim(h1_DOI, ylims1_DOI);
                xlim(h2_DOI, xlims2_DOI);
                ylim(h2_DOI, ylims2_DOI);

                cmax = 1;
                caxis(h1_DOI, [-cmax, cmax]);
                caxis(h2_DOI, [-cmax, cmax]);
            end
            update = @(r)update_fun(r,obj.Ecomp,xlims1,ylims1,xlims2,ylims2);
        end

        % Compare two functions u and v in the discrete l2 norm.
        function e = compareSolutions(obj, u, v)
            H = obj.diffOps.narrow.H;
            I_dim = speye(2,2);
            H = kron(H, I_dim);

            evec = u - v;
            e = sqrt( evec'*H*evec );
        end

        % Compare the grid function u to the analytical function g in the discrete l2 norm.
        function e = compareSolutionsAnalytical(obj, u, g, t)
            gr = obj.grid;
            v = grid.evalOn(gr, @(x,y)g(t,x,y));
            e = obj.compareSolutions(u, v);
        end

    end

    methods(Static)

    end
end