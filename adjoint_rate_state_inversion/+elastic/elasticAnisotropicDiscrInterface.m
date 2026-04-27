classdef elasticAnisotropicDiscrInterface < noname.Discretization
    properties
        name         = 'elastic 2D multiblock with point source'
        description  = 'rho*u_{i,tt} = d_j a d_i u_j + d_j a d_j u_i + F(x,y,t)'
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
        diffOp       %Various operators
        grid         %Multiblock grid in general

        RHO          %Density, diagonal matrix
        RHO_kron     %Kroneckered density
        C            %Elastic stiffness tensor

        Div          %Divergence operator
        Ecomp        % Ecomp{1}' picks out first component, Ecomp{2}' picks out 2nd comp.
    end

    methods

        % boundary_data should be a cell array of structs, one struct per BC.
        % Each struct contains the fields
        %       - boundary: e.g. {1, 'w'} for west boundary, 1st block
        %       - type:     e.g. {2, 'd'} for displacement condition on the 2nd component
        %       - data:     boundary data function handle

        % point_sources should be a struct containing cell arrays with one element per source
        %       - x:    cell array of coordinate vectors specifying source locations,
        %               e.g. x = {[x1, y1], [x2, y2]}
        %       - g:    cell array of cell arrays of scalar source time function handles. An array {gx, gy}
        %               specifies source strength gx in x-dir and gy in y-dir. Example for two sources:
        %               g = {{g1x, g1y}, {g2x, g2y}}
        function obj = elasticAnisotropicDiscrInterface(m, order, xlim, ylim, C, rho, F,...
                                             boundary_data, point_sources, opSet, optFlag)
            default_arg('optFlag',[]);
            default_arg('opSet',[]);
            default_arg('point_sources',[]);
            default_arg('F',[]);
            default_arg('boundary_data',[]);
            default_arg('rho',[]);
            default_arg('C',[]);
            default_arg('xlim',[-1,0,1]);
            default_arg('ylim',[1,0]);
            default_arg('order',4);
            dim = 2;
            obj.dim = dim;

            % Domain definition, grid, diffOp
            domain = multiblock.domain.Rectangle(xlim, ylim);
            g = domain.getGrid(m);
            obj.grid = g;

            nBlocks = g.nBlocks();
            doParam = cell(nBlocks, 1);

            if isa(rho, 'double') && ~isempty(rho);
                rho = g.splitFunc(rho);
            end

            for i = 1:nBlocks
                if iscell(rho)
                    doParam{i} = {rho{i}, C{i}, opSet, optFlag};
                elseif isa(rho, 'function_handle')
                    doParam{i} = {rho, C, opSet, optFlag};
                else
                    error('Inconsistent format for C and rho.');
                end
            end

            diffOp = multiblock.DiffOp(@scheme.Elastic2dVariableAnisotropic, g, order, doParam);
            D = diffOp.D;

            H = diffOp.H;
            obj.H = kron(H,speye(2));

            RHO = cell(nBlocks, nBlocks);
            RHO_kron = cell(nBlocks, nBlocks);
            E1 = cell(nBlocks, nBlocks);
            E2 = cell(nBlocks, nBlocks);

            C = cell(dim, dim, dim, dim);
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            C{i,j,k,l} = cell(nBlocks, nBlocks);
                        end
                    end
                end
            end

            for b = 1:nBlocks
                RHO_kron{b,b} = inv( diffOp.diffOps{b}.RHOi_kron );
                RHO{b,b} = diffOp.diffOps{b}.RHO;
                E1{b,b} = diffOp.diffOps{b}.E{1};
                E2{b,b} = diffOp.diffOps{b}.E{2};
                for i = 1:dim
                    for j = 1:dim
                        for k = 1:dim
                            for l = 1:dim
                                C{i,j,k,l}{b,b} = diffOp.diffOps{b}.C{i,j,k,l};
                            end
                        end
                    end
                end
            end
            obj.RHO_kron = blockmatrix.toMatrix(RHO_kron);
            obj.RHO = blockmatrix.toMatrix(RHO);

            obj.C = cell(dim, dim, dim, dim);
            for i = 1:dim
                for j = 1:dim
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

            % First derivatives
            D1 = cell(dim, 1);
            for j = 1:dim
                D1{j} = cell(nBlocks, nBlocks);
                for i = 1:nBlocks
                    D1{j}{i,i} = diffOp.diffOps{i}.D1{j};
                end
                D1{j} = blockmatrix.toMatrix(D1{j});
            end

            % Divergence operator
            obj.Div = D1{1}*Ec{1}' + D1{2}*Ec{2}';

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

            %----- Point sources -------
            if ~isempty(point_sources)

                x_s = point_sources.x;
                n_s = length(x_s);

                % Initialize
                source_func_cont = @(t) 0*t;
                source_func_discr = {};

                for s = 1:n_s

                    % Check which grid block source is in
                    blockId = [];
                    for i = 1:nBlocks-1
                        xlim = g.grids{i}.lim{1};
                        ylim = g.grids{i}.lim{2};
                        if ( xlim{1} <= x_s{s}(1) && x_s{s}(1) < xlim{2} && ...
                             ylim{1} <= x_s{s}(2) && x_s{s}(2) < ylim{2} )
                            blockId = i;
                        end
                    end
                    if isempty(blockId)
                        blockId = nBlocks;
                    end

                    % Operators local to block number blockId
                    H_local = diffOp.diffOps{blockId}.H_1D;

                    % Delta function corresponding to local grid
                    delta_fun_local = diracDiscr(x_s{s}, obj.grid.grids{blockId}.x, order, 0, H_local);

                    % Extend delta fun to global grid.
                    delta_fun = sparse(obj.grid.expandFunc(delta_fun_local, blockId));

                    % Loop over components
                    for i = 1:dim
                        % Project source for component i to full grid size and accumulate different sources.
                        if isa(point_sources.g{s}{i}, 'function_handle')
                            source_func_cont = @(t) source_func_cont(t) + ...
                                               RHO_kron\(point_sources.g{s}{i}(t)*obj.Ecomp{i}*delta_fun);
                        else
                            source_vec = sparse(point_sources.g{s}{i});
                            if iscolumn(source_vec)
                                source_vec = transpose(source_vec);
                            end
                            source_func_discr{end+1} = RHO_kron\kron(source_vec, obj.Ecomp{i}*delta_fun);
                        end
                    end
                end
                % If there were no discrete source functions, make it empty
                if isempty(source_func_discr)
                    source_func_discr = [];
                else
                    sourceCell = source_func_discr;
                    source_func_discr = sourceCell{1};
                    for i = 2:numel(sourceCell)
                        source_func_discr = source_func_discr + sourceCell{i};
                    end
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
            %----------------------------

            % Set BC
            [closure, penalty] = scheme.bcSetup(diffOp, boundary_data);
            D = D + closure;

            % Create data function S(t)
            data_funcs = cell(3,1);
            data_funcs{1} = Ft;
            data_funcs{2} = penalty;
            data_funcs{3} = source_func_cont;

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

            % Misc
            diffOp.D = D;
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
        function [ts, N] = getTimestepper(obj, method, time_align, k)
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
                nBlocks = length(obj.grid.grids);
                k = zeros(nBlocks, 1);
                for i = 1:nBlocks

                    % Very crude estimates of wave speeds
                    rho = diag(obj.diffOp.diffOps{i}.RHO);
                    C = obj.diffOp.diffOps{i}.C;
                    C1111 = diag(C{1,1,1,1});
                    C2222 = diag(C{2,2,2,2});

                    v_1 = sqrt(max( abs(C1111./rho) ));
                    v_2 = sqrt(max( abs(C2222./rho) ));
                    v = max(v_1,v_2);

                    h = min(obj.grid.grids{i}.h);
                    k(i) = cfl/v*h;
                end
                k = min(k);
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
            v0_1 = obj.Ecomp{1}'*v0;
            v0_2 = obj.Ecomp{2}'*v0;

            h1 = subplot(2,1,1);
            Sur1 = multiblock.Surface(obj.grid, v0_1);
            xlabel('x')
            ylabel('y')
            shading interp
            colorbar
            axis equal
            xlims1 = xlim;
            ylims1 = ylim;


            h2 = subplot(2,1,2);
            Sur2 = multiblock.Surface(obj.grid, v0_2);
            xlabel('x')
            ylabel('y')
            shading interp
            colorbar
            axis equal
            xlims2 = xlim;
            ylims2 = ylim;

            a = gca;

            function update_fun(r,E,xlims1,ylims1,xlims2,ylims2)
                t = r.t;
                v1 = E{1}'*r.v;
                v2 = E{2}'*r.v;
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
                % cmax = max(abs(Sur1.ZData));
                % caxis(h1, [-cmax, cmax]);
                % cmax = max(abs(Sur2.ZData));
                % caxis(h2, [-cmax, cmax]);
            end
            update = @(r)update_fun(r,obj.Ecomp,xlims1,ylims1,xlims2,ylims2);
        end

        % Generates delta functions.
        % Source coordinates should be a cell array of coordinate vectors.
        % Returns a (nSources x dim) cell array of delta functions.
        function deltaFunctions = generateDeltaFunctions(obj, sourceCoord, nMoment, nSmooth)
            default_arg('nMoment', obj.order);
            default_arg('nSmooth', 0);

            nSources = numel(sourceCoord);
            deltaFunctions = cell(nSources, obj.dim);

            nBlocks = obj.grid.nBlocks();
            E = obj.Ecomp;

            for s = 1:nSources
                x_s = sourceCoord{s};

                % Check which grid block source is in
                blockId = [];
                for i = 1:nBlocks-1
                    xlim = obj.grid.grids{i}.lim{1};
                    ylim = obj.grid.grids{i}.lim{2};
                    if ( xlim{1} <= x_s(1) && x_s(1) < xlim{2} && ...
                         ylim{1} <= x_s(2) && x_s(2) < ylim{2} )
                    blockId = i;
                    end
                end
                if isempty(blockId)
                    blockId = nBlocks;
                end

                % Operators local to block number blockId
                H_local = obj.diffOp.diffOps{blockId}.H_1D;

                % Delta function corresponding to local grid
                deltaFunLocal = diracDiscr(x_s, obj.grid.grids{blockId}.x, nMoment, nSmooth, H_local);

                % Extend delta fun to global grid.
                deltaFun = obj.grid.expandFunc(deltaFunLocal, blockId);

                % Extend to all components
                for d = 1:obj.dim
                    deltaFunctions{s,d} = E{d}*deltaFun;
                end
            end
        end

        % Compare two functions u and v in the discrete l2 norm.
        function e = compareSolutions(obj, u, v)
            g = obj.grid;
            nBlocks = length(g.grids);

            evec = u - v;
            ecell = g.splitFunc(evec);

            e = 0;
            for i = 1:nBlocks
                subgrid = g.grids{i};
                h = subgrid.h;
                e = e + sqrt(prod(h))*norm(ecell{i});
            end
        end

        % Compare the grid function u to the analytical function g in the discrete l2 norm.
        function e = compareSolutionsAnalytical(obj, u, g, t)
            gr = obj.grid;
            % v = multiblock.evalOn(gr, g, t);
            v = grid.evalOn(gr, @(x,y)g(t,x,y));
            e = obj.compareSolutions(u, v);
        end

    end

    methods(Static)

    end
end