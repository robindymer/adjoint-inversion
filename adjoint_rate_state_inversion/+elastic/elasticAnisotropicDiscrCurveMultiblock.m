classdef elasticAnisotropicDiscrCurveMultiblock < noname.Discretization
    properties
        name         = 'Anisotropic elastic 2D multiblock curvilinear'
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
        diffOp       %Various operators
        grid         %Multiblock grid in general

        RHO          %Density, diagonal matrix
        RHO_kron     %Kroneckered density
        C            %Elastic stiffness tensor

        Ecomp        % Ecomp{1}' picks out first component, Ecomp{2}' picks out 2nd comp.
    end

    methods

        % bc should be a cell array of structs, one struct per BC.
        % Each struct contains the fields
        %       - boundary: e.g. {1, 'w'} for west boundary, 1st block
        %       - type:     e.g. {2, 'd'} for displacement condition on the 2nd component
        %       - data:     boundary data function handle
        function obj = elasticAnisotropicDiscrCurveMultiblock(g, order, C, rho, F,...
                                             bc, pointSources, opSet, optFlag, intfTypes, hollow)

            default_arg('hollow', false);
            default_arg('intfTypes', []);
            default_arg('optFlag',[]);
            default_arg('opSet',[]);
            default_arg('F',[]);
            default_arg('pointSources',[]);
            default_arg('bc',[]);
            default_arg('rho',[]);
            default_arg('C',[]);
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

            for i = 1:nBlocks
                if iscell(rho)
                    doParam{i} = {rho{i}, C{i}, opSet, optFlag, hollow};
                elseif isa(rho, 'function_handle')
                    doParam{i} = {rho, C, opSet, optFlag, hollow};
                else
                    error('Inconsistent format for C and rho.');
                end
            end

            diffOp = multiblock.DiffOp(@scheme.Elastic2dCurvilinearAnisotropic, g, order, doParam, intfTypes);
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

            I_dim = speye(2,2);
            for b = 1:nBlocks
                RHO_kron{b,b} = kron(diffOp.diffOps{b}.RHO, I_dim);
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

            %---- Forcing function F ----
            if ~isempty(F)
                % Put the components of F in a vector, for each block
                F_comb = cell(nBlocks,1);
                for i = 1:nBlocks
                    F_comb{i} = @(t,x,y) [F{i}{1}(t,x,y); F{i}{2}(t,x,y)];
                end

                % Evaluate F
                Ft = @(t) RHO_kron\multiblock.evalOn(g, F_comb, t);
            else
                Ft = [];
            end
            %---------------------------

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
            data_funcs = {Ft, penalty, ps_cont};
            S = elastic.addFunctionHandles(data_funcs);
            obj.S_cont = S;
            obj.S_discr = ps_discr;

            % Zero initial data, two components of displacement
            v0_fun = @(x,y) [0*x; 0*y];
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
                    cfl = 0.5;
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
            default_arg('cfl',0.5);

            nBlocks = length(obj.grid.grids);
            dt = zeros(nBlocks, 1);
            for b = 1:nBlocks

                % Get transformed density and stiffness tensor
                rho = diag(obj.diffOp.diffOps{b}.refObj.RHO);
                C = obj.diffOp.diffOps{b}.refObj.C;

                dim = 2;
                for i = 1:dim
                    for j = 1:dim
                        for k = 1:dim
                            for l = 1:dim
                                C{i,j,k,l} = diag(C{i,j,k,l});
                            end
                        end
                    end
                end

                % Compute speed in a few different directions
                theta_range = linspace(0,pi,5);
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

                    m = obj.grid.grids{i}.logic.m;
                    hxi = 1/(m(1)-1);
                    heta = 1/(m(2)-1);
                    h = min(hxi,heta);

                    dt_local(th) = cfl/v_max*h;
                end
                dt(b) = min(dt_local);
            end
            k = min(dt);
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
            xlims1 = xlim;
            ylims1 = ylim;


            h2 = subplot(2,1,2);
            Sur2 = multiblock.Surface(obj.grid, v0_2);
            xlabel('x')
            ylabel('y')
            shading interp
            colorbar
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
                % caxis(h1, [-1,1]);
                % caxis(h2, [-1,1]);
                % cmax = max(abs(Sur1.ZData));
                % caxis(h1, [-cmax, cmax]);
                % cmax = max(abs(Sur2.ZData));
                % caxis(h2, [-cmax, cmax]);
            end
            update = @(r)update_fun(r,obj.Ecomp,xlims1,ylims1,xlims2,ylims2);
        end

        % Compare two functions u and v in the discrete l2 norm.
        function e = compareSolutions(obj, u, v)
            H = obj.diffOp.H;
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