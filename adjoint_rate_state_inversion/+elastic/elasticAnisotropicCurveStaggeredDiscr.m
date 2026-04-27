classdef elasticAnisotropicCurveStaggeredDiscr < noname.Discretization
    properties
        name         = 'Staggered Anisotropic elastic 2D multiblock curvilinear'
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

        %RHO          %Density, diagonal matrix
        RHO_kron     %Kroneckered density
        %C            %Elastic stiffness tensor
        %sigma        %Elastic stress operator

        Ecomp        % Ecomp{1}' picks out first component, Ecomp{2}' picks out 2nd comp.
        superGrid
        stencil

        U
        G
    end

    methods

        % bc should be a cell array of structs, one struct per BC.
        % Each struct contains the fields
        %       - boundary: e.g. {1, 'w'} for west boundary, 1st block
        %       - type:     e.g. {2, 'd'} for displacement condition on the 2nd component
        %       - data:     boundary data function handle
        function obj = elasticAnisotropicCurveStaggeredDiscr(g, order, C, rho, F,...
                                                            bc, pointSources)
            default_arg('pointSources',[]);
            default_arg('F',[]);
            default_arg('bc',[]);
            default_arg('rho',[]);
            default_arg('C',[]);
            default_arg('order',4);

            nGrids = 2;
            dim = 2;
            obj.dim = dim;

            % Domain definition, grid, diffOp
            obj.grid = g;

            nBlocks = g.nBlocks();
            doParam = cell(nBlocks, 1);

            for i = 1:nBlocks
                if iscell(rho)
                    doParam{i} = {rho{i}, C{i}};
                elseif isa(rho, 'function_handle')
                    doParam{i} = {rho, C};
                else
                    error('Inconsistent format for C and rho.');
                end
            end

            diffOp = multiblock.DiffOp(@scheme.Elastic2dStaggeredCurvilinearAnisotropic, g, order, doParam);

            D = diffOp.D;
            H = diffOp.H;

            RHO_kron = cell(nBlocks, nBlocks);

            U11 = cell(nBlocks, nBlocks);
            U12 = cell(nBlocks, nBlocks);
            U21 = cell(nBlocks, nBlocks);
            U22 = cell(nBlocks, nBlocks);
            G1 = cell(nBlocks, nBlocks);
            G2 = cell(nBlocks, nBlocks);

            I_dim = speye(2,2);
            for b = 1:nBlocks
                RHO1 = diffOp.diffOps{b}.RHO{1};
                RHO2 = diffOp.diffOps{b}.RHO{2};

                RHO1 = kron(RHO1, I_dim);
                RHO2 = kron(RHO2, I_dim);

                RHO_temp = {RHO1, [];...
                            [], RHO2};
                RHO_kron{b,b} = blockmatrix.toMatrix(RHO_temp);

                G1{b,b} = diffOp.diffOps{b}.G{1};
                G2{b,b} = diffOp.diffOps{b}.G{2};
                U11{b,b} = diffOp.diffOps{b}.U{1}{1};
                U21{b,b} = diffOp.diffOps{b}.U{2}{1};
                U12{b,b} = diffOp.diffOps{b}.U{1}{2};
                U22{b,b} = diffOp.diffOps{b}.U{2}{2};
            end
            RHO_kron = blockmatrix.toMatrix(RHO_kron);
            obj.RHO_kron = RHO_kron;

            G1 = blockmatrix.toMatrix(G1);
            G2 = blockmatrix.toMatrix(G2);
            obj.G = {G1, G2};

            U11 = blockmatrix.toMatrix(U11);
            U12 = blockmatrix.toMatrix(U12);
            U21 = blockmatrix.toMatrix(U21);
            U22 = blockmatrix.toMatrix(U22);

            obj.U = cell(nGrids, 1);
            for a = 1:nGrids
                obj.U{a} = cell(dim, 1);
            end
            obj.U{1}{1} = U11;
            obj.U{2}{1} = U21;
            obj.U{1}{2} = U12;
            obj.U{2}{2} = U22;

            %---- Forcing function F ----
            % Assume same forcing everywhere

            % Evaluate F
            if ~isempty(F)
                F_comb = @(t,x,y) [F{1}(t,x,y); F{2}(t,x,y)];
                Ft = @(t) RHO_kron\grid.evalOnStaggered(g, @(x,y)F_comb(t,x,y));
            else
                Ft = [];
            end
            %---------------------------

            % Set BC
            [closure, penalty] = scheme.bcSetupStaggered(diffOp, bc);
            D = D + closure;

            %----- Point sources ----------------
            % To be implemented
            ps_discr = [];
            ps_cont = [];
            % -----------------------------------

            % Create data function S(t)
            data_funcs = {Ft, penalty, ps_cont};
            S = elastic.addFunctionHandles(data_funcs);
            obj.S_cont = S;
            obj.S_discr = ps_discr;

            % Zero initial data, two components of displacement
            v0_fun = @(x,y) [0*x; 0*y];
            v0t_fun = @(x,y) [0*x; 0*y];

            obj.v0 = grid.evalOnStaggered(g, v0_fun);
            obj.v0t = grid.evalOnStaggered(g, v0t_fun);

            % Misc
            diffOp.D = D;
            obj.diffOp = diffOp;
            obj.D = D;
            obj.H = H;
            obj.order = order;
            obj.E = [];


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
                    cfl = 0.6;
                case 4
                    cfl = 0.6;
                case 6
                    cfl = 0.6;
                end
            end

            nBlocks = obj.grid.nBlocks();
            dt = zeros(nBlocks, 1);

            for b = 1:nBlocks

                % Get transformed density and stiffness tensor
                rho = diag(obj.diffOp.diffOps{b}.refObj.RHO{1});
                C = obj.diffOp.diffOps{b}.refObj.C{1};

                % Interpolate rho to stress grid
                g_u = obj.grid.grids{b}.gridGroups{1}{1};
                g_s = obj.grid.grids{b}.gridGroups{2}{1};

                rhoFun = scatteredInterpolant(g_u.points(), full(rho), 'linear', 'linear');
                rho = rhoFun(g_s.points());

                % Get grid sizes
                m = g_s.logic.m;
                hxi = 1/(m(1)-1);
                heta = 1/(m(2)-1);
                h = min(hxi,heta);

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
                theta_range = linspace(0,pi/2,5);
                dt_local = zeros(size(theta_range));
                for th = 1:length(theta_range)

                    theta = theta_range(th);
                    v_max = elastic.anisotropicWaveSpeed(rho, C, theta);
                    v_max = max(v_max);

                    dt_local(th) = cfl/v_max*h;
                end
                dt(b) = min(dt_local);
            end
            dt = min(dt);
            k = dt;
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

            subGrid = 1;
            G = obj.G{subGrid};
            U1 = obj.U{subGrid}{1};
            U2 = obj.U{subGrid}{2};

            E1 = (U1'*G')';
            E2 = (U2'*G')';

            g = obj.grid;
            figure_handle = figure();
            v0 = obj.v0;
            v0_1 = E1'*v0;
            v0_2 = E2'*v0;

            h1 = subplot(2,1,1);
            Sur1 = multiblock.StaggeredSurface(g, v0_1);
            xlabel('x')
            ylabel('y')
            shading interp
            colorbar
            xlims1 = xlim;
            ylims1 = ylim;


            h2 = subplot(2,1,2);
            Sur2 = multiblock.StaggeredSurface(g, v0_2);
            xlabel('x')
            ylabel('y')
            shading interp
            colorbar
            xlims2 = xlim;
            ylims2 = ylim;

            a = gca;

            function update_fun(r,E1,E2,xlims1,ylims1,xlims2,ylims2)
                t = r.t;
                v1 = E1'*r.v;
                v2 = E2'*r.v;
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
            update = @(r)update_fun(r,E1,E2,xlims1,ylims1,xlims2,ylims2);
        end

        % Compare two functions u and v in the discrete l2 norm.
        function e = compareSolutions(obj, u, v)
            H = obj.diffOp.H;

            evec = u - v;
            e = sqrt( evec'*H*evec );
        end

        % Compare the grid function u to the analytical function g in the discrete l2 norm.
        function e = compareSolutionsAnalytical(obj, u, g, t)
            gr = obj.grid;
            v = grid.evalOnStaggered(gr, @(x,y)g(t,x,y));
            e = obj.compareSolutions(u, v);
        end

    end

    methods(Static)

    end
end