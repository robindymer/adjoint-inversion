classdef elasticAnisotropicStaggeredDiscr < noname.Discretization
    properties
        name         = 'Anisotropic elastic 2D staggered Cartesian grid'
        description  = 'rho*u_{j,tt} = d_i C_ijkl d_k u_l + F(x,y,t)'
        order        %Order of accuracy
        dim          %Number of spatial dimensions

        D            %Discretization matrix including BC
        E            %Matrix for u_t
        H            % Quadrature for both comp of u, both grids

        % v_tt = D*v + E*v_t + S_cont(t) + S_discr
        S_cont       %Function handle, v_tt = D*v + E*v_t + S_cont(t) + S_discr
        S_discr      %Matrix of time-dependent data, one column per time-stepping stage.

        v0           %Initial data
        v0t          %Initial data for time derivative
        diffOp       %Various operators
        grid         %Multiblock grid in general

        rho          %Density, function handle
        RHO          %Density, diagonal matrix
        RHO_kron     %Density, both components, both grids
        C            %Elastic stiffness tensor

        U
        G
    end

    methods

        % bc should be a cell array of structs, one struct per BC.
        % Each struct contains the fields
        %       - boundary: e.g. {1, 'w'} for west boundary, 1st block
        %       - type:     e.g. {2, 'd'} for displacement condition on the 2nd component
        %       - data:     boundary data function handle
        function obj = elasticAnisotropicStaggeredDiscr(g, order, C, rho, F, bc, opSet)
            default_arg('opSet',[]);
            default_arg('F',[]);
            default_arg('bc',[]);
            default_arg('rho',[]);
            default_arg('C',[]);
            default_arg('order',4);
            dim = 2;
            nGrids = 2;
            obj.dim = dim;

            % Domain definition, grid, diffOp
            obj.grid = g;
            g_u = g.gridGroups{1};
            diffOp = scheme.Elastic2dStaggeredAnisotropic(g, order, rho, C);

            D = diffOp.D;
            obj.C = diffOp.C;
            G = diffOp.G;
            U = diffOp.U;
            RHO = diffOp.RHO;
            obj.RHO = RHO;
            obj.rho = rho;
            obj.U = U;
            obj.G = G;
            obj.H = diffOp.H;

            obj.RHO_kron = G{1}*(U{1}{1}*RHO{1}*U{1}{1}' + U{1}{2}*RHO{1}*U{1}{2}')*G{1}'...
                         + G{2}*(U{2}{1}*RHO{2}*U{2}{1}' + U{2}{2}*RHO{2}*U{2}{2}')*G{2}';

            %---- Forcing function F ----
            % Put the components of F in a vector
            % F_comb = @(t,x,y) [F{1}(t,x,y); F{2}(t,x,y)];

            % Evaluate F
            if ~isempty(F)
                Ft = @(t) G{1}*U{1}{1}*(RHO{1}\ grid.evalOn(g_u{1}, @(x,y)F{1}(t,x,y))) ...
                         +G{1}*U{1}{2}*(RHO{1}\ grid.evalOn(g_u{1}, @(x,y)F{2}(t,x,y))) ...
                         +G{2}*U{2}{1}*(RHO{2}\ grid.evalOn(g_u{2}, @(x,y)F{1}(t,x,y))) ...
                         +G{2}*U{2}{2}*(RHO{2}\ grid.evalOn(g_u{2}, @(x,y)F{2}(t,x,y)));
            else
                Ft = [];
            end
            %---------------------------

            % Set BC
            [closure, penalty] = scheme.bcSetupStaggered(diffOp, bc);
            D = D + closure;

            % -----------------------------------

            % Create data function S(t)
            data_funcs = {Ft, penalty};
            S = elastic.addFunctionHandles(data_funcs);
            obj.S_cont = S;
            obj.S_discr = [];

            % Zero initial data, two components of displacement
            v0_fun = @(x,y) [0*x; 0*y];
            v0t_fun = @(x,y) [0*x; 0*y];

            obj.v0 = G{1}*grid.evalOn(g_u{1}, v0_fun) + G{2}*grid.evalOn(g_u{2}, v0_fun);
            obj.v0t = G{1}*grid.evalOn(g_u{1}, v0t_fun) + G{2}*grid.evalOn(g_u{2}, v0t_fun);

            % Misc
            diffOp.D = D;
            obj.diffOp = diffOp;
            obj.D = D;
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
                    cfl = 0.2;
                case 4
                    cfl = 0.2;
                case 6
                    cfl = 0.2;
                end
            end

            dim = 2;

            % Get density and stiffness tensor
            rho = obj.rho;
            C = obj.C{1};

            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            C{i,j,k,l} = diag(C{i,j,k,l});
                        end
                    end
                end
            end

            % Get grid
            g = obj.grid.gridGroups{2}{1};
            RHO = grid.evalOn(g, rho);

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
                v_max_2 = 1./RHO .* (p/2 + sqrt( 1/4*p.^2 - q ) );
                v_max = sqrt(max(v_max_2));

                m = g.size();
                X = g.points();
                Lx = max(X(:,1)) - min(X(:,1));
                Ly = max(X(:,2)) - min(X(:,2));
                hx = Lx/(m(1) - 1);
                hy = Ly/(m(1) - 1);
                h = min(hx, hy);

                dt_local(th) = cfl/v_max*h;
            end
            dt = min(dt_local);
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

            % Get grid
            g = obj.grid.gridGroups{1}{1};

            figure_handle = figure();
            v = obj.G{1}'*obj.v0;
            v1 = obj.U{1}{1}'*v;
            v2 = obj.U{1}{2}'*v;

            v1 = grid.funcToPlotMatrix(g, v1);
            v2 = grid.funcToPlotMatrix(g, v2);
            X = g.matrices();

            % u1
            h1 = subplot(2,1,1);
            Sur1 = surf(X{1}', X{2}', v1);
            xlabel('x')
            ylabel('y')
            view(0,90)
            shading interp
            colorbar
            xlims1 = xlim;
            ylims1 = ylim;
            axis equal
            a = gca;

            % u2
            h2 = subplot(2,1,2);
            Sur2 = surf(X{1}', X{2}', v2);
            xlabel('x')
            ylabel('y')
            view(0,90)
            shading interp
            colorbar
            xlims2 = xlim;
            ylims2 = ylim;
            axis equal;

            function update_fun(r,G,U,xlims1,ylims1,xlims2,ylims2)
                t = r.t;
                v = G{1}'*r.v;
                v1 = U{1}{1}'*v;
                v2 = U{1}{2}'*v;
                v1 = grid.funcToPlotMatrix(g, v1);
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

                cmax = 1;
                caxis(h1, [-cmax, cmax]);
                caxis(h2, [-cmax, cmax]);

            end
            update = @(r)update_fun(r,obj.G,obj.U,xlims1,ylims1,xlims2,ylims2);
        end

        % Compare two functions u and v in the discrete l2 norm.
        function e = compareSolutions(obj, u, v)
            H = obj.diffOp.H_u{1};
            I_dim = speye(2,2);
            H = kron(H, I_dim);

            evec = u - v;
            e = sqrt( evec'*H*evec );
        end

        % Compare the grid function u to the analytical function g in the discrete l2 norm.
        function e = compareSolutionsAnalytical(obj, u, g, t)
            u = obj.G{1}'*u;
            gr = obj.grid.gridGroups{1}{1};
            v = grid.evalOn(gr, @(x,y)g(t,x,y));
            e = obj.compareSolutions(u, v);
        end

    end

    methods(Static)

    end
end