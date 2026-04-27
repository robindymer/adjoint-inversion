classdef viscoElasticQuasiStatic < noname.Discretization
    properties
        name         = '2d curvilinear viscoelastic quasi-static'
        description  = 'Momentum balance and viscous strain rate equations'
        order        %Order of accuracy
        dim          %Number of spatial dimensions
        diffOp       %The fully dynamic diffOp
        H_u

        v0           %Initial data
        grid         %Multiblock grid in general

        sigma_u, sigma_gamma %Stress operators

        Mech_u, Mech_gamma % Matrices in mechancical equilibrium equations
        Mech_u_factorized

        Dgu, Dgg    % Matrices in flow law

        mechForcing
        flowBodyForcing
        flowBoundaryForcing

        Eu, Egamma  % Pick out u/gamma from vector that contains both of them

        eta         % Struct for effective viscosity in flow law

    end

    methods

        % bc should be a cell array of structs, one struct per BC.
        % Each struct contains the fields
        %       - boundary: e.g. {1, 'w'} for west boundary, 1st block
        %       - type:     e.g. {2, 'd'} for displacement condition on the 2nd component
        %       - data:     boundary data function handle (or empty)
        %
        % Forcing: F.u and F.gamma
        function obj = viscoElasticQuasiStatic(g, order, C, eta, F,...
                                             bc, pointSources, intfTypes, computeLU)
            default_arg('computeLU', true);
            default_arg('intfTypes', []);
            default_arg('pointSources',[]);
            default_arg('F',[]);
            default_arg('bc',[]);
            default_arg('C',[]);
            default_arg('order',4);
            dim = 2;
            obj.dim = dim;

            rho = @(x,y) 0*x +1;

            % Default struct for eta (if a function handle, it is assumed to be independent of the solution)
            eta_default = struct;
            eta_default.type = 'linear';
            eta_default.fun = @(x,y) 0*x + 1;
            eta_default.nonlinFun = [];
            if isa(eta, 'function_handle')
                eta_default.fun = eta;
                eta = [];
            end
            default_struct('eta', eta_default);
            obj.eta = eta;

            % Domain definition, grid, diffOp
            obj.grid = g;

            %---- Build the dynamic discr and extract operators from it ----------
            dynamicDiscr = elastic.discrs.viscoElastic(g, order, C, rho, eta, F, bc, pointSources, intfTypes);
            dD = dynamicDiscr;

            % Mechanical equilibrium equation
            obj.Mech_u = -dD.Duu;
            obj.Mech_gamma = dD.Dug;

            % Factorize matrix
            Mech_u_factorized = struct;
            if computeLU
                tic
                [L, U, P, Q] = lu(obj.Mech_u);
                t = toc;
                fprintf('Time to compute LU factorization: %4.3es \n', t);

                Mech_u_factorized.L = L;
                Mech_u_factorized.U = U;
                Mech_u_factorized.P = P;
                Mech_u_factorized.Q = Q;
            end
            obj.Mech_u_factorized = Mech_u_factorized;

            % Flow law
            obj.Dgg = dD.Dgg;
            obj.Dgu = dD.Dgu;

            % Stress operators
            obj.sigma_u = dD.sigma_u;
            obj.sigma_gamma = dD.sigma_gamma;

            % Forcing
            obj.mechForcing = @(t) (dD.Ev)'*(dD.S_body(t) + dD.S_boundary(t));
            obj.flowBodyForcing = @(t) (dD.Egamma)' * dD.S_body(t);
            obj.flowBoundaryForcing = @(t) (dD.Egamma)' * dD.S_boundary(t);

            obj.Eu = dD.Eu_small;
            obj.Egamma = dD.Egamma_small;

            % Misc
            obj.order = order;
            obj.diffOp = dD.diffOp;
            obj.H_u = dD.H_u;

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

                    mms = struct;
                    if ~isempty(obj.mechForcing)
                        mms.mechForcing = obj.mechForcing;
                        mms.flowBoundaryForcing = obj.flowBoundaryForcing;
                        mms.flowBodyForcing = obj.flowBodyForcing;
                    else
                        mms = [];
                    end

                    ops = struct;
                    ops.Dgu = obj.Dgu;
                    ops.Dgg = obj.Dgg;

                    ops.Mech_u = obj.Mech_u;
                    ops.Mech_gamma = obj.Mech_gamma;
                    ops.Mech_u_factorized = obj.Mech_u_factorized;

                    ops.sigma_u = obj.sigma_u;
                    ops.sigma_gamma = obj.sigma_gamma;

                    F = @(U, t) elastic.visco.rhs.quasiStatic(t, U, ops, obj.eta.nonlinFun, mms);
                    ts = time.Rungekutta4proper(F, k, 0, obj.v0);

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
                case 3
                    cfl = 0.6;
                case 4
                    cfl = 0.6;
                case 5
                    cfl = 0.4;
                case 6
                    cfl = 0.4;
                end
            end
            DOI_IDs = 1:obj.grid.nBlocks();
            nBlocks = length(DOI_IDs);
            dt = zeros(nBlocks, 1);
            for ii = 1:length(DOI_IDs)
                b = DOI_IDs(ii);

                % Get transformed density and stiffness tensor
                rho = diag(obj.diffOp.diffOps{b}.elasticObj.refObj.RHO);
                C = obj.diffOp.diffOps{b}.elasticObj.refObj.C;
                for i = 1:numel(C)
                    C{i} = diag(C{i});
                end

                % Get grid spacings
                m = obj.grid.grids{b}.logic.m;
                hxi = 1/(m(1)-1);
                heta = 1/(m(2)-1);
                h = min(hxi,heta);

                % Compute speed in a few different directions
                theta_range = linspace(0, pi, 5);
                v_max = elastic.maxAnisotropicWaveSpeed(rho, C, theta_range);
                dt(ii) = cfl/v_max*h;
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

            % Get total grid and grid for DOI only
            g = obj.grid;
            DOI_IDs = 1:g.nBlocks();;
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

            % v0 = obj.v0;
            % v1 = v0(1:4:end-3);
            % v2 = v0(2:4:end-2);

            v0 = elastic.helpers.solveWithLU(obj.Mech_u_factorized, obj.Mech_gamma*obj.v0 + obj.mechForcing(0));
            v1 = v0(1:2:end-1);
            v2 = v0(2:2:end);

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

            function update_fun(r,xlims1,ylims1,xlims2,ylims2)
                t = r.t;

                % v = r.v;
                % v1 = v(1:4:end-3);
                % v2 = v(2:4:end-2);

                v0 = elastic.helpers.solveWithLU(obj.Mech_u_factorized, obj.Mech_gamma*r.v + obj.mechForcing(t));
                v1 = v0(1:2:end-1);
                v2 = v0(2:2:end);

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
                % caxis(h1, [-cmax, cmax]);
                % caxis(h2, [-cmax, cmax]);
            end
            update = @(r)update_fun(r,xlims1,ylims1,xlims2,ylims2);
        end

        % Compare two displacement functions u and v in the discrete l2 norm.
        function e = compareSolutions(obj, u, v)
            H = obj.H_u;

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