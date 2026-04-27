classdef viscoElastic < noname.Discretization
    properties
        name         = '2d curvilinear viscoelastic'
        description  = 'Momentum balance and viscous strain rate equations'
        order        %Order of accuracy
        dim          %Number of spatial dimensions

        D            %Discretization matrix including BC
        % E            %Matrix for u_t
        M            % First order in time: u_t = M*u + F
        H            %Total quadrature, for all components
        H_scalar     %Scalar quadrature
        H_u
        H_gamma

        % v_tt = D*v + E*v_t + S_cont(t) + S_discr
        S_cont       %Function handle, v_tt = D*v + E*v_t + S_cont(t) + S_discr
        S_discr      %Matrix of time-dependent data, one column per time-stepping stage.
        S_boundary
        S_body

        v0           %Initial data
        v0t          %Initial data for time derivative
        diffOp       %Various operators
        grid         %Multiblock grid in general

        % RHO          %Density, diagonal matrix
        % RHO_kron     %Kroneckered density
        C            %Elastic stiffness tensor
        sigma, sigma_u, sigma_gamma %Stress operators
        eta          %Struct for effective viscosity in the flow law

        Eu, Ev, Egamma   % Pick out all displacements/velocities/strains
        Eu_small, Egamma_small
        % eU, eGamma   % eU{i}' picks out displacement component i

        Duu, Dgu, Dug, Dgg % Blocks in total spatial operator
    end

    methods

        % bc should be a cell array of structs, one struct per BC.
        % Each struct contains the fields
        %       - boundary: e.g. {1, 'w'} for west boundary, 1st block
        %       - type:     e.g. {2, 'd'} for displacement condition on the 2nd component
        %       - data:     boundary data function handle (or empty)
        %
        % Forcing: F.u and F.gamma
        function obj = viscoElastic(g, order, C, rho, eta, F,...
                                             bc, pointSources, intfTypes)
            default_arg('intfTypes', []);
            default_arg('pointSources',[]);
            default_arg('F',[]);
            default_arg('bc',[]);
            default_arg('rho',[]);
            default_arg('C',[]);
            default_arg('order',4);
            dim = 2;
            obj.dim = dim;

            % Default struct for eta (if a function handle, it is assumed to be independent of the solution)
            eta_default = struct;
            eta_default.type = 'linear';
            eta_default.fun = @(x,y) 0*x + 1;
            eta_default.nonlinFun = [];
            if isa(eta, 'function_handle');
                eta_default.fun = eta;
                eta = [];
            end
            default_struct('eta', eta_default);
            obj.eta = eta;

            % Domain definition, grid, diffOp
            obj.grid = g;

            nBlocks = g.nBlocks();
            doParam = cell(nBlocks, 1);

            if isa(rho, 'double') && ~isempty(rho);
                rho = g.splitFunc(rho);
            end

            for i = 1:nBlocks
                if iscell(rho)
                    doParam{i} = {rho{i}, C{i}, eta.fun{i}};
                elseif isa(rho, 'function_handle')
                    doParam{i} = {rho, C, eta.fun};
                else
                    error('Inconsistent format for C and rho.');
                end
            end

            diffOp = multiblock.DiffOp(@scheme.ViscoElastic2d, g, order, doParam, intfTypes);
            D = diffOp.D;

            RHO_scalar = helpers.multiblockOperator(g, diffOp, 'RHO');
            Eu = helpers.multiblockOperator(g, diffOp, 'Eu');
            Egamma = helpers.multiblockOperator(g, diffOp, 'Egamma');
            obj.Eu_small = Eu;
            obj.Egamma_small = Egamma;
            % eU = helpers.multiblockTensor(g, diffOp, 'eU', dim);
            % eGamma = helpers.multiblockTwoTensor(g, diffOp, 'eGamma', dim);

            H_scalar = diffOp.H;
            I_u = speye(dim, dim);
            I_gamma = speye(dim^2, dim^2);

            H_u = kron(H_scalar, I_u);
            H_gamma = kron(H_scalar, I_gamma);

            obj.H = Eu*H_u*Eu' + Egamma*H_gamma*Egamma';

            RHO_u = kron(RHO_scalar, I_u);
            RHO_gamma = kron(RHO_scalar, I_gamma);

            C = helpers.multiblockStiffnessTensor(g, diffOp, dim);
            sigma = helpers.multiblockTwoTensor(g, diffOp, 'sigma', dim);

            sigma_u = cell(dim, dim);
            sigma_gamma = cell(dim, dim);
            for i = 1:dim
                for j = 1:dim
                    sigma_u{i,j} = sigma{i,j}*obj.Eu_small;
                    sigma_gamma{i,j} = sigma{i,j}*obj.Egamma_small;
                end
            end

            obj.C = C;
            obj.sigma = sigma;
            obj.sigma_u = sigma_u;
            obj.sigma_gamma = sigma_gamma;
            obj.H_u = H_u;
            obj.H_gamma = H_gamma;


            % Set BC
            [closure, penalty] = scheme.bcSetup(diffOp, bc);
            D = D + closure;

            % Write on first order form
            Duu = Eu'*D*Eu;
            Dug = Eu'*D*Egamma;
            Dgu = Egamma'*D*Eu;
            Dgg = Egamma'*D*Egamma;

            obj.Duu = Duu;
            obj.Dug = Dug;
            obj.Dgu = Dgu;
            obj.Dgg = Dgg;

            % Mu = H_u*RHO_u*Duu;
            % helpers.checkSymmetryAndEigenvalues(Mu);

            I = speye(dim*g.N(), dim*g.N());

            % Order: u-v-gamma
            M = {0*Duu,   I,   0*Dug;...
                   Duu,  0*Duu,  Dug;...
                   Dgu,  0*Dgu,  Dgg};
            M = blockmatrix.toMatrix(M);

            % Create corrected restriction operators
            [m, mu] = size(Eu);
            [~, mg] = size(Egamma);
            m = m + mu;

            Iu = speye(mu, mu);
            Ig = speye(mg, mg);
            Eu = cell2mat({Iu, 0*Iu, sparse(mu, mg)})';
            Ev = cell2mat({0*Iu, Iu, sparse(mu, mg)})';
            Egamma = cell2mat({sparse(mg, mu), sparse(mg, mu), Ig})';

            obj.Eu = Eu;
            obj.Ev = Ev;
            obj.Egamma = Egamma;

            %---- Forcing function F ----
            if ~isempty(F)
                if iscell(F.u)
                    Fu = @(t) multiblock.evalOn(g, F.u, t);
                    Fg = @(t) multiblock.evalOn(g, F.gamma, t);
                else
                    Fu = @(t) grid.evalOn(g, @(x,y)F.u(t,x,y));
                    Fg = @(t) grid.evalOn(g, @(x,y)F.gamma(t,x,y));
                end
                zeroVec = zeros(dim*g.N(), 1);
                Fv = @(t) zeroVec;

                Ft = @(t) [Fv(t); RHO_u\Fu(t); Fg(t)];
            else
                Ft = [];
            end
            %---------------------------

            % -- Penalty data -----
            if ~isempty(penalty)

                % Extract penalty parts
                penalty_u = @(t) (obj.Eu_small)'*penalty(t);
                penalty_gamma = @(t) (obj.Egamma_small)'*penalty(t);

                % Put penalty back together
                penalty = @(t) Ev*penalty_u(t) + Egamma*penalty_gamma(t);
            end
            %----------------------------

            % Create data function S(t)
            data_funcs = {Ft, penalty};
            S = elastic.addFunctionHandles(data_funcs);
            obj.S_cont = S;
            obj.S_discr = [];

            obj.S_body = Ft;
            obj.S_boundary = penalty;

            % Misc
            obj.M = M;
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

                    switch obj.eta.type
                    case 'linear'
                        ts = time.ExplicitRungeKuttaDiscreteData(...
                                            obj.M, obj.S_cont, obj.S_discr, k, t, obj.v0, 4);
                    case 'nonlinear'
                        if ~isempty(obj.S_discr)
                            error('Discrete forcing not implemented');
                        end

                        mms = struct;
                        if ~isempty(obj.S_body)
                            mms.bodyForcing = obj.S_body;
                        end
                        if ~isempty(obj.S_boundary)
                            mms.boundaryForcing = obj.S_boundary;
                        end
                        if isempty(obj.S_body) && isempty(obj.S_boundary)
                            mms = [];
                        end

                        ops = struct;
                        ops.Eu = obj.Eu;
                        ops.Ev = obj.Ev;
                        ops.Egamma = obj.Egamma;

                        ops.Duu = obj.Duu;
                        ops.Dug = obj.Dug;
                        ops.Dgu = obj.Dgu;
                        ops.Dgg = obj.Dgg;

                        ops.sigma_u = obj.sigma_u;
                        ops.sigma_gamma = obj.sigma_gamma;

                        F = @(U, t) elastic.visco.rhs.fullyDynamic(t, U, ops, obj.eta.nonlinFun, mms);
                        ts = time.Rungekutta4proper(F, k, 0, obj.v0);
                    end
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

            v0 = (obj.Eu)'*obj.v0;
            v1 = v0(1:2:end-1);
            v2 = v0(2:2:end);

            % v0 = (obj.Egamma)'*obj.v0;
            % v1 = v0(1:4:end-3);
            % v2 = v0(2:4:end-2);

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

            function update_fun(r,Eu,Egamma,xlims1,ylims1,xlims2,ylims2)
                t = r.t;
                v = Eu'*r.v;
                v1 = v(1:2:end-1);
                v2 = v(2:2:end);

                % v = Egamma'*r.v;
                % v1 = v(1:4:end-3);
                % v2 = v(2:4:end-2);
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
            update = @(r)update_fun(r,obj.Eu,obj.Egamma,xlims1,ylims1,xlims2,ylims2);
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