classdef elasticAcousticDiscr < noname.Discretization
    properties
        name         = 'Elastic-acoustic, Cartesian'
        description  = 'Solves for elastic displacement (u) and acoustic velocity potential (phi)'
        order        %Order of accuracy
        dim          %Number of spatial dimensions

        % A w_tt + B w_t + C w = G(t), w = [u; phi]
        A, B, C

        % w_tt + E w_t + D w = S(t), w = [u; phi]
        D, E

        S_cont       %Function handle, S(t) = S_cont(t) + S_discr
        S_discr      %Matrix of time-dependent data, one column per time-stepping stage.

        E_elastic    %E_elastic'*v picks out elastic unknowns
        E_acoustic

        D_elastic    %Not actually used. Separate discretization matrices without coupling terms.
        D_acoustic

        A_elastic    %Blocks of A
        A_acoustic

        H            %Quadratures, combined and separate
        H_elastic
        H_acoustic

        v0           %Initial data
        v0t          %Initial data for time derivative
        diffOp       %Various operators
        elasticDiscr

        g            %Elastic and acoustic grids
        domain       %Elastic and acoustic domains
        parameters   %Material parameters

        RHO_acoustic %Acoustic density
    end

    methods


        %       domain      --      struct with fields 'acoustic' and 'elastic'
        %                           that contain information about domain, see below.
        %       order       --      order of accuracy
        %       parameters  --      struct with fields 'acoustic' and 'elastic'
        %                           that contain material parameters
        %       bc          --      struct with fields 'acoustic' and 'elastic',
        %                           each of which is a cell array of bc structs
        %       F           --      struct with fields 'acoustic' and 'elastic'
        %                           Forcing functions
        %       pointSources--      struct with fields 'acoustic' and 'elastic'
        %
        function obj = elasticAcousticDiscr(domain, order, parameters, bc, F, pointSources, optFlag)

            % ---- Default values --------------
            default_arg('order',4);
            default_arg('optFlag',[]);

            % Domain
            domain_default.elastic.xlim = [0,1];
            domain_default.elastic.ylim = [0,-1];
            domain_default.elastic.m = [31, 31];
            domain_default.elastic.def = multiblock.domain.Rectangle([0,1], [0,-1]);
            domain_default.elastic.interfaceGroup = 'N';

            domain_default.acoustic.xlim = [0,1];
            domain_default.acoustic.ylim = [1,0];
            domain_default.acoustic.m = [31, 31];
            domain_default.acoustic.def = multiblock.domain.Rectangle([0,1], [1,0]);
            domain_default.acoustic.interfaceGroup = 'S';

            default_struct('domain', domain_default);

            % Parameters
            parameters_default.acoustic.c = @(x,y) 0*x + 1;
            parameters_default.acoustic.rho = @(x,y) 0*x + 1;

            parameters_default.elastic.lambda = [];
            parameters_default.elastic.mu = [];
            parameters_default.elastic.rho = [];

            default_struct('parameters', parameters_default);

            % Forcing
            F_default.elastic = [];
            F_default.acoustic = [];
            default_struct('F', F_default);

            % Point sources
            pointSources_default.elastic = [];
            pointSources_default.acoustic = [];
            default_struct('pointSources', pointSources_default);
            % -------------------------------

            dim = 2;
            obj.dim = dim;

            % Elastic discretization
            lambda = parameters.elastic.lambda;
            mu = parameters.elastic.mu;
            rho = parameters.elastic.rho;
            xlim = domain.elastic.xlim;
            ylim = domain.elastic.ylim;
            m = domain.elastic.m;

            elasticDiscr = elastic.elasticDiscrInterface(m, order, xlim, ylim, lambda, mu, rho, F.elastic,...
                                             bc.elastic, pointSources.elastic, [], optFlag);
            diffOp.elastic = elasticDiscr.diffOp;
            g.elastic = elasticDiscr.grid;
            D_elastic = elasticDiscr.D;

            % Acoustic parameters
            m = domain.acoustic.m;
            g.acoustic = domain.acoustic.def.getGrid(m);
            c = parameters.acoustic.c;
            rho_acoustic = parameters.acoustic.rho;

            if isa(rho_acoustic, 'function_handle')
                rho_acoustic = grid.evalOn(g.acoustic, rho_acoustic);
            end
            if isa(c, 'function_handle')
                c = grid.evalOn(g.acoustic, c);
            end
            RHO_acoustic = spdiag(rho_acoustic);

            % Acoustic diffOp
            a = c.^2 ./ rho_acoustic;
            b = rho_acoustic;

            % Split a and b to generate DiffOp
            a_cell = g.acoustic.splitFunc(a);
            b_cell = g.acoustic.splitFunc(b);
            doPar = cell(g.acoustic.nBlocks, 1);
            for i = 1:g.acoustic.nBlocks
                doPar{i} = {a_cell{i}, b_cell{i}};
            end
            diffOp.acoustic = multiblock.DiffOp(@scheme.LaplaceCurvilinearNewCorner, g.acoustic, order, doPar);
            D_acoustic = diffOp.acoustic.D;
            [closure, penalty_acoustic] = scheme.bcSetup(diffOp.acoustic, bc.acoustic);
            D_acoustic = D_acoustic + closure;

            % Energy coefficients
            A_acoustic = spdiag( rho_acoustic./c.^2 );
            A_elastic = elasticDiscr.RHO_kron;

            %----- Acoustic point sources ----------------
            discrAcousticSource = false;
            contAcousticSource = false;
            if isfield(pointSources, 'acoustic') && ~isempty(pointSources.acoustic)

                nBlocks = g.acoustic.nBlocks();
                x_s = pointSources.acoustic.x;
                n_s = length(x_s);

                % Initialize
                source_func_cont = @(t) 0*t;
                source_func_discr = 0;

                for s = 1:n_s

                    % Check which grid block source is in
                    blockId = [];
                    for i = 1:nBlocks-1
                        xlim = g.acoustic.grids{i}.lim{1};
                        ylim = g.acoustic.grids{i}.lim{2};
                        if ( xlim{1} <= x_s{s}(1) && x_s{s}(1) < xlim{2} && ...
                             ylim{1} <= x_s{s}(2) && x_s{s}(2) < ylim{2} )
                            blockId = i;
                        end
                    end
                    if isempty(blockId)
                        blockId = nBlocks;
                    end

                    % Operators local to block number blockId
                    H_local = diffOp.acoustic.diffOps{blockId}.H_1D;
                    xlim = g.acoustic.grids{blockId}.lim{1};
                    ylim = g.acoustic.grids{blockId}.lim{2};
                    Lx = xlim{2}-xlim{1};
                    Ly = ylim{2}-ylim{1};
                    H_local{1} = H_local{1}*Lx;
                    H_local{2} = H_local{2}*Ly;

                    % Delta function corresponding to local grid
                    delta_fun_local = diracDiscr(x_s{s}, g.acoustic.grids{blockId}.x, order, 0, H_local);

                    % Extend delta fun to global acoustic grid.
                    delta_fun = g.acoustic.expandFunc(delta_fun_local, blockId);

                    % Accumulate different sources.
                    if isa(pointSources.acoustic.g{s}, 'function_handle')
                        source_func_cont = @(t) source_func_cont(t) + ...
                                           A_acoustic \ (pointSources.acoustic.g{s}(t)*delta_fun);
                        contAcousticSource = true;
                    else
                        source_vec = pointSources.acoustic.g{s};
                        if iscolumn(source_vec)
                            source_vec = transpose(source_vec);
                        end
                        source_func_discr = source_func_discr + ...
                                           A_acoustic \ kron(source_vec, delta_fun);

                        discrAcousticSource = true;
                    end
                end
            end

            if discrAcousticSource
                S_discr_ac = source_func_discr;
            else
                S_discr_ac = [];
            end

            if contAcousticSource
                S_cont_ac = source_func_cont;
            else
                S_cont_ac = [];
            end
            % -----------------------------------

            % 'Adjoint' form of equations:
            % A w_tt + B w_t + C w = S, w = [u; phi].
            %
            % A -- Diagonal, positive 'energy coefficients'
            % B -- Skewsymmetric coupling terms (i.e. off-diagonal blocks only)
            % C -- Symmetric and positive semidefinite
            % S -- Forcing and penalties

            %----------- Interface ------------------------
            %
            %  Interface conditions:
            %  -- u_t dot n = dphi/dn
            %  -- tau = rho_acoustic*phi*n
            %
            elIG = domain.elastic.interfaceGroup;
            elasticBoundary = domain.elastic.def.boundaryGroups.(elIG);
            e_elastic = diffOp.elastic.getBoundaryOperator('e', elasticBoundary);
            E = elasticDiscr.Ecomp;
            e_elastic_small = transpose(E{1})*e_elastic;
            e_elastic_small = e_elastic_small(:, 1:obj.dim:end);
            e_elastic_normal = (e_elastic_small'* E{2}')';

            acIG = domain.acoustic.interfaceGroup;
            acousticBoundary = domain.acoustic.def.boundaryGroups.(acIG);
            e_acoustic = diffOp.acoustic.getBoundaryOperator('e', acousticBoundary);

            % -- Elastic side ---
            % Homogeneous traction BC for tangential component
            type = {'tangential', 't'};
            closure = diffOp.elastic.boundary_condition(elasticBoundary, type);
            D_elastic = D_elastic + closure;

            % Traction condition with acoustic data for normal component
            type = {'normal', 't'};
            [closure, boundaryToElastic] = diffOp.elastic.boundary_condition(elasticBoundary, type);
            D_elastic = D_elastic + closure;

            acousticToElastic_B = boundaryToElastic * e_acoustic' * RHO_acoustic;

            % -- Acoustic side --
            [closure, boundaryToAcoustic] = diffOp.acoustic.boundary_condition(acousticBoundary, 'n');
            D_acoustic = D_acoustic + closure;
            elasticToAcoustic_B = -1*boundaryToAcoustic * e_elastic_normal';
            %--------------------------------------------

            % ----- Build total matrices ----------------
            I_elastic = speye(size(D_elastic));
            I_acoustic = speye(size(D_acoustic));

            A = {A_elastic, [];...
                    [],     A_acoustic};

            B = {[]                 ,  -A_elastic*acousticToElastic_B;...
                 -A_acoustic*elasticToAcoustic_B,       []             };

            C = {-A_elastic*D_elastic,        [];...
                    [],            -A_acoustic*D_acoustic};

            A = blockmatrix.toMatrix(A);
            B = blockmatrix.toMatrix(B);
            C = blockmatrix.toMatrix(C);

            D = -A\C;
            E = -A\B;

            E_elastic = {I_elastic; sparse( length(I_acoustic), length(I_elastic) ) };
            E_acoustic = {sparse( length(I_elastic), length(I_acoustic) ); I_acoustic };

            E_elastic = blockmatrix.toMatrix(E_elastic);
            E_acoustic = blockmatrix.toMatrix(E_acoustic);
            % ----------------------------------------------

            % ----- Quadratures etc ------------------------
            H_elastic = elasticDiscr.H;
            H_acoustic = diffOp.acoustic.H;
            H = E_elastic*H_elastic*E_elastic' + E_acoustic*H_acoustic*E_acoustic';
            % ----------------------------------------------

            % ----- Point sources and forcing --------------
            S_ac_cell = {penalty_acoustic, S_cont_ac};
            S_ac = elastic.addFunctionHandles(S_ac_cell);
            if ~isempty(F.acoustic)
                S_ac = @(t) S_ac(t) + grid.evalOn(g.acoustic, @(x,y)F.acoustic(t,x,y));
            end

            S_el = elasticDiscr.S_cont;
            if ~isempty(S_ac) && ~isempty(S_el)
                S = @(t) [S_el(t); S_ac(t)];
            end
            if isempty(S_ac) && isempty(S_el)
                S = [];
            end
            if ~isempty(S_el) && isempty(S_ac)
                zero_ac = sparse(length(D_acoustic), 1);
                S = @(t) [S_el(t); zero_ac];
            end
            if isempty(S_el) && ~isempty(S_ac)
                zero_el = sparse(length(D_elastic), 1);
                S = @(t) [zero_el; S_ac(t)];
            end

            if ~isempty(elasticDiscr.S_discr) && isempty(S_discr_ac)
                S_discr_el = elasticDiscr.S_discr;
                S_discr_ac = sparse(size(D_acoustic, 1), size(S_discr_el, 2));
                S_discr = [S_discr_el; S_discr_ac];
            elseif isempty(elasticDiscr.S_discr) && ~isempty(S_discr_ac)
                S_discr_el = sparse(size(D_elastic, 1), size(S_discr_ac, 2));
                S_discr = [S_discr_el; S_discr_ac];
            elseif ~isempty(elasticDiscr.S_discr) && ~isempty(S_discr_ac)
                S_discr_el = elasticDiscr.S_discr;
                S_discr = [S_discr_el; S_discr_ac];
            else
                S_discr = [];
            end
            % ----------------------------------------------


            % ---- Set properties -----
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.D = D;
            obj.E = E;

            obj.A_acoustic = A_acoustic;
            obj.A_elastic = A_elastic;

            obj.D_acoustic = D_acoustic;
            obj.D_elastic = D_elastic;

            obj.E_acoustic = E_acoustic;
            obj.E_elastic = E_elastic;

            obj.H = H;
            obj.H_acoustic = H_acoustic;
            obj.H_elastic = H_elastic;

            obj.S_cont = S;
            obj.S_discr = S_discr;

            obj.v0 = [];
            obj.v0t = [];

            obj.RHO_acoustic = RHO_acoustic;

            obj.diffOp = diffOp;
            obj.elasticDiscr = elasticDiscr;
            obj.g = g;
            obj.domain = domain;
            obj.parameters = parameters;
            obj.order = order;

        end
        % Prints some info about the discretisation
        function printInfo(obj)
            fprintf('Name: %s\n',obj.name);
            fprintf('Size: %d\n',obj.size());
        end

        % Return the number of DOF
        function n = size(obj)
            n = length(obj.C);
        end

        % Generates delta functions on elastic side
        function deltaFunctions = generateElasticDeltaFunctions(obj, sourceCoord, nMoment, nSmooth)
            default_arg('nMoment', []);
            default_arg('nSmooth', []);

            nSources = numel(sourceCoord);
            d = obj.dim;

            % Use elasticDiscr to generate
            deltaFunctions = obj.elasticDiscr.generateDeltaFunctions(sourceCoord, nMoment, nSmooth);

            % Loop through sources and extend to elastic-acoustic grid
            for s = 1:nSources
                for d = 1:obj.dim
                    deltaFunctions{s,d} = obj.E_elastic*deltaFunctions{s,d};
                end
            end
        end

        % Generates delta functions on acoustic side
        % Source coordinates should be a cell array of coordinate vectors.
        % Returns a (nSources x dim) cell array of delta functions.
        function deltaFunctions = generateAcousticDeltaFunctions(obj, sourceCoord, nMoment, nSmooth)
            default_arg('nMoment', obj.order);
            default_arg('nSmooth', 0);

            nSources = numel(sourceCoord);
            deltaFunctions = cell(nSources, 1);

            g = obj.g.acoustic;
            nBlocks = g.nBlocks();

            for s = 1:nSources
                x_s = sourceCoord{s};

                % Check which grid block source is in
                blockId = [];
                for i = 1:nBlocks-1
                    xlim = g.grids{i}.lim{1};
                    ylim = g.grids{i}.lim{2};
                    if ( xlim{1} <= x_s(1) && x_s(1) < xlim{2} && ...
                         ylim{1} <= x_s(2) && x_s(2) < ylim{2} )
                    blockId = i;
                    end
                end
                if isempty(blockId)
                    blockId = nBlocks;
                end

                % Operators local to block number blockId
                H_local = obj.diffOp.acoustic.diffOps{blockId}.H_1D;
                xlim = g.grids{blockId}.lim{1};
                ylim = g.grids{blockId}.lim{2};
                Lx = xlim{2}-xlim{1};
                Ly = ylim{2}-ylim{1};
                H_local{1} = H_local{1}*Lx;
                H_local{2} = H_local{2}*Ly;

                % Delta function corresponding to local grid
                deltaFunLocal = diracDiscr(x_s, g.grids{blockId}.x, nMoment, nSmooth, H_local);

                % Extend delta fun to global acoustic grid.
                deltaFunctions{s} = g.expandFunc(deltaFunLocal, blockId);

                % Extend delta fun to acoustic-elastic grid
                deltaFunctions{s} = obj.E_acoustic*deltaFunctions{s};

            end
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
                error('Specify CFL')
            end

            % Get elastic time-step
            discr = obj.elasticDiscr;
            nBlocks = length(discr.grid.grids);
            g = obj.g;
            k = zeros(nBlocks, 1);

            for i = 1:nBlocks
                rho = diag(discr.diffOp.diffOps{i}.RHO);
                mu = diag(discr.diffOp.diffOps{i}.MU);
                lambda = diag(discr.diffOp.diffOps{i}.LAMBDA);

                v_s = sqrt(max( mu./rho ));
                v_p = sqrt(max( (lambda + 2*mu)./rho ));
                v = max(v_p,v_s);

                h = min(discr.grid.grids{i}.h);
                k(i) = cfl/v*h;
            end
            k_el = min(k);

            % Get acoustic time-step
            nBlocks = g.acoustic.nBlocks;
            k = zeros(nBlocks, 1);
            for i = 1:nBlocks
                a = diag(obj.diffOp.acoustic.diffOps{i}.a);
                b = diag(obj.diffOp.acoustic.diffOps{i}.b);
                c = sqrt(max(a.*b));

                h = min(obj.g.acoustic.grids{i}.h);
                k(i) = cfl/c*h;
            end
            k_ac = min(k);

            k = min(k_el, k_ac);
        end

        function r = getTimeSnapshot(obj, ts)
            if ts == 0
                r.t = 0;
                r.v = obj.v0;
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

            g = obj.g;
            v0 = obj.v0;

            v0_el = obj.E_elastic'*v0;
            v0_ac = obj.E_acoustic'*v0;

            E = obj.elasticDiscr.Ecomp;

            figure_handle = figure();
            % h_ac = subplot(2,1,1);
            Sur_ac = multiblock.Surface(obj.g.acoustic, v0_ac);
            hold on
            ylabel('y')
            shading interp
            colorbar
            xlims_ac = xlim;
            ylims_ac = ylim;
            % axis equal
            a = gca;

            % h_el = subplot(2,1,2);
            Sur_el = multiblock.Surface(obj.g.elastic, E{1}'*v0_el);
            xlabel('x')
            ylabel('y')
            shading interp
            colorbar
            xlims_el = xlim;
            ylims_el = ylim;
            % axis equal

            function update_fun(r,E,E_el,E_ac,xlims_el,ylims_el,xlims_ac,ylims_ac)
                t = r.t;
                v_el = E_el'*r.v;
                v_ac = E_ac'*r.v;
                if ishandle(a)
                    title(a,sprintf('T = %.3f',t))
                end
                Sur_el.ZData = E{1}'*v_el;
                Sur_el.CData = E{1}'*v_el;
                Sur_ac.ZData = v_ac;
                Sur_ac.CData = v_ac;
                cmax = max(abs([Sur_el.CData; Sur_ac.CData]));
                caxis(a, [-cmax, cmax])

            end
            update = @(r)update_fun(r,E,obj.E_elastic,obj.E_acoustic,xlims_el,ylims_el,xlims_ac,ylims_ac);
        end

        % Compare two functions u and v in the discrete l2 norm.
        function e = compareSolutions(obj, u, v)
            evec = u - v;
            e = 0;

            % Elastic error
            evec = (obj.E_elastic)'*(u - v);
            g = obj.g.elastic;
            nBlocks = length(g.grids);
            ecell = g.splitFunc(evec);
            for i = 1:nBlocks
                subgrid = g.grids{i};
                h = subgrid.h;
                e = e + sqrt(prod(h))*norm(ecell{i});
            end

            % Acoustic error
            evec = (obj.E_acoustic)'*(u - v);
            g = obj.g.acoustic;
            nBlocks = length(g.grids);
            ecell = g.splitFunc(evec);
            for i = 1:nBlocks
                subgrid = g.grids{i};
                h = subgrid.h;
                e = e + sqrt(prod(h))*norm(ecell{i});
            end
        end

        % Compare the grid function u to the analytical function g in the discrete l2 norm.
        function e = compareSolutionsAnalytical(obj, w, u_exact, phi_exact, t)
            v_el = grid.evalOn(obj.g.elastic, @(x,y)u_exact(t,x,y));
            v_ac = grid.evalOn(obj.g.acoustic, @(x,y)phi_exact(t,x,y));
            v = [v_el; v_ac];
            e = obj.compareSolutions(w, v);
        end

    end

    methods(Static)

    end
end