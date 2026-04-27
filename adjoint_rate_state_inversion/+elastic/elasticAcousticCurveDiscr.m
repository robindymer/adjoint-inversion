classdef elasticAcousticCurveDiscr < noname.Discretization
    properties
        name         = 'Elastic-acoustic, Curvilinear'
        description  = 'Solves for elastic displacement (u) and acoustic momentum potential (phi)'
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
        acousticDiscr

        g            %Elastic and acoustic grids
        domain       %Elastic and acoustic domains
        parameters   %Material parameters

        RHO_acoustic %Acoustic density
    end

    methods


        %       g           --      struct with fields 'acoustic' and 'elastic'
        %       order       --      order of accuracy
        %       parameters  --      struct with fields 'acoustic' and 'elastic'
        %                           that contain material parameters
        %       bc          --      struct with fields 'acoustic' and 'elastic',
        %                           each of which is a cell array of bc structs
        %       F           --      struct with fields 'acoustic' and 'elastic'
        %                           Forcing functions
        %       pointSources--      struct with fields 'acoustic' and 'elastic'
        %       intfForcing --      For MMS only. Struct with forcing functions at elastic-acoustic interface.
        function obj = elasticAcousticCurveDiscr(domain, order, parameters, bc, F, pointSources, intfForcing, optFlag)

            % ---- Default values --------------
            default_arg('order',4);
            default_arg('optFlag',[]);

            % Parameters
            parameters_default.acoustic.c = @(x,y) 0*x + 1;
            parameters_default.acoustic.rho = @(x,y) 0*x + 1;

            parameters_default.elastic.C = [];
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

            % Interface forcing (for MMS)
            default_arg('intfForcing', []);
            % -------------------------------

            dim = 2;
            obj.dim = dim;
            g = domain.g;

            % Elastic discretization
            C = parameters.elastic.C;
            rho = parameters.elastic.rho;

            elasticDiscr = elastic.elasticAnisotropicDiscrCurveMultiblock(g.elastic, order, C, rho, F.elastic,...
                                                                        bc.elastic, [], optFlag);
            diffOp.elastic = elasticDiscr.diffOp;
            D_elastic = elasticDiscr.D;

            % Acoustic parameters
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
            K = rho_acoustic .* c.^2;
            a = K;
            b = rho_acoustic.^(-1);

            % Split a and b to generate DiffOp
            a_cell = g.acoustic.splitFunc(a);
            b_cell = g.acoustic.splitFunc(b);
            acousticDiscr = elastic.acousticDiscrCurve(g.acoustic, order, a_cell, b_cell, F.acoustic, bc.acoustic, pointSources.acoustic);
            diffOp.acoustic = acousticDiscr.diffOp;
            D_acoustic = acousticDiscr.D;

            % Energy coefficients
            A_acoustic = spdiag(1./K);
            A_elastic = elasticDiscr.RHO_kron;

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
            %  -- u_t dot n = rho_acoustic^{-1} dphi/dn
            %  -- tau_i = n_i phi_t
            %
            elasticBoundary = domain.elastic.interfaceGroup;
            e1_elastic = diffOp.elastic.getBoundaryOperator('e1', elasticBoundary);
            e2_elastic = diffOp.elastic.getBoundaryOperator('e2', elasticBoundary);
            en_elastic = diffOp.elastic.getBoundaryOperator('en', elasticBoundary);

            acousticBoundary = domain.acoustic.interfaceGroup;
            e_acoustic = diffOp.acoustic.getBoundaryOperator('e', acousticBoundary);

            % -- Elastic side: Traction BC with acoustic data ---

            % Generate n_i phi
            acousticToBoundary_1 =  e1_elastic'* en_elastic * e_acoustic';
            acousticToBoundary_2 =  e2_elastic'* en_elastic * e_acoustic';

            % Traction comp 1
            type = {1, 't'};
            [closure, boundaryToElastic_1] = diffOp.elastic.boundary_condition(elasticBoundary, type);
            D_elastic = D_elastic + closure;

            % Tration comp 2
            type = {2, 't'};
            [closure, boundaryToElastic_2] = diffOp.elastic.boundary_condition(elasticBoundary, type);
            D_elastic = D_elastic + closure;

            % Coupling terms
            acousticToElastic_B = boundaryToElastic_1 * acousticToBoundary_1 ...
                                + boundaryToElastic_2 * acousticToBoundary_2;;

            % -- Acoustic side --
            [closure, boundaryToAcoustic] = diffOp.acoustic.boundary_condition(acousticBoundary, 'n');
            D_acoustic = D_acoustic + closure;

            % Coupling
            elasticToAcoustic_B = -1*RHO_acoustic*boundaryToAcoustic * en_elastic';
            %--------------------------------------------

            %---- Interface forcing ---------
            if ~isempty(intfForcing)
                X = e_acoustic'*g.acoustic.points();
                x = X(:,1);
                y = X(:,2);
                intfEl1 = @(t) boundaryToElastic_1 * intfForcing.tau1_diff(t,x,y);
                intfEl2 = @(t) boundaryToElastic_2 * intfForcing.tau2_diff(t,x,y);
                intfEl = @(t) intfEl1(t) + intfEl2(t);

                intfAc = @(t) -boundaryToAcoustic * intfForcing.n_der_diff(t,x,y);
            else
                intfAc = [];
                intfEl = [];
            end
            % -------------------------------

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
            % S_ac_cell = {acousticDiscr.S_cont, S_cont_ac, intfAc};
            S_ac_cell = {acousticDiscr.S_cont, intfAc};
            S_ac = elastic.addFunctionHandles(S_ac_cell);

            S_el_cell = {elasticDiscr.S_cont, intfEl};
            S_el = elastic.addFunctionHandles(S_el_cell);

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

            % if ~isempty(elasticDiscr.S_discr) && isempty(S_discr_ac)
            S_discr_el = elasticDiscr.S_discr;
            S_discr_ac = acousticDiscr.S_discr;
            if ~isempty(S_discr_el) && isempty(S_discr_ac)
                S_discr_ac = sparse(size(D_acoustic, 1), size(S_discr_el, 2));
                S_discr = [S_discr_el; S_discr_ac];
            elseif isempty(S_discr_el) && ~isempty(S_discr_ac)
                S_discr_el = sparse(size(D_elastic, 1), size(S_discr_ac, 2));
                S_discr = [S_discr_el; S_discr_ac];
            elseif ~isempty(S_discr_el) && ~isempty(S_discr_ac)
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
            obj.acousticDiscr = acousticDiscr;
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
            error('Not implemented for Curvilinear')

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
        function deltaFunctions = generateAcousticDeltaFunctions(obj, sourceCoord, blockIds, nMoment, nSmooth)
            default_arg('nMoment', obj.order);
            default_arg('nSmooth', 0);

            nSources = numel(sourceCoord);
            deltaFunctions = cell(nSources, 1);

            g = obj.g.acoustic;

            for s = 1:nSources
                x_s = sourceCoord{s};

                % Check which grid block source is in
                blockId = blockIds(s);

                % Delta function corresponding to local grid
                deltaFunLocal = elastic.diracDiscrCurve(x_s{s}, g.acoustic.grids{blockId}, order, 0);

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
            default_arg('cfl',0.5);

            % Get elastic time-step
            k_el = obj.elasticDiscr.getTimestep('rk4', []);

            % Get acoustic time-step
            k_ac = obj.acousticDiscr.getTimestep('rk4', []);

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
                % caxis(a, [-cmax, cmax])
                caxis(a, [-1, 1]);

            end
            update = @(r)update_fun(r,E,obj.E_elastic,obj.E_acoustic,xlims_el,ylims_el,xlims_ac,ylims_ac);
        end

        % Compare two functions u and v in the discrete l2 norm.
        function e = compareSolutions(obj, u, v)
            evec = u - v;
            e = sqrt(evec'*obj.H*evec);
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