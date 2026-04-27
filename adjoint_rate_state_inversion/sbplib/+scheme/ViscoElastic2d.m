classdef ViscoElastic2d < scheme.Scheme

% Discretizes the visco-elastic wave equation in curvilinear coordinates.
% Assumes fully compatible operators.

    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid
        dim

        order % Order of accuracy for the approximation

        % Diagonal matrices for variable coefficients
        % J, Ji
        RHO % Density
        C   % Elastic stiffness tensor
        ETA % Effective viscosity, used in strain rate eq

        D % Total operator
        Delastic        % Elastic operator (momentum balance)
        Dviscous        % Acts on viscous strains in momentum balance
        DstrainRate     % Acts on u and gamma, returns strain rate gamma_t

        D1, D1Tilde % Physical derivatives
        sigma % Cell matrix of physical stress operators

        % Inner products
        H

        % Boundary inner products (for scalar field)
        H_w, H_e, H_s, H_n

        % Restriction operators
        Eu, Egamma  % Pick out all components of u/gamma
        eU, eGamma  % Pick out one specific component

        % Bundary restriction ops
        e_scalar_w, e_scalar_e, e_scalar_s, e_scalar_n

        n_w, n_e, n_s, n_n % Physical normals
        tangent_w, tangent_e, tangent_s, tangent_n % Physical tangents

        tau1_w, tau1_e, tau1_s, tau1_n  % Return scalar field
        tau2_w, tau2_e, tau2_s, tau2_n  % Return scalar field
        tau_n_w, tau_n_e, tau_n_s, tau_n_n % Return scalar field
        tau_t_w, tau_t_e, tau_t_s, tau_t_n % Return scalar field

        elasticObj

    end

    methods

        % The coefficients can either be function handles or grid functions
        function obj = ViscoElastic2d(g, order, rho, C, eta)
            default_arg('rho', @(x,y) 0*x+1);
            default_arg('eta', @(x,y) 0*x+1);
            dim = 2;

            C_default = cell(dim,dim,dim,dim);
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            C_default{i,j,k,l} = @(x,y) 0*x ;
                        end
                    end
                end
            end
            default_arg('C', C_default);

            assert(isa(g, 'grid.Curvilinear'));

            if isa(rho, 'function_handle')
                rho = grid.evalOn(g, rho);
            end

            if isa(eta, 'function_handle')
                eta = grid.evalOn(g, eta);
            end

            C_mat = cell(dim,dim,dim,dim);
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            if isa(C{i,j,k,l}, 'function_handle')
                                C{i,j,k,l} = grid.evalOn(g, C{i,j,k,l});
                            end
                            C_mat{i,j,k,l} = spdiag(C{i,j,k,l});
                        end
                    end
                end
            end
            obj.C = C_mat;

            elasticObj = scheme.Elastic2dCurvilinearAnisotropic(g, order, rho, C);

            % Construct a pair of first derivatives
            K = elasticObj.K;
            for i = 1:dim
                for j = 1:dim
                    K{i,j} = spdiag(K{i,j});
                end
            end
            J = elasticObj.J;
            Ji = elasticObj.Ji;
            D_ref = elasticObj.refObj.D1;
            D1 = cell(dim, 1);
            D1Tilde = cell(dim, 1);
            for i = 1:dim
                D1{i} = 0*D_ref{i};
                D1Tilde{i} = 0*D_ref{i};
                for j = 1:dim
                    D1{i} = D1{i} + K{i,j}*D_ref{j};
                    D1Tilde{i} = D1Tilde{i} + Ji*D_ref{j}*J*K{i,j};
                end
            end
            obj.D1 = D1;
            obj.D1Tilde = D1Tilde;

            eU = elasticObj.E;

            % Storage order for gamma: 11-12-21-22.
            I = speye(g.N(), g.N());
            eGamma = cell(dim, dim);
            e = cell(dim, dim);
            e{1,1} = [1;0;0;0];
            e{1,2} = [0;1;0;0];
            e{2,1} = [0;0;1;0];
            e{2,2} = [0;0;0;1];
            for i = 1:dim
                for j = 1:dim
                    eGamma{i,j} = kron(I, e{i,j});
                end
            end

            % Store u first, then gamma
            mU = dim*g.N();
            mGamma = dim^2*g.N();
            Iu = speye(mU, mU);
            Igamma = speye(mGamma, mGamma);

            Eu = cell2mat({Iu, sparse(mU, mGamma)})';
            Egamma = cell2mat({sparse(mGamma, mU), Igamma})';

            for i = 1:dim
                eU{i} = Eu*eU{i};
            end
            for i = 1:dim
                for j = 1:dim
                    eGamma{i,j} = Egamma*eGamma{i,j};
                end
            end

            obj.eGamma = eGamma;
            obj.eU = eU;
            obj.Egamma = Egamma;
            obj.Eu = Eu;

            % Build stress operator
            sigma = cell(dim, dim);
            C = obj.C;
            for i = 1:dim
                for j = 1:dim
                    sigma{i,j} = spalloc(g.N(), (dim^2 + dim)*g.N(), order^2*g.N());
                    for k = 1:dim
                        for l = 1:dim
                            sigma{i,j} = sigma{i,j} + C{i,j,k,l}*(D1{k}*eU{l}' - eGamma{k,l}');
                        end
                    end
                end
            end

            % Elastic operator
            Delastic = Eu*elasticObj.D*Eu';

            % Add viscous strains to momentum balance
            RHOi = spdiag(1./rho);
            Dviscous = spalloc((dim^2 + dim)*g.N(), (dim^2 + dim)*g.N(), order^2*(dim^2 + dim)*g.N());
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            Dviscous = Dviscous - eU{j}*RHOi*D1Tilde{i}*C{i,j,k,l}*eGamma{k,l}';
                        end
                    end
                end
            end

            ETA = spdiag(eta);
            DstrainRate = 0*Delastic;
            for i = 1:dim
                for j = 1:dim
                    DstrainRate = DstrainRate + eGamma{i,j}*(ETA\sigma{i,j});
                end
            end

            obj.D = Delastic + Dviscous + DstrainRate;
            obj.Delastic = Delastic;
            obj.Dviscous = Dviscous;
            obj.DstrainRate = DstrainRate;
            obj.sigma = sigma;

            %---- Set remaining object properties ------
            obj.RHO = elasticObj.RHO;
            obj.ETA = ETA;
            obj.H = elasticObj.H;

            obj.n_w = elasticObj.n_w;
            obj.n_e = elasticObj.n_e;
            obj.n_s = elasticObj.n_s;
            obj.n_n = elasticObj.n_n;

            obj.tangent_w = elasticObj.tangent_w;
            obj.tangent_e = elasticObj.tangent_e;
            obj.tangent_s = elasticObj.tangent_s;
            obj.tangent_n = elasticObj.tangent_n;

            obj.H_w = elasticObj.H_w;
            obj.H_e = elasticObj.H_e;
            obj.H_s = elasticObj.H_s;
            obj.H_n = elasticObj.H_n;

            obj.e_scalar_w = elasticObj.e_scalar_w;
            obj.e_scalar_e = elasticObj.e_scalar_e;
            obj.e_scalar_s = elasticObj.e_scalar_s;
            obj.e_scalar_n = elasticObj.e_scalar_n;

            % -- Create new traction operators including viscous strain contribution --
            tau1 = struct;
            tau2 = struct;
            tau_n = struct;
            tau_t = struct;
            boundaries = {'w', 'e', 's', 'n'};
            for bNumber = 1:numel(boundaries)
                b = boundaries{bNumber};

                n = elasticObj.getNormal(b);
                t = elasticObj.getTangent(b);
                e = elasticObj.getBoundaryOperatorForScalarField('e', b);
                tau1.(b) = (elasticObj.getBoundaryOperator('tau1', b)'*Eu')';
                tau2.(b) = (elasticObj.getBoundaryOperator('tau2', b)'*Eu')';

                % Add viscous contributions
                for i = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            tau1.(b) = tau1.(b) - (n{i}*e'*C{i,1,k,l}*eGamma{k,l}')';
                            tau2.(b) = tau2.(b) - (n{i}*e'*C{i,2,k,l}*eGamma{k,l}')';
                        end
                    end
                end

                % Compute normal and tangential components
                tau_n.(b) = tau1.(b)*n{1} + tau2.(b)*n{2};
                tau_t.(b) = tau1.(b)*t{1} + tau2.(b)*t{2};
            end

            obj.tau1_w = tau1.w;
            obj.tau1_e = tau1.e;
            obj.tau1_s = tau1.s;
            obj.tau1_n = tau1.n;

            obj.tau2_w = tau2.w;
            obj.tau2_e = tau2.e;
            obj.tau2_s = tau2.s;
            obj.tau2_n = tau2.n;

            obj.tau_n_w = tau_n.w;
            obj.tau_n_e = tau_n.e;
            obj.tau_n_s = tau_n.s;
            obj.tau_n_n = tau_n.n;

            obj.tau_t_w = tau_t.w;
            obj.tau_t_e = tau_t.e;
            obj.tau_t_s = tau_t.s;
            obj.tau_t_n = tau_t.n;
            %----------------------------------------

            % Misc.
            obj.elasticObj = elasticObj;
            obj.m = elasticObj.m;
            obj.h = elasticObj.h;

            obj.order = order;
            obj.grid = g;
            obj.dim = dim;

        end


        % Closure functions return the operators applied to the own domain to close the boundary
        % Penalty functions return the operators to force the solution. In the case of an interface it returns the operator applied to the other doamin.
        %       boundary            is a string specifying the boundary e.g. 'l','r' or 'e','w','n','s'.
        %       bc                  is a cell array of component and bc type, e.g. {1, 'd'} for Dirichlet condition
        %                           on the first component. Can also be e.g.
        %                           {'normal', 'd'} or {'tangential', 't'} for conditions on
        %                           tangential/normal component.
        %       data                is a function returning the data that should be applied at the boundary.
        %       neighbour_scheme    is an instance of Scheme that should be interfaced to.
        %       neighbour_boundary  is a string specifying which boundary to interface to.

        % For displacement bc:
        % bc = {comp, 'd', dComps},
        % where
        % dComps = vector of components with displacement BC. Default: 1:dim.
        % In this way, we can specify one BC at a time even though the SATs depend on all BC.
        function [closure, penalty] = boundary_condition(obj, boundary, bc, tuning)
            default_arg('tuning', 1.0);
            assert( iscell(bc), 'The BC type must be a 2x1 or 3x1 cell array' );

            component = bc{1};
            type = bc{2};
            dim = obj.dim;

            n       = obj.getNormal(boundary);
            H_gamma = obj.getBoundaryQuadratureForScalarField(boundary);
            e       = obj.getBoundaryOperatorForScalarField('e', boundary);

            H       = obj.H;
            RHO     = obj.RHO;
            ETA     = obj.ETA;
            C       = obj.C;
            Eu      = obj.Eu;
            eU      = obj.eU;
            eGamma  = obj.eGamma;

            % Get elastic closure and penalty
            [closure, penalty] = obj.elasticObj.boundary_condition(boundary, bc, tuning);
            closure = Eu*closure*Eu';
            penalty = Eu*penalty;

            switch component
            case 't'
                dir = obj.getTangent(boundary);
                tau = obj.getBoundaryOperator('tau_t', boundary);
            case 'n'
                dir = obj.getNormal(boundary);
                tau = obj.getBoundaryOperator('tau_n', boundary);
            case 1
                dir = {1, 0};
                tau = obj.getBoundaryOperator('tau1', boundary);
            case 2
                dir = {0, 1};
                tau = obj.getBoundaryOperator('tau2', boundary);
            end

            switch type
            case {'F','f','Free','free','traction','Traction','t','T'}

                % Set elastic closure to zero
                closure = 0*closure;

                for m = 1:dim
                    closure = closure - eU{m}*( (RHO*H)\(e*dir{m}*H_gamma*tau') );
                end

            case {'D','d','dirichlet','Dirichlet','displacement','Displacement'}

                % Add penalty to strain rate eq
                for i = 1:dim
                    for j = 1:dim
                        for k = 1:dim
                            for l = 1:dim
                                for m = 1:dim
                                    closure = closure - eGamma{i,j}*( (H*ETA)\(C{i,j,k,l}*e*H_gamma*dir{l}*dir{m}*n{k}*e'*eU{m}') );
                                end
                                penalty = penalty + eGamma{i,j}*( (H*ETA)\(C{i,j,k,l}*e*H_gamma*dir{l}*n{k}) );
                            end
                        end
                    end
                end

            end

        end

        % type     Struct that specifies the interface coupling.
        %          Fields:
        %          -- tuning:           penalty strength, defaults to 1.0
        %          -- interpolation:    type of interpolation, default 'none'
        function [closure, penalty, forcingPenalties] = interface(obj,boundary,neighbour_scheme,neighbour_boundary,type)

            defaultType.tuning = 1.0;
            defaultType.interpolation = 'none';
            defaultType.type = 'standard';
            default_struct('type', defaultType);

            forcingPenalties = [];

            switch type.type
            case 'standard'
                [closure, penalty] = obj.interfaceStandard(boundary,neighbour_scheme,neighbour_boundary,type);
            case 'frictionalFault'
                [closure, penalty] = obj.interfaceFrictionalFault(boundary,neighbour_scheme,neighbour_boundary,type);
            case 'normalTangential'
                [closure, penalty, forcingPenalties] = obj.interfaceNormalTangential(boundary,neighbour_scheme,neighbour_boundary,type);
            end

        end

        function [closure, penalty] = interfaceStandard(obj,boundary,neighbour_scheme,neighbour_boundary,type)

            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            u = obj;
            v = neighbour_scheme;

            dim = obj.dim;

            n       = u.getNormal(boundary);
            H_gamma = u.getBoundaryQuadratureForScalarField(boundary);
            e       = u.getBoundaryOperatorForScalarField('e', boundary);

            ev      = v.getBoundaryOperatorForScalarField('e', neighbour_boundary);

            H       = u.H;
            RHO     = u.RHO;
            ETA     = u.ETA;
            C       = u.C;
            Eu      = u.Eu;
            eU      = u.eU;
            eGamma  = u.eGamma;

            CV       = v.C;
            Ev      = v.Eu;
            eV      = v.eU;
            eGammaV = v.eGamma;
            nV      = v.getNormal(neighbour_boundary);


            % Get elastic closure and penalty
            [closure, penalty] = obj.elasticObj.interface(boundary, v.elasticObj, neighbour_boundary, type);
            closure = Eu*closure*Eu';
            penalty = Eu*penalty*Ev';

            % Add viscous part of traction coupling
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            closure = closure + 1/2*eU{j}*( (RHO*H)\(C{i,j,k,l}*e*H_gamma*n{i}*e'*eGamma{k,l}') );
                            penalty = penalty + 1/2*eU{j}*( (RHO*H)\(e*H_gamma*nV{i}*(ev'*CV{i,j,k,l}*ev)*ev'*eGammaV{k,l}') );
                        end
                    end
                end
            end

            % Add penalty to strain rate eq
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            closure = closure - 1/2*eGamma{i,j}*( (H*ETA)\(C{i,j,k,l}*e*H_gamma*n{k}*e'*eU{l}') );
                            penalty = penalty + 1/2*eGamma{i,j}*( (H*ETA)\(C{i,j,k,l}*e*H_gamma*n{k}*ev'*eV{l}') );
                        end
                    end
                end
            end


        end

        function [closure, penalty] = interfaceFrictionalFault(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            tuning = type.tuning;

            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            u = obj;
            v = neighbour_scheme;

            dim = obj.dim;

            n       = u.getNormal(boundary);
            H_gamma = u.getBoundaryQuadratureForScalarField(boundary);
            e       = u.getBoundaryOperatorForScalarField('e', boundary);

            ev      = v.getBoundaryOperatorForScalarField('e', neighbour_boundary);

            H       = u.H;
            RHO     = u.RHO;
            ETA     = u.ETA;
            C       = u.C;
            Eu      = u.Eu;
            eU      = u.eU;
            eGamma  = u.eGamma;
            Egamma  = u.Egamma;

            CV       = v.C;
            Ev      = v.Eu;
            eV      = v.eU;
            eGammaV = v.eGamma;
            nV      = v.getNormal(neighbour_boundary);

            % Reduce stiffness tensors to boundary size
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            C{i,j,k,l} = e'*C{i,j,k,l}*e;
                            CV{i,j,k,l} = ev'*CV{i,j,k,l}*ev;
                        end
                    end
                end
            end

            % Get elastic closure and penalty
            [closure, penalty] = obj.elasticObj.interface(boundary, v.elasticObj, neighbour_boundary, type);
            closure = Eu*closure*Eu';
            penalty = Eu*penalty*Ev';

            % ---- Tangential tractions are imposed just like traction BC ------
            % We only need the viscous part
            closure_tangential = obj.boundary_condition(boundary, {'t', 't'});
            closure = closure + closure_tangential*Egamma*Egamma';


            % ------ Coupling of normal component -----------
            % Add viscous part of traction coupling
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            for m = 1:dim
                                closure = closure + 1/2*eU{m}*( (RHO*H)\(e*n{m}*H_gamma*n{j}*n{i}*C{i,j,k,l}*e'*eGamma{k,l}') );
                                penalty = penalty - 1/2*eU{m}*( (RHO*H)\(e*n{m}*H_gamma*nV{j}*nV{i}*CV{i,j,k,l}*ev'*eGammaV{k,l}') );
                            end
                        end
                    end
                end
            end

            % Add penalty to strain rate eq
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            for m = 1:dim
                                closure = closure - 1/2*eGamma{i,j}*( (H*ETA)\(e*n{l}*n{k}*C{i,j,k,l}*H_gamma*n{m}*e'*eU{m}') );
                                penalty = penalty - 1/2*eGamma{i,j}*( (H*ETA)\(e*n{l}*n{k}*C{i,j,k,l}*H_gamma*nV{m}*ev'*eV{m}') );
                            end
                        end
                    end
                end
            end
            %-------------------------------------------------

        end

        function [closure, penalty, forcingPenalties] = interfaceNormalTangential(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            tuning = type.tuning;

            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            u = obj;
            v = neighbour_scheme;

            dim = obj.dim;

            n       = u.getNormal(boundary);
            t       = u.getTangent(boundary);
            H_gamma = u.getBoundaryQuadratureForScalarField(boundary);
            e       = u.getBoundaryOperatorForScalarField('e', boundary);

            ev      = v.getBoundaryOperatorForScalarField('e', neighbour_boundary);

            H       = u.H;
            RHO     = u.RHO;
            ETA     = u.ETA;
            C       = u.C;
            Eu      = u.Eu;
            eU      = u.eU;
            eGamma  = u.eGamma;
            Egamma  = u.Egamma;

            CV       = v.C;
            Ev      = v.Eu;
            eV      = v.eU;
            eGammaV = v.eGamma;
            nV      = v.getNormal(neighbour_boundary);
            tV      = v.getTangent(neighbour_boundary);

            % Reduce stiffness tensors to boundary size
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            C{i,j,k,l} = e'*C{i,j,k,l}*e;
                            CV{i,j,k,l} = ev'*CV{i,j,k,l}*ev;
                        end
                    end
                end
            end

            % Get elastic closure and penalty
            [closure, penalty, forcingPenalties] = obj.elasticObj.interface(boundary, v.elasticObj, neighbour_boundary, type);
            closure = Eu*closure*Eu';
            penalty = Eu*penalty*Ev';

            for i = 1:numel(forcingPenalties)
                forcingPenalties{i} = Eu*forcingPenalties{i};
            end
            forcing_u_n = forcingPenalties{1};
            forcing_u_t = forcingPenalties{2};

            % ------ Traction coupling, viscous part -----------
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            for m = 1:dim
                                % Normal component
                                closure = closure + 1/2*eU{m}*( (RHO*H)\(e*n{m}*H_gamma*n{j}*n{i}*C{i,j,k,l}*e'*eGamma{k,l}') );
                                penalty = penalty - 1/2*eU{m}*( (RHO*H)\(e*n{m}*H_gamma*nV{j}*nV{i}*CV{i,j,k,l}*ev'*eGammaV{k,l}') );

                                % Tangential component
                                closure = closure + 1/2*eU{m}*( (RHO*H)\(e*t{m}*H_gamma*t{j}*n{i}*C{i,j,k,l}*e'*eGamma{k,l}') );
                                penalty = penalty - 1/2*eU{m}*( (RHO*H)\(e*t{m}*H_gamma*tV{j}*nV{i}*CV{i,j,k,l}*ev'*eGammaV{k,l}') );
                            end
                        end
                    end
                end
            end
            %-------------------------------------------------

            % --- Displacement coupling ----------------------
            % Add penalty to strain rate eq
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            for m = 1:dim
                                % Normal component
                                closure = closure           - 1/2*eGamma{i,j}*( (H*ETA)\(e*n{l}*n{k}*C{i,j,k,l}*H_gamma*n{m}*e'*eU{m}') );
                                penalty = penalty           - 1/2*eGamma{i,j}*( (H*ETA)\(e*n{l}*n{k}*C{i,j,k,l}*H_gamma*nV{m}*ev'*eV{m}') );


                                % Tangential component
                                closure = closure           - 1/2*eGamma{i,j}*( (H*ETA)\(e*t{l}*n{k}*C{i,j,k,l}*H_gamma*t{m}*e'*eU{m}') );
                                penalty = penalty           - 1/2*eGamma{i,j}*( (H*ETA)\(e*t{l}*n{k}*C{i,j,k,l}*H_gamma*tV{m}*ev'*eV{m}') );
                            end
                            forcing_u_n = forcing_u_n   + 1/2*eGamma{i,j}*( (H*ETA)\(e*n{l}*n{k}*C{i,j,k,l}*H_gamma) );
                            forcing_u_t = forcing_u_t   + 1/2*eGamma{i,j}*( (H*ETA)\(e*t{l}*n{k}*C{i,j,k,l}*H_gamma) );
                        end
                    end
                end
            end
            %-------------------------------------------------

            forcingPenalties{1} = forcing_u_n;
            forcingPenalties{2} = forcing_u_t;

        end

        % Returns the outward unit normal vector for the boundary specified by the string boundary.
        % n is a cell of diagonal matrices for each normal component, n{1} = n_1, n{2} = n_2.
        function n = getNormal(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            n = obj.(['n_' boundary]);
        end

        % Returns the unit tangent vector for the boundary specified by the string boundary.
        % t is a cell of diagonal matrices for each normal component, t{1} = t_1, t{2} = t_2.
        function t = getTangent(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            t = obj.(['tangent_' boundary]);
        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op -- string
        function o = getBoundaryOperator(obj, op, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})
            assertIsMember(op, {'e', 'e1', 'e2', 'tau', 'tau1', 'tau2', 'en', 'et', 'tau_n', 'tau_t'})

            o = obj.([op, '_', boundary]);

        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op -- string
        function o = getBoundaryOperatorForScalarField(obj, op, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})
            assertIsMember(op, {'e'})

            switch op

                case 'e'
                    o = obj.(['e_scalar', '_', boundary]);
            end

        end

        % Returns the boundary operator T_ij (cell format) for the boundary specified by the string boundary.
        % Formula: tau_i = T_ij u_j
        % op -- string
        function T = getBoundaryTractionOperator(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            T = obj.(['T', '_', boundary]);
        end

        % Returns square boundary quadrature matrix, of dimension
        % corresponding to the number of boundary unknowns
        %
        % boundary -- string
        function H = getBoundaryQuadrature(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            H = obj.getBoundaryQuadratureForScalarField(boundary);
            I_dim = speye(obj.dim, obj.dim);
            H = kron(H, I_dim);
        end

        % Returns square boundary quadrature matrix, of dimension
        % corresponding to the number of boundary grid points
        %
        % boundary -- string
        function H_b = getBoundaryQuadratureForScalarField(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            H_b = obj.(['H_', boundary]);
        end

        function N = size(obj)
            N = (obj.dim + obj.dim^2)*prod(obj.m);
        end
    end
end
