classdef Elastic2dCurvilinearAnisotropic < scheme.Scheme

% Discretizes the elastic wave equation:
% rho u_{i,tt} = dj C_{ijkl} dk u_j
% in curvilinear coordinates.
% opSet should be cell array of opSets, one per dimension. This
% is useful if we have periodic BC in one direction.
% Assumes fully compatible operators.

    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid
        dim

        order % Order of accuracy for the approximation

        % Diagonal matrices for variable coefficients
        J, Ji
        RHO % Density
        C   % Elastic stiffness tensor

        D  % Total operator

        K % Transformation gradient
        Dx, Dy % Physical derivatives
        sigma % Cell matrix of physical stress operators
        n_w, n_e, n_s, n_n % Physical normals
        tangent_w, tangent_e, tangent_s, tangent_n % Physical tangents

        % Boundary operators in cell format, used for BC
        T_w, T_e, T_s, T_n

        % Traction operators
        tau_w, tau_e, tau_s, tau_n      % Return vector field
        tau1_w, tau1_e, tau1_s, tau1_n  % Return scalar field
        tau2_w, tau2_e, tau2_s, tau2_n  % Return scalar field
        tau_n_w, tau_n_e, tau_n_s, tau_n_n % Return scalar field
        tau_t_w, tau_t_e, tau_t_s, tau_t_n % Return scalar field

        % Inner products
        H

        % Boundary inner products (for scalar field)
        H_w, H_e, H_s, H_n

        % Surface Jacobian vectors
        s_w, s_e, s_s, s_n

        % Boundary restriction operators
        e_w, e_e, e_s, e_n      % Act on vector field, return vector field at boundary
        e1_w, e1_e, e1_s, e1_n  % Act on vector field, return scalar field at boundary
        e2_w, e2_e, e2_s, e2_n  % Act on vector field, return scalar field at boundary
        e_scalar_w, e_scalar_e, e_scalar_s, e_scalar_n; % Act on scalar field, return scalar field
        en_w, en_e, en_s, en_n  % Act on vector field, return normal component
        et_w, et_e, et_s, et_n  % Act on vector field, return tangential component

        % E{i}^T picks out component i
        E

        % Elastic2dVariableAnisotropic object for reference domain
        refObj
    end

    methods

        % The coefficients can either be function handles or grid functions
        % optFlag -- if true, extra computations are performed, which may be helpful for optimization.
        function obj = Elastic2dCurvilinearAnisotropic(g, order, rho, C, opSet, optFlag, hollow)
            default_arg('hollow', false);
            default_arg('rho', @(x,y) 0*x+1);
            default_arg('opSet',{@sbp.D2VariableCompatible, @sbp.D2VariableCompatible});
            default_arg('optFlag', false);
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

            m = g.size();
            m_tot = g.N();

            % 1D operators
            m_u = m(1);
            m_v = m(2);
            ops_u = opSet{1}(m_u, {0, 1}, order);
            ops_v = opSet{2}(m_v, {0, 1}, order);

            h_u = ops_u.h;
            h_v = ops_v.h;

            I_u = speye(m_u);
            I_v = speye(m_v);

            D1_u = ops_u.D1;
            H_u =  ops_u.H;
            Hi_u = ops_u.HI;
            e_l_u = ops_u.e_l;
            e_r_u = ops_u.e_r;
            d1_l_u = ops_u.d1_l;
            d1_r_u = ops_u.d1_r;

            D1_v = ops_v.D1;
            H_v =  ops_v.H;
            Hi_v = ops_v.HI;
            e_l_v = ops_v.e_l;
            e_r_v = ops_v.e_r;
            d1_l_v = ops_v.d1_l;
            d1_r_v = ops_v.d1_r;


            % Logical operators
            Du = kr(D1_u,I_v);
            Dv = kr(I_u,D1_v);

            e_w  = kr(e_l_u,I_v);
            e_e  = kr(e_r_u,I_v);
            e_s  = kr(I_u,e_l_v);
            e_n  = kr(I_u,e_r_v);

            % Metric coefficients
            coords = g.points();
            x = coords(:,1);
            y = coords(:,2);
            
            % When computing the metric derivatives, we can't use periodic operators
            if isequal(opSet{1},@sbp.D2VariablePeriodic)
                ops_metric = sbp.D2Standard(m_u, {0, 1-1/m_u}, order);
                Du_metric = kr(ops_metric.D1,I_v);
            else
                Du_metric = Du;
            end
            
            if isequal(opSet{2},@sbp.D2VariablePeriodic)
                ops_metric = sbp.D2Standard(m_v, {0, 1-1/m_v}, order);
                Dv_metric = kr(I_u,ops_metric.D1);
            else
                Dv_metric = Dv;
            end

            x_u = Du_metric*x;
            x_v = Dv_metric*x;
            y_u = Du_metric*y;
            y_v = Dv_metric*y;

            J = x_u.*y_v - x_v.*y_u;

            K = cell(dim, dim);
            K{1,1} = y_v./J;
            K{1,2} = -y_u./J;
            K{2,1} = -x_v./J;
            K{2,2} = x_u./J;
            obj.K = K;

            % Physical derivatives
            % TODO: Fix for periodic operators
            obj.Dx = spdiag( y_v./J)*Du + spdiag(-y_u./J)*Dv;
            obj.Dy = spdiag(-x_v./J)*Du + spdiag( x_u./J)*Dv;

            % Wrap around Aniosotropic Cartesian
            rho_tilde = J.*rho;

            PHI = cell(dim,dim,dim,dim);
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            PHI{i,j,k,l} = 0*C{i,j,k,l};
                            for m = 1:dim
                                for n = 1:dim
                                    PHI{i,j,k,l} = PHI{i,j,k,l} + J.*K{m,i}.*C{m,j,n,l}.*K{n,k};
                                end
                            end
                        end
                    end
                end
            end

            refObj = scheme.Elastic2dVariableAnisotropic(g.logic, order, rho_tilde, PHI, opSet, [], hollow);

            %---- Set object properties ------
            obj.RHO = spdiag(rho);

            % Volume quadrature
            obj.J = spdiag(J);
            obj.Ji = spdiag(1./J);
            obj.H = obj.J*kr(H_u,H_v);

            % Boundary quadratures
            s_w = sqrt((e_w'*x_v).^2 + (e_w'*y_v).^2);
            s_e = sqrt((e_e'*x_v).^2 + (e_e'*y_v).^2);
            s_s = sqrt((e_s'*x_u).^2 + (e_s'*y_u).^2);
            s_n = sqrt((e_n'*x_u).^2 + (e_n'*y_u).^2);
            obj.s_w = s_w;
            obj.s_e = s_e;
            obj.s_s = s_s;
            obj.s_n = s_n;

            obj.H_w = H_v*spdiag(s_w);
            obj.H_e = H_v*spdiag(s_e);
            obj.H_s = H_u*spdiag(s_s);
            obj.H_n = H_u*spdiag(s_n);

            % Restriction operators
            obj.e_w = refObj.e_w;
            obj.e_e = refObj.e_e;
            obj.e_s = refObj.e_s;
            obj.e_n = refObj.e_n;

            % Adapt things from reference object
            obj.D = refObj.D;
            obj.E = refObj.E;

            obj.e1_w = refObj.e1_w;
            obj.e1_e = refObj.e1_e;
            obj.e1_s = refObj.e1_s;
            obj.e1_n = refObj.e1_n;

            obj.e2_w = refObj.e2_w;
            obj.e2_e = refObj.e2_e;
            obj.e2_s = refObj.e2_s;
            obj.e2_n = refObj.e2_n;

            obj.e_scalar_w = refObj.e_scalar_w;
            obj.e_scalar_e = refObj.e_scalar_e;
            obj.e_scalar_s = refObj.e_scalar_s;
            obj.e_scalar_n = refObj.e_scalar_n;

            e1_w = obj.e1_w;
            e1_e = obj.e1_e;
            e1_s = obj.e1_s;
            e1_n = obj.e1_n;

            e2_w = obj.e2_w;
            e2_e = obj.e2_e;
            e2_s = obj.e2_s;
            e2_n = obj.e2_n;

            obj.tau1_w = (spdiag(1./s_w)*refObj.tau1_w')';
            obj.tau1_e = (spdiag(1./s_e)*refObj.tau1_e')';
            obj.tau1_s = (spdiag(1./s_s)*refObj.tau1_s')';
            obj.tau1_n = (spdiag(1./s_n)*refObj.tau1_n')';

            obj.tau2_w = (spdiag(1./s_w)*refObj.tau2_w')';
            obj.tau2_e = (spdiag(1./s_e)*refObj.tau2_e')';
            obj.tau2_s = (spdiag(1./s_s)*refObj.tau2_s')';
            obj.tau2_n = (spdiag(1./s_n)*refObj.tau2_n')';

            obj.tau_w = (refObj.e_w'*obj.e1_w*obj.tau1_w')' + (refObj.e_w'*obj.e2_w*obj.tau2_w')';
            obj.tau_e = (refObj.e_e'*obj.e1_e*obj.tau1_e')' + (refObj.e_e'*obj.e2_e*obj.tau2_e')';
            obj.tau_s = (refObj.e_s'*obj.e1_s*obj.tau1_s')' + (refObj.e_s'*obj.e2_s*obj.tau2_s')';
            obj.tau_n = (refObj.e_n'*obj.e1_n*obj.tau1_n')' + (refObj.e_n'*obj.e2_n*obj.tau2_n')';

            % Physical normals
            e_w = obj.e_scalar_w;
            e_e = obj.e_scalar_e;
            e_s = obj.e_scalar_s;
            e_n = obj.e_scalar_n;

            e_w_vec = obj.e_w;
            e_e_vec = obj.e_e;
            e_s_vec = obj.e_s;
            e_n_vec = obj.e_n;

            nu_w = [-1,0];
            nu_e = [1,0];
            nu_s = [0,-1];
            nu_n = [0,1];

            obj.n_w = cell(2,1);
            obj.n_e = cell(2,1);
            obj.n_s = cell(2,1);
            obj.n_n = cell(2,1);

            % Compute normal and rotate (exactly!) 90 degrees counter-clockwise to get tangent
            n_w_1 = (1./s_w).*e_w'*(J.*(K{1,1}*nu_w(1) + K{1,2}*nu_w(2)));
            n_w_2 = (1./s_w).*e_w'*(J.*(K{2,1}*nu_w(1) + K{2,2}*nu_w(2)));
            obj.n_w{1} = spdiag(n_w_1);
            obj.n_w{2} = spdiag(n_w_2);
            obj.tangent_w = {-obj.n_w{2}, obj.n_w{1}};

            n_e_1 = (1./s_e).*e_e'*(J.*(K{1,1}*nu_e(1) + K{1,2}*nu_e(2)));
            n_e_2 = (1./s_e).*e_e'*(J.*(K{2,1}*nu_e(1) + K{2,2}*nu_e(2)));
            obj.n_e{1} = spdiag(n_e_1);
            obj.n_e{2} = spdiag(n_e_2);
            obj.tangent_e = {-obj.n_e{2}, obj.n_e{1}};

            n_s_1 = (1./s_s).*e_s'*(J.*(K{1,1}*nu_s(1) + K{1,2}*nu_s(2)));
            n_s_2 = (1./s_s).*e_s'*(J.*(K{2,1}*nu_s(1) + K{2,2}*nu_s(2)));
            obj.n_s{1} = spdiag(n_s_1);
            obj.n_s{2} = spdiag(n_s_2);
            obj.tangent_s = {-obj.n_s{2}, obj.n_s{1}};

            n_n_1 = (1./s_n).*e_n'*(J.*(K{1,1}*nu_n(1) + K{1,2}*nu_n(2)));
            n_n_2 = (1./s_n).*e_n'*(J.*(K{2,1}*nu_n(1) + K{2,2}*nu_n(2)));
            obj.n_n{1} = spdiag(n_n_1);
            obj.n_n{2} = spdiag(n_n_2);
            obj.tangent_n = {-obj.n_n{2}, obj.n_n{1}};

            % Operators that extract the normal component
            obj.en_w = (obj.n_w{1}*obj.e1_w' + obj.n_w{2}*obj.e2_w')';
            obj.en_e = (obj.n_e{1}*obj.e1_e' + obj.n_e{2}*obj.e2_e')';
            obj.en_s = (obj.n_s{1}*obj.e1_s' + obj.n_s{2}*obj.e2_s')';
            obj.en_n = (obj.n_n{1}*obj.e1_n' + obj.n_n{2}*obj.e2_n')';

            % Operators that extract the tangential component
            obj.et_w = (obj.tangent_w{1}*obj.e1_w' + obj.tangent_w{2}*obj.e2_w')';
            obj.et_e = (obj.tangent_e{1}*obj.e1_e' + obj.tangent_e{2}*obj.e2_e')';
            obj.et_s = (obj.tangent_s{1}*obj.e1_s' + obj.tangent_s{2}*obj.e2_s')';
            obj.et_n = (obj.tangent_n{1}*obj.e1_n' + obj.tangent_n{2}*obj.e2_n')';

            obj.tau_n_w = (obj.n_w{1}*obj.tau1_w' + obj.n_w{2}*obj.tau2_w')';
            obj.tau_n_e = (obj.n_e{1}*obj.tau1_e' + obj.n_e{2}*obj.tau2_e')';
            obj.tau_n_s = (obj.n_s{1}*obj.tau1_s' + obj.n_s{2}*obj.tau2_s')';
            obj.tau_n_n = (obj.n_n{1}*obj.tau1_n' + obj.n_n{2}*obj.tau2_n')';

            obj.tau_t_w = (obj.tangent_w{1}*obj.tau1_w' + obj.tangent_w{2}*obj.tau2_w')';
            obj.tau_t_e = (obj.tangent_e{1}*obj.tau1_e' + obj.tangent_e{2}*obj.tau2_e')';
            obj.tau_t_s = (obj.tangent_s{1}*obj.tau1_s' + obj.tangent_s{2}*obj.tau2_s')';
            obj.tau_t_n = (obj.tangent_n{1}*obj.tau1_n' + obj.tangent_n{2}*obj.tau2_n')';

            % Stress operators
            sigma = cell(dim, dim);
            D1 = {obj.Dx, obj.Dy};
            E = obj.E;
            N = length(obj.RHO);
            for i = 1:dim
                for j = 1:dim
                    sigma{i,j} = sparse(N,2*N);
                    for k = 1:dim
                        for l = 1:dim
                            sigma{i,j} = sigma{i,j} + obj.C{i,j,k,l}*D1{k}*E{l}';
                        end
                    end
                end
            end
            obj.sigma = sigma;

            % Misc.
            obj.refObj = refObj;
            obj.m = refObj.m;
            obj.h = refObj.h;
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

            switch component

            % If conditions on Cartesian components
            case {1,2}
                [closure, penalty] = obj.refObj.boundary_condition(boundary, bc, tuning);

            % If conditions on normal or tangential components
            case {'n', 't'}

                switch component
                    case 'n'
                        en = obj.getBoundaryOperator('en', boundary);
                    case 't'
                        en = obj.getBoundaryOperator('et', boundary);
                end
                e1 = obj.getBoundaryOperator('e1', boundary);

                bc1 = {1, type};
                [c1, p1] = obj.refObj.boundary_condition(boundary, bc1, tuning);
                bc2 = {2, type};
                c2 = obj.refObj.boundary_condition(boundary, bc2, tuning);

                switch type
                case {'F','f','Free','free','traction','Traction','t','T'}
                    closure = en*en'*(c1+c2);
                    penalty = en*e1'*p1;
                case {'D','d','dirichlet','Dirichlet','displacement','Displacement'}
                    [closure, penalty] = obj.displacementBCNormalTangential(boundary, bc, tuning);
                end

            end

            switch type
            case {'F','f','Free','free','traction','Traction','t','T'}

                s = obj.(['s_' boundary]);
                s = spdiag(s);
                penalty = penalty*s;

            end
        end

        function [closure, penalty] = displacementBCNormalTangential(obj, boundary, bc, tuning)
            u = obj;

            component = bc{1};
            type = bc{2};

            switch component
            case 'n'
                en      = u.getBoundaryOperator('en', boundary);
                tau_n   = u.getBoundaryOperator('tau_n', boundary);
                N       = u.getNormal(boundary);
            case 't'
                en      = u.getBoundaryOperator('et', boundary);
                tau_n   = u.getBoundaryOperator('tau_t', boundary);
                N       = u.getTangent(boundary);
            end

            % Operators
            e       = u.getBoundaryOperatorForScalarField('e', boundary);
            h11     = u.getBorrowing(boundary);
            n      = u.getNormal(boundary);

            C = u.C;
            Ji = u.Ji;
            s = spdiag(u.(['s_' boundary]));
            m_tot = u.grid.N();

            Hi      = u.E{1}*inv(u.H)*u.E{1}' + u.E{2}*inv(u.H)*u.E{2}';
            RHOi    = u.E{1}*inv(u.RHO)*u.E{1}' + u.E{2}*inv(u.RHO)*u.E{2}';

            H_gamma         = u.getBoundaryQuadratureForScalarField(boundary);
            dim             = u.dim;

            % Preallocate
            [~, m_int] = size(H_gamma);
            closure = sparse(dim*m_tot, dim*m_tot);
            penalty = sparse(dim*m_tot, m_int);

            % Term 1: The symmetric term
            Z = sparse(m_int, m_int);
            for i = 1:dim
                for j = 1:dim
                    for l = 1:dim
                        for k = 1:dim
                            Z = Z + n{i}*N{j}*e'*Ji*C{i,j,k,l}*e*n{k}*N{l};
                        end
                    end
                end
            end

            Z = -tuning*dim*1/h11*s*Z;
            closure = closure + en*H_gamma*Z*en';
            penalty = penalty - en*H_gamma*Z;

            % Term 2: The symmetrizing term
            closure = closure + tau_n*H_gamma*en';
            penalty = penalty - tau_n*H_gamma;

            % Multiply all terms by inverse of density x quadraure
            closure = RHOi*Hi*closure;
            penalty = RHOi*Hi*penalty;
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
                [closure, penalty] = obj.refObj.interface(boundary,neighbour_scheme.refObj,neighbour_boundary,type);
            case 'normalTangential'
                [closure, penalty, forcingPenalties] = obj.interfaceNormalTangential(boundary,neighbour_scheme,neighbour_boundary,type);
            case 'frictionalFault'
                [closure, penalty, forcingPenalties] = obj.interfaceFrictionalFault(boundary,neighbour_scheme,neighbour_boundary,type);
            end

        end

        function [closure, penalty, forcingPenalties] = interfaceFrictionalFault(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            tuning = type.tuning;

            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain

            forcingPenalties = cell(1, 1);
            u = obj;
            v = neighbour_scheme;

            % Operators, u side
            e_u       = u.getBoundaryOperatorForScalarField('e', boundary);
            en_u       = u.getBoundaryOperator('en', boundary);
            tau_n_u     = u.getBoundaryOperator('tau_n', boundary);
            h11_u     = u.getBorrowing(boundary);
            n_u      = u.getNormal(boundary);

            C_u = u.C;
            Ji_u = u.Ji;
            s_u = spdiag(u.(['s_' boundary]));
            m_tot_u = u.grid.N();

            % Operators, v side
            e_v       = v.getBoundaryOperatorForScalarField('e', neighbour_boundary);
            en_v       = v.getBoundaryOperator('en', neighbour_boundary);
            tau_n_v     = v.getBoundaryOperator('tau_n', neighbour_boundary);
            h11_v     = v.getBorrowing(neighbour_boundary);
            n_v      = v.getNormal(neighbour_boundary);

            C_v = v.C;
            Ji_v = v.Ji;
            s_v = spdiag(v.(['s_' neighbour_boundary]));
            m_tot_v = v.grid.N();

            % Operators that are only required for own domain
            Hi      = u.E{1}*inv(u.H)*u.E{1}' + u.E{2}*inv(u.H)*u.E{2}';
            RHOi    = u.E{1}*inv(u.RHO)*u.E{1}' + u.E{2}*inv(u.RHO)*u.E{2}';

            % Shared operators
            H_gamma         = u.getBoundaryQuadratureForScalarField(boundary);
            dim             = u.dim;

            % Preallocate
            [~, m_int] = size(H_gamma);
            closure = sparse(dim*m_tot_u, dim*m_tot_u);
            penalty = sparse(dim*m_tot_u, dim*m_tot_v);

            % Continuity of normal displacement, term 1: The symmetric term
            Z_u = sparse(m_int, m_int);
            Z_v = sparse(m_int, m_int);
            for i = 1:dim
                for j = 1:dim
                    for l = 1:dim
                        for k = 1:dim
                            Z_u = Z_u + n_u{i}*n_u{j}*e_u'*Ji_u*C_u{i,j,k,l}*e_u*n_u{k}*n_u{l};
                            Z_v = Z_v + n_v{i}*n_v{j}*e_v'*Ji_v*C_v{i,j,k,l}*e_v*n_v{k}*n_v{l};
                        end
                    end
                end
            end

            Z = -tuning*dim*( 1/(4*h11_u)*s_u*Z_u + 1/(4*h11_v)*s_v*Z_v );
            closure = closure + en_u*H_gamma*Z*en_u';
            penalty = penalty + en_u*H_gamma*Z*en_v';

            % Continuity of normal displacement, term 2: The symmetrizing term
            closure = closure + 1/2*tau_n_u*H_gamma*en_u';
            penalty = penalty + 1/2*tau_n_u*H_gamma*en_v';

            % Continuity of normal traction
            closure = closure - 1/2*en_u*H_gamma*tau_n_u';
            penalty = penalty + 1/2*en_u*H_gamma*tau_n_v';
            forcing_tau_n = 1/2*en_u*H_gamma;

            % Multiply all normal component terms by inverse of density x quadraure
            closure = RHOi*Hi*closure;
            penalty = RHOi*Hi*penalty;
            forcing_tau_n = RHOi*Hi*forcing_tau_n;

            % ---- Tangential tractions are imposed just like traction BC ------
            closure = closure + obj.boundary_condition(boundary, {'t', 't'});

            forcingPenalties{1} = forcing_tau_n;

        end

        % Same interface conditions as in interfaceStandard, but imposed in the normal-tangential
        % coordinate system. For the isotropic case, the components decouple nicely.
        % The resulting scheme is not identical to that of interfaceStandard. This appears to be better.
        function [closure, penalty, forcingPenalties, Zt] = interfaceNormalTangential(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            tuning = type.tuning;

            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain

            forcingPenalties = cell(2, 1);
            u = obj;
            v = neighbour_scheme;

            % Operators, u side
            e_u         = u.getBoundaryOperatorForScalarField('e', boundary);
            en_u        = u.getBoundaryOperator('en', boundary);
            et_u        = u.getBoundaryOperator('et', boundary);
            tau_n_u     = u.getBoundaryOperator('tau_n', boundary);
            tau_t_u     = u.getBoundaryOperator('tau_t', boundary);
            h11_u       = u.getBorrowing(boundary);
            n_u         = u.getNormal(boundary);
            t_u         = u.getTangent(boundary);

            C_u = u.C;
            Ji_u = u.Ji;
            s_u = spdiag(u.(['s_' boundary]));
            m_tot_u = u.grid.N();

            % Operators, v side
            e_v         = v.getBoundaryOperatorForScalarField('e', neighbour_boundary);
            en_v        = v.getBoundaryOperator('en', neighbour_boundary);
            et_v        = v.getBoundaryOperator('et', neighbour_boundary);
            tau_n_v     = v.getBoundaryOperator('tau_n', neighbour_boundary);
            tau_t_v     = v.getBoundaryOperator('tau_t', neighbour_boundary);
            h11_v       = v.getBorrowing(neighbour_boundary);
            n_v         = v.getNormal(neighbour_boundary);
            t_v         = v.getTangent(neighbour_boundary);

            C_v = v.C;
            Ji_v = v.Ji;
            s_v = spdiag(v.(['s_' neighbour_boundary]));
            m_tot_v = v.grid.N();

            % Operators that are only required for own domain
            Hi      = u.E{1}*inv(u.H)*u.E{1}' + u.E{2}*inv(u.H)*u.E{2}';
            RHOi    = u.E{1}*inv(u.RHO)*u.E{1}' + u.E{2}*inv(u.RHO)*u.E{2}';

            % Shared operators
            H_gamma         = u.getBoundaryQuadratureForScalarField(boundary);
            dim             = u.dim;

            % Preallocate
            [~, m_int] = size(H_gamma);
            closure = sparse(dim*m_tot_u, dim*m_tot_u);
            penalty = sparse(dim*m_tot_u, dim*m_tot_v);

            % -- Continuity of displacement, term 1: The symmetric term ---
            Zn_u = sparse(m_int, m_int);
            Zn_v = sparse(m_int, m_int);
            Zt_u = sparse(m_int, m_int);
            Zt_v = sparse(m_int, m_int);
            for i = 1:dim
                for j = 1:dim
                    for l = 1:dim
                        for k = 1:dim
                            % Penalty strength for normal component
                            Zn_u = Zn_u + n_u{i}*n_u{j}*e_u'*Ji_u*C_u{i,j,k,l}*e_u*n_u{k}*n_u{l};
                            Zn_v = Zn_v + n_v{i}*n_v{j}*e_v'*Ji_v*C_v{i,j,k,l}*e_v*n_v{k}*n_v{l};

                            % Penalty strength for tangential component
                            Zt_u = Zt_u + n_u{i}*t_u{j}*e_u'*Ji_u*C_u{i,j,k,l}*e_u*n_u{k}*t_u{l};
                            Zt_v = Zt_v + n_v{i}*t_v{j}*e_v'*Ji_v*C_v{i,j,k,l}*e_v*n_v{k}*t_v{l};
                        end
                    end
                end
            end

            Zn = -tuning*dim*( 1/(4*h11_u)*s_u*Zn_u + 1/(4*h11_v)*s_v*Zn_v );
            Zt = -tuning*dim*( 1/(4*h11_u)*s_u*Zt_u + 1/(4*h11_v)*s_v*Zt_v );

            % Continuity of normal component
            closure = closure + en_u*H_gamma*Zn*en_u';
            penalty = penalty + en_u*H_gamma*Zn*en_v';
            forcing_u_n = -en_u*H_gamma*Zn;

            % Continuity of tangential component
            closure = closure + et_u*H_gamma*Zt*et_u';
            penalty = penalty + et_u*H_gamma*Zt*et_v';
            forcing_u_t = -et_u*H_gamma*Zt;
            %------------------------------------------------------------------

            % --- Continuity of displacement, term 2: The symmetrizing term

            % Continuity of normal displacement
            closure = closure + 1/2*tau_n_u*H_gamma*en_u';
            penalty = penalty + 1/2*tau_n_u*H_gamma*en_v';
            forcing_u_n = forcing_u_n - 1/2*tau_n_u*H_gamma;

            % Continuity of tangential displacement
            closure = closure + 1/2*tau_t_u*H_gamma*et_u';
            penalty = penalty + 1/2*tau_t_u*H_gamma*et_v';
            forcing_u_t = forcing_u_t - 1/2*tau_t_u*H_gamma;
            % ------------------------------------------------------------------

            % --- Continuity of tractions -----------------------------

            % Continuity of normal traction
            closure = closure - 1/2*en_u*H_gamma*tau_n_u';
            penalty = penalty + 1/2*en_u*H_gamma*tau_n_v';

            % Continuity of tangential traction
            closure = closure - 1/2*et_u*H_gamma*tau_t_u';
            penalty = penalty + 1/2*et_u*H_gamma*tau_t_v';
            %--------------------------------------------------------------------

            % Multiply all terms by inverse of density x quadraure
            closure = RHOi*Hi*closure;
            penalty = RHOi*Hi*penalty;
            forcing_u_n = RHOi*Hi*forcing_u_n;
            forcing_u_t = RHOi*Hi*forcing_u_t;

            forcingPenalties{1} = forcing_u_n;
            forcingPenalties{2} = forcing_u_t;

        end


        % Returns h11 for the boundary specified by the string boundary.
        % op -- string
        function h11 = getBorrowing(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            switch boundary
            case {'w','e'}
                h11 = obj.refObj.h11{1};
            case {'s', 'n'}
                h11 = obj.refObj.h11{2};
            end
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
            N = obj.dim*prod(obj.m);
        end
    end
end
