classdef Elastic2dVariable < scheme.Scheme

% Discretizes the elastic wave equation:
% rho u_{i,tt} = di lambda dj u_j + dj mu di u_j + dj mu dj u_i
% opSet should be cell array of opSets, one per dimension. This
% is useful if we have periodic BC in one direction.

    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid
        dim

        order % Order of accuracy for the approximation

        % Diagonal matrices for variable coefficients
        LAMBDA % Lame's first parameter, related to dilation
        MU     % Shear modulus
        RHO, RHOi, RHOi_kron % Density

        D % Total operator
        D1 % First derivatives

        % Second derivatives
        D2_lambda
        D2_mu

        % Boundary operators in cell format, used for BC
        T_w, T_e, T_s, T_n

        % Traction operators
        tau_w, tau_e, tau_s, tau_n      % Return vector field
        tau1_w, tau1_e, tau1_s, tau1_n  % Return scalar field
        tau2_w, tau2_e, tau2_s, tau2_n  % Return scalar field

        % Inner products
        H, Hi, Hi_kron, H_1D

        % Boundary inner products (for scalar field)
        H_w, H_e, H_s, H_n

        % Boundary restriction operators
        e_w, e_e, e_s, e_n      % Act on vector field, return vector field at boundary
        e1_w, e1_e, e1_s, e1_n  % Act on vector field, return scalar field at boundary
        e2_w, e2_e, e2_s, e2_n  % Act on vector field, return scalar field at boundary
        e_scalar_w, e_scalar_e, e_scalar_s, e_scalar_n; % Act on scalar field, return scalar field

        % E{i}^T picks out component i
        E

        % Borrowing constants of the form gamma*h, where gamma is a dimensionless constant.
        theta_R % Borrowing (d1- D1)^2 from R
        theta_H % First entry in norm matrix
        theta_M % Borrowing d1^2 from M.

        % Structures used for adjoint optimization
        B
    end

    methods

        % The coefficients can either be function handles or grid functions
        % optFlag -- if true, extra computations are performed, which may be helpful for optimization.
        function obj = Elastic2dVariable(g ,order, lambda, mu, rho, opSet, optFlag)
            default_arg('opSet',{@sbp.D2Variable, @sbp.D2Variable});
            default_arg('lambda', @(x,y) 0*x+1);
            default_arg('mu', @(x,y) 0*x+1);
            default_arg('rho', @(x,y) 0*x+1);
            default_arg('optFlag', false);
            dim = 2;

            assert(isa(g, 'grid.Cartesian'))

            if isa(lambda, 'function_handle')
                lambda = grid.evalOn(g, lambda);
            end
            if isa(mu, 'function_handle')
                mu = grid.evalOn(g, mu);
            end
            if isa(rho, 'function_handle')
                rho = grid.evalOn(g, rho);
            end

            m = g.size();
            m_tot = g.N();

            h = g.scaling();
            lim = g.lim;
            if isempty(lim)
                x = g.x;
                lim = cell(length(x),1);
                for i = 1:length(x)
                    lim{i} = {min(x{i}), max(x{i})};
                end
            end

            % 1D operators
            ops = cell(dim,1);
            for i = 1:dim
                ops{i} = opSet{i}(m(i), lim{i}, order);
            end

            % Borrowing constants
            for i = 1:dim
                obj.theta_R{i} = h(i)*ops{i}.borrowing.R.delta_D;
                obj.theta_H{i} = h(i)*ops{i}.borrowing.H11;
                obj.theta_M{i} = h(i)*ops{i}.borrowing.M.d1;
            end

            I = cell(dim,1);
            D1 = cell(dim,1);
            D2 = cell(dim,1);
            H = cell(dim,1);
            Hi = cell(dim,1);
            e_0 = cell(dim,1);
            e_m = cell(dim,1);
            d1_0 = cell(dim,1);
            d1_m = cell(dim,1);

            for i = 1:dim
                I{i} = speye(m(i));
                D1{i} = ops{i}.D1;
                D2{i} = ops{i}.D2;
                H{i} =  ops{i}.H;
                Hi{i} = ops{i}.HI;
                e_0{i} = ops{i}.e_l;
                e_m{i} = ops{i}.e_r;
                d1_0{i} = ops{i}.d1_l;
                d1_m{i} = ops{i}.d1_r;
            end

            %====== Assemble full operators ========
            LAMBDA = spdiag(lambda);
            obj.LAMBDA = LAMBDA;
            MU = spdiag(mu);
            obj.MU = MU;
            RHO = spdiag(rho);
            obj.RHO = RHO;
            obj.RHOi = inv(RHO);

            obj.D1 = cell(dim,1);
            obj.D2_lambda = cell(dim,1);
            obj.D2_mu = cell(dim,1);

            % D1
            obj.D1{1} = kron(D1{1},I{2});
            obj.D1{2} = kron(I{1},D1{2});

            % Boundary restriction operators
            e_l = cell(dim,1);
            e_r = cell(dim,1);
            e_l{1} = kron(e_0{1}, I{2});
            e_l{2} = kron(I{1}, e_0{2});
            e_r{1} = kron(e_m{1}, I{2});
            e_r{2} = kron(I{1}, e_m{2});

            e_scalar_w = e_l{1};
            e_scalar_e = e_r{1};
            e_scalar_s = e_l{2};
            e_scalar_n = e_r{2};

            I_dim = speye(dim, dim);
            e_w = kron(e_scalar_w, I_dim);
            e_e = kron(e_scalar_e, I_dim);
            e_s = kron(e_scalar_s, I_dim);
            e_n = kron(e_scalar_n, I_dim);

            % Boundary derivatives
            d1_l = cell(dim,1);
            d1_r = cell(dim,1);
            d1_l{1} = kron(d1_0{1}, I{2});
            d1_l{2} = kron(I{1}, d1_0{2});
            d1_r{1} = kron(d1_m{1}, I{2});
            d1_r{2} = kron(I{1}, d1_m{2});


            % E{i}^T picks out component i.
            E = cell(dim,1);
            I = speye(m_tot,m_tot);
            for i = 1:dim
                e = sparse(dim,1);
                e(i) = 1;
                E{i} = kron(I,e);
            end
            obj.E = E;

            e1_w = (e_scalar_w'*E{1}')';
            e1_e = (e_scalar_e'*E{1}')';
            e1_s = (e_scalar_s'*E{1}')';
            e1_n = (e_scalar_n'*E{1}')';

            e2_w = (e_scalar_w'*E{2}')';
            e2_e = (e_scalar_e'*E{2}')';
            e2_s = (e_scalar_s'*E{2}')';
            e2_n = (e_scalar_n'*E{2}')';


            % D2
            for i = 1:dim
                obj.D2_lambda{i} = sparse(m_tot, m_tot);
                obj.D2_mu{i} = sparse(m_tot, m_tot);
            end
            ind = grid.funcToMatrix(g, 1:m_tot);

            for i = 1:m(2)
                D_lambda = D2{1}(lambda(ind(:,i)));
                D_mu = D2{1}(mu(ind(:,i)));

                p = ind(:,i);
                obj.D2_lambda{1}(p,p) = D_lambda;
                obj.D2_mu{1}(p,p) = D_mu;
            end

            for i = 1:m(1)
                D_lambda = D2{2}(lambda(ind(i,:)));
                D_mu = D2{2}(mu(ind(i,:)));

                p = ind(i,:);
                obj.D2_lambda{2}(p,p) = D_lambda;
                obj.D2_mu{2}(p,p) = D_mu;
            end

            % Quadratures
            obj.H = kron(H{1},H{2});
            obj.Hi = inv(obj.H);
            obj.H_w = H{2};
            obj.H_e = H{2};
            obj.H_s = H{1};
            obj.H_n = H{1};
            obj.H_1D = {H{1}, H{2}};

            % Differentiation matrix D (without SAT)
            D2_lambda = obj.D2_lambda;
            D2_mu = obj.D2_mu;
            D1 = obj.D1;
            D = sparse(dim*m_tot,dim*m_tot);
            d = @kroneckerDelta;    % Kronecker delta
            db = @(i,j) 1-d(i,j); % Logical not of Kronecker delta
            for i = 1:dim
                for j = 1:dim
                    D = D + E{i}*inv(RHO)*( d(i,j)*D2_lambda{i}*E{j}' +...
                                            db(i,j)*D1{i}*LAMBDA*D1{j}*E{j}' ...
                                          );
                    D = D + E{i}*inv(RHO)*( d(i,j)*D2_mu{i}*E{j}' +...
                                            db(i,j)*D1{j}*MU*D1{i}*E{j}' + ...
                                            D2_mu{j}*E{i}' ...
                                          );
                end
            end
            obj.D = D;
            %=========================================%'

            % Numerical traction operators for BC.
            % Because d1 =/= e0^T*D1, the numerical tractions are different
            % at every boundary.
            %
            % Formula at boundary j: % tau^{j}_i = sum_k T^{j}_{ik} u_k
            %
            T_l = cell(dim,1);
            T_r = cell(dim,1);
            tau_l = cell(dim,1);
            tau_r = cell(dim,1);

            D1 = obj.D1;

            % Loop over boundaries
            for j = 1:dim
                T_l{j} = cell(dim,dim);
                T_r{j} = cell(dim,dim);
                tau_l{j} = cell(dim,1);
                tau_r{j} = cell(dim,1);

                LAMBDA_l = e_l{j}'*LAMBDA*e_l{j};
                LAMBDA_r = e_r{j}'*LAMBDA*e_r{j};
                MU_l = e_l{j}'*MU*e_l{j};
                MU_r = e_r{j}'*MU*e_r{j};

                [~, n_l] = size(e_l{j});
                [~, n_r] = size(e_r{j});

                % Loop over components
                for i = 1:dim
                    tau_l{j}{i} = sparse(dim*m_tot, n_l);
                    tau_r{j}{i} = sparse(dim*m_tot, n_r);
                    for k = 1:dim
                        T_l{j}{i,k} = ...
                        (-d(i,j)*LAMBDA_l*(d(i,k)*d1_l{j}' + db(i,k)*e_l{j}'*D1{k})...
                         -d(j,k)*MU_l*(d(i,j)*d1_l{j}' + db(i,j)*e_l{j}'*D1{i})...
                         -d(i,k)*MU_l*d1_l{j}')';

                        T_r{j}{i,k} = ...
                        (d(i,j)*LAMBDA_r*(d(i,k)*d1_r{j}' + db(i,k)*e_r{j}'*D1{k})...
                        +d(j,k)*MU_r*(d(i,j)*d1_r{j}' + db(i,j)*e_r{j}'*D1{i})...
                        +d(i,k)*MU_r*d1_r{j}')';

                        tau_l{j}{i} = tau_l{j}{i} + (T_l{j}{i,k}'*E{k}')';
                        tau_r{j}{i} = tau_r{j}{i} + (T_r{j}{i,k}'*E{k}')';
                    end

                end
            end

            % Traction tensors, T_ij
            obj.T_w = T_l{1};
            obj.T_e = T_r{1};
            obj.T_s = T_l{2};
            obj.T_n = T_r{2};

            % Restriction operators
            obj.e_w = e_w;
            obj.e_e = e_e;
            obj.e_s = e_s;
            obj.e_n = e_n;

            obj.e1_w = e1_w;
            obj.e1_e = e1_e;
            obj.e1_s = e1_s;
            obj.e1_n = e1_n;

            obj.e2_w = e2_w;
            obj.e2_e = e2_e;
            obj.e2_s = e2_s;
            obj.e2_n = e2_n;

            obj.e_scalar_w = e_scalar_w;
            obj.e_scalar_e = e_scalar_e;
            obj.e_scalar_s = e_scalar_s;
            obj.e_scalar_n = e_scalar_n;

            % First component of traction
            obj.tau1_w = tau_l{1}{1};
            obj.tau1_e = tau_r{1}{1};
            obj.tau1_s = tau_l{2}{1};
            obj.tau1_n = tau_r{2}{1};

            % Second component of traction
            obj.tau2_w = tau_l{1}{2};
            obj.tau2_e = tau_r{1}{2};
            obj.tau2_s = tau_l{2}{2};
            obj.tau2_n = tau_r{2}{2};

            % Traction vectors
            obj.tau_w = (e_w'*e1_w*obj.tau1_w')' + (e_w'*e2_w*obj.tau2_w')';
            obj.tau_e = (e_e'*e1_e*obj.tau1_e')' + (e_e'*e2_e*obj.tau2_e')';
            obj.tau_s = (e_s'*e1_s*obj.tau1_s')' + (e_s'*e2_s*obj.tau2_s')';
            obj.tau_n = (e_n'*e1_n*obj.tau1_n')' + (e_n'*e2_n*obj.tau2_n')';

            % Kroneckered norms and coefficients
            obj.RHOi_kron = kron(obj.RHOi, I_dim);
            obj.Hi_kron = kron(obj.Hi, I_dim);

            % Misc.
            obj.m = m;
            obj.h = h;
            obj.order = order;
            obj.grid = g;
            obj.dim = dim;

            % B, used for adjoint optimization
            B = [];
            if optFlag
                B = cell(dim, 1);
                for i = 1:dim
                    B{i} = cell(m_tot, 1);
                end

                B0 = sparse(m_tot, m_tot);
                for i = 1:dim
                    for j = 1:m_tot
                        B{i}{j} = B0;
                    end
                end

                ind = grid.funcToMatrix(g, 1:m_tot);

                % Direction 1
                for k = 1:m(1)
                    c = sparse(m(1),1);
                    c(k) = 1;
                    [~, B_1D] = ops{1}.D2(c);
                    for l = 1:m(2)
                        p = ind(:,l);
                        B{1}{(k-1)*m(2) + l}(p, p) = B_1D;
                    end
                end

                % Direction 2
                for k = 1:m(2)
                    c = sparse(m(2),1);
                    c(k) = 1;
                    [~, B_1D] = ops{2}.D2(c);
                    for l = 1:m(1)
                        p = ind(l,:);
                        B{2}{(l-1)*m(2) + k}(p, p) = B_1D;
                    end
                end
            end
            obj.B = B;


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
        function [closure, penalty] = boundary_condition(obj, boundary, bc, tuning)
            default_arg('tuning', 1.2);

            assert( iscell(bc), 'The BC type must be a 2x1 cell array' );
            comp = bc{1};
            type = bc{2};
            if ischar(comp)
                comp = obj.getComponent(comp, boundary);
            end

            e       = obj.getBoundaryOperatorForScalarField('e', boundary);
            tau     = obj.getBoundaryOperator(['tau' num2str(comp)], boundary);
            T       = obj.getBoundaryTractionOperator(boundary);
            alpha   = obj.getBoundaryOperatorForScalarField('alpha', boundary);
            H_gamma = obj.getBoundaryQuadratureForScalarField(boundary);

            E = obj.E;
            Hi = obj.Hi;
            LAMBDA = obj.LAMBDA;
            MU = obj.MU;
            RHOi = obj.RHOi;

            dim = obj.dim;
            m_tot = obj.grid.N();

            % Preallocate
            [~, col] = size(tau);
            closure = sparse(dim*m_tot, dim*m_tot);
            penalty = sparse(dim*m_tot, col);

            k = comp;
            switch type

            % Dirichlet boundary condition
            case {'D','d','dirichlet','Dirichlet'}

                % Loop over components that Dirichlet penalties end up on
                for i = 1:dim
                    C = transpose(T{k,i});
                    A = -tuning*e*transpose(alpha{i,k});
                    B = A + e*C;
                    closure = closure + E{i}*RHOi*Hi*B'*e*H_gamma*(e'*E{k}' );
                    penalty = penalty - E{i}*RHOi*Hi*B'*e*H_gamma;
                end

            % Free boundary condition
            case {'F','f','Free','free','traction','Traction','t','T'}
                    closure = closure - E{k}*RHOi*Hi*e*H_gamma*tau';
                    penalty = penalty + E{k}*RHOi*Hi*e*H_gamma;

            % Unknown boundary condition
            otherwise
                error('No such boundary condition: type = %s',type);
            end
        end

        % type     Struct that specifies the interface coupling.
        %          Fields:
        %          -- tuning:           penalty strength, defaults to 1.2
        %          -- interpolation:    type of interpolation, default 'none'
        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary,type)

            defaultType.tuning = 1.2;
            defaultType.interpolation = 'none';
            default_struct('type', defaultType);

            switch type.interpolation
            case {'none', ''}
                [closure, penalty] = interfaceStandard(obj,boundary,neighbour_scheme,neighbour_boundary,type);
            case {'op','OP'}
                [closure, penalty] = interfaceNonConforming(obj,boundary,neighbour_scheme,neighbour_boundary,type);
            otherwise
                error('Unknown type of interpolation: %s ', type.interpolation);
            end
        end

        function [closure, penalty] = interfaceStandard(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            tuning = type.tuning;

            % u denotes the solution in the own domain
            % v denotes the solution in the neighbour domain
            % Operators without subscripts are from the own domain.

            % Get boundary operators
            e   = obj.getBoundaryOperator('e', boundary);
            tau = obj.getBoundaryOperator('tau', boundary);

            e_v   = neighbour_scheme.getBoundaryOperator('e', neighbour_boundary);
            tau_v = neighbour_scheme.getBoundaryOperator('tau', neighbour_boundary);

            H_gamma = obj.getBoundaryQuadrature(boundary);

            % Operators and quantities that correspond to the own domain only
            Hi = obj.Hi_kron;
            RHOi = obj.RHOi_kron;

            % Penalty strength operators
            alpha_u = 1/4*tuning*obj.getBoundaryOperator('alpha', boundary);
            alpha_v = 1/4*tuning*neighbour_scheme.getBoundaryOperator('alpha', neighbour_boundary);

            closure = -RHOi*Hi*e*H_gamma*(alpha_u' + alpha_v'*e_v*e');
            penalty = RHOi*Hi*e*H_gamma*(alpha_u'*e*e_v' + alpha_v');

            closure = closure - 1/2*RHOi*Hi*e*H_gamma*tau';
            penalty = penalty - 1/2*RHOi*Hi*e*H_gamma*tau_v';

            closure = closure + 1/2*RHOi*Hi*tau*H_gamma*e';
            penalty = penalty - 1/2*RHOi*Hi*tau*H_gamma*e_v';

        end

        function [closure, penalty] = interfaceNonConforming(obj,boundary,neighbour_scheme,neighbour_boundary,type)
            error('Non-conforming interfaces not implemented yet.');
        end

        % Returns the component number that is the tangential/normal component
        % at the specified boundary
        function comp = getComponent(obj, comp_str, boundary)
            assertIsMember(comp_str, {'normal', 'tangential'});
            assertIsMember(boundary, {'w', 'e', 's', 'n'});

            switch boundary
            case {'w', 'e'}
                switch comp_str
                case 'normal'
                    comp = 1;
                case 'tangential'
                    comp = 2;
                end
            case {'s', 'n'}
                switch comp_str
                case 'normal'
                    comp = 2;
                case 'tangential'
                    comp = 1;
                end
            end
        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op -- string
        function o = getBoundaryOperator(obj, op, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})
            assertIsMember(op, {'e', 'e1', 'e2', 'tau', 'tau1', 'tau2', 'alpha', 'alpha1', 'alpha2'})

            switch op

                case {'e', 'e1', 'e2', 'tau', 'tau1', 'tau2'}
                    o = obj.([op, '_', boundary]);

                % Yields vector-valued penalty strength given displacement BC on all components
                case 'alpha'
                    e               = obj.getBoundaryOperator('e', boundary);
                    e_scalar        = obj.getBoundaryOperatorForScalarField('e', boundary);
                    alpha_scalar    = obj.getBoundaryOperatorForScalarField('alpha', boundary);
                    E = obj.E;
                    [m, n] = size(alpha_scalar{1,1});
                    alpha = sparse(m*obj.dim, n*obj.dim);
                    for i = 1:obj.dim
                        for l = 1:obj.dim
                            alpha = alpha + (e'*E{i}*e_scalar*alpha_scalar{i,l}'*E{l}')';
                        end
                    end
                    o = alpha;

                % Yields penalty strength for component 1 given displacement BC on all components
                case 'alpha1'
                    alpha   = obj.getBoundaryOperator('alpha', boundary);
                    e       = obj.getBoundaryOperator('e', boundary);
                    e1      = obj.getBoundaryOperator('e1', boundary);

                    alpha1 = (e1'*e*alpha')';
                    o = alpha1;

                % Yields penalty strength for component 2 given displacement BC on all components
                case 'alpha2'
                    alpha   = obj.getBoundaryOperator('alpha', boundary);
                    e       = obj.getBoundaryOperator('e', boundary);
                    e2      = obj.getBoundaryOperator('e2', boundary);

                    alpha2 = (e2'*e*alpha')';
                    o = alpha2;
            end

        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op -- string
        function o = getBoundaryOperatorForScalarField(obj, op, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})
            assertIsMember(op, {'e', 'alpha'})

            switch op

                case 'e'
                    o = obj.(['e_scalar', '_', boundary]);

                case 'alpha'

                    % alpha{i,j} is the penalty strength on component i due to
                    % displacement BC for component j.
                    e = obj.getBoundaryOperatorForScalarField('e', boundary);

                    LAMBDA = obj.LAMBDA;
                    MU = obj.MU;
                    dim = obj.dim;

                    switch boundary
                        case {'w', 'e'}
                            k = 1;
                        case {'s', 'n'}
                            k = 2;
                    end

                    theta_R = obj.theta_R{k};
                    theta_H = obj.theta_H{k};
                    theta_M = obj.theta_M{k};

                    a_lambda = dim/theta_H + 1/theta_R;
                    a_mu_i = 2/theta_M;
                    a_mu_ij = 2/theta_H + 1/theta_R;

                    d = @kroneckerDelta;  % Kronecker delta
                    db = @(i,j) 1-d(i,j); % Logical not of Kronecker delta

                    alpha_func = @(i,j) d(i,j)* a_lambda*LAMBDA ...
                                        + d(i,j)* a_mu_i*MU ...
                                        + db(i,j)*a_mu_ij*MU;

                    alpha = cell(obj.dim, obj.dim);
                    for i = 1:obj.dim
                        for j = 1:obj.dim
                            alpha{i,j} = d(i,j)*alpha_func(i,k)*e;
                        end
                    end
                    o = alpha;
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
