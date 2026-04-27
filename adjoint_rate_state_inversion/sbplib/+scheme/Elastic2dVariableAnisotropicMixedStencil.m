classdef Elastic2dVariableAnisotropicMixedStencil < scheme.Scheme

% Discretizes the elastic wave equation:
% rho u_{i,tt} = dj C_{ijkl} dk u_j
% opSet should be cell array of opSets, one per dimension. This
% is useful if we have periodic BC in one direction.
% Assumes fully compatible operators

    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid
        dim

        order % Order of accuracy for the approximation

        % Diagonal matrices for variable coefficients
        RHO, RHOi, RHOi_kron % Density
        C, C_D1, C_D2           % Elastic stiffness tensor, C = C_D1 + C_D2.

        D  % Total operator
        D1 % First derivatives
        % D2 % Second derivatives

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
        h11 % First entry in norm matrix

    end

    methods

        % Uses D1*D1 for the C_D1 part of the stiffness tensor C
        % Uses narrow D2 whenever possible for the C_D2 part of C
        % The coefficients can either be function handles or grid functions
        function obj = Elastic2dVariableAnisotropicMixedStencil(g, order, rho, C_D1, C_D2, opSet)
            default_arg('rho', @(x,y) 0*x+1);
            default_arg('opSet',{@sbp.D2VariableCompatible, @sbp.D2VariableCompatible});
            default_arg('optFlag', false);
            dim = 2;

            C_D1_default = cell(dim,dim,dim,dim);
            C_D2_default = cell(dim,dim,dim,dim);
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            C_D1_default{i,j,k,l} = @(x,y) 0*x + 0;
                            C_D2_default{i,j,k,l} = @(x,y) 0*x + 1;
                        end
                    end
                end
            end
            default_arg('C_D1', C_D1_default);
            default_arg('C_D2', C_D2_default);
            assert(isa(g, 'grid.Cartesian'))

            if isa(rho, 'function_handle')
                rho = grid.evalOn(g, rho);
            end

            C_mat = cell(dim,dim,dim,dim);
            C_D1_mat = cell(dim,dim,dim,dim);
            C_D2_mat = cell(dim,dim,dim,dim);
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            if isa(C_D1{i,j,k,l}, 'function_handle')
                                C_D1{i,j,k,l} = grid.evalOn(g, C_D1{i,j,k,l});
                            end
                            if isa(C_D2{i,j,k,l}, 'function_handle')
                                C_D2{i,j,k,l} = grid.evalOn(g, C_D2{i,j,k,l});
                            end
                            C_D1_mat{i,j,k,l} = spdiag(C_D1{i,j,k,l});
                            C_D2_mat{i,j,k,l} = spdiag(C_D2{i,j,k,l});
                            C_mat{i,j,k,l} = C_D1_mat{i,j,k,l} + C_D2_mat{i,j,k,l};
                        end
                    end
                end
            end
            obj.C = C_mat;
            obj.C_D1 = C_D1_mat;
            obj.C_D2 = C_D2_mat;

            m = g.size();
            m_tot = g.N();
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
            h = zeros(dim,1);
            for i = 1:dim
                ops{i} = opSet{i}(m(i), lim{i}, order);
                h(i) = ops{i}.h;
            end

            % Borrowing constants
            for i = 1:dim
                obj.h11{i} = h(i)*ops{i}.borrowing.H11;
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
            I_dim = speye(dim, dim);
            RHO = spdiag(rho);
            obj.RHO = RHO;
            obj.RHOi = inv(RHO);
            obj.RHOi_kron = kron(obj.RHOi, I_dim);

            obj.D1 = cell(dim,1);
            D2_temp = cell(dim,dim,dim);

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

            e_w = kron(e_scalar_w, I_dim);
            e_e = kron(e_scalar_e, I_dim);
            e_s = kron(e_scalar_s, I_dim);
            e_n = kron(e_scalar_n, I_dim);

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
            switch order
            case 2
                width = 3;
            case 4
                width = 5;
            case 6
                width = 7;
            end
            for j = 1:dim
                for k = 1:dim
                    for l = 1:dim
                        D2_temp{j,k,l} = spalloc(m_tot, m_tot, width*m_tot);
                    end
                end
            end
            ind = grid.funcToMatrix(g, 1:m_tot);

            k = 1;
            for r = 1:m(2)
                p = ind(:,r);
                for j = 1:dim
                    for l = 1:dim
                        coeff = C_D2{k,j,k,l};
                        D_kk = D2{1}(coeff(p));
                        D2_temp{j,k,l}(p,p) = D_kk;
                    end
                end
            end

            k = 2;
            for r = 1:m(1)
                p = ind(r,:);
                for j = 1:dim
                    for l = 1:dim
                        coeff = C_D2{k,j,k,l};
                        D_kk = D2{2}(coeff(p));
                        D2_temp{j,k,l}(p,p) = D_kk;
                    end
                end
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
            D1 = obj.D1;
            D = sparse(dim*m_tot,dim*m_tot);
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            if i == k
                                D = D + E{j}*D2_temp{j,k,l}*E{l}';
                                D2_temp{j,k,l} = [];
                            else
                                D = D + E{j}*D1{i}*C_D2_mat{i,j,k,l}*D1{k}*E{l}';
                            end
                            D = D + E{j}*D1{i}*C_D1_mat{i,j,k,l}*D1{k}*E{l}';
                        end
                    end
                end
            end
            clear D2_temp;
            D = obj.RHOi_kron*D;
            obj.D = D;
            clear D;
            %=========================================%'

            % Numerical traction operators for BC.
            %
            % Formula at boundary j: % tau^{j}_i = sum_l T^{j}_{il} u_l
            %
            T_l = cell(dim,1);
            T_r = cell(dim,1);
            tau_l = cell(dim,1);
            tau_r = cell(dim,1);

            D1 = obj.D1;

            % Boundary j
            for j = 1:dim
                T_l{j} = cell(dim,dim);
                T_r{j} = cell(dim,dim);
                tau_l{j} = cell(dim,1);
                tau_r{j} = cell(dim,1);

                [~, n_l] = size(e_l{j});
                [~, n_r] = size(e_r{j});

                % Traction component i
                for i = 1:dim
                    tau_l{j}{i} = sparse(dim*m_tot, n_l);
                    tau_r{j}{i} = sparse(dim*m_tot, n_r);

                    % Displacement component l
                    for l = 1:dim
                        T_l{j}{i,l} = sparse(m_tot, n_l);
                        T_r{j}{i,l} = sparse(m_tot, n_r);

                        % Derivative direction k
                        for k = 1:dim
                            T_l{j}{i,l} = T_l{j}{i,l} ...
                                        - (e_l{j}'*C_mat{j,i,k,l}*D1{k})';
                            T_r{j}{i,l} = T_r{j}{i,l} ...
                                        + (e_r{j}'*C_mat{j,i,k,l}*D1{k})';
                        end
                        tau_l{j}{i} = tau_l{j}{i} + (T_l{j}{i,l}'*E{l}')';
                        tau_r{j}{i} = tau_r{j}{i} + (T_r{j}{i,l}'*E{l}')';
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
            obj.Hi_kron = kron(obj.Hi, I_dim);

            % Misc.
            obj.m = m;
            obj.h = h;
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
            comp = bc{1};
            type = bc{2};
            if ischar(comp)
                comp = obj.getComponent(comp, boundary);
            end

            e       = obj.getBoundaryOperatorForScalarField('e', boundary);
            tau     = obj.getBoundaryOperator(['tau' num2str(comp)], boundary);
            T       = obj.getBoundaryTractionOperator(boundary);
            h11     = obj.getBorrowing(boundary);
            H_gamma = obj.getBoundaryQuadratureForScalarField(boundary);
            nu      = obj.getNormal(boundary);

            E = obj.E;
            Hi = obj.Hi;
            RHOi = obj.RHOi;
            C = obj.C;

            dim = obj.dim;
            m_tot = obj.grid.N();

            % Preallocate
            [~, col] = size(tau);
            closure = sparse(dim*m_tot, dim*m_tot);
            penalty = sparse(dim*m_tot, col);

            j = comp;
            switch type

            % Dirichlet boundary condition
            case {'D','d','dirichlet','Dirichlet','displacement','Displacement'}

                if numel(bc) >= 3
                    dComps = bc{3};
                else
                    dComps = 1:dim;
                end

                % Loops over components that Dirichlet penalties end up on
                % Y: symmetrizing part of penalty
                % Z: symmetric part of penalty
                % X = Y + Z.

                % Nonsymmetric part goes on all components to
                % yield traction in discrete energy rate
                for i = 1:dim
                    Y = T{j,i}';
                    X = e*Y;
                    closure = closure + E{i}*RHOi*Hi*X'*e*H_gamma*(e'*E{j}' );
                    penalty = penalty - E{i}*RHOi*Hi*X'*e*H_gamma;
                end

                % Symmetric part only required on components with displacement BC.
                % (Otherwise it's not symmetric.)
                for i = dComps
                    Z = sparse(m_tot, m_tot);
                    for l = 1:dim
                        for k = 1:dim
                            Z = Z + nu(l)*C{l,i,k,j}*nu(k);
                        end
                    end
                    Z = -tuning*dim/h11*Z;
                    X = Z;
                    closure = closure + E{i}*RHOi*Hi*X'*e*H_gamma*(e'*E{j}' );
                    penalty = penalty - E{i}*RHOi*Hi*X'*e*H_gamma;
                end

            % Free boundary condition
            case {'F','f','Free','free','traction','Traction','t','T'}
                    closure = closure - E{j}*RHOi*Hi*e*H_gamma*tau';
                    penalty = penalty + E{j}*RHOi*Hi*e*H_gamma;

            % Unknown boundary condition
            otherwise
                error('No such boundary condition: type = %s',type);
            end
        end

        % type     Struct that specifies the interface coupling.
        %          Fields:
        %          -- tuning:           penalty strength, defaults to 1.0
        %          -- interpolation:    type of interpolation, default 'none'
        function [closure, penalty] = interface(obj,boundary,neighbour_scheme,neighbour_boundary,type)

            defaultType.tuning = 1.0;
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

            u = obj;
            v = neighbour_scheme;

            % Operators, u side
            e_u       = u.getBoundaryOperatorForScalarField('e', boundary);
            tau_u     = u.getBoundaryOperator('tau', boundary);
            h11_u     = u.getBorrowing(boundary);
            nu_u      = u.getNormal(boundary);

            E_u = u.E;
            C_u = u.C;
            m_tot_u = u.grid.N();

            % Operators, v side
            e_v       = v.getBoundaryOperatorForScalarField('e', neighbour_boundary);
            tau_v     = v.getBoundaryOperator('tau', neighbour_boundary);
            h11_v     = v.getBorrowing(neighbour_boundary);
            nu_v      = v.getNormal(neighbour_boundary);

            E_v = v.E;
            C_v = v.C;
            m_tot_v = v.grid.N();

            % Fix {'e', 's'}, {'w', 'n'}, and {'x','x'} couplings
            flipFlag = false;
            e_v_flip = e_v;
            if (strcmp(boundary,'s') && strcmp(neighbour_boundary,'e')) || ...
               (strcmp(boundary,'e') && strcmp(neighbour_boundary,'s')) || ...
               (strcmp(boundary,'w') && strcmp(neighbour_boundary,'n')) || ...
               (strcmp(boundary,'n') && strcmp(neighbour_boundary,'w')) || ...
               (strcmp(boundary,'s') && strcmp(neighbour_boundary,'s')) || ...
               (strcmp(boundary,'n') && strcmp(neighbour_boundary,'n')) || ...
               (strcmp(boundary,'w') && strcmp(neighbour_boundary,'w')) || ...
               (strcmp(boundary,'e') && strcmp(neighbour_boundary,'e'))

                flipFlag = true;
                e_v_flip = fliplr(e_v);

                t1 = tau_v(:,1:2:end-1);
                t2 = tau_v(:,2:2:end);

                t1 = fliplr(t1);
                t2 = fliplr(t2);

                tau_v(:,1:2:end-1) = t1;
                tau_v(:,2:2:end) = t2;
            end

            % Operators that are only required for own domain
            Hi      = u.Hi_kron;
            RHOi    = u.RHOi_kron;
            e_kron  = u.getBoundaryOperator('e', boundary);
            T_u     = u.getBoundaryTractionOperator(boundary);

            % Shared operators
            H_gamma         = u.getBoundaryQuadratureForScalarField(boundary);
            H_gamma_kron    = u.getBoundaryQuadrature(boundary);
            dim             = u.dim;

            % Preallocate
            [~, m_int] = size(H_gamma);
            closure = sparse(dim*m_tot_u, dim*m_tot_u);
            penalty = sparse(dim*m_tot_u, dim*m_tot_v);

            % ---- Continuity of displacement ------

            % Y: symmetrizing part of penalty
            % Z: symmetric part of penalty
            % X = Y + Z.

            % Loop over components to couple across interface
            for j = 1:dim

                % Loop over components that penalties end up on
                for i = 1:dim
                    Y = 1/2*T_u{j,i}';
                    Z_u = sparse(m_int, m_int);
                    Z_v = sparse(m_int, m_int);
                    for l = 1:dim
                        for k = 1:dim
                            Z_u = Z_u + e_u'*nu_u(l)*C_u{l,i,k,j}*nu_u(k)*e_u;
                            Z_v = Z_v + e_v'*nu_v(l)*C_v{l,i,k,j}*nu_v(k)*e_v;
                        end
                    end

                    if flipFlag
                        Z_v = rot90(Z_v,2);
                    end

                    Z = -tuning*dim*( 1/(4*h11_u)*Z_u + 1/(4*h11_v)*Z_v );
                    X = Y + Z*e_u';
                    closure = closure + E_u{i}*X'*H_gamma*e_u'*E_u{j}';
                    penalty = penalty - E_u{i}*X'*H_gamma*e_v_flip'*E_v{j}';

                end
            end

            % ---- Continuity of traction ------
            closure = closure - 1/2*e_kron*H_gamma_kron*tau_u';
            penalty = penalty - 1/2*e_kron*H_gamma_kron*tau_v';

            % ---- Multiply by inverse of density x quadraure ----
            closure = RHOi*Hi*closure;
            penalty = RHOi*Hi*penalty;

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

        % Returns h11 for the boundary specified by the string boundary.
        % op -- string
        function h11 = getBorrowing(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            switch boundary
            case {'w','e'}
                h11 = obj.h11{1};
            case {'s', 'n'}
                h11 = obj.h11{2};
            end
        end

        % Returns the outward unit normal vector for the boundary specified by the string boundary.
        function nu = getNormal(obj, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            switch boundary
            case 'w'
                nu = [-1,0];
            case 'e'
                nu = [1,0];
            case 's'
                nu = [0,-1];
            case 'n'
                nu = [0,1];
            end
        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op -- string
        function o = getBoundaryOperator(obj, op, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})
            assertIsMember(op, {'e', 'e1', 'e2', 'tau', 'tau1', 'tau2'})

            switch op
                case {'e', 'e1', 'e2', 'tau', 'tau1', 'tau2'}
                    o = obj.([op, '_', boundary]);
            end

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
