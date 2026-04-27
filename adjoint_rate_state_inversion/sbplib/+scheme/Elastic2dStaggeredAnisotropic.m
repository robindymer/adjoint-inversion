classdef Elastic2dStaggeredAnisotropic < scheme.Scheme

% Discretizes the elastic wave equation:
% rho u_{i,tt} = dj C_{ijkl} dk u_j
% Uses a staggered Lebedev grid
% The solution (displacement) is stored on g_u
% Stresses (and hance tractions) appear on g_s
% Density is evaluated on g_u
% The stiffness tensor is evaluated on g_s

    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid
        dim
        nGrids
        N       % Total number of unknowns stored (2 displacement components on 2 grids)

        order % Order of accuracy for the approximation

        % Diagonal matrices for variable coefficients
        RHO  % Density
        C    % Elastic stiffness tensor

        D  % Total operator

        % Boundary operators in cell format, used for BC
        % T_w, T_e, T_s, T_n

        % Traction operators
        tau_w, tau_e, tau_s, tau_n      % Return scalar field
        T_w, T_e, T_s, T_n              % Act on scalar, return scalar

        % Inner products
        H, H_u, H_s
        % , Hi, Hi_kron, H_1D

        % Boundary inner products (for scalar field)
        H_w_u, H_e_u, H_s_u, H_n_u
        H_w_s, H_e_s, H_s_s, H_n_s

        % Boundary restriction operators
        e_w_u, e_e_u, e_s_u, e_n_u      % Act on scalar field, return scalar field at boundary
        e_w_s, e_e_s, e_s_s, e_n_s      % Act on scalar field, return scalar field at boundary

        % U{i}^T picks out component i
        U

        % G{i}^T picks out displacement grid i
        G

        % Borrowing constants of the form gamma*h, where gamma is a dimensionless constant.
        h11 % First entry in norm matrix

    end

    methods

        % The coefficients can either be function handles or grid functions
        function obj = Elastic2dStaggeredAnisotropic(g, order, rho, C)
            default_arg('rho', @(x,y) 0*x+1);
            dim = 2;
            nGrids = 2;

            C_default = cell(dim,dim,dim,dim);
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            C_default{i,j,k,l} = @(x,y) 0*x + 1;
                        end
                    end
                end
            end
            default_arg('C', C_default);
            assert(isa(g, 'grid.Staggered'))

            g_u = g.gridGroups{1};
            g_s = g.gridGroups{2};

            m_u = {g_u{1}.size(), g_u{2}.size()};
            m_s = {g_s{1}.size(), g_s{2}.size()};

            if isa(rho, 'function_handle')
                rho_vec = cell(nGrids, 1);
                for i = 1:nGrids
                    rho_vec{i} = grid.evalOn(g_u{i}, rho);
                end
                rho = rho_vec;
            end
            for i = 1:nGrids
                RHO{i} = spdiag(rho{i});
            end
            obj.RHO = RHO;

            C_mat = cell(nGrids, 1);
            for a = 1:nGrids
                C_mat{a} = cell(dim,dim,dim,dim);
            end
            for a = 1:nGrids
                for i = 1:dim
                    for j = 1:dim
                        for k = 1:dim
                            for l = 1:dim
                                if numel(C) == nGrids
                                    C_mat{a}{i,j,k,l} = spdiag(C{a}{i,j,k,l});
                                else
                                    C_mat{a}{i,j,k,l} = spdiag(grid.evalOn(g_s{a}, C{i,j,k,l}));
                                end
                            end
                        end
                    end
                end
            end
            C = C_mat;
            obj.C = C;

            % Reference m for primal grid
            m = g_u{1}.size();
            X = g_u{1}.points;
            lim = cell(dim, 1);
            for i = 1:dim
                lim{i} = {min(X(:,i)), max(X(:,i))};
            end

            % 1D operators
            ops = cell(dim,1);
            D1p = cell(dim, 1);
            D1d = cell(dim, 1);
            mp = cell(dim, 1);
            md = cell(dim, 1);
            Ip = cell(dim, 1);
            Id = cell(dim, 1);
            Hp = cell(dim, 1);
            Hd = cell(dim, 1);

            opSet = @sbp.D1StaggeredUpwind;
            for i = 1:dim
                ops{i} = opSet(m(i), lim{i}, order);
                D1p{i} = ops{i}.D1_dual;
                D1d{i} = ops{i}.D1_primal;
                mp{i} = length(ops{i}.x_primal);
                md{i} = length(ops{i}.x_dual);
                Ip{i} = speye(mp{i}, mp{i});
                Id{i} = speye(md{i}, md{i});
                Hp{i} = ops{i}.H_primal;
                Hd{i} = ops{i}.H_dual;
                ep_l{i} = ops{i}.e_primal_l;
                ep_r{i} = ops{i}.e_primal_r;
                ed_l{i} = ops{i}.e_dual_l;
                ed_r{i} = ops{i}.e_dual_r;
            end

            % Borrowing constants
            % for i = 1:dim
            %     obj.h11{i} = ops{i}.H_dual(1,1);
            % end
            obj.h11{1}{1} = ops{1}.H_dual(1,1);
            obj.h11{1}{2} = ops{2}.H_primal(1,1);
            obj.h11{2}{1} = ops{1}.H_primal(1,1);
            obj.h11{2}{2} = ops{2}.H_dual(1,1);

            %---- Grid layout -------
            % gu1 = xp o yp;
            % gu2 = xd o yd;
            % gs1 = xd o yp;
            % gs2 = xp o yd;
            %------------------------

            % Quadratures
            obj.H_u = cell(nGrids, 1);
            obj.H_s = cell(nGrids, 1);
            obj.H_u{1} = kron(Hp{1}, Hp{2});
            obj.H_u{2} = kron(Hd{1}, Hd{2});
            obj.H_s{1} = kron(Hd{1}, Hp{2});
            obj.H_s{2} = kron(Hp{1}, Hd{2});

            obj.H_w_s = cell(nGrids, 1);
            obj.H_e_s = cell(nGrids, 1);
            obj.H_s_s = cell(nGrids, 1);
            obj.H_n_s = cell(nGrids, 1);

            obj.H_w_s{1} = Hp{2};
            obj.H_w_s{2} = Hd{2};
            obj.H_e_s{1} = Hp{2};
            obj.H_e_s{2} = Hd{2};

            obj.H_s_s{1} = Hd{1};
            obj.H_s_s{2} = Hp{1};
            obj.H_n_s{1} = Hd{1};
            obj.H_n_s{2} = Hp{1};

            % Boundary restriction ops
            e_w_u = cell(nGrids, 1);
            e_s_u = cell(nGrids, 1);
            e_e_u = cell(nGrids, 1);
            e_n_u = cell(nGrids, 1);

            e_w_s = cell(nGrids, 1);
            e_s_s = cell(nGrids, 1);
            e_e_s = cell(nGrids, 1);
            e_n_s = cell(nGrids, 1);

            e_w_u{1} = kron(ep_l{1}, Ip{2});
            e_e_u{1} = kron(ep_r{1}, Ip{2});
            e_s_u{1} = kron(Ip{1}, ep_l{2});
            e_n_u{1} = kron(Ip{1}, ep_r{2});

            e_w_u{2} = kron(ed_l{1}, Id{2});
            e_e_u{2} = kron(ed_r{1}, Id{2});
            e_s_u{2} = kron(Id{1}, ed_l{2});
            e_n_u{2} = kron(Id{1}, ed_r{2});

            e_w_s{1} = kron(ed_l{1}, Ip{2});
            e_e_s{1} = kron(ed_r{1}, Ip{2});
            e_s_s{1} = kron(Id{1}, ep_l{2});
            e_n_s{1} = kron(Id{1}, ep_r{2});

            e_w_s{2} = kron(ep_l{1}, Id{2});
            e_e_s{2} = kron(ep_r{1}, Id{2});
            e_s_s{2} = kron(Ip{1}, ed_l{2});
            e_n_s{2} = kron(Ip{1}, ed_r{2});

            obj.e_w_u = e_w_u;
            obj.e_e_u = e_e_u;
            obj.e_s_u = e_s_u;
            obj.e_n_u = e_n_u;

            obj.e_w_s = e_w_s;
            obj.e_e_s = e_e_s;
            obj.e_s_s = e_s_s;
            obj.e_n_s = e_n_s;


            % D1_u2s{a, b}{i} approximates ddi and
            % takes from u grid number b to s grid number a
            % Some of D1_x2y{a, b} are 0.
            D1_u2s = cell(nGrids, nGrids);
            D1_s2u = cell(nGrids, nGrids);

            N_u = cell(nGrids, 1);
            N_s = cell(nGrids, 1);
            for a = 1:nGrids
                N_u{a} = g_u{a}.N();
                N_s{a} = g_s{a}.N();
            end

            %---- Grid layout -------
            % gu1 = xp o yp;
            % gu2 = xd o yd;
            % gs1 = xd o yp;
            % gs2 = xp o yd;
            %------------------------

            D1_u2s{1,1}{1} = kron(D1p{1}, Ip{2});
            D1_s2u{1,1}{1} = kron(D1d{1}, Ip{2});

            D1_u2s{1,2}{2} = kron(Id{1}, D1d{2});
            D1_u2s{2,1}{2} = kron(Ip{1}, D1p{2});

            D1_s2u{1,2}{2} = kron(Ip{1}, D1d{2});
            D1_s2u{2,1}{2} = kron(Id{1}, D1p{2});

            D1_u2s{2,2}{1} = kron(D1d{1}, Id{2});
            D1_s2u{2,2}{1} = kron(D1p{1}, Id{2});

            D1_u2s{1,1}{2} = sparse(N_s{1}, N_u{1});
            D1_s2u{1,1}{2} = sparse(N_u{1}, N_s{1});

            D1_u2s{2,2}{2} = sparse(N_s{2}, N_u{2});
            D1_s2u{2,2}{2} = sparse(N_u{2}, N_s{2});

            D1_u2s{1,2}{1} = sparse(N_s{1}, N_u{2});
            D1_s2u{1,2}{1} = sparse(N_u{1}, N_s{2});

            D1_u2s{2,1}{1} = sparse(N_s{2}, N_u{1});
            D1_s2u{2,1}{1} = sparse(N_u{2}, N_s{1});


            %---- Combine grids and components -----

            % U{a}{i}^T picks out u component i on grid a
            U = cell(nGrids, 1);
            for a = 1:2
                U{a} = cell(dim, 1);
                I = speye(N_u{a}, N_u{a});
                for i = 1:dim
                    E = sparse(dim,1);
                    E(i) = 1;
                    U{a}{i} = kron(I, E);
                end
            end
            obj.U = U;

            % Order grids
            % u1, u2
            Iu1 = speye(dim*N_u{1}, dim*N_u{1});
            Iu2 = speye(dim*N_u{2}, dim*N_u{2});

            Gu1 = cell2mat( {Iu1; sparse(dim*N_u{2}, dim*N_u{1})} );
            Gu2 = cell2mat( {sparse(dim*N_u{1}, dim*N_u{2}); Iu2} );

            G = {Gu1; Gu2};
            obj.G = G;

            obj.H = G{1}*(U{1}{1}*obj.H_u{1}*U{1}{1}' + U{1}{2}*obj.H_u{1}*U{1}{2}')*G{1}'...
                  + G{2}*(U{2}{1}*obj.H_u{2}*U{2}{1}' + U{2}{2}*obj.H_u{2}*U{2}{2}')*G{2}';

            % e1_w = (e_scalar_w'*E{1}')';
            % e1_e = (e_scalar_e'*E{1}')';
            % e1_s = (e_scalar_s'*E{1}')';
            % e1_n = (e_scalar_n'*E{1}')';

            % e2_w = (e_scalar_w'*E{2}')';
            % e2_e = (e_scalar_e'*E{2}')';
            % e2_s = (e_scalar_s'*E{2}')';
            % e2_n = (e_scalar_n'*E{2}')';

            stencilWidth = order;
            % Differentiation matrix D (without SAT)
            N = dim*(N_u{1} + N_u{2});
            D = spalloc(N, N, stencilWidth^2*N);
            for a = 1:nGrids
                for b = 1:nGrids
                    for c = 1:nGrids
                        for i = 1:dim
                            for j = 1:dim
                                for k = 1:dim
                                    for l = 1:dim
                                        D = D + (G{a}*U{a}{j})*(RHO{a}\(D1_s2u{a,b}{i}*C{b}{i,j,k,l}*D1_u2s{b,c}{k}*U{c}{l}'*G{c}'));
                                    end
                                end
                            end
                        end
                    end
                end
            end
            obj.D = D;
            clear D;
            obj.N = N;
            %=========================================%'

            % Numerical traction operators for BC.
            %
            % Formula at boundary j: % tau^{j}_i = sum_l T^{j}_{il} u_l
            %

            n_w = obj.getNormal('w');
            n_e = obj.getNormal('e');
            n_s = obj.getNormal('s');
            n_n = obj.getNormal('n');

            tau_w = cell(nGrids, 1);
            tau_e = cell(nGrids, 1);
            tau_s = cell(nGrids, 1);
            tau_n = cell(nGrids, 1);

            T_w = cell(nGrids, nGrids);
            T_e = cell(nGrids, nGrids);
            T_s = cell(nGrids, nGrids);
            T_n = cell(nGrids, nGrids);
            for b = 1:nGrids
                [~, m_we] = size(e_w_s{b});
                [~, m_sn] = size(e_s_s{b});
                for c = 1:nGrids
                    T_w{b,c} = cell(dim, dim);
                    T_e{b,c} = cell(dim, dim);
                    T_s{b,c} = cell(dim, dim);
                    T_n{b,c} = cell(dim, dim);

                    for i = 1:dim
                        for j = 1:dim
                            T_w{b,c}{i,j} = sparse(N_u{c}, m_we);
                            T_e{b,c}{i,j} = sparse(N_u{c}, m_we);
                            T_s{b,c}{i,j} = sparse(N_u{c}, m_sn);
                            T_n{b,c}{i,j} = sparse(N_u{c}, m_sn);
                        end
                    end
                end
            end

            for b = 1:nGrids
                tau_w{b} = cell(dim, 1);
                tau_e{b} = cell(dim, 1);
                tau_s{b} = cell(dim, 1);
                tau_n{b} = cell(dim, 1);

                for j = 1:dim
                    tau_w{b}{j} = sparse(N, m_s{b}(2));
                    tau_e{b}{j} = sparse(N, m_s{b}(2));
                    tau_s{b}{j} = sparse(N, m_s{b}(1));
                    tau_n{b}{j} = sparse(N, m_s{b}(1));
                end

                for c = 1:nGrids
                    for i = 1:dim
                        for j = 1:dim
                            for k = 1:dim
                                for l = 1:dim
                                    sigma_b_ij = C{b}{i,j,k,l}*D1_u2s{b,c}{k}*U{c}{l}'*G{c}';

                                    tau_w{b}{j} = tau_w{b}{j} + (e_w_s{b}'*n_w(i)*sigma_b_ij)';
                                    tau_e{b}{j} = tau_e{b}{j} + (e_e_s{b}'*n_e(i)*sigma_b_ij)';
                                    tau_s{b}{j} = tau_s{b}{j} + (e_s_s{b}'*n_s(i)*sigma_b_ij)';
                                    tau_n{b}{j} = tau_n{b}{j} + (e_n_s{b}'*n_n(i)*sigma_b_ij)';

                                    S_bc_ijl = C{b}{i,j,k,l}*D1_u2s{b,c}{k};

                                    T_w{b,c}{j,l} = T_w{b,c}{j,l} + (e_w_s{b}'*n_w(i)*S_bc_ijl)';
                                    T_e{b,c}{j,l} = T_e{b,c}{j,l} + (e_e_s{b}'*n_e(i)*S_bc_ijl)';
                                    T_s{b,c}{j,l} = T_s{b,c}{j,l} + (e_s_s{b}'*n_s(i)*S_bc_ijl)';
                                    T_n{b,c}{j,l} = T_n{b,c}{j,l} + (e_n_s{b}'*n_n(i)*S_bc_ijl)';
                                end
                            end
                        end
                    end
                end
            end

            obj.tau_w = tau_w;
            obj.tau_e = tau_e;
            obj.tau_s = tau_s;
            obj.tau_n = tau_n;

            obj.T_w = T_w;
            obj.T_e = T_e;
            obj.T_s = T_s;
            obj.T_n = T_n;

            % D1 = obj.D1;

            % Traction tensors, T_ij
            % obj.T_w = T_l{1};
            % obj.T_e = T_r{1};
            % obj.T_s = T_l{2};
            % obj.T_n = T_r{2};

            % Restriction operators
            % obj.e_w = e_w;
            % obj.e_e = e_e;
            % obj.e_s = e_s;
            % obj.e_n = e_n;

            % obj.e1_w = e1_w;
            % obj.e1_e = e1_e;
            % obj.e1_s = e1_s;
            % obj.e1_n = e1_n;

            % obj.e2_w = e2_w;
            % obj.e2_e = e2_e;
            % obj.e2_s = e2_s;
            % obj.e2_n = e2_n;

            % obj.e_scalar_w = e_scalar_w;
            % obj.e_scalar_e = e_scalar_e;
            % obj.e_scalar_s = e_scalar_s;
            % obj.e_scalar_n = e_scalar_n;

            % % First component of traction
            % obj.tau1_w = tau_l{1}{1};
            % obj.tau1_e = tau_r{1}{1};
            % obj.tau1_s = tau_l{2}{1};
            % obj.tau1_n = tau_r{2}{1};

            % % Second component of traction
            % obj.tau2_w = tau_l{1}{2};
            % obj.tau2_e = tau_r{1}{2};
            % obj.tau2_s = tau_l{2}{2};
            % obj.tau2_n = tau_r{2}{2};

            % % Traction vectors
            % obj.tau_w = (e_w'*e1_w*obj.tau1_w')' + (e_w'*e2_w*obj.tau2_w')';
            % obj.tau_e = (e_e'*e1_e*obj.tau1_e')' + (e_e'*e2_e*obj.tau2_e')';
            % obj.tau_s = (e_s'*e1_s*obj.tau1_s')' + (e_s'*e2_s*obj.tau2_s')';
            % obj.tau_n = (e_n'*e1_n*obj.tau1_n')' + (e_n'*e2_n*obj.tau2_n')';

            % Misc.
            obj.m = m;
            obj.h = [];
            obj.order = order;
            obj.grid = g;
            obj.dim = dim;
            obj.nGrids = nGrids;

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

            e_u       = obj.getBoundaryOperatorForScalarField('e_u', boundary);
            e_s       = obj.getBoundaryOperatorForScalarField('e_s', boundary);
            tau     = obj.getBoundaryOperator('tau', boundary);
            T       = obj.getBoundaryTractionOperator(boundary);
            H_gamma = obj.getBoundaryQuadratureForScalarField(boundary);
            nu      = obj.getNormal(boundary);

            U = obj.U;
            G = obj.G;
            H = obj.H_u;
            RHO = obj.RHO;
            C = obj.C;

            %---- Grid layout -------
            % gu1 = xp o yp;
            % gu2 = xd o yd;
            % gs1 = xd o yp;
            % gs2 = xp o yd;
            %------------------------

            switch boundary
                case {'w', 'e'}
                    gridCombos = {{1,1}, {2,2}};
                case {'s', 'n'}
                    gridCombos = {{2,1}, {1,2}};
            end

            dim = obj.dim;
            nGrids = obj.nGrids;

            m_tot = obj.N;

            % Preallocate
            [~, col] = size(tau{1}{1});
            closure = sparse(m_tot, m_tot);
            penalty = cell(1, nGrids);
            for a = 1:nGrids
                [~, col] = size(e_u{a});
                penalty{a} = sparse(m_tot, col);
            end

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
                for c = 1:nGrids
                    for m = 1:numel(gridCombos)
                        gc = gridCombos{m};
                        a = gc{1};
                        b = gc{2};

                        for i = 1:dim
                            Y = T{a,c}{j,i}';
                            closure = closure + G{c}*U{c}{i}*((RHO{c}*H{c})\(Y'*H_gamma{a}*(e_u{b}'*U{b}{j}'*G{b}') ));
                            penalty{b} = penalty{b} - G{c}*U{c}{i}*((RHO{c}*H{c})\(Y'*H_gamma{a}) );
                        end
                    end
                end

                % Symmetric part only required on components with displacement BC.
                % (Otherwise it's not symmetric.)
                for m = 1:numel(gridCombos)
                    gc = gridCombos{m};
                    a = gc{1};
                    b = gc{2};

                    h11 = obj.getBorrowing(b, boundary);

                    for i = dComps
                        Z = 0*C{b}{1,1,1,1};
                        for l = 1:dim
                            for k = 1:dim
                                Z = Z + nu(l)*C{b}{l,i,k,j}*nu(k);
                            end
                        end
                        Z = -tuning*dim/h11*Z;
                        X = e_s{b}'*Z*e_s{b};
                        closure = closure + G{a}*U{a}{i}*((RHO{a}*H{a})\(e_u{a}*X'*H_gamma{b}*(e_u{a}'*U{a}{j}'*G{a}' ) ));
                        penalty{a} = penalty{a} - G{a}*U{a}{i}*((RHO{a}*H{a})\(e_u{a}*X'*H_gamma{b} ));
                    end
                end

            % Free boundary condition
            case {'F','f','Free','free','traction','Traction','t','T'}
                for m = 1:numel(gridCombos)
                    gc = gridCombos{m};
                    a = gc{1};
                    b = gc{2};
                    closure = closure - G{a}*U{a}{j}*(RHO{a}\(H{a}\(e_u{a}*H_gamma{b}*tau{b}{j}')));
                    penalty{b} = G{a}*U{a}{j}*(RHO{a}\(H{a}\(e_u{a}*H_gamma{b})));
                end

            % Unknown boundary condition
            otherwise
                error('No such boundary condition: type = %s',type);
            end

            penalty = cell2mat(penalty);
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
            eu_u       = u.getBoundaryOperatorForScalarField('e_u', boundary);
            es_u       = u.getBoundaryOperatorForScalarField('e_s', boundary);
            tau_u     = u.getBoundaryOperator('tau', boundary);
            nu_u      = u.getNormal(boundary);

            G_u = u.G;
            U_u = u.U;
            C_u = u.C;
            m_tot_u = u.N;

            % Operators, v side
            eu_v       = v.getBoundaryOperatorForScalarField('e_u', neighbour_boundary);
            es_v       = v.getBoundaryOperatorForScalarField('e_s', neighbour_boundary);
            tau_v     = v.getBoundaryOperator('tau', neighbour_boundary);
            nu_v      = v.getNormal(neighbour_boundary);

            G_v = v.G;
            U_v = v.U;
            C_v = v.C;
            m_tot_v = v.N;

            % Operators that are only required for own domain
            % Hi      = u.Hi_kron;
            % RHOi    = u.RHOi_kron;
            % e_kron  = u.getBoundaryOperator('e', boundary);
            H       = u.H_u;
            RHO     = u.RHO;
            T_u     = u.getBoundaryTractionOperator(boundary);

            % Shared operators
            H_gamma         = u.getBoundaryQuadratureForScalarField(boundary);
            % H_gamma_kron    = u.getBoundaryQuadrature(boundary);
            dim             = u.dim;
            nGrids          = obj.nGrids;

            % Preallocate
            % [~, m_int] = size(H_gamma);
            closure = sparse(m_tot_u, m_tot_u);
            penalty = sparse(m_tot_u, m_tot_v);

            %---- Grid layout -------
            % gu1 = xp o yp;
            % gu2 = xd o yd;
            % gs1 = xd o yp;
            % gs2 = xp o yd;
            %------------------------

            switch boundary
                case {'w', 'e'}
                    switch neighbour_boundary
                    case {'w', 'e'}
                        gridCombos = {{1,1,1}, {2,2,2}};
                    case {'s', 'n'}
                        gridCombos = {{1,1,2}, {2,2,1}};
                    end
                case {'s', 'n'}
                    switch neighbour_boundary
                    case {'s', 'n'}
                        gridCombos = {{2,1,1}, {1,2,2}};
                    case {'w', 'e'}
                        gridCombos = {{2,1,2}, {1,2,1}};
                    end
            end

            % Symmetrizing part
            for c = 1:nGrids
                for m = 1:numel(gridCombos)
                    gc = gridCombos{m};
                    a = gc{1};
                    b = gc{2};

                    for i = 1:dim
                        for j = 1:dim
                            Y = 1/2*T_u{a,c}{j,i}';
                            closure = closure + G_u{c}*U_u{c}{i}*((RHO{c}*H{c})\(Y'*H_gamma{a}*(eu_u{b}'*U_u{b}{j}'*G_u{b}') ));
                            penalty = penalty - G_u{c}*U_u{c}{i}*((RHO{c}*H{c})\(Y'*H_gamma{a}*(eu_v{b}'*U_v{b}{j}'*G_v{b}') ));
                        end
                    end
                end
            end

            % Symmetric part
            for m = 1:numel(gridCombos)
                gc = gridCombos{m};
                a = gc{1};
                b = gc{2};
                bv = gc{3};

                h11_u = u.getBorrowing(b, boundary);
                h11_v = v.getBorrowing(bv, neighbour_boundary);

                for i = 1:dim
                    for j = 1:dim
                        Z_u = 0*es_u{b}'*es_u{b};
                        Z_v = 0*es_v{bv}'*es_v{bv};
                        for l = 1:dim
                            for k = 1:dim
                                Z_u = Z_u + es_u{b}'*nu_u(l)*C_u{b}{l,i,k,j}*nu_u(k)*es_u{b};
                                Z_v = Z_v + es_v{bv}'*nu_v(l)*C_v{bv}{l,i,k,j}*nu_v(k)*es_v{bv};
                            end
                        end
                        Z = -tuning*dim*( 1/(4*h11_u)*Z_u + 1/(4*h11_v)*Z_v );
                        % X = es_u{b}'*Z*es_u{b};
                        % X = Z;
                        closure = closure + G_u{a}*U_u{a}{i}*((RHO{a}*H{a})\(eu_u{a}*Z'*H_gamma{b}*(eu_u{a}'*U_u{a}{j}'*G_u{a}' ) ));
                        penalty = penalty - G_u{a}*U_u{a}{i}*((RHO{a}*H{a})\(eu_u{a}*Z'*H_gamma{b}*(eu_v{a}'*U_v{a}{j}'*G_v{a}' ) ));
                    end
                end
            end

            % Continuity of traction
            for j = 1:dim
                for m = 1:numel(gridCombos)
                    gc = gridCombos{m};
                    a = gc{1};
                    b = gc{2};
                    bv = gc{3};
                    closure = closure - 1/2*G_u{a}*U_u{a}{j}*(RHO{a}\(H{a}\(eu_u{a}*H_gamma{b}*tau_u{b}{j}')));
                    penalty = penalty - 1/2*G_u{a}*U_u{a}{j}*(RHO{a}\(H{a}\(eu_u{a}*H_gamma{b}*tau_v{bv}{j}')));
                end
            end

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
        function h11 = getBorrowing(obj, stressGrid, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            switch boundary
            case {'w','e'}
                h11 = obj.h11{stressGrid}{1};
            case {'s', 'n'}
                h11 = obj.h11{stressGrid}{2};
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
            assertIsMember(op, {'e_u', 'e_s'})

            switch op
                case 'e_u'
                    o = obj.(['e_', boundary, '_u']);
                case 'e_s'
                    o = obj.(['e_', boundary, '_s']);
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

            H_b = obj.(['H_', boundary, '_s']);
        end

        function N = size(obj)
            N = length(obj.D);
        end
    end
end
