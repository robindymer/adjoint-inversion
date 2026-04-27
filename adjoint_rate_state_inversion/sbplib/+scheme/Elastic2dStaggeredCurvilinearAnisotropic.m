classdef Elastic2dStaggeredCurvilinearAnisotropic < scheme.Scheme

% Discretizes the elastic wave equation:
% rho u_{i,tt} = dj C_{ijkl} dk u_j
% in curvilinear coordinates, using Lebedev staggered grids

    properties
        m % Number of points in each direction, possibly a vector
        h % Grid spacing

        grid
        dim
        nGrids

        order % Order of accuracy for the approximation

        % Diagonal matrices for variable coefficients
        J, Ji
        RHO % Density
        C   % Elastic stiffness tensor

        D  % Total operator

        % Dx, Dy % Physical derivatives
        n_w, n_e, n_s, n_n % Physical normals

        % Boundary operators in cell format, used for BC
        T_w, T_e, T_s, T_n

        % Traction operators
        % tau_w, tau_e, tau_s, tau_n      % Return vector field
        % tau1_w, tau1_e, tau1_s, tau1_n  % Return scalar field
        % tau2_w, tau2_e, tau2_s, tau2_n  % Return scalar field

        % Inner products
        H

        % Boundary inner products (for scalar field)
        H_w, H_e, H_s, H_n

        % Surface Jacobian vectors
        s_w, s_e, s_s, s_n

        % Boundary restriction operators
        e_w_u, e_e_u, e_s_u, e_n_u      % Act on scalar field, return scalar field at boundary
        e_w_s, e_e_s, e_s_s, e_n_s      % Act on scalar field, return scalar field at boundary
        % e1_w, e1_e, e1_s, e1_n  % Act on vector field, return scalar field at boundary
        % e2_w, e2_e, e2_s, e2_n  % Act on vector field, return scalar field at boundary
        % e_scalar_w, e_scalar_e, e_scalar_s, e_scalar_n; % Act on scalar field, return scalar field
        % en_w, en_e, en_s, en_n  % Act on vector field, return normal component

        % U{i}^T picks out component i
        U

        % G{i}^T picks out displacement grid i
        G

        % Elastic2dVariableAnisotropic object for reference domain
        refObj
    end

    methods

        % The coefficients can either be function handles or grid functions
        % optFlag -- if true, extra computations are performed, which may be helpful for optimization.
        function obj = Elastic2dStaggeredCurvilinearAnisotropic(g, order, rho, C)
            default_arg('rho', @(x,y) 0*x+1);

            opSet = @sbp.D1StaggeredUpwind;
            dim = 2;
            nGrids = 2;

            C_default = cell(dim,dim,dim,dim);
            for i = 1:dim
                for j = 1:dim
                    for k = 1:dim
                        for l = 1:dim
                            C_default{i,j,k,l} = @(x,y) 0*x;
                        end
                    end
                end
            end
            default_arg('C', C_default);
            assert(isa(g, 'grid.Staggered'));

            g_u = g.gridGroups{1};
            g_s = g.gridGroups{2};

            m_u = {g_u{1}.size(), g_u{2}.size()};
            m_s = {g_s{1}.size(), g_s{2}.size()};

            if isa(rho, 'function_handle')
                rho_vec = cell(nGrids, 1);
                for a = 1:nGrids
                    rho_vec{a} = grid.evalOn(g_u{a}, rho);
                end
                rho = rho_vec;
            end
            for a = 1:nGrids
                RHO{a} = spdiag(rho{a});
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
                                if numel(C) == dim
                                    C_mat{a}{i,j,k,l} = spdiag(C{a}{i,j,k,l});
                                else
                                    C_mat{a}{i,j,k,l} = spdiag(grid.evalOn(g_s{a}, C{i,j,k,l}));
                                end
                            end
                        end
                    end
                end
            end
            obj.C = C_mat;

            C = cell(nGrids, 1);
            for a = 1:nGrids
                C{a} = cell(dim,dim,dim,dim);
                for i = 1:dim
                    for j = 1:dim
                        for k = 1:dim
                            for l = 1:dim
                                C{a}{i,j,k,l} = diag(C_mat{a}{i,j,k,l});
                            end
                        end
                    end
                end
            end

            % Reference m for primal grid
            m = g_u{1}.size();

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
                ops{i} = opSet(m(i), {0,1}, order);
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

            % Logical operators
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
            for a = 1:nGrids
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

            %---- Grid layout -------
            % gu1 = xp o yp;
            % gu2 = xd o yd;
            % gs1 = xd o yp;
            % gs2 = xp o yd;
            %------------------------

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

            %  --- Metric coefficients on stress grids -------
            x = cell(nGrids, 1);
            y = cell(nGrids, 1);
            J = cell(nGrids, 1);
            x_xi = cell(nGrids, 1);
            x_eta = cell(nGrids, 1);
            y_xi = cell(nGrids, 1);
            y_eta = cell(nGrids, 1);

            for a = 1:nGrids
                coords = g_u{a}.points();
                x{a} = coords(:,1);
                y{a} = coords(:,2);
            end

            for a = 1:nGrids
                x_xi{a} = zeros(N_s{a}, 1);
                y_xi{a} = zeros(N_s{a}, 1);
                x_eta{a} = zeros(N_s{a}, 1);
                y_eta{a} = zeros(N_s{a}, 1);

                for b = 1:nGrids
                    x_xi{a} = x_xi{a} + D1_u2s{a,b}{1}*x{b};
                    y_xi{a} = y_xi{a} + D1_u2s{a,b}{1}*y{b};
                    x_eta{a} = x_eta{a} + D1_u2s{a,b}{2}*x{b};
                    y_eta{a} = y_eta{a} + D1_u2s{a,b}{2}*y{b};
                end
            end

            for a = 1:nGrids
                J{a} = x_xi{a}.*y_eta{a} - x_eta{a}.*y_xi{a};
            end

            K = cell(nGrids, 1);
            for a = 1:nGrids
                K{a} = cell(dim, dim);
                K{a}{1,1} = y_eta{a}./J{a};
                K{a}{1,2} = -y_xi{a}./J{a};
                K{a}{2,1} = -x_eta{a}./J{a};
                K{a}{2,2} = x_xi{a}./J{a};
            end
            % ----------------------------------------------

            %  --- Metric coefficients on displacement grids -------
            x_s = cell(nGrids, 1);
            y_s = cell(nGrids, 1);
            J_u = cell(nGrids, 1);
            x_xi_u = cell(nGrids, 1);
            x_eta_u = cell(nGrids, 1);
            y_xi_u = cell(nGrids, 1);
            y_eta_u = cell(nGrids, 1);

            for a = 1:nGrids
                coords = g_s{a}.points();
                x_s{a} = coords(:,1);
                y_s{a} = coords(:,2);
            end

            for a = 1:nGrids
                x_xi_u{a} = zeros(N_u{a}, 1);
                y_xi_u{a} = zeros(N_u{a}, 1);
                x_eta_u{a} = zeros(N_u{a}, 1);
                y_eta_u{a} = zeros(N_u{a}, 1);

                for b = 1:nGrids
                    x_xi_u{a} = x_xi_u{a} + D1_s2u{a,b}{1}*x_s{b};
                    y_xi_u{a} = y_xi_u{a} + D1_s2u{a,b}{1}*y_s{b};
                    x_eta_u{a} = x_eta_u{a} + D1_s2u{a,b}{2}*x_s{b};
                    y_eta_u{a} = y_eta_u{a} + D1_s2u{a,b}{2}*y_s{b};
                end
            end

            for a = 1:nGrids
                J_u{a} = x_xi_u{a}.*y_eta_u{a} - x_eta_u{a}.*y_xi_u{a};
            end
            % ----------------------------------------------

            % x_u = Du*x;
            % x_v = Dv*x;
            % y_u = Du*y;
            % y_v = Dv*y;

            % J = x_u.*y_v - x_v.*y_u;

            % K = cell(dim, dim);
            % K{1,1} = y_v./J;
            % K{1,2} = -y_u./J;
            % K{2,1} = -x_v./J;
            % K{2,2} = x_u./J;

            % Physical derivatives
            % obj.Dx = spdiag( y_v./J)*Du + spdiag(-y_u./J)*Dv;
            % obj.Dy = spdiag(-x_v./J)*Du + spdiag( x_u./J)*Dv;

            % Wrap around Aniosotropic Cartesian. Transformed density and stiffness
            rho_tilde = cell(nGrids, 1);
            PHI = cell(nGrids, 1);

            for a = 1:nGrids
                rho_tilde{a} = J_u{a}.*rho{a};
            end

            for a = 1:nGrids
                PHI{a} = cell(dim,dim,dim,dim);
                for i = 1:dim
                    for j = 1:dim
                        for k = 1:dim
                            for l = 1:dim
                                PHI{a}{i,j,k,l} = 0*C{a}{i,j,k,l};
                                for m = 1:dim
                                    for n = 1:dim
                                        PHI{a}{i,j,k,l} = PHI{a}{i,j,k,l} + J{a}.*K{a}{m,i}.*C{a}{m,j,n,l}.*K{a}{n,k};
                                    end
                                end
                            end
                        end
                    end
                end
            end

            refObj = scheme.Elastic2dStaggeredAnisotropic(g.logic, order, rho_tilde, PHI);

            G = refObj.G;
            U = refObj.U;
            H_u = refObj.H_u;

            % Volume quadrature
            [m, n] = size(refObj.H);
            obj.H = sparse(m, n);
            obj.J = sparse(m, n);
            for a = 1:nGrids
                for i = 1:dim
                    obj.H = obj.H + G{a}*U{a}{i}*spdiag(J_u{a})*refObj.H_u{a}*U{a}{i}'*G{a}';
                    obj.J = obj.J + G{a}*U{a}{i}*spdiag(J_u{a})*U{a}{i}'*G{a}';
                end
            end
            obj.Ji = inv(obj.J);

            % Boundary quadratures on stress grids
            s_w = cell(nGrids, 1);
            s_e = cell(nGrids, 1);
            s_s = cell(nGrids, 1);
            s_n = cell(nGrids, 1);

            % e_w_u = refObj.e_w_u;
            % e_e_u = refObj.e_e_u;
            % e_s_u = refObj.e_s_u;
            % e_n_u = refObj.e_n_u;

            e_w_s = refObj.e_w_s;
            e_e_s = refObj.e_e_s;
            e_s_s = refObj.e_s_s;
            e_n_s = refObj.e_n_s;

            for a = 1:nGrids
                s_w{a} = sqrt((e_w_s{a}'*x_eta{a}).^2 + (e_w_s{a}'*y_eta{a}).^2);
                s_e{a} = sqrt((e_e_s{a}'*x_eta{a}).^2 + (e_e_s{a}'*y_eta{a}).^2);
                s_s{a} = sqrt((e_s_s{a}'*x_xi{a}).^2 + (e_s_s{a}'*y_xi{a}).^2);
                s_n{a} = sqrt((e_n_s{a}'*x_xi{a}).^2 + (e_n_s{a}'*y_xi{a}).^2);
            end

            obj.s_w = s_w;
            obj.s_e = s_e;
            obj.s_s = s_s;
            obj.s_n = s_n;

            % for a = 1:nGrids
                % obj.H_w_s{a} = refObj.H_w_s{a}*spdiag(s_w{a});
                % obj.H_e_s{a} = refObj.H_e_s{a}*spdiag(s_e{a});
                % obj.H_s_s{a} = refObj.H_s_s{a}*spdiag(s_s{a});
                % obj.H_n_s{a} = refObj.H_n_s{a}*spdiag(s_n{a});
            % end

            % Restriction operators
            obj.e_w_u = refObj.e_w_u;
            obj.e_e_u = refObj.e_e_u;
            obj.e_s_u = refObj.e_s_u;
            obj.e_n_u = refObj.e_n_u;

            obj.e_w_s = refObj.e_w_s;
            obj.e_e_s = refObj.e_e_s;
            obj.e_s_s = refObj.e_s_s;
            obj.e_n_s = refObj.e_n_s;

            % Adapt things from reference object
            obj.D = refObj.D;
            obj.U = refObj.U;
            obj.G = refObj.G;

            % obj.e1_w = refObj.e1_w;
            % obj.e1_e = refObj.e1_e;
            % obj.e1_s = refObj.e1_s;
            % obj.e1_n = refObj.e1_n;

            % obj.e2_w = refObj.e2_w;
            % obj.e2_e = refObj.e2_e;
            % obj.e2_s = refObj.e2_s;
            % obj.e2_n = refObj.e2_n;

            % obj.e_scalar_w = refObj.e_scalar_w;
            % obj.e_scalar_e = refObj.e_scalar_e;
            % obj.e_scalar_s = refObj.e_scalar_s;
            % obj.e_scalar_n = refObj.e_scalar_n;

            % e1_w = obj.e1_w;
            % e1_e = obj.e1_e;
            % e1_s = obj.e1_s;
            % e1_n = obj.e1_n;

            % e2_w = obj.e2_w;
            % e2_e = obj.e2_e;
            % e2_s = obj.e2_s;
            % e2_n = obj.e2_n;

            % obj.tau1_w = (spdiag(1./s_w)*refObj.tau1_w')';
            % obj.tau1_e = (spdiag(1./s_e)*refObj.tau1_e')';
            % obj.tau1_s = (spdiag(1./s_s)*refObj.tau1_s')';
            % obj.tau1_n = (spdiag(1./s_n)*refObj.tau1_n')';

            % obj.tau2_w = (spdiag(1./s_w)*refObj.tau2_w')';
            % obj.tau2_e = (spdiag(1./s_e)*refObj.tau2_e')';
            % obj.tau2_s = (spdiag(1./s_s)*refObj.tau2_s')';
            % obj.tau2_n = (spdiag(1./s_n)*refObj.tau2_n')';

            % obj.tau_w = (refObj.e_w'*obj.e1_w*obj.tau1_w')' + (refObj.e_w'*obj.e2_w*obj.tau2_w')';
            % obj.tau_e = (refObj.e_e'*obj.e1_e*obj.tau1_e')' + (refObj.e_e'*obj.e2_e*obj.tau2_e')';
            % obj.tau_s = (refObj.e_s'*obj.e1_s*obj.tau1_s')' + (refObj.e_s'*obj.e2_s*obj.tau2_s')';
            % obj.tau_n = (refObj.e_n'*obj.e1_n*obj.tau1_n')' + (refObj.e_n'*obj.e2_n*obj.tau2_n')';

            % % Physical normals
            % e_w = obj.e_scalar_w;
            % e_e = obj.e_scalar_e;
            % e_s = obj.e_scalar_s;
            % e_n = obj.e_scalar_n;

            % e_w_vec = obj.e_w;
            % e_e_vec = obj.e_e;
            % e_s_vec = obj.e_s;
            % e_n_vec = obj.e_n;

            % nu_w = [-1,0];
            % nu_e = [1,0];
            % nu_s = [0,-1];
            % nu_n = [0,1];

            % obj.n_w = cell(2,1);
            % obj.n_e = cell(2,1);
            % obj.n_s = cell(2,1);
            % obj.n_n = cell(2,1);

            % n_w_1 = (1./s_w).*e_w'*(J.*(K{1,1}*nu_w(1) + K{1,2}*nu_w(2)));
            % n_w_2 = (1./s_w).*e_w'*(J.*(K{2,1}*nu_w(1) + K{2,2}*nu_w(2)));
            % obj.n_w{1} = spdiag(n_w_1);
            % obj.n_w{2} = spdiag(n_w_2);

            % n_e_1 = (1./s_e).*e_e'*(J.*(K{1,1}*nu_e(1) + K{1,2}*nu_e(2)));
            % n_e_2 = (1./s_e).*e_e'*(J.*(K{2,1}*nu_e(1) + K{2,2}*nu_e(2)));
            % obj.n_e{1} = spdiag(n_e_1);
            % obj.n_e{2} = spdiag(n_e_2);

            % n_s_1 = (1./s_s).*e_s'*(J.*(K{1,1}*nu_s(1) + K{1,2}*nu_s(2)));
            % n_s_2 = (1./s_s).*e_s'*(J.*(K{2,1}*nu_s(1) + K{2,2}*nu_s(2)));
            % obj.n_s{1} = spdiag(n_s_1);
            % obj.n_s{2} = spdiag(n_s_2);

            % n_n_1 = (1./s_n).*e_n'*(J.*(K{1,1}*nu_n(1) + K{1,2}*nu_n(2)));
            % n_n_2 = (1./s_n).*e_n'*(J.*(K{2,1}*nu_n(1) + K{2,2}*nu_n(2)));
            % obj.n_n{1} = spdiag(n_n_1);
            % obj.n_n{2} = spdiag(n_n_2);

            % % Operators that extract the normal component
            % obj.en_w = (obj.n_w{1}*obj.e1_w' + obj.n_w{2}*obj.e2_w')';
            % obj.en_e = (obj.n_e{1}*obj.e1_e' + obj.n_e{2}*obj.e2_e')';
            % obj.en_s = (obj.n_s{1}*obj.e1_s' + obj.n_s{2}*obj.e2_s')';
            % obj.en_n = (obj.n_n{1}*obj.e1_n' + obj.n_n{2}*obj.e2_n')';

            % Misc.
            obj.refObj = refObj;
            obj.m = refObj.m;
            obj.h = refObj.h;
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

            [closure, penalty] = obj.refObj.boundary_condition(boundary, bc, tuning);

            type = bc{2};

            switch type
            case {'F','f','Free','free','traction','Traction','t','T'}
                s = obj.(['s_' boundary]);
                s = spdiag(cell2mat(s));
                penalty = penalty*s;
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

            [closure, penalty] = obj.refObj.interface(boundary,neighbour_scheme.refObj,neighbour_boundary,type);
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
            error('Not implemented');
            assertIsMember(boundary, {'w', 'e', 's', 'n'})

            n = obj.(['n_' boundary]);
        end

        % Returns the boundary operator op for the boundary specified by the string boundary.
        % op -- string
        function o = getBoundaryOperator(obj, op, boundary)
            assertIsMember(boundary, {'w', 'e', 's', 'n'})
            assertIsMember(op, {'e', 'e1', 'e2', 'tau', 'tau1', 'tau2', 'en'})

            o = obj.([op, '_', boundary]);

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

            H_b = obj.(['H_', boundary]);
        end

        function N = size(obj)
            N = length(obj.D);
        end
    end
end
