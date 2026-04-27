function d = diracDiscrCurve(x_s, g, m_order, s_order, order, scheme)
    % 2-dimensional delta function for single-block curvilinear grid
    % x_s:      source point coordinate vector, e.g. [x, y] or [x, y, z].
    % g:        single-block grid containing the source
    % m_order:  Number of moment conditions
    % s_order:  Number of smoothness conditions
    % order:    Order of SBP derivative approximations
    % opSet:    Cell array of function handle to opSet generator

    default_arg('scheme', 'standard')
    default_arg('order', m_order);
    % default_arg('opSet', {@sbp.D2Variable, @sbp.D2Variable});

    dim = length(x_s);
    assert(dim == 2, 'diracDiscrCurve: Only implemented for 2d.');
    % assert(isa(g, 'grid.Curvilinear'));

    switch scheme
    case {'standard', 'narrow'}
        opSetMetric = {@sbp.D2Variable, @sbp.D2Variable};
        orderMetric = order;
        opSet = {@sbp.D2Variable, @sbp.D2Variable};
    case 'upwind'
        opSetMetric = {@sbp.D2Variable, @sbp.D2Variable};
        orderMetric = ceil(order/2)*2;
        opSet = {@sbp.D1Upwind, @sbp.D1Upwind};
    end

    switch scheme
    case {'standard', 'narrow', 'upwind'}

        % Compute metric terms
        m = g.size();
        m_u = m(1);
        m_v = m(2);

        ops_u = opSetMetric{1}(m_u, {0, 1}, orderMetric);
        ops_v = opSetMetric{2}(m_v, {0, 1}, orderMetric);
        I_u = speye(m_u);
        I_v = speye(m_v);

        D1_u = ops_u.D1;
        D1_v = ops_v.D1;

        Du = kr(D1_u,I_v);
        Dv = kr(I_u,D1_v);

        u = ops_u.x;
        v = ops_v.x;

        % Compute Jacobian
        coords = g.points();
        x = coords(:,1);
        y = coords(:,2);

        x_u = Du*x;
        x_v = Dv*x;
        y_u = Du*y;
        y_v = Dv*y;

        J = x_u.*y_v - x_v.*y_u;

        % Find approximate logical coordinates of point source
        [U, V] = meshgrid(u, v);
        U_interp = scatteredInterpolant(coords, U(:));
        V_interp = scatteredInterpolant(coords, V(:));
        uS = U_interp(x_s);
        vS = V_interp(x_s);

        % Make sure that we don't accidentally end up outside domain
        tol = 1e-12;
        if abs(uS) < tol;
            uS = 0;
        end
        if abs(uS-1) < tol;
            uS = 1;
        end
        if abs(vS) < tol;
            vS = 0;
        end
        if abs(vS-1) < tol;
            vS = 1;
        end

        % Get quadratures
        ops_u = opSet{1}(m_u, {0, 1}, order);
        ops_v = opSet{2}(m_v, {0, 1}, order);
        u = ops_u.x;
        v = ops_v.x;
        H_u = ops_u.H;
        H_v = ops_v.H;

        d = (1./J).*diracDiscr([uS, vS], {u, v}, m_order, s_order, {H_u, H_v});

    case 'staggered'

        g_u = g.gridGroups{1};
        g_s = g.gridGroups{2};
        nGrids = 2;

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
        D1_s2u{1,1}{1} = kron(D1d{1}, Ip{2});
        D1_s2u{1,2}{2} = kron(Ip{1}, D1d{2});
        D1_s2u{2,1}{2} = kron(Id{1}, D1p{2});
        D1_s2u{2,2}{1} = kron(D1p{1}, Id{2});

        D1_s2u{1,1}{2} = sparse(N_u{1}, N_s{1});
        D1_s2u{2,2}{2} = sparse(N_u{2}, N_s{2});
        D1_s2u{1,2}{1} = sparse(N_u{1}, N_s{2});
        D1_s2u{2,1}{1} = sparse(N_u{2}, N_s{1});

        %  --- Metric coefficients on displacement grids -------
        x_sg = cell(nGrids, 1);
        y_sg = cell(nGrids, 1);
        J = cell(nGrids, 1);
        x_xi = cell(nGrids, 1);
        x_eta = cell(nGrids, 1);
        y_xi = cell(nGrids, 1);
        y_eta = cell(nGrids, 1);

        for a = 1:nGrids
            coords = g_s{a}.points();
            x_sg{a} = coords(:,1);
            y_sg{a} = coords(:,2);
        end

        for a = 1:nGrids
            x_xi{a} = zeros(N_u{a}, 1);
            y_xi{a} = zeros(N_u{a}, 1);
            x_eta{a} = zeros(N_u{a}, 1);
            y_eta{a} = zeros(N_u{a}, 1);

            for b = 1:nGrids
                x_xi{a} = x_xi{a} + D1_s2u{a,b}{1}*x_sg{b};
                y_xi{a} = y_xi{a} + D1_s2u{a,b}{1}*y_sg{b};
                x_eta{a} = x_eta{a} + D1_s2u{a,b}{2}*x_sg{b};
                y_eta{a} = y_eta{a} + D1_s2u{a,b}{2}*y_sg{b};
            end
        end

        for a = 1:nGrids
            J{a} = x_xi{a}.*y_eta{a} - x_eta{a}.*y_xi{a};
        end
        % ----------------------------------------------

        u1 = ops{1}.x_primal;
        v1 = ops{2}.x_primal;
        u2 = ops{1}.x_dual;
        v2 = ops{2}.x_dual;

        % Find approximate logical coordinates of point source
        [U, V] = meshgrid(u1, v1);
        coords = g_u{1}.points();
        U_interp = scatteredInterpolant(coords, U(:));
        V_interp = scatteredInterpolant(coords, V(:));
        uS = U_interp(x_s);
        vS = V_interp(x_s);

        % Make sure that we don't accidentally end up outside domain
        tol = 1e-12;
        if abs(uS) < tol;
            uS = 0;
        end
        if abs(uS-1) < tol;
            uS = 1;
        end
        if abs(vS) < tol;
            vS = 0;
        end
        if abs(vS-1) < tol;
            vS = 1;
        end

        d.d1 = (1./J{1}).*diracDiscr([uS, vS], {u1, v1}, m_order, s_order, {Hp{1}, Hp{2}});
        d.d2 = (1./J{2}).*diracDiscr([uS, vS], {u2, v2}, m_order, s_order, {Hd{1}, Hd{2}});

    end

end