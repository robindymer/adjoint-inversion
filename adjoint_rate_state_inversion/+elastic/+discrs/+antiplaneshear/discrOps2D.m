function [D, E, bc_struct, ic_struct] = discrOps2D(domain, mbGrid, mbDiffOp, material, bc)
    E = componentOps(domain, mbGrid, mbDiffOp, bc);
    [closure_bc, bc_struct] = bcOps(mbGrid, mbDiffOp, material, bc, E);
    [closure_ic, ic_struct] = faultInterfaceOps(domain, mbGrid, mbDiffOp, material, E);
    D = diffOp(domain, mbGrid, mbDiffOp, material, E, closure_bc, closure_ic, bc_struct, ic_struct);
end

function D = diffOp(domain, mbGrid, mbDiffOp, material, E, closure_bc, closure_ic, bc_struct, ic_struct)
     % ---- 1st order system -----
     % u = [u_m; u_p];
     % v = u_t
 
     % -- Standard:
     % u = [u_m; u_p];
     % v = u_t
     % w = [u; v; Psi]
     % D = [0         I 0;
     %      D_laplace 0 0;
     %      0         0 0];
     % F = [0; S + penalty_bc*bc_data + penalty_fault*fault_fun; g]
     % w_t = D*w + F(w)
 
     % -- Erickson2022:
     % us: boundary fluxes. The ODE for us is linear and can be built into D.
     % usim, usip: interface fluxes (minus and plus)
     % The ODEs for the interface fluxes have both linear and nonlinear parts.
     %
     % w = [u; v; Psi; us; usim; usip]
     % F = [0; S + penalty_bc*bc_data + penalty_fault*fault_fun; g; 0; nonlinear part of evolution for interface u*]
     % w_t = D*w + F(w)

    % ---- Construct scheme for multiblock Laplace including the closures -----
    D_laplace = mbDiffOp.D + closure_bc + closure_ic;

    %---- Build system matrix-----
    D = E.v'*D_laplace*E.u;
    D = D + E.u'*E.v;   
      
     %--- Add contributions from erickson2022 BCs to the difference operator
     if ~isempty(bc_struct.R)
        %---- Unpack erickson struct ------
        % Block-diagonal matrices
        R = elastic.helpers.cell_row_to_diag_blockmatrix(bc_struct.R);
        Z = elastic.helpers.cell_row_to_diag_blockmatrix(bc_struct.Z);
        gamma = elastic.helpers.cell_row_to_diag_blockmatrix(bc_struct.gamma);

        % Tall and skinny matrices, horizontally stacked
        tau_star_penalty = blockmatrix.toMatrix(bc_struct.tau_star_penalty);
        u_star_penalty = blockmatrix.toMatrix(bc_struct.u_star_penalty);
        tau_op = blockmatrix.toMatrix(bc_struct.tau_op);
        e = blockmatrix.toMatrix(bc_struct.e);

        [ns, ~] = size(E.us);
        IR = speye(ns);
        %------------------------------------
 
        % ----Build operators for fluxes ----
 
        % Contributions to tau* from u, u_t and u*.
        taus_u = (IR - R)/2*tau_op' + gamma/2*(R - IR)*e';
        taus_v = Z*(R - IR)/2*e';
        taus_us = -gamma/2*(R - IR);

        % Contributions to u*_t from u, u_t and u_*
        ust_u = -inv(Z)*((IR + R)/2*tau_op' - gamma/2*(IR + R)*e');
        ust_v = (R + IR)/2*e';
        ust_us = -inv(Z)*gamma/2*(IR + R);
        
        %- Add contributions to system matrix-         
        % SATs with tau*
        D = D + E.v'*tau_star_penalty*(taus_u*E.u + taus_v*E.v + taus_us*E.us);

        % SATs with u*
        D = D + E.v'*u_star_penalty*E.us;

        % Evolution of u*
        D = D + E.us'*(ust_u*E.u + ust_v*E.v + ust_us*E.us);
     end

    em = ic_struct.em;
    ep = ic_struct.ep;
    Zm = ic_struct.Zm;
    Zp = ic_struct.Zp;
    Tm = ic_struct.Tm;
    Tp = ic_struct.Tp;
    gammam = ic_struct.gammam;
    gammap = ic_struct.gammap;

    %---- ODEs for us_p and us_m --------------------------------
    Zm_mat = spdiag(Zm);
    Zp_mat = spdiag(Zp);
    D = D + E.usim'*(em'*E.v - Zm_mat\(Tm'*E.u + gammam*(E.usim - em'*E.u)));
    D = D + E.usip'*(ep'*E.v - Zp_mat\(Tp'*E.u + gammap*(E.usip - ep'*E.u)));

    %---- Penalties ------
    H = mbDiffOp.H;
    bidm = domain.boundaryGroups.fault_minus;
    H_b = mbDiffOp.getBoundaryQuadrature(bidm);
    rho_mat = spdiag(multiblock.evalOn(mbGrid, material.rho));

    D = D - E.v' * inv(rho_mat*H) * Tm*H_b*E.usim;
    D = D - E.v' * inv(rho_mat*H) * Tp*H_b*E.usip;
end

function E = componentOps(domain, mbGrid, mbDiffOp, bc)
    E = struct;
    [n, ns, nsip, nsim] = elastic.discrs.antiplaneshear.nunknowns2D(domain, mbGrid, mbDiffOp, bc);
    nf = nsim;
    % ---- Component operators -----   
    E.u = [speye(n), sparse(n,n), sparse(n,nf), sparse(n,ns), sparse(n,nsim), sparse(n,nsip)];%, sparse(n,nf)];
    E.v = [sparse(n,n), speye(n), sparse(n,nf), sparse(n,ns), sparse(n,nsim), sparse(n,nsip)];%, sparse(n,nf)];
    E.Psi = [sparse(nf,n), sparse(nf,n), speye(nf), sparse(nf,ns), sparse(nf,nsim), sparse(nf,nsip)];%, sparse(nf,nf)];
    E.us = [sparse(ns,n), sparse(ns,n), sparse(ns,nf), speye(ns), sparse(ns,nsim), sparse(ns, nsip)];%, sparse(ns,nf)];
    E.usim = [sparse(nsim,n), sparse(nsim,n), sparse(nsim,nf), sparse(nsim,ns), speye(nsim), sparse(nsim,nsip)];%, sparse(nsim,nf)];
    E.usip = [sparse(nsip,n), sparse(nsip,n), sparse(nsip,nf), sparse(nsip,ns), sparse(nsip,nsim), speye(nsip)];%, sparse(nsip,nf)];
    %E.Vs = [sparse(nsim,n), sparse(nsim,n), sparse(nsim,nf), sparse(nsim,ns), sparse(nsim,nsim), sparse(nsim,nsip), speye(nf)];
end

function [closure_tot, bc_struct] = bcOps(mbGrid, mbDiffOp, material, bc, E)

    closure_tot = 0*mbDiffOp.D;
    bc_struct = struct;
    bc_struct.R = {};
    bc_struct.Z = {};
    bc_struct.gamma = {};
    bc_struct.tau_star_penalty = {};
    bc_struct.u_star_penalty = {};
    bc_struct.tau_op = {};
    bc_struct.e = {};
    bc_struct.penalized_data = [];

    % Loop over boundaries
    for i = 1:numel(bc.ids)
        bid = bc.ids{i};
        type = bc.type{i};
        data = bc.data{i};
        if ~isempty(data)
            [closure, penalty] = mbDiffOp.boundary_condition(bid,type);
            penalty = E.v'*penalty;
            if ~isempty(bc_struct.penalized_data)
                bc_struct.penalized_data = @(t) bc_struct.penalized_data(t) + penalty*data(t);
            else
                bc_struct.penalized_data = @(t) penalty*data(t);
            end
        else
            e = mbDiffOp.getBoundaryOperator('e',bid);
            % Get closure only 
            [closure, ~] = mbDiffOp.boundary_condition(bid,'erickson2022');

            % Get operators
            rho_vec = multiblock.evalOn(mbGrid, material.rho);
            rho_mat = spdiag(rho_vec);
            mu_vec = multiblock.evalOn(mbGrid, material.mu);
            mu_boundary = e'*mu_vec;
            rho_boundary = e'*rho_vec;
            Z = spdiag(sqrt(mu_boundary.*rho_boundary));
            mu_boundary = diag(mu_boundary);
            rho_boundary = diag(rho_boundary);
            H = mbDiffOp.H;

            tau_op = mbDiffOp.getBoundaryOperator('traction', bid);
            H_b = mbDiffOp.getBoundaryQuadrature(bid);

            % Set up penalties in front of fluxes
            tau_star_penalty = (rho_mat*H)\(e*H_b);
            u_star_penalty = -(rho_mat*H)\(tau_op*H_b);

            I = eye(size(rho_boundary));
            switch type
            case 'outflow'
                R = 0*I;
            case 'dirichlet'
                R = -1*I;
            case {'neumann', 'traction'}
                R = I;
            end
            R = sparse(R);

            % Penalty strength gamma
            gamma = mbDiffOp.getBoundaryField('gamma', bid);

            bc_struct.R{end+1} = R;
            bc_struct.Z{end+1} = Z;
            bc_struct.gamma{end+1} = gamma;
            bc_struct.tau_star_penalty{end+1} = tau_star_penalty;
            bc_struct.u_star_penalty{end+1} = u_star_penalty;
            bc_struct.tau_op{end+1} = tau_op;
            bc_struct.e{end+1} = e;
        end
        closure_tot = closure_tot + closure;
    end
end

function [closure, ic_struct] = faultInterfaceOps(domain, mbGrid, mbDiffOp, material, E)
    ic_struct = struct;
    % ---- Friction  -----
    bid = domain.boundaryGroups.fault;
    e = mbDiffOp.getBoundaryOperator('e',bid);
    % Penalty is not used, is a zero matrix.
    [closure, ~] =  mbDiffOp.boundary_condition(bid, 'erickson2022');

    bidm = domain.boundaryGroups.fault_minus;
    bidp = domain.boundaryGroups.fault_plus;
    em = mbDiffOp.getBoundaryOperator('e', bidm);
    ep = mbDiffOp.getBoundaryOperator('e', bidp);

    % Compute impedances ( Z = sqrt(rho*mu) ) at interface
    rho_vec = multiblock.evalOn(mbGrid, material.rho);
    mu_vec = multiblock.evalOn(mbGrid, material.mu);
    rhom = em'*rho_vec;
    mum = em'*mu_vec;
    Zm = sqrt(rhom.*mum);
    rhop = ep'*rho_vec;
    mup = ep'*mu_vec;
    Zp = sqrt(rhop.*mup);
    Zm_mat = spdiag(Zm);
    Zp_mat = spdiag(Zp);

    %---- Operators that compute characteristic variables ----
    % Traction operators, with sign correction
    Tm = mbDiffOp.getBoundaryOperator('traction', bidm);
    Tp = mbDiffOp.getBoundaryOperator('traction', bidp);

    % Penalty strengths
    gammam = mbDiffOp.getBoundaryField('gamma', bidm);
    gammap = mbDiffOp.getBoundaryField('gamma', bidp);

    Wm = Zm_mat*em'*E.v - Tm'*E.u - gammam*(E.usim - em'*E.u);
    Wp = Zp_mat*ep'*E.v - Tp'*E.u - gammap*(E.usip - ep'*E.u);
    
    %---- Penalties ------
    H = mbDiffOp.H;
    H_b = mbDiffOp.getBoundaryQuadrature(bidm);
    rho_mat = spdiag(rho_vec);
    interface_ode = -(E.usip'*inv(Zp_mat) - E.usim'*inv(Zm_mat));
    interface_penalties = -E.v'*inv(rho_mat*H)*(ep-em)*H_b;

    ic_struct.e = e;
    ic_struct.em = em;
    ic_struct.ep = ep;
    ic_struct.Zm = Zm;
    ic_struct.Zp = Zp;
    ic_struct.Tm = Tm;
    ic_struct.Tp = Tp;
    ic_struct.Wm = Wm;
    ic_struct.Wp = Wp;
    ic_struct.gammam = gammam;
    ic_struct.gammap = gammap;
    ic_struct.penalty = interface_ode+interface_penalties;
end

