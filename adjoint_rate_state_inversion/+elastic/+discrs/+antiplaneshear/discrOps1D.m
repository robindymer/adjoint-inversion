function [D, E, bc_struct, ic_struct] = discrOps1D(domain, mbGrid, mbDiffOp, material, bc, ic_method)
    E = componentOps(domain, mbGrid, mbDiffOp, bc.method, ic_method);
    [closure_bc, bc_struct] = bcOps(domain, mbGrid, mbDiffOp, material, bc);
    [closure_ic, ic_struct] = faultInterfaceOps(domain, mbGrid, mbDiffOp, material, E, ic_method);
    D = diffOp(mbGrid, mbDiffOp, material, E, closure_bc, closure_ic, bc.method, bc_struct, ic_method, ic_struct);
end

function D = diffOp(mbGrid, mbDiffOp, material, E, closure_bc, closure_ic, bc_method, bc_struct, ic_method, ic_struct)
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
      
     %--- Add contributions from standrd outflow BC:s or erickson2022 bcs
     switch bc_method
     case 'standard'
         if ~isempty(bc_struct.outflow_penalty) % Hacky : (
             penalty = blockmatrix.toMatrix(bc_struct.outflow_penalty);
             c_inv = elastic.helpers.cell_row_to_diag_blockmatrix(bc_struct.c_inv);
             e = blockmatrix.toMatrix(bc_struct.e);
             %- Add contributions to system matrix-
             D = D - E.v'*penalty*c_inv*e'*E.v;
         end
     case 'erickson2022'
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
 
         ns = length(R);
         IR = eye(size(R));
         %------------------------------------
 
         % ----Build operators for fluxes ----
 
         % Contributions to tau* from u, u_t and u*.
         taus_u = (IR - R)/2*tau_op' + gamma/2*(R - IR)*e';
         taus_v = Z*(R - IR)/2*e';
         taus_us = -gamma/2*(R - IR);
 
         % Contributions to u*_t from u, u_t and u_*
         ust_u = -Z\((IR + R)/2*tau_op' - gamma/2*(IR + R)*e');
         ust_v = (R + IR)/2*e';
         ust_us = -Z\(gamma/2*(IR + R));
         
         %- Add contributions to system matrix-         
         % SATs with tau*
         D = D + E.v'*tau_star_penalty*(taus_u*E.u + taus_v*E.v + taus_us*E.us);
 
         % SATs with u*
         D = D + E.v'*u_star_penalty*E.us;
 
         % Evolution of u*
         D = D + E.us'*(ust_u*E.u + ust_v*E.v + ust_us*E.us);
     end

     switch ic_method
     case 'erickson2022'
        em = ic_struct.em;
        ep = ic_struct.ep;
        Zm = ic_struct.Zm;
        Zp = ic_struct.Zp;
        Tm = ic_struct.Tm;
        Tp = ic_struct.Tp;
        gammam = ic_struct.gammam;
        gammap = ic_struct.gammap;


         %---- ODEs for us_p and us_m --------------------------------
         D = D + E.usim'*(em'*E.v - Zm\(Tm'*E.u + gammam*(E.usim - em'*E.u)));
         D = D + E.usip'*(ep'*E.v - Zp\(Tp'*E.u + gammap*(E.usip - ep'*E.u)));

         %---- Penalties ------
         H = mbDiffOp.H;
         rho_mat = spdiag(multiblock.evalOn(mbGrid, material.rho));

         % Terms with u* are linear
         D = D - E.v' * inv(rho_mat*H) * Tm*E.usim;
         D = D - E.v' * inv(rho_mat*H) * Tp*E.usip;
     end
end

function E = componentOps(domain, mbGrid, mbDiffOp, bc_method, ic_method)
    E = struct;
    [n, ns, nsip, nsim] = elastic.discrs.antiplaneshear.nunknowns(domain, mbGrid, mbDiffOp, bc_method, ic_method);
    % ---- Component operators -----   
    E.u = [speye(n), sparse(n,n), sparse(n,1), sparse(n,ns), sparse(n,nsim), sparse(n,nsip)];
    E.v = [sparse(n,n), speye(n), sparse(n,1), sparse(n,ns), sparse(n,nsim), sparse(n,nsip)];
    E.Psi = [sparse(1,n), sparse(1,n), 1, sparse(1,ns), sparse(1,nsim), sparse(1,nsip)];
    % erickson2022 component operators. Are all zero sparse matrices if standard bc/ic is used.
    E.us = [sparse(ns,n), sparse(ns,n), sparse(ns,1), speye(ns), sparse(ns,nsim), sparse(ns, nsip)];
    E.usim = [sparse(nsim,n), sparse(nsim,n), sparse(nsim,1), sparse(nsim,ns), speye(nsim), sparse(nsim,nsip)];
    E.usip = [sparse(nsip,n), sparse(nsip,n), sparse(nsip,1), sparse(nsip,ns), sparse(nsip,nsim), speye(nsip)];
end

function [closure_tot, bc_struct] = bcOps(domain, mbGrid, mbDiffOp, material, bc)

    closure_tot = 0*mbDiffOp.D;
    bc_struct = struct;
    bc_struct.boundary_data_fun = {};
    switch bc.method 
    case 'standard'
        bc_struct.c_inv = {};
        bc_struct.e = {};
        bc_struct.outflow_penalty = {};
        bc_struct.e = {};
    case 'erickson2022'
        bc_struct.R = {};
        bc_struct.Z = {};
        bc_struct.gamma = {};
        bc_struct.tau_star_penalty = {};
        bc_struct.u_star_penalty = {};
        bc_struct.tau_op = {};
        bc_struct.e = {};
    end

    % Initialize before looping over bc
    bc.ids = {domain.boundaryGroups.left, domain.boundaryGroups.right};
    % Loop over boundaries to impose bc on
    for i = 1:numel(bc.ids)
        bid = bc.ids{i};
        type = bc.type{i};
        e = mbDiffOp.getBoundaryOperator('e',bid);

        switch bc.method
        case 'standard'
            switch type
            case 'outflow'
               % Create operator for 1/c on boundary
               mu_boundary = e'*multiblock.evalOn(mbGrid, material.mu);
               rho_boundary = e'*multiblock.evalOn(mbGrid, material.rho);
               c_boundary_inv = diag(1./(sqrt(mu_boundary./rho_boundary)));

               % Setup closures and penalties as a neumann condition
               % (where the penalty will be included in the system matrix to
               % act on velocity)
               % Store in struct for setting up system later
               [closure, penalty] = mbDiffOp.boundary_condition(bid,'neumann');
               bc_struct.c_inv{end+1} = c_boundary_inv;
               bc_struct.outflow_penalty{end+1} = penalty;
               bc_struct.e{end+1} = e;
            otherwise % Dirichlet or Neumann conditions
                [closure, penalty] = mbDiffOp.boundary_condition(bid,type);
                data = bc.data{i}
                if ~isempty(data)
                    if ~isempty(boundary_data_fun)
                        bc_struct.boundary_data_fun = @(t,u,ut) bc_struct.boundary_data_fun(t,u,ut) + penalty*data(t,u,ut);
                    else
                        bc_struct.boundary_data_fun = @(t,u,ut) penalty*data(t,u,ut);
                    end
                end
            end
        case 'erickson2022'
            % Get closure only 
            [closure, ~] = mbDiffOp.boundary_condition(bid,'erickson2022');

            % Get operators
            rho_vec = multiblock.evalOn(mbGrid, material.rho);
            mu_vec = multiblock.evalOn(mbGrid, material.mu);
            mu_boundary = e'*mu_vec;
            rho_boundary = e'*rho_vec;
            Z = spdiag(sqrt(mu_boundary.*rho_boundary));
            mu_boundary = diag(mu_boundary);
            rho_boundary = diag(rho_boundary);
            H = mbDiffOp.H;

            tau_op = mbDiffOp.getBoundaryOperator('d', bid)*mu_boundary;
            % Sign for normal derivate. Hacky solution :(
            if i == 1
                tau_op = -tau_op;
            end

            % Set up penalties in front of fluxes
            H = mbDiffOp.H;
            rho_mat = spdiag(rho_vec);
            tau_star_penalty = (rho_mat*H)\e;
            u_star_penalty = -(rho_mat*H)\tau_op;

            I = eye(size(rho_boundary));
            switch type
            case 'outflow'
                R = 0*I;
            case 'dirichlet'
                R = -1*I;
            case 'neumann'
                R = I;
            end
            R = sparse(R);

            % Penalty strength gamma
            gamma = 1.1/mbDiffOp.diffOps{1}.gamm*mu_boundary;

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

function [closure, ic_struct] = faultInterfaceOps(domain, mbGrid, mbDiffOp, material, E, ic_method)
    ic_struct = struct;
    % ---- Friction  -----
    bid = domain.boundaryGroups.interface;
    ic_struct.e = mbDiffOp.getBoundaryOperator('e',bid);
    switch ic_method
    case 'standard'
        [closure, ic_struct.penalty] = mbDiffOp.boundary_condition(bid, 'traction');
    case 'erickson2022'
        % Penalty is not used, is a zero matrix.
        [closure, ic_struct.penalty] =  mbDiffOp.boundary_condition(bid, 'erickson2022');

        bidm = domain.boundaryGroups.interface{1};
        bidp = domain.boundaryGroups.interface{2};
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

        %---- Operators that compute characteristic variables ----
        % Traction operators, with sign correction
        Tm = (-1 + 2*(bidm{2}=='r')) * mbDiffOp.getBoundaryOperator('d', bidm)*mum;
        Tp = (-1 + 2*(bidp{2}=='r')) * mbDiffOp.getBoundaryOperator('d', bidp)*mup;

        % Penalty strengths
        gammam = 1.1/mbDiffOp.diffOps{1}.gamm*mum;
        gammap = 1.1/mbDiffOp.diffOps{1}.gamm*mup;

        Wm = Zm*em'*E.v - Tm'*E.u - gammam*(E.usim - em'*E.u);
        Wp = Zp*ep'*E.v - Tp'*E.u - gammap*(E.usip - ep'*E.u);
       
        %---- ODEs for us_p and us_m --------------------------------
        H = mbDiffOp.H;
        rho_mat = spdiag(rho_vec);
        interface_ode = (E.usip'*inv(Zp) - E.usim'*inv(Zm)); % To be multiplied by tau*_+(t, U)
        interface_penalties = E.v'* inv(rho_mat*H)*(ep-em); % To be multiplied by tau*_+(t, U)

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
end

