function [Pf2g,Pg2f,H,M] = projOpsGridToAlignedGlue(optype,x,h,order,acc)
    [stencil_g2f,BCU_g2f,HU,M] = get_stencils(optype,order,acc);
      
    %%% Setup dimensions of operators %%%
    m_fd = length(x);
    basis_size = length(M);
    n_intervals = length(x)-1;
    m_gg = basis_size*n_intervals;

    %%% Norm matrices %%%
    % FD norm matrix
    H = speye(m_fd,m_fd);
    HUm = length(HU); 
    H(1:HUm,1:HUm) = spdiags(HU,0,HUm,HUm);
    H(m_fd-HUm+1:m_fd,m_fd-HUm+1:m_fd) = spdiags(rot90(HU,2),0,HUm,HUm);
    H = H*h;

    % Glue grid mass matrix
    M = spdiags(M,0,basis_size,basis_size); % Mass matrix for basis
    intervals = x(2:end)-x(1:end-1);
    R = 1/2*spdiags(intervals,0,n_intervals,n_intervals); % Rate of change in integration limits going from basis to glue
    M = kron(R,M); % Mass matrix for glue grid
    %%%%%%%%%%%%%%%%%%%%%%

    %%% Create Pg2f from stencil and BC %%%
    stencil_width = length(stencil_g2f);
    [BC_rows,BC_cols] = size(BCU_g2f);
    n_intervals_int = stencil_width/basis_size;
    n_intervals_bnd = BC_cols/basis_size;

    Pg2f = sparse(m_fd,m_gg);
    % Interior
    for row = BC_rows+1 : m_fd-BC_rows
        col_start = basis_size*(row-1-n_intervals_int/2)+1; 
        cols = col_start:col_start+stencil_width-1;
        Pg2f(row,cols) = stencil_g2f;  %#ok<SPRIX>
    end

    % Boundary closures
    % Insert upper block
    Pg2f(1:BC_rows,1:BC_cols) = BCU_g2f;
    % The lower block is symmetric w.r.t to grid points and glue interval, but
    % keeps the ordering of the glue grid modes. Also apply symmetry conditions
    % to even/odd modes:
    % 1. Reorder the rows such that the last row corresponds to the
    % final point
    BCL_g2f = flipud(BCU_g2f);
    BCL_by_glue_interval = reshape(BCL_g2f,BC_rows,basis_size,n_intervals_bnd);
    % 3. Reorder the glue intervals, such that the 0th interval corresponds to 
    % the final point.
    BCL_by_glue_interval = flip(BCL_by_glue_interval,3); 
    % 4. Reshape back into a 2D array containing all the glue intervals
    BCL_g2f = reshape(BCL_by_glue_interval,BC_rows,BC_cols);
    % 5. Apply symmetry conditions to odd modes
    even_odd_modes = 2*mod(1:basis_size,2)-1;
    symmetry_cond_bc_r = repmat(even_odd_modes,1,n_intervals_bnd);
    BCL_g2f = symmetry_cond_bc_r.*BCL_g2f;
    % Insert into operator
    Pg2f(end-BC_rows+1:end,end-BC_cols+1:end) = BCL_g2f;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Create Pf2g using stability (adjoint) relation %%%%
    Pf2g = M\(Pg2f'*H);
end

function [stencil_g2f,BCU_g2f,HU,M] = get_stencils(optype,order,good_dir)
    switch optype
    case {'Traditional','trad','Standard','standard'}
        [stencil_g2f,BCU_g2f,HU,M] = get_eq_stencils(order,good_dir);
    case {'boundary-optimized','bopt','noneq'}
        [stencil_g2f,BCU_g2f,HU,M] = get_noneq_stencils(order,good_dir);
    end
end

% TODO: Constructed by minimizing the spectrum rather than the projection error.
%       Consider chnaging to operators with minimzed error instead.
function [stencil_g2f,BCU_g2f,HU,M] = get_eq_stencils(order,good_dir)
    switch order
    case 2
        switch good_dir
        case 'none'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOp_2_accGlue2Fd1_accFd2Glue1();
        case 'glue2fd'
            error('Not implemented');
        case 'fd2glue'
            error('Not implemented');
        end
    case 4
        switch good_dir
        case 'glue2fd'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_4_accGlue2Fd3_accFd2Glue2();
        case 'fd2glue'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_4_accGlue2Fd2_accFd2Glue3();                        
        end
    case 6
        switch good_dir
        case 'glue2fd'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_6_accGlue2Fd4_accFd2Glue3();
        case 'fd2glue'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_6_accGlue2Fd3_accFd2Glue4();
        end
    case 8
        switch good_dir
        case 'glue2fd'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_8_accGlue2Fd5_accFd2Glue4();
        case 'fd2glue'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_8_accGlue2Fd4_accFd2Glue5();
        end
    otherwise
        error('Order %d not implemented', order);
    end
end

function [stencil_g2f,BCU_g2f,HU,M] = get_noneq_stencils(order,good_dir)
    switch order
    case 4
        switch good_dir
        case 'glue2fd'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_4_noneq_accGlue2Fd3_accFd2Glue2();
        case 'fd2glue'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_4_noneq_accGlue2Fd2_accFd2Glue3();
        end
    case 6
        switch good_dir
        case 'glue2fd'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_6_noneq_accGlue2Fd4_accFd2Glue3();
        case 'fd2glue'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_6_noneq_accGlue2Fd3_accFd2Glue4();
        end
    case 8
        switch good_dir
        case 'none'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOp_8_noneq_accGlue2Fd4_accFd2Glue4(); % TODO: Recreate with minimized boundary error!
        case 'glue2fd'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_8_noneq_accGlue2Fd5_accFd2Glue4();
        case 'fd2glue'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_8_noneq_accGlue2Fd4_accFd2Glue5();
        end
    case 10
        switch good_dir
        case 'glue2fd'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_10_noneq_accGlue2Fd6_accFd2Glue5();
        case 'fd2glue'
            [stencil_g2f,BCU_g2f,HU,M] = sbp.implementations.projOpOp_10_noneq_accGlue2Fd5_accFd2Glue6();    
        end
    otherwise
        error('Order %d not implemented', order);
    end
end