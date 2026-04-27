function [Px2b,Pb2x] = projOpsGlueToBaseGlue(qb, xb, x, tol)
  % [Px2b,Pb2x] = projOpsGlueToBaseGlue(qb, xb, x, tol)
  %
  % Generate projection operators to go from the left, right grids
  % xl, xr to the iterior base glue grid xb and back
  % using the same order of accuracy qb:
  %
  %    x   xb
  %    o   o 
  %    |   | 
  %    o   o 
  %    |   | 
  %    |   o 
  %    |   | 
  %    |   | 
  %    o   o 
  %    |   | 
  %    |   o 
  %    |   | 
  %    |   | 
  %    |   | 
  %    o   o 
  %
  % Here
  %
  %  Vxi_l(i,j) is the right mesh of xi cell edges
  %  Vxi_r(i,j) is the left mesh of xi cell edges
  %
  %  Vxb_lr(i,j) if the glue mesh xi cell edges

  % Dimensional parameters
  basis_size = qb+1; % Basis size on the base glue grid xb
  n_intervals_b = length(xb) - 1;
  n_intervals_x = length(x) - 1;

  % Store the intervals of the base glue grid in a 2xn matrix, where 
  % the first row holds the left edges and the second row holds the
  % right edges of the interval. 
  xb_l = xb(1:end-1)';
  xb_r = xb(2:end)';
  intervals = [xb_l; xb_r];
  
  m_b = basis_size*n_intervals_b; % dimension of base glue grid
  m_x = basis_size*n_intervals_x; % dimension of glue grid x
  Pb2x = sparse(m_x, m_b);
  Px2b = sparse(m_b, m_x);

  % Construct polynomial basis functions and mass matrix for the 
  % basis space [-1,1], evaluated at the Gaussian quadrate points.
  [r, w] = sbp.implementations.gaussianQuad(qb);
  M_basis = massMatrixEntries(r, w, qb); % Mass matrix. Note: Uses non-normalized basis
  V_basis  = sbp.implementations.legendrePolynomials(r, qb, true); % Generalized vandermonde matrix. Note: Uses normalized basis
  
  %% Loop through the elements
  for k = 1:n_intervals_x
    x_l = x(k); % Left node value on glue grid x for interval k
    x_r = x(k+1); % Right node value on glue grid x for interval k

    % Get the intervals xb_k of xb which lie within [x_l, x_r]
    ind_b_k = find(all([xb' < x_r - tol;  xb' >= x_l - tol]));
    intervals_k = intervals(:, ind_b_k);
    xb_k = [intervals_k(1,:), intervals_k(2,end)];

    % Get the projection from [x_l, x_r] to xb_k
    [Px2b_k, Pb2x_k] = projOpsGlueIntervalToBaseGlueInterval(V_basis, M_basis, r, xb_k);

    ind_x = (k-1)*basis_size+1;
    ind_b = (ind_b_k(1)-1)*basis_size+1;

    Pb2x(ind_x:ind_x+qb, ind_b:ind_b+basis_size*length(ind_b_k)-1) = Pb2x_k; %#ok<SPRIX>
    Px2b(ind_b:ind_b+basis_size*length(ind_b_k)-1, ind_x:ind_x+qb) = Px2b_k;  %#ok<SPRIX>
    
  end

end

% Compute the entries of the mass matrix for the basis on the basis interval [-1,1]
function M = massMatrixEntries(r,w,qb)
    Ps  = sbp.implementations.legendrePolynomials(r, qb, false);
    % Since the Legendre polynomials are orthogonal, M is diagonal
    basis_size = qb+1;
    M = zeros(basis_size,1);
    for i = 1:basis_size
        % Compute the norm-squared of the i:th polynomial on [-1,1]
        M(i) = sum(w.*Ps(:,i).^2);
    end
end


function [Px2b_k, Pb2x_k] = projOpsGlueIntervalToBaseGlueInterval(V_basis, M_basis, rq, xb_k)
  % [Px2b_k, Pb2x_k] = projOpsGlueIntervalToBaseGlueInterval(V_basis, M_basis, rq, xb_k)
  %
  % Generate projection operators to go from the f grid to the c grid both with
  % the same order of accuracy N:
  %
  %    x_k    xb_k
  %    o      o
  %    |      |
  %    |      |
  %    |      |
  %    |      o
  %    |      |
  %    |      |
  %    |      |
  %    |      o
  %    |      |
  %    |      |
  %    |      |
  %    o      o
  %
  %  xb_k   are the glue intervals on the base grid contained which lie
  %         with in x_k.
  %
  %  It is assumed that x_k(1) == xb_k(1), x_k(end) == xb_k(end)
  basis_size = length(M_basis);
  qb = basis_size-1; %Polynomial order on the grids
  n_intervals_k = length(xb_k)-1;
  xl = xb_k(1);
  xr = xb_k(end);

  h_k = xr-xl;
  r_k = 2/h_k * xb_k + 1 - 2*xr/h_k; % Shift elements of xb_k to r_k in [-1,1]

  % Evaluate all intervals of r_k in gaussian quadrature points rq
  ind_left = (1:n_intervals_k)'; 
  ind_right = (2:n_intervals_k+1)';
  rq_k = ones(basis_size,1)*r_k(ind_left) + 0.5*(rq+1)*(r_k(ind_right)-r_k(ind_left));

  % Construct projection from xb_k to x_k
  Px2b_k = zeros(length(rq_k(:)), basis_size);
  sqrtMinv_Vinv = diag(sqrt(1./M_basis)) / V_basis;
  sqrtM = diag(sqrt(M_basis));
  for i = 1:n_intervals_k
    Vi = sbp.implementations.legendrePolynomials(rq_k(:,i),qb,true); %Note: Uses normalized basis
    Px2b_k((i-1)*(basis_size)+1:i*(basis_size), :) = sqrtMinv_Vinv*Vi*sqrtM;
  end
  % Metric Jacobians for scaling mass matrices to glue grid interval
  Jx_k = (r_k(ind_right) - r_k(ind_left))/2;
  Jb_k = 1;
  % Construct projection from x_k to xb_k using stability relation 
  Pb2x_k = (1/Jb_k) * diag(1./M_basis) * Px2b_k' * kron(diag(Jx_k), diag(M_basis));
end