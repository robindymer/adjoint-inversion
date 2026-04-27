% projOpsGridToGrid Create glue interpolation operators for two SBP blocks which have
% different numbers of grid points
% [Iu2v,Iv2u] = projOpsGridToGrid(optype_u, x_u, h_u, q_u, acc_u, optype_v, x_v, h_v, q_v, acc_v)
%
% inputs:
%   x_u:   finite difference grid on side u
%   h_u:   grid spacing on side u
%   q_u:   finite difference order on side u
%   acc_u: finite difference order on side u
%   x_v:   finite difference grid on side v
%   h_v:   grid spacing on side v
%   q_v:   finite difference order on side v
%		
%
% outputs:
%   Iu2v: glue interpolation from side u to side v
%   Iv2u: glue interpolation from side v to side u
function [Iu2v,Iv2u] = projOpsGridToGrid(optype_u, x_u, h_u, q_u, acc_u, optype_v, x_v, h_v, q_v, acc_v)
	% Turn x_u and x_v into column vectors in case they are not
	x_u = x_u(:);
	x_v = x_v(:);
    
	% Get the projection operators for each side to the aligned FD-glue grids.
	[Px2g_u,Pg2x_u] = sbp.implementations.projOpsGridToAlignedGlue(optype_u,x_u,h_u,q_u,acc_u); % fd and glue grid for u
	[Px2g_v,Pg2x_v] = sbp.implementations.projOpsGridToAlignedGlue(optype_v,x_v,h_v,q_v,acc_v); % fd and glue grid for v

	% Align polynomial order of the projection operators to max(q_u,q_v)-1
	q_b = max(q_u,q_v)-1; % Degree of polynomials on the base glue grid
	[Px2g_u, Pg2x_u] = alignPolynomialOrder(x_u,q_u,q_b,Px2g_u,Pg2x_u);
	[Px2g_v, Pg2x_v] = alignPolynomialOrder(x_v,q_v,q_b,Px2g_v,Pg2x_v);

	% Create the base glue grid as the union of x_u and x_v
	tol = 100*eps;
	x_b = createBaseGlueGrid(x_u,x_v,tol);

	% Create the projection operators from glue spaces to 
	% the base glue space
	[Pg2b_u, Pb2g_u] = sbp.implementations.projOpsGlueToBaseGlue(q_b, x_b, x_u, tol);
	[Pg2b_v, Pb2g_v] = sbp.implementations.projOpsGlueToBaseGlue(q_b, x_b, x_v, tol);

	% Create interpolation operators from x_u to x_v, by 
	% stacking the projection operators together
	Iu2v = Pg2x_v * Pb2g_v * Pg2b_u * Px2g_u;
	Iv2u = Pg2x_u * Pb2g_u * Pg2b_v * Px2g_v;
end

function [Pf2g, Pg2f] = alignPolynomialOrder(x,q,q_b,Pf2g,Pg2f)
	ngi = length(x)-1; %number of glue intervals
	I = speye(ngi);
	Pf2g = kron(I,speye(q_b+1,q))*Pf2g;
	Pg2f = Pg2f*kron(I,speye(q,q_b+1));
end

% Compares wether x,y are equal, up to tolerance tol
function cmp = isequalRel(x,y,tol)
    cmp = (abs(x-y) <= (tol*max(abs(x),abs(y)) + eps));
end 

  %    x_u   x_b   x_v
  %    o    o    o
  %    |    |    |
  %    o    o    |
  %    |    |    |
  %    |    o    o
  %    |    |    |
  %    |    |    |
  %    o    o    |
  %    |    |    |
  %    |    o    o
  %    |    |    |
  %    |    |    |
  %    |    |    |
  %    o    o    o
function x_b = createBaseGlueGrid(x_u,x_v,tol)
	  % Check to make sure the first and last grid points lineup
	  assert(isequalRel(x_u(1),   x_v(1),   tol));
	  assert(isequalRel(x_u(end), x_v(end), tol));

	  x_b = union(x_u, x_v);
	  ind_equal = arrayfun(@isequalRel, x_b(1:end-1), x_b(2:end), ...
	                	   tol*ones(length(x_b)-1,1));
	  x_b(ind_equal) = [];
end


