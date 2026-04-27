% C: cell array of vectors
% rhoJI: Vector of (rho*J)^{-1}
function out = elasticOperatorOMPOpt(U, rhoJI, C, ops, numThreads)

	u1 = U(1:2:end-1);
	u2 = U(2:2:end);

	%------- Compute d_i C_{ijkl} d_k u_l ---------

	% Compute all strains
	[d1u1, d2u1] = elastic.mex.D1x_and_D1y_OMP(u1, ops.dx, ops.dy, ops.Nx, ops.Ny, ops.order, numThreads);
	[d1u2, d2u2] = elastic.mex.D1x_and_D1y_OMP(u2, ops.dx, ops.dy, ops.Nx, ops.Ny, ops.order, numThreads);

	% Some of the narrow-stencil second derivatives
	Du1 = elastic.mex.D2xVariablePlusD2yVariableOMP(u1, C{1,1,1,1}, C{2,1,2,1}, ops.dx, ops.dy, ops.Nx, ops.Ny, ops.order, numThreads);
	Du2 = elastic.mex.D2xVariablePlusD2yVariableOMP(u2, C{1,2,1,2}, C{2,2,2,2}, ops.dx, ops.dy, ops.Nx, ops.Ny, ops.order, numThreads);

	% Apply remaining D1 operator for all mixed derivatives
	[du1_1, du2_1] = elastic.mex.D1LeftOMP(d2u1, d2u2, d1u1, d1u2, C{1,1,2,1}, C{1,1,2,2}, C{1,2,2,1}, C{1,2,2,2}, ops.dx, ops.dy, ops.Nx, ops.Ny, ops.order, numThreads);

	% Remaining narrow-stencil applications
	[du1, du2] = elastic.mex.D2VariableCombinedXPlusYOMP(u2, u1, C{1,1,1,2}, C{2,1,2,2}, ops.dx, ops.dy, ops.Nx, ops.Ny, ops.order, numThreads);
	%------------------------------------------------

	dim = 2;
	out = zeros(dim*ops.Nx*ops.Ny, 1);

	out(1:2:end-1) = rhoJI.*(Du1 + du1 + du1_1);
	out(2:2:end) = rhoJI.*(Du2 + du2 + du2_1);

end