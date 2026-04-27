function [RHOJi, PHI, ops] = setupFromDiscr(discr)

	ops = struct;
	ops.order = discr.order;
	ops.E = discr.Ecomp;
	N = discr.grid.grids{1}.size();
	ops.Nx = N(1);
	ops.Ny = N(2);
	ops.dx = 1/(N(1)-1);
	ops.dy = 1/(N(2)-1);

	PHI = discr.diffOp.diffOps{1}.refObj.C;
	RHO = full(diag(discr.diffOp.diffOps{1}.RHO));
	J = full(diag(discr.diffOp.diffOps{1}.J));

	RHOi = 1./RHO;
	Ji = 1./J;
	RHOJi = RHOi.*Ji;
	dim = 2;

	for i = 1:dim
		for j = 1:dim
			for k = 1:dim
				for l = 1:dim
					PHI{i,j,k,l} = full(diag(PHI{i,j,k,l}));
				end
			end
		end
	end

end