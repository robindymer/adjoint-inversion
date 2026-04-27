function [preFactor, C, ops] = setupFromMultiblockDiscr(discr, symmetrize)

	default_arg('symmetrize', false);

	g = discr.grid;
	nBlocks = g.nBlocks();

	ops = cell(nBlocks, 1);
	preFactor = cell(nBlocks, 1);
	C = cell(nBlocks, 1);


	for b = 1:nBlocks

		ops{b} = struct;
		ops{b}.order = discr.order;
		ops{b}.E = discr.diffOp.diffOps{b}.E;
		N = g.grids{b}.size();
		ops{b}.Nx = N(1);
		ops{b}.Ny = N(2);
		ops{b}.dx = 1/(N(1)-1);
		ops{b}.dy = 1/(N(2)-1);

		PHI = discr.diffOp.diffOps{b}.refObj.C;
		RHO = full(diag(discr.diffOp.diffOps{b}.RHO));
		J = full(diag(discr.diffOp.diffOps{b}.J));
		H = full(diag(discr.diffOp.diffOps{b}.H));

		if symmetrize
			% Multiply by H for reference domain
			preFactor{b} = H./J;
		else
			RHOi = 1./RHO;
			Ji = 1./J;
			preFactor{b} = RHOi.*Ji;
		end

		dim = 2;
		C{b} = cell(dim, dim, dim, dim);
		for i = 1:dim
			for j = 1:dim
				for k = 1:dim
					for l = 1:dim
						C{b}{i,j,k,l} = full(diag(PHI{i,j,k,l}));
					end
				end
			end
		end
	end

end