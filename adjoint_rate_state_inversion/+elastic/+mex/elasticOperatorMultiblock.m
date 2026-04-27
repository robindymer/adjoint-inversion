% C: cell array of vectors
% rhoJI: Vector of (rho*J)^{-1}
function out = elasticOperatorMultiblock(U, g, rhoJI, C, ops)

	U = g.splitFunc(U);
	out = cell(g.nBlocks(), 1);
	for i = 1:g.nBlocks()
		U{i} = full(U{i});
		out{i} = elastic.mex.elasticOperator(U{i}, rhoJI{i}, C{i}, ops{i});
	end
	out = cell2mat(out);

end