% C: cell array of vectors
% rhoJI: Vector of (rho*J)^{-1}
function out = elasticOperatorMultiblockOMP(U, g, rhoJI, C, ops, numThreads)

	U = g.splitFunc(U);
	out = cell(g.nBlocks(), 1);
	for i = 1:g.nBlocks()
		U{i} = full(U{i});
		% out{i} = elastic.mex.elasticOperatorOMP(U{i}, rhoJI{i}, C{i}, ops{i}, numThreads);
		out{i} = elastic.mex.elasticOperatorOMPOpt(U{i}, rhoJI{i}, C{i}, ops{i}, numThreads);
	end
	out = cell2mat(out);

end