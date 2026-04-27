function globalOp = multiblockOperator(g, diffOp, name)

nBlocks = g.nBlocks;

globalOp = cell(nBlocks, nBlocks);
for i = 1:nBlocks
	localOp = eval(['diffOp.diffOps{i}.', name]);
	globalOp{i,i} = localOp;
end

globalOp = blockmatrix.toMatrix(globalOp);
