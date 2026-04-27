function gf = evalOnStaggered(g, f, gridGroup, subGrids)

	default_arg('gridGroup', 1);

	gf = [];

	if isa(g, 'multiblock.Grid')

		nSubGrids = numel(g.grids{1}.gridGroups{gridGroup});
		default_arg('subGrids', 1:nSubGrids );

		for i = 1:g.nBlocks()
			for j = subGrids
				gf = [gf; grid.evalOn(g.grids{i}.gridGroups{gridGroup}{j}, f)];
			end
		end

	else

		nSubGrids = numel(g.gridGroups{gridGroup});
		default_arg('subGrids', 1:nSubGrids );

		for j = subGrids
			gf = [gf; grid.evalOn(g.grids{i}.gridGroups{gridGroup}{j}, f)];
		end

	end

end