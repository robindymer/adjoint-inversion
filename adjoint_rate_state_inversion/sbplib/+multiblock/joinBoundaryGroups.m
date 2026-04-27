% bg = joinBoundaryGroups(bg1, bg2, bg3, ...)
function boundaryGroup = joinBoundaryGroups(varargin)

	boundaryIdentifiers = [];
	k = 1;
	for i = 1:nargin
		bg = varargin{i};
		for j = 1:length(bg)
			boundaryIdentifiers{k} = bg{j};
			k = k+1;
		end
	end

	boundaryGroup = multiblock.BoundaryGroup(boundaryIdentifiers);
end