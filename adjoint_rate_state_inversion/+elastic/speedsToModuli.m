function [lambda, mu] = speedsToModuli(vp, vs, rho, dim)
	default_arg('dim', 2);

	% Double
	if isa(vp, 'double') || isa(vp, 'symfun')
		mu = rho.*vs.^2;
		lambda = rho.*vp.^2 - 2*mu;
	end

	% Single-block, function handle
	if isa(vp, 'function_handle')
		mu = @(x,y) rho(x,y).*vs(x,y).^2;
		lambda = @(x,y) rho(x,y).*vp(x,y).^2 - 2*mu(x,y);
	end

	% Multi-block, function handle or double
	if isa(vp, 'cell')
		nBlocks = numel(vp);
		lambda = cell(nBlocks, 1);
		mu = cell(nBlocks, 1);

		for i = 1:nBlocks
			if isa(vp{i}, 'function_handle')
				mu{i} = @(x,y) rho(x,y).*vs(x,y).^2;
				lambda{i} = @(x,y) rho(x,y).*vp(x,y).^2 - 2*mu(x,y);
			else
				mu{i} = rho.*vs.^2;
				lambda{i} = rho.*vp.^2 - 2*mu;
			end
		end
	end

end