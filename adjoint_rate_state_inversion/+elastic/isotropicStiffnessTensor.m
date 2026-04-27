function [C, C_lambda, C_mu] = isotropicStiffnessTensor(lambda, mu, dim)

	default_arg('dim', 2);

	d = @kroneckerDelta;

	% Double
	if isa(lambda, 'double') || isa(lambda, 'symfun')
		C = cell(dim,dim,dim,dim);
		C_lambda = cell(dim,dim,dim,dim);
		C_mu = cell(dim,dim,dim,dim);
		for i = 1:dim
			for j = 1:dim
				for k = 1:dim
					for l = 1:dim
						C{i,j,k,l} = lambda*d(i,j)*d(k,l) + ...
									+ mu*(d(i,k)*d(j,l) + d(i,l)*d(j,k));

						C_lambda{i,j,k,l} = lambda*d(i,j)*d(k,l);

						C_mu{i,j,k,l} = mu*(d(i,k)*d(j,l) + d(i,l)*d(j,k));
					end
				end
			end
		end
	end

	% Single-block, function handle
	if isa(lambda, 'function_handle')
		C = cell(dim,dim,dim,dim);
		C_lambda = cell(dim,dim,dim,dim);
		C_mu = cell(dim,dim,dim,dim);
		for i = 1:dim
			for j = 1:dim
				for k = 1:dim
					for l = 1:dim
						C{i,j,k,l} = @(x,y) lambda(x,y)*d(i,j)*d(k,l) + ...
									+ mu(x,y)*(d(i,k)*d(j,l) + d(i,l)*d(j,k));

						C_lambda{i,j,k,l} = @(x,y) lambda(x,y)*d(i,j)*d(k,l);

						C_mu{i,j,k,l} = @(x,y) mu(x,y)*(d(i,k)*d(j,l) + d(i,l)*d(j,k));
					end
				end
			end
		end
	end

	% Multiblock, cell array of function handles or doubles
	if isa(lambda, 'cell')
		nBlocks = numel(lambda);
		C = cell(nBlocks, 1);
		C_lambda = cell(nBlocks, 1);
		C_mu = cell(nBlocks, 1);
		for b = 1:nBlocks
			C{b} = cell(dim,dim,dim,dim);
			C_lambda{b} = cell(dim,dim,dim,dim);
			C_mu{b} = cell(dim,dim,dim,dim);
			for i = 1:dim
				for j = 1:dim
					for k = 1:dim
						for l = 1:dim

							if isa(lambda{b}, 'function_handle')
								C{b}{i,j,k,l} = @(x,y) lambda{b}(x,y)*d(i,j)*d(k,l) + ...
										+ mu{b}(x,y)*(d(i,k)*d(j,l) + d(i,l)*d(j,k));

								C_lambda{b}{i,j,k,l} = @(x,y) lambda{b}(x,y)*d(i,j)*d(k,l);

								C_mu{b}{i,j,k,l} = @(x,y) mu{b}(x,y)*(d(i,k)*d(j,l) + d(i,l)*d(j,k));
							else
								C{b}{i,j,k,l} = lambda{b}*d(i,j)*d(k,l) + ...
										+ mu{b}*(d(i,k)*d(j,l) + d(i,l)*d(j,k));

								C_lambda{b}{i,j,k,l} = lambda{b}*d(i,j)*d(k,l);

								C_mu{b}{i,j,k,l} = mu{b}*(d(i,k)*d(j,l) + d(i,l)*d(j,k));
							end
						end
					end
				end
			end
		end
	end


end