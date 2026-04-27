% C: cell array of vectors
% rhoJI: Vector of (rho*J)^{-1}
function out = elasticOperatorOMP(U, rhoJI, C, ops, numThreads)

	u = cell(2,1);
	u{1} = U(1:2:end-1);
	u{2} = U(2:2:end);

	dd = cell(2,1);
	dd{1} = @(u) elastic.mex.D1xOMP(u, ops.dx, ops.Nx, ops.Ny, ops.order, numThreads);
	dd{2} = @(u) elastic.mex.D1yOMP(u, ops.dy, ops.Nx, ops.Ny, ops.order, numThreads);

	dd2 = cell(2,1);
	dd2{1} = @(u, b) elastic.mex.D2xVariableOMP(u, b, ops.dx, ops.Nx, ops.Ny, ops.order, numThreads);
	dd2{2} = @(u, b) elastic.mex.D2yVariableOMP(u, b, ops.dy, ops.Nx, ops.Ny, ops.order, numThreads);

	dd2Combined = cell(2,1);
	dd2Combined{1} = @(u, v, b) elastic.mex.D2xVariableDoubleOMP(u, v, b, ops.dx, ops.Nx, ops.Ny, ops.order, numThreads);
	dd2Combined{2} = @(u, v, b) elastic.mex.D2yVariableDoubleOMP(u, v, b, ops.dy, ops.Nx, ops.Ny, ops.order, numThreads);

	d = {struct, struct};
	d{1}.u{1} = dd{1}(u{1});
	d{2}.u{1} = dd{2}(u{1});
	d{1}.u{2} = dd{1}(u{2});
	d{2}.u{2} = dd{2}(u{2});

	dim = 2;
	Du = cell(dim, 1);
	Du{1} = zeros(ops.Nx*ops.Ny, 1);
	Du{2} = zeros(ops.Nx*ops.Ny, 1);

	% Compute d_i C_{ijkl} d_k u_l
	for i = 1:dim
		for j = 1:dim
			for k = 1:dim
				if i == k
					for l = j
						Du{j} = Du{j} + dd2{i}(u{l}, C{i,j,k,l});
					end
				else
					Du{j} = Du{j} + dd{i}(C{i,j,k,1}.*d{k}.u{1} + C{i,j,k,2}.*d{k}.u{2});
				end
			end
		end
	end
	for i = 1:dim
		[du1, du2] = dd2Combined{i}(u{2}, u{1}, C{i,1,i,2});
		Du{1} = Du{1} + du1;
		Du{2} = Du{2} + du2;
	end

	out = zeros(dim*ops.Nx*ops.Ny, 1);
	out(1:2:end-1) = rhoJI.*Du{1};
	out(2:2:end) = rhoJI.*Du{2};

end