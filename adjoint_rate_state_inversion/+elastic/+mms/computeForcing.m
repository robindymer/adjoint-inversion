% Accepts exact solution and returns body forcing, time derivatives and tractions at boundaries
function funcs = computeForcing(u1, u2, C, rho, x, y, t, dim)

	default_arg('dim', 2);

	u1t = diff(u1, t);
	u2t = diff(u2, t);

	u1tt = diff(u1t, t);
	u2tt = diff(u2t, t);

	sigma = cell(dim,dim);
	force = cell(dim,1);
	u = {u1, u2};
	ut = {u1t, u2t};

	for j = 1:dim
		force{j} = 0;
		for i = 1:dim
			sigma{i,j} = 0;
			for k = 1:dim
				for l = 1:dim
					sigma{i,j} = sigma{i,j} + C{i,j,k,l}*dd(k, u{l}, x, y);
				end
			end
			force{j} = force{j} + dd(i, sigma{i,j}, x, y);
		end

	end

	f1 = rho*u1tt - force{1};
	f2 = rho*u2tt - force{2};

	% Normal. OBS! Only for circular domain
	nx = symfun(x/sqrt(x^2 + y^2), [t, x, y]);
	ny = symfun(y/sqrt(x^2 + y^2), [t, x, y]);
	n = {nx, ny};
	nt = {-ny, nx};

	tau = cell(dim,1);
	for j = 1:dim
		tau{j} = 0;
		for i = 1:dim
			tau{j} = tau{j} + sigma{i,j}*n{i};
		end
	end

	% Normal-tangential
	tau_n = 0;
	tau_t = 0;
	u_n = 0;
	u_t = 0;
	ut_n = 0;
	ut_t = 0;
	u = {u1, u2};
	for i = 1:dim
		tau_n = tau_n + n{i}*tau{i};
		tau_t = tau_t + nt{i}*tau{i};
		u_n = u_n + n{i}*u{i};
		u_t = u_t + nt{i}*u{i};
		ut_n = ut_n + n{i}*ut{i};
		ut_t = ut_t + nt{i}*ut{i};
	end

	rho = symfun(rho(0,x,y), [x, y]);

	funcs = struct;
	for i = 1:dim
		for j = 1:dim
			for k = 1:dim
				for l = 1:dim
					C{i,j,k,l} = symfun(C{i,j,k,l}(0,x,y), [x, y]);
					funcs.C{i,j,k,l} = matlabFunctionSizePreserving(C{i,j,k,l});
				end
			end
		end
	end

	funcs.rho = matlabFunctionSizePreserving(rho);
	funcs.u1 = matlabFunctionSizePreserving(u1);
	funcs.u2 = matlabFunctionSizePreserving(u2);
	funcs.u1t = matlabFunctionSizePreserving(u1t);
	funcs.u2t = matlabFunctionSizePreserving(u2t);
	funcs.f1 = matlabFunctionSizePreserving(f1);
	funcs.f2 = matlabFunctionSizePreserving(f2);
	funcs.tau1 = matlabFunctionSizePreserving(tau{1});
	funcs.tau2 = matlabFunctionSizePreserving(tau{2});

	funcs.tau_n = matlabFunctionSizePreserving(tau_n);
	funcs.tau_t = matlabFunctionSizePreserving(tau_t);

	funcs.u_n = matlabFunctionSizePreserving(u_n);
	funcs.u_t = matlabFunctionSizePreserving(u_t);
	funcs.ut_t = matlabFunctionSizePreserving(ut_t);
	funcs.ut_n = matlabFunctionSizePreserving(ut_n);

	funcs.sym = struct;
	funcs.sym.rho = rho;
	funcs.sym.u1 = u1;
	funcs.sym.u2 = u2;
	funcs.sym.u1t = u1t;
	funcs.sym.u2t = u2t;
	funcs.sym.f1 = f1;
	funcs.sym.f2 = f2;
	funcs.sym.tau1 = tau{1};
	funcs.sym.tau2 = tau{2};

	funcs.sym.tau_n = tau_n;
	funcs.sym.tau_t = tau_t;

	funcs.sym.u_n = u_n;
	funcs.sym.u_t = u_t;

	funcs.sym.ut_t = ut_t;
	funcs.sym.ut_n = ut_n;

end


function u_i = dd(i, u, x, y)
	switch i
	case 1
		u_i = diff(u, x);
	case 2
		u_i = diff(u, y);
	end
end


