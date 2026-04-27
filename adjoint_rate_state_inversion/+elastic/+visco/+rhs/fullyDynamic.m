function U_t = fullyDynamic(t, U, ops, etaFun, mms)

	% u: 		Displacements
	% v: 		Particle velocities, du/dt.
	% gamma:  	Viscous strains

	u = (ops.Eu)'*U;
	v = (ops.Ev)'*U;
	gamma = (ops.Egamma)'*U;

	% Compute stress and effective viscosity eta
	dim = 2;
	sigma = cell(dim, dim);
	for i = 1:dim
		for j = 1:dim
			sigma{i,j} = ops.sigma_u{i,j}*u + ops.sigma_gamma{i,j}*gamma;
		end
	end
	etaInvScalar = etaFun(sigma);
	etaInv = kron(etaInvScalar, ones(4,1));

	% Compute rates
	u_t = v;
	v_t = ops.Duu*u + ops.Dug*gamma;
	gamma_t = etaInv .* (ops.Dgu*u + ops.Dgg*gamma);

	U_t = ops.Eu*u_t + ops.Ev*v_t + ops.Egamma*gamma_t;

	if ~isempty(mms)
		U_t = U_t + mms.bodyForcing(t);

		factor = [ones(size(etaInv)); etaInv];
		U_t = U_t + factor.*mms.boundaryForcing(t);
	end
end