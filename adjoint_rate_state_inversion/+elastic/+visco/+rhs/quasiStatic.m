function gamma_t = quasiStatic(t, gamma, ops, etaFun, mms)

	% u: 		Displacements
	% gamma:  	Viscous strains
	RHS = ops.Mech_gamma*gamma;
	if ~isempty(mms)
		RHS = RHS + mms.mechForcing(t);
	end
	u = elastic.helpers.solveWithLU(ops.Mech_u_factorized, RHS);

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
	gamma_t = etaInv .* (ops.Dgu*u + ops.Dgg*gamma);

	if ~isempty(mms)
		gamma_t = gamma_t + mms.flowBodyForcing(t);
		gamma_t = gamma_t + etaInv .* mms.flowBoundaryForcing(t);
	end
end