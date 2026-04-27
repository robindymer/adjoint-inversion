function [U_t, state] = quasiStatic(t, U, ops, friction, preStress, flowLaw, mms, tectonic, output, state)

	%---- Components of U --
	% gamma:  	Viscous strains
	% psi:  	Fault state variable
	% delta:	Slip on the fault

	% --- Other variables ---
	% u: 		Displacements
	% f:		Friction coefficient

	gamma = (ops.Egamma)'*U;
	psi = (ops.Epsi)'*U;
	delta = (ops.Edelta)'*U;

	%--- Solve mechanical equilibrium for displacement ---

	% Account for viscous strains
	RHS = ops.Mech_gamma*gamma;

	% Account for slip on fault
	RHS = RHS + ops.slipInserter_u*delta;

	% Add forcing if doing mms
	if ~isempty(mms)
		% Shear and normal traction on the fault may be included here
		RHS = RHS + mms.mechForcing(t);
	end

	% Add tectonic plate movement
	if ~isempty(tectonic)
		for i = 1:numel(tectonic)
			displ = tectonic{i}.u0 + tectonic{i}.v*t;
			RHS = RHS + ops.tectonicInserter_u{i} * displ;
		end
	end

	%--- Solve linear system ------
	switch ops.solver

	% Use LU factorization
	case 'lu'
		u = elastic.helpers.solveWithLU(ops.Mech_u_factorized, RHS);

	% Preconditioned conjugate gradient
	case 'pcg'

		pcgTol = 1e-10;

		% Initial guess. Extrapolate from previous time levels if data available.
		if ~isempty(state) && isfield(state, 'unm1') && ~isempty(state.unm1)
			un = state.un;
			unm1 = state.unm1;
			tnp1 = t;
			tn = state.tn;
			tnm1 = state.tnm1;
			u0 = un + (tnp1-tn)/(tn-tnm1)*(un-unm1);
		elseif ~isempty(state) && isfield(state, 'u')
			u0 = state.u;
		else
			u0 = [];
		end

		% Modified RHS, to match symmetrized matrix.
		b = ops.H_u*RHS;

		% Solve
		switch ops.elasticOperator
		case 'mex'
			[u, flag] = pcg(ops.mexFun, b, pcgTol, length(ops.D), ops.precond, ops.precond', u0);
		otherwise
			[u, flag] = pcg(ops.D, b, pcgTol, length(ops.D), ops.precond, ops.precond', u0);
		end
		if flag ~=0
			disp('Some problem with pcg');
		end

		% Update state
		if ~isempty(state)
			state.u = u;
		end
	end

	deltaAfterSolve = ops.slipExtractor*u;
	%-----------------------------------------------------

	%--- Compute slip velocity by solving force balance on fault -------

	% Extract shear stress perturbation on the fault
	tau1 = (ops.tauShearFault1_u)'*u + (ops.tauShearFault1_gamma)'*gamma;
	tau2 = (ops.tauShearFault2_u)'*u + (ops.tauShearFault2_gamma)'*gamma;

	% Average shear stress from the two sides of the fault. Minus sign by convention.
	tau = -1/2*(tau1 + tau2);
	if ~isempty(mms)
		tau = tau + mms.frictionTractionForcing(t);
	end

	% Add prestress
	shearStress = preStress.shear + tau;

	% SAT contribution to effective shear stress
	effShearStress = shearStress + ops.Zt*(deltaAfterSolve - delta);

	% Extract normal stress perturbation on the fault
	sigma1 = (ops.tauNormalFault1_u)'*u + (ops.tauNormalFault1_gamma)'*gamma;
	sigma2 = (ops.tauNormalFault2_u)'*u + (ops.tauNormalFault2_gamma)'*gamma;

	% Average normal stress (positive in compression) from the two sides of the fault
	sigma = -1/2*(sigma1 + sigma2);

	% Add prestress and perturbation and compute compressive normal stress
	normalStress = max(preStress.sigma0 + sigma, 0);

	% ----- Solve force balance for slip velocity ------------
	if ~isempty(mms)
		mms.t = t;
	end
	V = elastic.viscofriction.slipVelocityBisection(friction, effShearStress, normalStress, psi, mms);
	delta_t = V;

	% Compute resulting friction
	f = friction.coefficient(V, psi, friction);
	if ~isempty(mms)
		f = f + mms.frictionForcing(t);
	end
	%---------------------------------

	%--- Compute stresses everywhere ----
	dim = 2;
	sigma = cell(dim, dim);
	for i = 1:dim
		for j = 1:dim
			sigma{i,j} = ops.sigma_u{i,j}*u + ops.sigma_gamma{i,j}*gamma;
		end
	end
	%----------------------------------------------------


	%---- Evolve viscous strains -----------------------

	% Compute effective viscosity eta
	etaInvScalar = flowLaw.eta.nonlinFun(sigma);
	etaInv = kron(etaInvScalar, ones(4,1));

	% Flow law
	gamma_t = etaInv .* (ops.Dgu*u + ops.Dgg*gamma);

	% Penalty contribution from slip on fault
	gamma_t = gamma_t + etaInv .* ops.slipInserter_gamma*delta;

	% Add tectonic plate movement
	if ~isempty(tectonic)
		for i = 1:numel(tectonic)
			displ = tectonic{i}.u0 + tectonic{i}.v*t;
			gamma_t = gamma_t + etaInv .* ops.tectonicInserter_gamma{i} * displ;
	end

	% Add forcing if doing mms
	if ~isempty(mms)
		gamma_t = gamma_t + mms.flowBodyForcing(t);
		gamma_t = gamma_t + etaInv .* mms.flowBoundaryForcing(t);
	end
	%----------------------------------------------------

	% Evolve the state variable
	psi_t = elastic.friction.stateEvolution(V, psi, f, friction);
	if ~isempty(mms)
		psi_t = psi_t + mms.stateForcing(t);
	end

	% Reassemble rates
	U_t = ops.Egamma*gamma_t + ops.Epsi*psi_t + ops.Edelta*delta_t;

	%------ Write fields to disk -------------
	if ~isempty(output)

		% Load MAT-file with count if it exists
		if helpers.fileExists([output.dir, 'countArray.mat'])
			s = load([output.dir, 'countArray.mat']);
			countArray = s.countArray;

		% If this is the first time, initialize saving structures
		else
			countArray = ones(1, numel(output.fields));
			for i = 1:numel(output.fields)
				dat = struct;
				dat.t = [];
				dat.field = [];
				fileName = [output.dir, output.fields{i}];
				save(fileName, 'dat');
			end
		end

		for i = 1:numel(output.fields)
			fieldName = output.fields{i};
			count = countArray(i);

			% If it is time to write field i, then write that field to disk
			if output.writeFlags{i}

				% Load files
				fileName = [output.dir, fieldName];
				s = load([fileName, '.mat']);
				dat = s.dat;

				switch fieldName
				case 'slipVelocity'
					field = max(abs(delta_t));
				case 'shearStressMax'
					field = max(abs(effShearStress));
				case 'shearStressMean'
					field = mean(abs(effShearStress));
				case 'psiMin'
					field = min(psi);
				case 'psiMax'
					field = max(psi);
				case 'slipProfile'
					field = delta;
				case 'slipVelocityProfile'
					field = V;
				case 'shearStressProfile'
					field = shearStress;
				case 'normalStressProfile'
					field = normalStress;
				case 'psiProfile'
					field = psi;
				case 'surfaceDisplacementX'
					field = (ops.surfaceDisplX)'*u;
				case 'surfaceDisplacementY'
					field = (ops.surfaceDisplY)'*u;
				end

				% Store fault coordinates if required
				if contains(fieldName, 'Profile')
					dat.xFault = ops.xFault;
				end

				% Store surface coordinates if required
				if contains(fileName, 'surface')
					dat.xSurface = ops.xSurface;
				end

				% Append current time and field
				dat.t = [dat.t, t];
				dat.field = [dat.field, field];

				% Write to file
				save(fileName, 'dat');

				% Update count and write to file
				countArray(i) = countArray(i) + 1;
				save([output.dir, 'countArray.mat'], 'countArray');
			end
		end
	end
	%-----------------------------------------

end