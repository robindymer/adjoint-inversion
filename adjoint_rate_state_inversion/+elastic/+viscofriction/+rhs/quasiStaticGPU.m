function [U_t, state] = quasiStaticGPU(t, U, ops, friction, preStress, flowLaw, mms, tectonic, output, state)

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
	deltaGPU = gpuArray(delta);

	% Account for slip on fault
	RHS = ops.slipInserter_u*deltaGPU;

	% Add forcing if doing mms
	if ~isempty(mms)
		% Shear and normal traction on the fault may be included here
		RHS = RHS + mms.mechForcing(t);
	end

	% Add tectonic plate movement
	if ~isempty(tectonic)
		displ = tectonic.u0 + tectonic.v*t;
		RHS = RHS + ops.tectonicInserter_u * gpuArray(displ);
	end

	%--- Solve linear system ------
	switch ops.solver

	% Use LU factorization
	case 'lu'
		u = elastic.helpers.solveWithLU(ops.Mech_u_factorized, RHS);

	% Preconditioned conjugate gradient
	case 'pcg'

		% Initial guess. Extrapolate from previous time levels if data available.
		if ~isempty(state) && isfield(state, 'unm1') && ~isempty(state.unm1)
			un = state.un;
			unm1 = state.unm1;
			tnp1 = t;
			tn = state.tn;
			tnm1 = state.tnm1;
			u0 = un + gpuArray((tnp1-tn)/(tn-tnm1))*(un-unm1);
		elseif ~isempty(state) && isfield(state, 'u')
			u0 = state.u;
		else
			u0 = [];
		end

		% Modified RHS, to match symmetrized matrix.
		b = ops.H_u*RHS;

		% Solve
		% tic
		% precond = spdiags(diag(ops.precond),0,length(ops.D),length(ops.D));

		% DOES NOT WORK AT ALL ?!?
		% precond = ops.precond;
		% precondFun = @(x) (precond')\(precond\x);

		% [u, flag] = pcg(ops.D, b, [], length(ops.D), ops.precond, [], u0);
		% [u, flag] = pcg(ops.D, b, [], length(ops.D), precondFun, [], u0);
		[u, flag] = pcg(ops.D, b, [], length(ops.D), [], [], u0);
		if flag ~=0
			disp('Some problem with pcg');
		end

		% tcomp = toc;
    	% fprintf('Time to solve using pcg/gpu: %4.3es \n', tcomp);

		% Update state
		if ~isempty(state)
			state.u = u;
		end
	end

	deltaAfterSolve = ops.slipExtractor*u;

	% Extract shear stress perturbation on the fault
	tau1 = (ops.tauShearFault1_u)'*u;
	tau2 = (ops.tauShearFault2_u)'*u;

	% Average shear stress from the two sides of the fault. Minus sign by convention.
	tau = -1/2*(tau1 + tau2);

	% Extract normal stress perturbation on the fault
	sigma1 = (ops.tauNormalFault1_u)'*u;
	sigma2 = (ops.tauNormalFault2_u)'*u;

	% Average normal stress (positive in compression) from the two sides of the fault
	sigma = -1/2*(sigma1 + sigma2);

	% Transfer back to CPU
	% tic
	deltaAfterSolve = gather(deltaAfterSolve);
	tau = gather(tau);
	sigma = gather(sigma);
	% tcomp = toc;
 %    fprintf('Time to transfer from GPU: %4.3es \n', tcomp);
	%-----------------------------------------------------

	%--- Compute slip velocity by solving force balance on fault -------
	if ~isempty(mms)
		tau = tau + mms.frictionTractionForcing(t);
	end

	% Add prestress
	shearStress = preStress.shear + tau;

	% SAT contribution to effective shear stress
	effShearStress = shearStress + ops.Zt*(deltaAfterSolve - delta);

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
	% dim = 2;
	% sigma = cell(dim, dim);
	% for i = 1:dim
	% 	for j = 1:dim
	% 		sigma{i,j} = ops.sigma_u{i,j}*u + ops.sigma_gamma{i,j}*gamma;
	% 	end
	% end
	%----------------------------------------------------


	%---- Evolve viscous strains -----------------------

	% Compute effective viscosity eta
	% etaInvScalar = flowLaw.eta.nonlinFun(sigma);
	% etaInv = kron(etaInvScalar, ones(4,1));

	% Flow law
	% gamma_t = etaInv .* (ops.Dgu*u + ops.Dgg*gamma);
	gamma_t = 0*gamma;

	% Penalty contribution from slip on fault
	% gamma_t = gamma_t + etaInv .* ops.slipInserter_gamma*delta;

	% Add tectonic plate movement
	if ~isempty(tectonic)
		displ = tectonic.u0 + tectonic.v*t;
		% gamma_t = gamma_t + etaInv .* ops.tectonicInserter_gamma * displ;
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
					field = max(abs(shearStress));
				case 'shearStressMean'
					field = mean(abs(shearStress));
				case 'psiMin'
					field = min(psi);
				case 'psiMax'
					field = max(psi);
				case 'slipProfile'
					field = delta;
					dat.xFault = ops.xFault;
				case 'shearStressProfile'
					field = shearStress;
					dat.xFault = ops.xFault;
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