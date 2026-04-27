function [ops, plotting] = setupOpsGPU(discr, faultBoundaryGroups, intfType, plotting, tectonic, material, timeStepper)
	plotting_default = struct;
	plotting_default.faultBoundaryGroup = [];
	plotting_default.surfaceBoundaryGroup = [];
	default_struct('plotting', plotting_default);

	ops = struct;
	ops.solver = timeStepper.solver;

	diffOp = discr.diffOp;

	helpDiffOp = multiblock.DiffOp(@scheme.Elastic2dCurvilinearAnisotropic, discr.grid, discr.order, {[], material.CFun});
	et_1 = helpDiffOp.getBoundaryOperator('et', faultBoundaryGroups{1});
	et_2 = helpDiffOp.getBoundaryOperator('et', faultBoundaryGroups{2});
	ops.slipExtractor = et_1' + et_2';

	% --- Penalty factors for traction forcing on the fault --------
	bc = {'t', 't'};
	[~, shearTractionInserter] 	= diffOp.boundary_condition(faultBoundaryGroups{1}, bc);
	ops.shearTractionInserter1 	= 1/2*(discr.Eu)'*shearTractionInserter;

	[~, shearTractionInserter] 	= diffOp.boundary_condition(faultBoundaryGroups{2}, bc);
	ops.shearTractionInserter2 	= 1/2*(discr.Eu)'*shearTractionInserter;

	bc = {'n', 't'};
	[~, normalTractionInserter] = diffOp.boundary_condition(faultBoundaryGroups{1}, bc);
	ops.normalTractionInserter1 = 1/2*(discr.Eu)'*normalTractionInserter;

	[~, normalTractionInserter] = diffOp.boundary_condition(faultBoundaryGroups{2}, bc);
	ops.normalTractionInserter2 = 1/2*(discr.Eu)'*normalTractionInserter;
	%----------------------------------------------------------------------

	% --- Penalty factors for slip on the fault --------
	intfPenalty = discr.diffOp.interfaceForcing(faultBoundaryGroups{1}, faultBoundaryGroups{2}, intfType);
	slipInserter1 = intfPenalty{2};

	intfPenalty = discr.diffOp.interfaceForcing(faultBoundaryGroups{2}, faultBoundaryGroups{1}, intfType);
	slipInserter2 = intfPenalty{2};

	ops.slipInserter_u 		= (discr.Eu)'*(slipInserter1 + slipInserter2);
	ops.slipInserter_gamma 	= (discr.Egamma)'*(slipInserter1 + slipInserter2);
	%----------------------------------------------------------------------

	% --- Get factor in SAT appearing in force balance, yielding "effective shear stress"
	nBlocks = length(faultBoundaryGroups{1});
	Zt = cell(nBlocks, nBlocks);
	for i = 1:nBlocks
		B1 = faultBoundaryGroups{1}{i}{1};
		b1 = faultBoundaryGroups{1}{i}{2};

		B2 = faultBoundaryGroups{2}{i}{1};
		b2 = faultBoundaryGroups{2}{i}{2};

		dO1 = helpDiffOp.diffOps{B1};
		dO2 = helpDiffOp.diffOps{B2};
		[~, ~, ~, Zt_local] = dO1.interfaceNormalTangential(b1, dO2, b2, intfType);
		Zt{i,i} = Zt_local;
	end
	Zt = blockmatrix.toMatrix(Zt);
	ops.Zt = -Zt;
	clear helpDiffOp;

	% --- Penalty factors for tectonic plate movement ---
	if ~isempty(tectonic)

		% Tectonic forcing on boundaries
		if isfield(tectonic, 'bc')
			[~, tectonicInserter] = diffOp.boundary_condition(tectonic.bc.boundary, tectonic.bc.type);

			% Both boundaries are sliding with half the plate rate
			tectonicInserter = 1/2*tectonicInserter;
		end

		% Tectonic forcing on interface
		if isfield(tectonic, 'interface')
			intfPenalty = discr.diffOp.interfaceForcing(tectonic.interface.boundaryGroups{1}, tectonic.interface.boundaryGroups{2}, intfType);
			tectonicInserter1 = intfPenalty{2};

			intfPenalty = discr.diffOp.interfaceForcing(tectonic.interface.boundaryGroups{2}, tectonic.interface.boundaryGroups{1}, intfType);
			tectonicInserter2 = intfPenalty{2};

			% Check if there are tectonic boundaries too.
			if isfield(tectonic, 'bc')
				tectonicInserter = tectonicInserter + tectonicInserter1 + tectonicInserter2;
			else
				tectonicInserter = tectonicInserter1 + tectonicInserter2;
			end
		end

		% Use that the boundary data is homogenous in space
		tectonicInserter = sum(tectonicInserter, 2);

		ops.tectonicInserter_u 		= (discr.Eu)'*tectonicInserter;
		ops.tectonicInserter_gamma 	= (discr.Egamma)'*tectonicInserter;
	end

	% --- Operators for displacement and tractions on the fault -----
	tauShearFault1 = diffOp.getBoundaryOperator('tau_t', faultBoundaryGroups{1});
	ops.tauShearFault1_u = ( tauShearFault1'*discr.Eu )';
	ops.tauShearFault1_gamma = ( tauShearFault1'*discr.Egamma )';

	tauShearFault2 = diffOp.getBoundaryOperator('tau_t', faultBoundaryGroups{2});
	ops.tauShearFault2_u = ( tauShearFault2'*discr.Eu )';
	ops.tauShearFault2_gamma = ( tauShearFault2'*discr.Egamma )';

	tauNormalFault1 = diffOp.getBoundaryOperator('tau_n', faultBoundaryGroups{1});
	ops.tauNormalFault1_u = ( tauNormalFault1'*discr.Eu )';
	ops.tauNormalFault1_gamma = ( tauNormalFault1'*discr.Egamma )';

	tauNormalFault2 = diffOp.getBoundaryOperator('tau_n', faultBoundaryGroups{2});
	ops.tauNormalFault2_u = ( tauNormalFault2'*discr.Eu )';
	ops.tauNormalFault2_gamma = ( tauNormalFault2'*discr.Egamma )';
	%------------------------------------------------------------------

	%--- Storage order: gamma - psi - delta -----
	dim = 2;
	mGamma = dim^2*discr.grid.N();
	[mFault, ~] = size(tauShearFault1');

	Igamma = speye(mGamma, mGamma);
	Ifault = speye(mFault, mFault);
	gammaZero = sparse(mFault, mGamma);
	faultZero = sparse(mGamma, mFault);

	ops.Egamma 	= cell2mat({Igamma, faultZero, faultZero})';
	ops.Epsi 	= cell2mat({gammaZero, Ifault, 0*Ifault})';
	ops.Edelta 	= cell2mat({gammaZero, 0*Ifault, Ifault})';

	ops.mGamma = mGamma;
	ops.mFault = mFault;
	%---------------------------------------------

	% --- Copy matrices from discr ----
	ops.Dgu = discr.Dgu;
    ops.Dgg = discr.Dgg;

    ops.Mech_gamma = discr.Mech_gamma;
    ops.sigma_u = discr.sigma_u;
    ops.sigma_gamma = discr.sigma_gamma;

    switch ops.solver
    case 'lu'
    	ops.Mech_u_factorized = discr.Mech_u_factorized;

    case 'pcg'
    	ops.D = discr.H_u * discr.Mech_u;
    	ops.H_u = discr.H_u;

	    tic
	    ops.precond = ichol(ops.D);
	    t = toc;
	    fprintf('Time to compute preconditioner: %4.3es \n', t);
	end
    %-----------------------------------

	% --- Plotting ---------------------------
	ops.xFault = discr.grid.getBoundary(faultBoundaryGroups{1});
	if ~isempty(plotting.surfaceBoundaryGroup)
		plotting.xSurface = discr.grid.getBoundary(plotting.surfaceBoundaryGroup);
	else
		plotting.xSurface = [];
	end
	if ~isempty(plotting.faultBoundaryGroup)
		plotting.xFault = discr.grid.getBoundary(plotting.faultBoundaryGroup);
	else
		plotting.xFault = ops.xFault;
	end
	%-------------------------------------------

	%--- Transfer ops to GPU ---------
	tic
	% precond = spdiags(diag(ops.D),0,length(ops.D),length(ops.D));
	% ops.precond = gpuArray(precond);

	% precondDiagInv = 1./diag(ops.D);
	% precondDiagInv = gpuArray(precondDiagInv);
	% ops.precondFun = @(x) precondDiagInv.*x;

	% [L, U] = ilu(ops.D);
	% ops.precond = gpuArray(L*U);

	% ops.precond = gpuArray(ops.precond);

	% precond = ops.precond*(ops.precond)';
	% ops.precond = gpuArray(precond);

	ops.slipInserter_u = gpuArray(ops.slipInserter_u);
	ops.tectonicInserter_u = gpuArray(ops.tectonicInserter_u);
	ops.H_u = gpuArray(ops.H_u);
	ops.D = gpuArray(ops.D);
	ops.slipExtractor = gpuArray(ops.slipExtractor);
	ops.tauShearFault1_u = gpuArray(ops.tauShearFault1_u);
	ops.tauShearFault2_u = gpuArray(ops.tauShearFault2_u);
	ops.tauNormalFault1_u = gpuArray(ops.tauNormalFault1_u);
	ops.tauNormalFault2_u = gpuArray(ops.tauNormalFault2_u);
	t = toc;
    fprintf('Time to transfer matrices to GPU: %4.3es \n', t);

end