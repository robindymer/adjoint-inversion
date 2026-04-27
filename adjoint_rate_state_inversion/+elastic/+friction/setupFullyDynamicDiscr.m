function [iS, ops, discr] = setupFullyDynamicDiscr(inputFile, m, T, order)
	% iS: inputStruct
	%--------------------------------------------------

	% Call input file to get input struct
	iS = inputFile();
	rho = iS.material.rhoFun;
	C = iS.material.CFun;

	if isempty(iS.mms)
		F = [];
	else
		F = iS.mms.F;
	end

	% ---- Domain ----------------
	g = iS.domain.def.getGrid(m);
	bc = iS.domain.bc;
	faultBoundaryGroups = iS.domain.faultBoundaryGroups;

	if isempty(iS.domain.intfTypes)
		type = struct;
		type.type = 'frictionalFault';
		iS.domain.intfTypes = multiblock.setAllInterfaceTypes(g, type);
	end

	% Create discretization object
	if iS.solver.mex
		hollow = true;
		discr = elastic.elasticAnisotropicDiscrCurveMultiblock(g, order, C, rho, F, bc,[],[],[], iS.domain.intfTypes, hollow);
		% Setup operators needed to run the mex stencil code
		[RHOJi, PHI, ops] = elastic.mex.setupFromMultiblockDiscr(discr);
		mexOps = struct;
		mexOps.RHOJi = RHOJi;
		mexOps.PHI = PHI;
		mexOps.ops = ops;
		mexOps.g = g;
	else
		discr = elastic.elasticAnisotropicDiscrCurveMultiblock(g, order, C, rho, F, bc,[],[],[], iS.domain.intfTypes);
		mexOps = [];
	end

	% Create structures for fault data, plotting, etc.
	[ops, iS.plotting] = elastic.friction.setupOps(discr, faultBoundaryGroups, iS.plotting);
	ops.mex = mexOps;


	% --------- Evaluate fields on the fault ---------
	% Evaluate initial stresses on fault
	iS.preStress.shear = elastic.helpers.evalOnLine(ops.xFault, iS.preStress.shear);
	iS.preStress.sigma0 = elastic.helpers.evalOnLine(ops.xFault, iS.preStress.sigma0);

	iS.friction.a = elastic.helpers.evalOnLine(ops.xFault, iS.friction.a);
	iS.friction.b = elastic.helpers.evalOnLine(ops.xFault, iS.friction.b);
	% -------------------------------------------------

	% MMS forcing for friction coefficient and state evolution
	if ~isempty(iS.mms)
		iS.mms.forcing = discr.S_cont;

		iS.mms.frictionForcing = @(t) elastic.helpers.evalOnLine(ops.xFault, @(x,y) iS.mms.frictionForcing(t,x,y));
		iS.mms.stateForcing = @(t) elastic.helpers.evalOnLine(ops.xFault, @(x,y) iS.mms.stateForcing(t,x,y));

		iS.mms.tractionForcing1 = @(t) elastic.helpers.evalOnLine(ops.xFault, @(x,y) iS.mms.tractionForcing1(t,x,y));
		iS.mms.tractionForcing2 = @(t) elastic.helpers.evalOnLine(ops.xFault, @(x,y) iS.mms.tractionForcing2(t,x,y));

		iS.mms.normalTractionForcing1 = @(t) elastic.helpers.evalOnLine(ops.xFault, @(x,y) iS.mms.normalTractionForcing1(t,x,y));
		iS.mms.normalTractionForcing2 = @(t) elastic.helpers.evalOnLine(ops.xFault, @(x,y) iS.mms.normalTractionForcing2(t,x,y));
	end

end