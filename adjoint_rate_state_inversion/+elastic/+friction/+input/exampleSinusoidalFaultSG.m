function iS = exampleSinusoidalFaultSG()

	% Start with default options
	iS = elastic.friction.input.default();

	% ----- Domain and boundary conditions -------
	% Stretching to cluster points near fault
	stretchDist = 1/4;
	stretchRatio = 2;

	def = domains.SinusoidalFaultSG(stretchDist, stretchRatio);
	g = def.getGrid(11);
	bgs = def.getBoundaryGroups(g);

	% --- Traction BC on surface ----
	bcFree1.boundary = bgs.surface;
	bcFree1.type = {1, 't'};

	bcFree2.boundary = bgs.surface;
	bcFree2.type = {2, 't'};

	% Displacement BC on remaining boundaries
	bcD1.boundary = bgs.farfield;
	bcD1.type = {1, 'd'};

	bcD2.boundary = bgs.farfield;
	bcD2.type = {2, 'd'};

	bcBottom1.boundary = bgs.bottom;
	bcBottom1.type = {1, 'd'};

	bcBottom2.boundary = bgs.bottom;
	bcBottom2.type = {2, 'd'};

	bc = {bcFree1, bcFree2, bcD1, bcD2, bcBottom1, bcBottom2};

	% Fault boundary groups
	faultBoundaryGroups = {bgs.faultBottom, bgs.faultTop};

	type.type = 'standard';
	intfTypes = multiblock.setAllInterfaceTypes(g, type);
	type.type = 'frictionalFault';

	% DOI fault
	intfTypes{1,2} = type;
	intfTypes{2,1} = type;

	% Supergrid fault
	intfTypes{3,4} = type;
	intfTypes{4,3} = type;

	intfTypes{5,6} = type;
	intfTypes{6,5} = type;

	iS.domain.intfTypes = intfTypes;
	iS.domain = elastic.friction.input.domainStruct(iS.domain, def, bc, faultBoundaryGroups);
	%------------------------------------------------

	% ---- Plotting options ----------------------
	iS.plotting.xlims = 1.2*[-def.DOIdef.L/2, def.DOIdef.L/2];
	iS.plotting.surfaceBoundaryGroup = bgs.surface;

	% ----- Friction parameters ------------------
	iS.friction.a = 0.01; 			% Direct effect [ ]
	iS.friction.b = 0.02;			% State evolution effect [ ]
	iS.friction.L = 0.25;			% State evolution distance [ m ]
	iS.friction.V0 = 1e-6;			% Ref velocity [ m/s ]
	iS.friction.f0 = 0.6;			% Ref friction coeff for steady sliding [ ]
	iS.friction.psi0 = 0.4367;		% Initial state variable

	iS.friction.coefficient = 'standard';
	iS.friction.coefficientSteadyState = [];

	iS.friction.law = 'aging';
	%----------------------------------------------

	% ---- Pre stress ----------------------------
	x0 = 0;
	s = 1/10*1e3;
	A = 100*1e6;
	bump = @(x,y) exp(-(x-x0).^2/(2*s^2) );

	iS.preStress.sigma0 = @(x,y) 126*1e6 + 0*x;
	iS.preStress.shear = @(x,y) 0.3*iS.preStress.sigma0(x,y) + A*bump(x,y);
	% ---------------------------------------------

	% ------ Material parameters -----------------
	rho = 2e3;
	vs = 2e3;
	vp = 6e3;

	iS.material = elastic.friction.input.materialStruct(rho, vs, vp);

	% ------ Time stepper ----------------------
	iS.timeStepper.type = 'RK4';

	% ------ Spatial operator ------------------
	iS.solver.mex = false;
	iS.solver.check_ev_flag = false;
end