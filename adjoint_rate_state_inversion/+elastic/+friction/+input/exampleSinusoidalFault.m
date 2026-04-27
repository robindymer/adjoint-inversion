function iS = exampleSinusoidalFault()
	% iS: inputStruct
	%--------------------------------------------------

	% Start with default options
	iS = elastic.friction.input.default();

	% ----- Domain and boundary conditions -------
	stretchingDist = 1/4;
	stretchingRatio = 2;
	def = domains.SinusoidalFault(stretchingDist, stretchingRatio);

	% --- Traction BC on surface ----
	bcFree1.boundary = def.boundaryGroups.surface;
	bcFree1.type = {1, 't'};

	bcFree2.boundary = def.boundaryGroups.surface;
	bcFree2.type = {2, 't'};

	% Displacement BC on remaining boundaries
	bcD1.boundary = def.boundaryGroups.outer;
	bcD1.type = {1, 'd'};

	bcD2.boundary = def.boundaryGroups.outer;
	bcD2.type = {2, 'd'};

	bc = {bcFree1, bcFree2, bcD1, bcD2};

	% Fault boundary groups
	faultBoundaryGroups = {def.boundaryGroups.faultBottom, def.boundaryGroups.faultTop};

	iS.domain = elastic.friction.input.domainStruct(iS.domain, def, bc, faultBoundaryGroups);
	%------------------------------------------------

	% ----- Friction parameters ------------------
	iS.friction.a = 0.03; 			% Direct effect [ ]
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
	A = 1000*1e6;
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

	% ------ Plotting ----------------------
end