function [domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper] = sinusoidalFault()

	% Start with default options
	[domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper] = elastic.viscofriction.input.default();

	% ------ Material parameters -----------------
	rho = @(x,y) 2 + 0*x;
	vs = @(x,y) 2 + 0*x;
	vp = @(x,y) 6 + 0*x;
	[lambda, mu] = elastic.speedsToModuli(vp, vs, rho);
	C = elastic.isotropicStiffnessTensor(lambda, mu);

	material = elastic.viscofriction.input.materialStruct(rho, vs, vp);
	% ---------------------------------------------

	% ---- Pre stress ----------------------------
	preStress.sigma0 = @(x,y) 1e8 + 0*x;
	preStress.shear = @(x,y) 0.3*preStress.sigma0(x,y);
	% ---------------------------------------------

	% ----- Friction parameters ------------------
	friction.a = 0.03; 			% Direct effect [ ]
	friction.b = 0.02;			% State evolution effect [ ]
	friction.dc = 0.25;			% State evolution distance [ m ]
	friction.V0 = 1e-6;			% Ref velocity [ m/s ]
	friction.f0 = 0.6;			% Ref friction coeff for steady sliding [ ]
	friction.psi0 = 0.4367;		% Initial state variable

	friction.coefficient = 'Rice2001';
	friction.coefficientSteadyState = 'Rice2001';
	friction.law = 'slip';
	%----------------------------------------------

	% ----- Domain and boundary conditions -------
	L = 2;
	Bottom = -2;
	faultPos = -1;

	stretchingDist = 1/4;
	stretchingRatio = 1;
	def = domains.SinusoidalFault(stretchingDist, stretchingRatio, L, faultPos, Bottom);

	% --- Traction BC on west and east sides ----
	bcFree1.boundary = def.boundaryGroups.sides;
	bcFree1.type = {1, 't'};

	bcFree2.boundary = def.boundaryGroups.sides;
	bcFree2.type = {2, 't'};

	% Tectonic plate BC on north and south boundaries
	bcTectonic1.boundary = def.boundaryGroups.topAndBottom;
	bcTectonic1.type = {'t', 'd'};

	bcTectonic2.boundary = def.boundaryGroups.topAndBottom;
	bcTectonic2.type = {'n', 't'};

	bc = {bcFree1, bcFree2, bcTectonic1, bcTectonic2};

	% Fault boundary groups
	faultBoundaryGroups = {def.boundaryGroups.faultBottom, def.boundaryGroups.faultTop};

	domain = elastic.viscofriction.input.domainStruct(domain, def, bc, faultBoundaryGroups);
	%------------------------------------------------

	% ------ Tectonic plate movement -------
	tectonic = struct;
	tectonic.v = 1e-9;
	tectonic.bc = struct;
	tectonic.bc.boundary = def.boundaryGroups.topAndBottom;
	tectonic.bc.type = {'t', 'd'};

	% ------ Plotting ----------------------

end
