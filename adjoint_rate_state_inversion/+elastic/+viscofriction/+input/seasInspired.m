function [domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = seasInspired(outputDir, MODE)
	% - outputDir: must include / at the end.
	% - MODE: 'startSimulation' or 'plot'

	default_arg('outputDir', 'seasInspired/');
	default_arg('MODE', 'startSimulation');

	% Start with default options
	[domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = elastic.viscofriction.input.default();

	% Time-stepping
	timeStepper.relTol = 1e-3;

	% ------ Material parameters -----------------
	rho = @(x,y) 2670 + 0*x;
	vs = @(x,y) 3464 + 0*x;
	vp = @(x,y) sqrt(3)*vs(x,y);
	[lambda, mu] = elastic.speedsToModuli(vp, vs, rho);
	C = elastic.isotropicStiffnessTensor(lambda, mu);

	material = elastic.viscofriction.input.materialStruct(rho, vs, vp);
	% ---------------------------------------------

	% ---- Pre stress ----------------------------
	preStress.sigma0 = @(x,y) 50*1e6 + 0*x;
	preStress.shear = @(x,y) 0*preStress.sigma0(x,y);
	% ---------------------------------------------

	% ----- Friction parameters ------------------
	% VW patch, b is kept constant but a varies.
	amin = 0.01;
	amax = 0.025;
	transitionLength = 3*1e3;
	VWLength = 44*1e3;
	aFun = @(x,y) helpers.linearBoxCar(x, 0, VWLength, transitionLength, amax, amin);
	% aFun = @(x,y) helpers.smoothedBoxCar(x, 0, VWLength, transitionLength, amax, amin);

	friction.a = aFun; 			% Direct effect [ ]
	friction.b = 0.015;			% State evolution effect [ ]
	friction.dc = [];			% State evolution distance [ m ]
	friction.V0 = 1e-6;			% Ref velocity [ m/s ]
	friction.f0 = 0.6;			% Ref friction coeff for steady sliding [ ]
	friction.L = 4*0.008;			% Critical slip distance (used in reg. friction coeff)

	% Set state variable to steady-state
	Vp = 1e-9;
	theta0 = friction.L/Vp;
	friction.psi0 = friction.f0 + friction.b.*log(friction.V0*theta0/friction.L);

	friction.coefficient = 'standard';
	friction.coefficientSteadyState = [];
	friction.law = 'aging';

	rho = 2670;
	vs = 3464;
	mu = rho*vs^2;
	friction.eta = mu/(2*vs);	% Radiation damping
	%----------------------------------------------

	% ----- Domain and boundary conditions -------
	L = 2*3*VWLength;
	H = 2*VWLength;

	% Relevant length scales
	qsProcessZone = mu*friction.L/(friction.b*preStress.sigma0(0,0));
	nucleationZoneSize = 2/pi*mu*friction.b*friction.L/((friction.b-amin)^2*preStress.sigma0(0,0) );
	fprintf('Velocity-weakening zone: %3.1f km \n', VWLength/1e3);
	fprintf('Quasi-static process zone: %3.1f km \n', qsProcessZone/1e3);
	fprintf('Nucleation zone size: %3.1f km \n', nucleationZoneSize/1e3);
	fprintf('Grid spacing with m = 100: %3.1f km \n', L/1e3/100);

	stretchingDist = 1/2;
	stretchingRatio = 3;

	def = domains.MapViewFault(stretchingDist, stretchingRatio, L, H, 0);
	boundaryGroups = def.boundaryGroups;

	% --- Traction BC on west and east sides ----
	bcFree1.boundary = boundaryGroups.sides;
	bcFree1.type = {'n', 't'};

	bcFree2.boundary = boundaryGroups.sides;
	bcFree2.type = {'t', 't'};

	% Tectonic plate BC on north and south boundaries
	bcTectonic1.boundary = boundaryGroups.topAndBottom;
	bcTectonic1.type = {'t', 'd'};

	bcTectonic2.boundary = boundaryGroups.topAndBottom;
	bcTectonic2.type = {'n', 'd'};

	bc = {bcFree1, bcFree2, bcTectonic1, bcTectonic2};

	% Fault boundary groups
	faultBoundaryGroups = {boundaryGroups.faultBottom, boundaryGroups.faultTop};

	domain = elastic.viscofriction.input.domainStruct(domain, def, bc, faultBoundaryGroups);
	%------------------------------------------------

	% ------ Tectonic plate movement -------
	tectonic = struct;
	tectonic.v = 1e-9;
	tectonic.bc = struct;
	tectonic.bc.boundary = boundaryGroups.topAndBottom;
	tectonic.bc.type = {'t', 'd'};

	tectonic.u0 = 0*1*H/mu*friction.f0*preStress.sigma0(0,0);

	% ------ Flow law --------
	flowLaw.eta.nonlinFun = @flowLawFun;


	% ------ Plotting ----------------------
	plotting.setupFun = @elastic.viscofriction.plotting.stresses;

	% ------ Output ------------------------

	% strides: number of time steps between writes
	scalarStride = 5;
	faultFieldStride = 20;

	% Fields to output
	output.fields = {'slipVelocity', 'shearStressMax', 'shearStressMean', 'slipProfile', 'shearStressProfile'};
	output.strides = {scalarStride, scalarStride, scalarStride, faultFieldStride, faultFieldStride};


	% ---- Create directories etc. --------
	output.dir = outputDir;
	switch MODE
	case 'startSimulation'
		helpers.createDirectory(output.dir);
		if helpers.fileExists([output.dir, 'countArray.mat'])
			delete([output.dir, 'countArray.mat']);
		end

	case 'plot'

	end
	%---------------------------------------


end

function etaInv = flowLawFun(sigma)
	% Compute an invariant
	tauBar = sqrt(sigma{1,2}.^2 + sigma{2,1}.^2);

	n = 2.4;
	A = 1.3*1e-3 * (1e-6)^n;
	Q = 219*1e3;
	R = 8.31;
	T = 273.15 + 225;

	etaInv = 0*A*exp(-Q/(R*T))*tauBar.^(n-1);
end

