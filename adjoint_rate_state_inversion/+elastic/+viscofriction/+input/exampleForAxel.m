function [domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = exampleForAxel(outputDir, MODE, dipAngle, faultMotion)
	% Example for Axel and Eric that solves the SEAS BP3-QD benchmark
	%-------------------------------------------------------------
	% - outputDir: must include / at the end.
	% - MODE: 'startSimulation' or 'plot'
	% - dipAngle: Angle in radians
	% - faultMotion: 'normal' or 'thrust'

	default_arg('outputDir', 'SEAS_BP3QS_LU/');
	default_arg('MODE', 'startSimulation');

	% Start with default options
	[domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = elastic.viscofriction.input.default();

	% Time-stepping
	timeStepper.relTol = 1e-3;
	timeStepper.solver = 'lu';

	% ------ Material parameters -----------------
	rho = @(x,y) 2670 + 0*x;
	vs = @(x,y) 3464 + 0*x;
	% nu = 0.25 => lambda = mu.
	vp = @(x,y) sqrt(3)*vs(x,y);
	[lambda, mu] = elastic.speedsToModuli(vp, vs, rho);
	C = elastic.isotropicStiffnessTensor(lambda, mu);

	material = elastic.viscofriction.input.materialStruct(rho, vs, vp);
	% ---------------------------------------------

	% ----- Friction parameters ------------------
	% VW patch, b is kept constant but a varies.
	amin = 0.01;
	amax = 0.025;
	transitionLength = 3*1e3;
	VWLength = 15*1e3;

	r = @(x,y) sqrt(x.^2 + y.^2);
	aFun = @(x,y) helpers.linearBoxCar(r(x,y), 0, 2*VWLength, transitionLength, amax, amin);

	friction.a = aFun; 			% Direct effect [ ]
	friction.b = 0.015;			% State evolution effect [ ]
	friction.dc = [];			% State evolution distance [ m ]
	friction.V0 = 1e-6;			% Ref velocity [ m/s ]
	friction.f0 = 0.6;			% Ref friction coeff for steady sliding [ ]
	friction.L = 0.008;			% Critical slip distance (used in reg. friction coeff)

	rho = 2670;
	vs = 3464;
	mu = rho*vs^2;
	friction.eta = mu/(2*vs);	% Radiation damping

	friction.coefficient = 'standard';
	friction.law = 'aging';
	%----------------------------------------------

	% ---- Other parameters ----
	Vp = 1e-9;
	Vinit = 1e-9;
	V0 = friction.V0;
	L = friction.L;
	a = friction.a;
	b = friction.b;
	f0 = friction.f0;
	eta = friction.eta;

	% ---- Pre stress ----------------------------
	preStress.sigma0 = @(x,y) 50*1e6 + 0*x;
	sigma0 = preStress.sigma0(0,0);

	preStress.shear = @(x,y) preStress.sigma0(x,y)*amax.*asinh(Vinit/(2*V0)*exp((f0+b*log(V0/Vinit))/amax)) + eta*Vinit;
	% ---------------------------------------------

	%---- Set state variable to steady-state ---
	tau = preStress.shear;

	theta0 = @(x,y) L/V0*exp(a(x,y)/b.*log(2*V0/Vinit*sinh( (tau(x,y)-eta*Vinit)./(a(x,y)*sigma0) ) ) - f0/b);
	friction.psi0 = @(x,y) friction.f0 + friction.b.*log(friction.V0*theta0(x,y)/friction.L);
	%----------------------------------------------


	% ----- Domain and boundary conditions -------
	default_arg('dipAngle', pi/2);
	Wf = 40*1e3;
	Lx = 150*1e3;
	Ly = 100*1e3;

	% Relevant length scales
	% qsProcessZone = mu*friction.L/(friction.b*preStress.sigma0(0,0));
	% nucleationZoneSize = 2/pi*mu*friction.b*friction.L/((friction.b-amin)^2*preStress.sigma0(0,0) );
	% fprintf('Velocity-weakening zone: %3.1f km \n', VWLength/1e3);
	% fprintf('Quasi-static process zone: %3.1f km \n', qsProcessZone/1e3);
	% fprintf('Nucleation zone size: %3.1f km \n', nucleationZoneSize/1e3);
	% fprintf('Grid spacing with m = 100: %3.1f km \n', L/1e3/100);

	stretchingDist = [];
	stretchingRatio = 20;

	def = domains.DippingFault(dipAngle, Wf, Lx, Ly, stretchingDist, stretchingDist, stretchingRatio, stretchingRatio);
	boundaryGroups = def.boundaryGroups;

	% {'t', 'd'} means specifying tangential ('t') component of displacement ('d')
	% {'n', 't'} means specifying normal ('n') componenet of traction ('t')
	% etc.

	% Displacement BC on west and east sides
	bcSidesN.boundary = boundaryGroups.sides;
	bcSidesN.type = {'n', 'd'};

	bcSidesT.boundary = boundaryGroups.sides;
	bcSidesT.type = {'t', 'd'};

	% Traction BC on bottom boundary
	bcBottomN.boundary = boundaryGroups.bottom;
	bcBottomN.type = {'n', 't'};

	bcBottomT.boundary = boundaryGroups.bottom;
	bcBottomT.type = {'t', 't'};

	% Free BC on surface
	bcSurfaceN.boundary = boundaryGroups.surface;
	bcSurfaceN.type = {'n', 't'};

	bcSurfaceT.boundary = boundaryGroups.surface;
	bcSurfaceT.type = {'t', 't'};

	bc = {bcSidesN, bcSidesT, bcBottomT, bcBottomN, bcSurfaceT, bcSurfaceN};

	% Fault boundary groups
	faultBoundaryGroups = {boundaryGroups.faultLeft, boundaryGroups.faultRight};

	domain = elastic.viscofriction.input.domainStruct(domain, def, bc, faultBoundaryGroups);

	% Add surface boundary groups (necessary for surface output)
	domain.surfaceBoundaryGroups = boundaryGroups.surface;
	%------------------------------------------------

	% ------ Tectonic plate movement -------
	default_arg('faultMotion', 'normal');

	tectonic = cell(2, 1);
	for i = 1:numel(tectonic)
		tectonic{i} = struct;
	end

	% Fault sliding
	tectonic{1}.interface = struct;
	tectonic{1}.interface.boundaryGroups = boundaryGroups.tectonic;
	tectonic{1}.u0 = 0;

	% Side boundaries
	tectonic{2}.bc = struct;
	tectonic{2}.u0 = 0;
	tectonic{2}.bc.boundary = boundaryGroups.sides;
	tectonic{2}.bc.type = {'t', 'd'};

	switch faultMotion
	case {'reverse', 'thrust'}
		tectonic{1}.v = -Vp;
		preStress.shear = @(x,y) -preStress.shear(x,y);

	case 'normal'
		tectonic{1}.v = Vp;
	end

	% Switch sign compared to fault sliding.
	tectonic{2}.v = -tectonic{1}.v;

	% ------ Flow law --------
	% 0 gives no viscous flow
	flowLaw.eta.nonlinFun = @(sigma) 0*sigma{1,1};


	% ------ Plotting ----------------------
	plotting.setupFun = @elastic.viscofriction.plotting.stresses;
	plotting.faultCoordinate = 'y';

	% ------ Output ------------------------

	% Fields to output
	output.fields = {'slipVelocity', 'slipVelocityProfile', 'slipProfile',...
					 'shearStressProfile', 'normalStressProfile', 'psiProfile',...
					 'surfaceDisplacementX', 'surfaceDisplacementY'};

	% strides: number of time steps between writes
	scalarStride = 5;
	faultFieldStride = 5;
	output.strides = cell(1, numel(output.fields));
	for i = 1:length(output.fields)
		if contains(output.fields{i}, 'Profile')
			output.strides{i} = faultFieldStride;
		elseif contains(output.fields{i}, 'surface')
			output.strides{i} = faultFieldStride;
		else
			output.strides{i} = scalarStride;
		end
	end

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

