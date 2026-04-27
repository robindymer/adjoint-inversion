function [domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = mapViewFault()

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

	friction.a = aFun; 			% Direct effect [ ]
	friction.b = 0.015;			% State evolution effect [ ]
	friction.dc = [];			% State evolution distance [ m ]
	friction.V0 = 1e-6;			% Ref velocity [ m/s ]
	friction.f0 = 0.6;			% Ref friction coeff for steady sliding [ ]

	friction.L = 4*0.008;			% Critical slip distance (used in reg. friction coeff)

	friction.coefficient = 'standard';
	friction.coefficientSteadyState = [];
	friction.law = 'aging';

	% Set state variable to steady-state
	Vp = 1e-9;
	theta0 = friction.L/Vp;
	friction.psi0 = friction.f0 + friction.b.*log(friction.V0*theta0/friction.L);

	rho = rho(0,0);
	vs = vs(0,0);
	mu = rho*vs^2;
	friction.eta = mu/(2*vs);	% Radiation damping
	%----------------------------------------------

	% ----- Domain and boundary conditions -------
	L = 2*VWLength;
	H = VWLength;

	% Relevant length scales
	qsProcessZone = mu*friction.L/(friction.b*preStress.sigma0(0,0));
	nucleationZoneSize = 2/pi*mu*friction.b*friction.L/((friction.b-amin)^2*preStress.sigma0(0,0) );
	fprintf('Velocity-weakening zone: %3.1f km \n', VWLength/1e3);
	fprintf('Quasi-static process zone: %3.1f km \n', qsProcessZone/1e3);
	fprintf('Nucleation zone size: %3.1f km \n', nucleationZoneSize/1e3);
	fprintf('Grid spacing with m = 100: %3.1f km \n', L/1e3/100);

	stretchingDist = 1/2;
	stretchingRatio = 3;

	mSG = 41;
	ratioSG = 2;

	def = domains.MapViewFaultSG(stretchingDist, stretchingRatio, L, H, mSG, ratioSG, 0);
	g = def.getGrid(11); % Get dummy grid to obtain boundary groups
	boundaryGroups = def.getBoundaryGroups(g);

	% --- Traction BC on west and east sides ----
	bcFree1.boundary = boundaryGroups.WE;
	bcFree1.type = {1, 't'};

	bcFree2.boundary = boundaryGroups.WE;
	bcFree2.type = {2, 't'};

	% Tectonic plate BC on north and south boundaries
	bcTectonic1.boundary = boundaryGroups.SN;
	bcTectonic1.type = {'t', 'd'};

	bcTectonic2.boundary = boundaryGroups.SN;
	bcTectonic2.type = {'n', 'd'};

	bc = {bcFree1, bcFree2, bcTectonic1, bcTectonic2};

	% Fault boundary groups
	faultBoundaryGroups = {boundaryGroups.faultBottom, boundaryGroups.faultTop};

	domain = elastic.viscofriction.input.domainStruct(domain, def, bc, faultBoundaryGroups);
	%------------------------------------------------

	% ------ Tectonic plate movement -------
	tectonic.v = 1e-9;
	tectonic.bc = struct;
	tectonic.bc.boundary = boundaryGroups.SN;
	tectonic.bc.type = {'t', 'd'};

	% ------ Flow law --------
	flowLaw.eta.nonlinFun = @flowLawFun;


	% ------ Plotting ----------------------
	plotting.setupFun = @elastic.viscofriction.plotting.stresses;

	% ------ Output ------------------------
	% output.dir = 'output/test/';
	% helpers.createDirectory(output.dir);

	% if helpers.fileExists([output.dir, 'countArray.mat'])
	% 	delete([output.dir, 'countArray.mat']);
	% end

	% output.fields = {'slipVelocity', 'shearStressMax', 'shearStressMean', 'psiMin', 'psiMax', 'slipProfile', 'shearStressProfile'};

	% oneYear = 365.25*24*3600;
	% times = (0:0.00005:100)*oneYear;

	% fieldTimes = (0:0.05:100)*oneYear;
	% output.targetTimes = {times, times, times, times, times, fieldTimes, fieldTimes};


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

