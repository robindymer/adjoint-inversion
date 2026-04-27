function [domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = mms(outputDir, MODE)

	% - outputDir: must include / at the end.
	% - MODE: 'startSimulation' or 'plot'

	default_arg('outputDir', 'mms/');
	default_arg('MODE', 'startSimulation');

	% Load MAT-file with MMS solution if it exists
	if (exist('+elastic/+viscofriction/+input/mms.mat','file') == 2 )
		load('+elastic/+viscofriction/+input/mms.mat');
		return;
	end

	% Start with default options
	[domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = elastic.viscofriction.input.default();

	% --------- Time stepping -----------------
	timeStepper.relTol = 1e-5;

	% --------- Exact solution -----------------
	syms x y t

	% Polar coordinates
	theta = symfun( atan2(y,x), [t, x, y]);
	r = symfun( sqrt(x^2 + y^2), [t, x, y] );
	rHat = {x/r, y/r};
	thetaHat = {-rHat{2}, rHat{1}};

	% ------- Solution block 1 ------------------------
	% Radial component of solution
	wvl = 1e3;
	dir = [1,1];
	c = 2e3;
	k = elastic.wavelengthToWavenumber(wvl, dir);
	w = elastic.angularFrequency(c, wvl);

	ur = symfun( sin(k(1)*x + k(2)*y - w*t), [t, x, y]);

	% Angular component of solution
	wvl = 1e3;
	dir = [1,-1];
	c = 2e3;
	k = elastic.wavelengthToWavenumber(wvl, dir);
	w = elastic.angularFrequency(c, wvl);

	utheta = symfun( sin(k(1)*x + k(2)*y - w*t), [t, x, y]);

	% Get Cartesian components
	u1_1 = rHat{1}*ur + thetaHat{1}*utheta;
	u2_1 = rHat{2}*ur + thetaHat{2}*utheta;

	% Viscous strains
	gamma11_1 = symfun( sin(k(1)/2*x + k(2)*y - w*t), [t, x, y]);
	gamma12_1 = symfun( sin(k(1)*x + k(2)/2*y - w*t), [t, x, y]);
	gamma21_1 = symfun( sin(k(1)/3*x + k(2)*y - w*t), [t, x, y]);
	gamma22_1 = symfun( sin(k(1)*x + k(2)/3*y - w*t), [t, x, y]);
	gamma1 = {gamma11_1, gamma12_1; gamma21_1, gamma22_1};
	% -------------------------------------------------

	% ------- Solution block 2 ------------------------
	% Use same radial component

	% Angular component of solution
	wvl = 2e3;
	dir = [-1,1];
	c = 6e3;
	k = elastic.wavelengthToWavenumber(wvl, dir);
	w = elastic.angularFrequency(c, wvl);

	utheta = symfun( sin(k(1)*x + k(2)*y - w*t), [t, x, y]);

	% Get Cartesian components
	u1_2 = rHat{1}*ur + thetaHat{1}*utheta;
	u2_2 = rHat{2}*ur + thetaHat{2}*utheta;

	% Viscous strains
	gamma11_2 = symfun( sin(k(1)*x + k(2)/2*y - w*t), [t, x, y]);
	gamma12_2 = symfun( sin(k(1)/2*x + k(2)*y - w*t), [t, x, y]);
	gamma21_2 = symfun( sin(k(1)*x + k(2)/3*y - w*t), [t, x, y]);
	gamma22_2 = symfun( sin(k(1)/3*x + k(2)*y - w*t), [t, x, y]);
	gamma2 = {gamma11_2, gamma12_2; gamma21_2, gamma22_2};
	% -------------------------------------------------

	% ------ Material parameters -----------------
	rho = symfun( 2e0 + 0*x, [t,x,y] );
	vs = symfun( 2e0*(1 + 0*x), [t,x,y] );
	vp = symfun( 6e0*(1 + 0*x), [t,x,y] );
	[lambda, mu] = elastic.speedsToModuli(vp, vs, rho);
	C = elastic.isotropicStiffnessTensor(lambda, mu);

	rhoFun = matlabFunctionSizePreserving( symfun(rho(0,x,y), [x, y]) );
	vsFun = matlabFunctionSizePreserving( symfun(vs(0,x,y), [x, y]) );
	vpFun = matlabFunctionSizePreserving( symfun(vp(0,x,y), [x, y]) );

	material = elastic.viscofriction.input.materialStruct(rhoFun, vsFun, vpFun);
	% ---------------------------------------------

	% ---- Pre stress ----------------------------
	sigma0 = symfun( 1e11 + 0*x, [t,x,y] );
	tau0 = 0.3*sigma0;

	preStress.sigma0 = matlabFunctionSizePreserving( symfun(sigma0(0,x,y), [x, y]) );
	preStress.shear = matlabFunctionSizePreserving( symfun(tau0(0,x,y), [x, y]) );
	% ---------------------------------------------

	% ----- Friction parameters ------------------
	friction.a = symfun(0.03+0.01*sin(2*k(1)*x), [t, x, y]); 			% Direct effect [ ]
	friction.b = symfun(0.02+0.01*cos(k(1)*x+2*k(2)*y), [t, x, y]);			% State evolution effect [ ]
	friction.dc = 0.25;			% State evolution distance [ m ]
	friction.V0 = 1e-6;			% Ref velocity [ m/s ]
	friction.f0 = 0.6;			% Ref friction coeff for steady sliding [ ]
	% friction.psi0 = 0.4367;		% Initial state variable

	friction.coefficient = 'standard';
	friction.coefficientSteadyState = [];

	friction.law = 'slip';
	% friction.law = 'aging';

	friction.eta = mu/(2*vs);

	% Choose state variable psi
	wvl = 1e3;
	c = 2e3;
	k = elastic.wavelengthToWavenumber(wvl, dir);
	w = elastic.angularFrequency(c, wvl);
	psi = symfun( 0.4*(1+1/2*sin(k(1)*x + k(2)*y - w*t)) , [t,x,y]);

	friction.psi0 = matlabFunctionSizePreserving( symfun(psi(0,x,y), [x, y]) );
	%----------------------------------------------

	%-------- Compute PDE forcing ------------
	rhoEff = symfun( 0 + 0*x, [t,x,y] );
	forcingFuncs1 = elastic.mms.computeViscoElasticForcing({u1_1, u2_1}, gamma1, C, rhoEff, flowLaw.eta.nonlinFun, x, y, t);
	forcingFuncs2 = elastic.mms.computeViscoElasticForcing({u1_2, u2_2}, gamma2, C, rhoEff, flowLaw.eta.nonlinFun, x, y, t);

	% ----- Fault slip -----------------------
	delta = forcingFuncs1.sym.u_t - forcingFuncs2.sym.u_t;

	% ------- Compute friction coeff forcing ---
	V = forcingFuncs1.sym.ut_t - forcingFuncs2.sym.ut_t;

	% Prescribe f
	f = symfun (0.5 + 0.1*cos(k(1)*x - k(2)*y), [t,x,y]);

	% Compute (compressive) normal stress on fault as average of block 1 and 2.
	normalStressPerturbation = -1/2*(forcingFuncs1.sym.tau_n + forcingFuncs2.sym.tau_n);
	normalStress = sigma0 + normalStressPerturbation;
	normalStress = 1/2*(1+sign(normalStress))*normalStress;

	% Traction corresponding to the friction coefficient
	frictionTraction = normalStress*f*sign(V) + friction.eta*V;
	frictionTractionPerturbation = frictionTraction - tau0;

	% The traction on the fault is computed as mean of traction on the two sides
	estimatedFrictionTraction = -1/2*(forcingFuncs1.sym.tau_t + forcingFuncs2.sym.tau_t);
	frictionTractionForcing = frictionTractionPerturbation - estimatedFrictionTraction;

	% Tractions should be continuous across fault -> difference is added as forcing
	shearTractionForcing1 = forcingFuncs1.sym.tau_t - forcingFuncs2.sym.tau_t;
	shearTractionForcing2 = forcingFuncs2.sym.tau_t - forcingFuncs1.sym.tau_t;

	normalTractionForcing1 = forcingFuncs1.sym.tau_n - forcingFuncs2.sym.tau_n;
	normalTractionForcing2 = forcingFuncs2.sym.tau_n - forcingFuncs1.sym.tau_n;

	frictionForcing = elastic.friction.mms.computeFrictionForcing(f, V, psi, friction);

	% ------- Compute state evolution forcing ---
	stateForcing = elastic.friction.mms.computeStateForcing(psi, V, f, friction, t);

	% ---- Turn friction symfuns into function handles ------
	friction.a = matlabFunctionSizePreserving( symfun(friction.a(0,x,y), [x, y]) );
	friction.b = matlabFunctionSizePreserving( symfun(friction.b(0,x,y), [x, y]) );
	friction.eta = matlabFunctionSizePreserving( symfun(friction.eta(0,x,y), [x, y]) );

	% -------- Fill mms struct --------------------
	mms.delta0 = matlabFunctionSizePreserving( symfun(delta(0,x,y), [x,y]) );
	mms.delta = matlabFunctionSizePreserving(delta);
	mms.V = matlabFunctionSizePreserving(V);

	u_1 = @(t,x,y)[forcingFuncs1.u1(t,x,y); forcingFuncs1.u2(t,x,y)];
	u_2 = @(t,x,y)[forcingFuncs2.u1(t,x,y); forcingFuncs2.u2(t,x,y)];
	mms.u = {u_1, u_2};

	gamma_1 = @(t,x,y)[forcingFuncs1.g11(t,x,y); forcingFuncs1.g12(t,x,y);...
					   forcingFuncs1.g21(t,x,y); forcingFuncs1.g22(t,x,y)];

    gamma_2 = @(t,x,y)[forcingFuncs2.g11(t,x,y); forcingFuncs2.g12(t,x,y);...
					   forcingFuncs2.g21(t,x,y); forcingFuncs2.g22(t,x,y)];
    mms.gamma = {gamma_1, gamma_2};


	mms.F = struct;

	Fu1 = @(t,x,y) [forcingFuncs1.f1(t,x,y); forcingFuncs1.f2(t,x,y)];
	Fu2 = @(t,x,y) [forcingFuncs2.f1(t,x,y); forcingFuncs2.f2(t,x,y)];
	mms.F.u = {Fu1, Fu2};

	Fg1 = @(t,x,y) [forcingFuncs1.fg11(t,x,y); forcingFuncs1.fg12(t,x,y);...
					forcingFuncs1.fg21(t,x,y); forcingFuncs1.fg22(t,x,y)];
	Fg2 = @(t,x,y) [forcingFuncs2.fg11(t,x,y); forcingFuncs2.fg12(t,x,y);...
					forcingFuncs2.fg21(t,x,y); forcingFuncs2.fg22(t,x,y)];
	mms.F.gamma = {Fg1, Fg2};

	mms.psi = matlabFunctionSizePreserving(psi);
	mms.frictionForcing = frictionForcing;
	mms.stateForcing = stateForcing;

	mms.shearTractionForcing1 = matlabFunctionSizePreserving(shearTractionForcing1);
	mms.shearTractionForcing2 = matlabFunctionSizePreserving(shearTractionForcing2);

	mms.normalTractionForcing1 = matlabFunctionSizePreserving(normalTractionForcing1);
	mms.normalTractionForcing2 = matlabFunctionSizePreserving(normalTractionForcing2);

	mms.frictionTractionForcing = matlabFunctionSizePreserving(frictionTractionForcing);

	% ----- Domain and boundary conditions -------
	def = domains.CircularFault();

	% -- Displacement BC everywhere on block 1 --
	bc11.boundary = def.boundaryGroups.block1;
	bc11.type = {1, 'd'};
	bc11.data = forcingFuncs1.u1;

	bc12.boundary = def.boundaryGroups.block1;
	bc12.type = {2, 'd'};
	bc12.data = forcingFuncs1.u2;

	% -- Displacement BC for sides of block 2
	bc21.boundary = def.boundaryGroups.linesBlock2;
	bc21.type = {1, 'd'};
	bc21.data = forcingFuncs2.u1;

	bc22.boundary = def.boundaryGroups.linesBlock2;
	bc22.type = {2, 'd'};
	bc22.data = forcingFuncs2.u2;

	% -- Traction BC for free surface
	bcFree1.boundary = def.boundaryGroups.outerArc;
	bcFree1.type = {1, 't'};
	bcFree1.data = forcingFuncs2.tau1;

	bcFree2.boundary = def.boundaryGroups.outerArc;
	bcFree2.type = {2, 't'};
	bcFree2.data = forcingFuncs2.tau2;

	bc = {bc11, bc12, bc21, bc22, bcFree1, bcFree2};

	% Fault boundary groups
	faultBoundaryGroups = {def.boundaryGroups.innerInterface, def.boundaryGroups.outerInterface};

	domain = elastic.viscofriction.input.domainStruct(domain, def, bc, faultBoundaryGroups);
	%------------------------------------------------

	% ------ Plotting ----------------------

	save('+elastic/+viscofriction/+input/mms.mat', 'domain', 'friction', 'preStress', ...
					 'material', 'mms', 'flowLaw', 'plotting', 'tectonic', 'timeStepper', 'output');

end
