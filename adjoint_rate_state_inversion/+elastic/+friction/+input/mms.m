function iS = mms()
	% iS: inputStruct
	%--------------------------------------------------

	% Load MAT-file with MMS solution if it exists
	if (exist('+elastic/+friction/+input/mms.mat','file') == 2 )
		load('+elastic/+friction/+input/mms.mat');
		return;
	end

	% Start with default options
	iS = elastic.friction.input.default();

	% --------- Time stepper ------------------
	iS.timeStepper.type = 'RK4';

	% ------ Spatial operator ------------------
	iS.solver.mex = false;
	iS.solver.check_ev_flag = false;

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

	utheta = symfun( -w*t + sin(k(1)*x + k(2)*y - w*t), [t, x, y]);

	% Get Cartesian components
	u1_1 = rHat{1}*ur + thetaHat{1}*utheta;
	u2_1 = rHat{2}*ur + thetaHat{2}*utheta;
	% -------------------------------------------------

	% ------- Solution block 2 ------------------------
	% Use same radial component

	% Angular component of solution
	wvl = 2e3;
	dir = [-1,1];
	c = 6e3;
	k = elastic.wavelengthToWavenumber(wvl, dir);
	w = elastic.angularFrequency(c, wvl);

	utheta = symfun( w*t + sin(k(1)*x + k(2)*y - w*t), [t, x, y]);

	% Get Cartesian components
	u1_2 = rHat{1}*ur + thetaHat{1}*utheta;
	u2_2 = rHat{2}*ur + thetaHat{2}*utheta;
	% -------------------------------------------------

	% ------ Material parameters -----------------
	rho = symfun( 2e3*(1 + 0*x), [t,x,y] );
	vs = symfun( 2e3*(1 + 0*x), [t,x,y] );
	vp = symfun( 6e3*(1 + 0*x), [t,x,y] );
	[lambda, mu] = elastic.speedsToModuli(vp, vs, rho);
	C = elastic.isotropicStiffnessTensor(lambda, mu);

	rhoFun = matlabFunctionSizePreserving( symfun(rho(0,x,y), [x, y]) );
	vsFun = matlabFunctionSizePreserving( symfun(vs(0,x,y), [x, y]) );
	vpFun = matlabFunctionSizePreserving( symfun(vp(0,x,y), [x, y]) );

	iS.material = elastic.friction.input.materialStruct(rhoFun, vsFun, vpFun);
	% ---------------------------------------------

	% ---- Pre stress ----------------------------
	sigma0 = symfun( 1/1000*126*1e6 + 0*x, [t,x,y] );
	tau0 = 0.3*sigma0;

	iS.preStress.sigma0 = matlabFunctionSizePreserving( symfun(sigma0(0,x,y), [x, y]) );
	iS.preStress.shear = matlabFunctionSizePreserving( symfun(tau0(0,x,y), [x, y]) );
	% ---------------------------------------------

	% ----- Friction parameters ------------------
	iS.friction.a = symfun(0.03+0.01*sin(2*k(1)*x), [t, x, y]); 			% Direct effect [ ]
	iS.friction.b = symfun(0.02+0.01*cos(k(1)*x+2*k(2)*y), [t, x, y]);			% State evolution effect [ ]
	iS.friction.L = 0.25;			% State evolution distance [ m ]
	iS.friction.V0 = 1e-6;			% Ref velocity [ m/s ]
	iS.friction.f0 = 0.6;			% Ref friction coeff for steady sliding [ ]
	% iS.friction.psi0 = 0.4367;		% Initial state variable

	iS.friction.coefficient = 'standard';
	iS.friction.coefficientSteadyState = [];

	% iS.friction.law = 'slip';
	% iS.timeStepper.CFL = 1;

	iS.friction.law = 'aging';
	iS.timeStepper.CFL = 0.1;

	% Choose state variable psi
	wvl = 1e3;
	c = 2e3;
	k = elastic.wavelengthToWavenumber(wvl, dir);
	w = elastic.angularFrequency(c, wvl);
	psi = symfun( 0.4*(1+1/2*sin(k(1)*x + k(2)*y - w*t)) , [t,x,y]);

	iS.friction.psi0 = matlabFunctionSizePreserving( symfun(psi(0,x,y), [x, y]) );
	%----------------------------------------------

	%-------- Compute PDE forcing ------------
	forcingFuncs1 = elastic.mms.computeForcing(u1_1, u2_1, C, rho, x, y, t);
	forcingFuncs2 = elastic.mms.computeForcing(u1_2, u2_2, C, rho, x, y, t);

	% ------- Compute friction coeff forcing ---
	V = forcingFuncs1.sym.ut_t - forcingFuncs2.sym.ut_t;
	% f = abs( (forcingFuncs1.sym.tau_t + tau0)/sigma0 );

	% Prescribe f
	f = symfun (0.5 + 0.1*cos(k(1)*x - k(2)*y), [t,x,y]);

	% Compute (compressive) normal stress on fault as average of block 1 and 2.
	normalStressPerturbation = -1/2*(forcingFuncs1.sym.tau_n + forcingFuncs2.sym.tau_n);
	normalStress = sigma0 + normalStressPerturbation;
	normalStress = 1/2*(1+sign(normalStress))*normalStress;

	frictionTraction = -normalStress*f*sign(V);
	tractionForcing1 = forcingFuncs1.sym.tau_t - (frictionTraction - tau0);
	tractionForcing2 = forcingFuncs2.sym.tau_t - (frictionTraction - tau0);

	normalTractionForcing1 = forcingFuncs1.sym.tau_n - forcingFuncs2.sym.tau_n;
	normalTractionForcing2 = forcingFuncs2.sym.tau_n - forcingFuncs1.sym.tau_n;

	frictionForcing = elastic.friction.mms.computeFrictionForcing(f, V, psi, iS.friction);

	% ------- Compute state evolution forcing ---
	stateForcing = elastic.friction.mms.computeStateForcing(psi, V, f, iS.friction, t);

	% ---- Turn friction symfuns into function handles ------
	iS.friction.a = matlabFunctionSizePreserving( symfun(iS.friction.a(0,x,y), [x, y]) );
	iS.friction.b = matlabFunctionSizePreserving( symfun(iS.friction.b(0,x,y), [x, y]) );

	% -------- Fill mms struct --------------------
	mms = struct;

	u_1 = @(t,x,y)[forcingFuncs1.u1(t,x,y); forcingFuncs1.u2(t,x,y)];
	u_2 = @(t,x,y)[forcingFuncs2.u1(t,x,y); forcingFuncs2.u2(t,x,y)];
	mms.u = {u_1, u_2};

	ut_1 = @(t,x,y)[forcingFuncs1.u1t(t,x,y); forcingFuncs1.u2t(t,x,y)];
	ut_2 = @(t,x,y)[forcingFuncs2.u1t(t,x,y); forcingFuncs2.u2t(t,x,y)];
	mms.ut = {ut_1, ut_2};

	F1 = {forcingFuncs1.f1, forcingFuncs1.f2};
	F2 = {forcingFuncs2.f1, forcingFuncs2.f2};
	mms.F = {F1, F2};

	mms.psi = matlabFunctionSizePreserving(psi);
	mms.frictionForcing = frictionForcing;
	mms.stateForcing = stateForcing;

	mms.tractionForcing1 = matlabFunctionSizePreserving(tractionForcing1);
	mms.tractionForcing2 = matlabFunctionSizePreserving(tractionForcing2);

	mms.normalTractionForcing1 = matlabFunctionSizePreserving(normalTractionForcing1);
	mms.normalTractionForcing2 = matlabFunctionSizePreserving(normalTractionForcing2);

	iS.mms = mms;

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

	iS.domain = elastic.friction.input.domainStruct(...
										iS.domain, def, bc, faultBoundaryGroups);
	%------------------------------------------------

	% ------ Plotting ----------------------
	save('+elastic/+friction/+input/mms.mat', 'iS');

end
