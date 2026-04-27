function defaultInput = default()


	mms = [];

	% ----- Friction parameters ------------------
	friction = struct;
	friction.a = 0.03; 			% Direct effect [ ]
	friction.b = 0.02;			% State evolution effect [ ]
	friction.L = 0.25;			% State evolution distance [ m ]
	friction.V0 = 1e-6;			% Ref velocity [ m/s ]
	friction.f0 = 0.6;			% Ref friction coeff for steady sliding [ ]

	friction.coefficient = 'Rice2001';
	friction.coefficientSteadyState = 'Rice2001';
	friction.law = 'slip';
	%----------------------------------------------

	% --- Domain -------
	domain = struct;
	domain.intfTypes = [];

	% --- Material parameters -------
	material = struct;

	% --- Prestress -------
	preStress = struct;

	% --- Time-stepping ---
	timeStepper = struct;
	timeStepper.type = 'EmbeddedRK';
	timeStepper.order = 3;
	timeStepper.relTol = 1e-3;
	timeStepper.dt = 1e-2;
	timeStepper.CFL = 1;

	% --- Plotting -------
	plotting = struct;

	% --- Create return value ---
	defaultInput = struct;
	defaultInput.domain = domain;
	defaultInput.friction = friction;
	defaultInput.preStress = preStress;
	defaultInput.material = material;
	defaultInput.mms = mms;
	defaultInput.plotting = plotting;
	defaultInput.timeStepper = timeStepper;

	% --- Solver type ---
	defaultInput.solver = struct;
	defaultInput.solver.mex = false;
	defaultInput.solver.check_ev_flag = false;

end