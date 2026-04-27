function [domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = default()


	mms = [];

	% ----- Friction parameters ------------------
	friction = struct;
	friction.a = 0.03; 			% Direct effect [ ]
	friction.b = 0.02;			% State evolution effect [ ]
	friction.L = 0.25;			% State evolution distance [ m ]
	friction.V0 = 1e-6;			% Ref velocity [ m/s ]
	friction.f0 = 0.6;			% Ref friction coeff for steady sliding [ ]

	friction.coefficient = 'standard';
	friction.coefficientSteadyState = 'standard';
	friction.law = 'slip';

	% Radiation damping
	friction.eta = 0;
	%----------------------------------------------

	% ---- Flow law ------
	flowLaw = struct;
	eta = struct;
    eta.type = 'nonlinear';
    eta.nonlinFun = @elastic.visco.flowLawFuns.basic;
    flowLaw.eta = eta;

	% --- Domain -------
	domain = struct;
	domain.intfTypes = [];

	% --- Material parameters -------
	material = struct;

	% --- Prestress -------
	preStress = struct;

	% --- Plotting -------
	plotting = struct;
	plotting.setupFun = @elastic.viscofriction.plotting.basic;
	plotting.faultCoordinate = 'x';

	% --- Tectonic plate movement -------
	tectonic = [];

	% --- Time-stepping ---
	timeStepper = struct;
	timeStepper.type = @timeSteppers.EmbeddedRungeKutta;
	timeStepper.order = 3;
	timeStepper.relTol = 1e-3;
	timeStepper.dt = 1e-2;
	timeStepper.solver = 'lu';
	timeStepper.hardware = 'cpu';
	timeStepper.elasticOperator = 'sparse';
	timeStepper.numThreads = 1;

	% --- Quantities to output ----
	output = [];

end