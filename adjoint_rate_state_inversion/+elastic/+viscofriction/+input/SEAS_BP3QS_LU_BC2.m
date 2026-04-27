function [domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = SEAS_BP3QS_LU_BC2(outputDir, MODE, dipAngle, faultMotion)
	% - outputDir: must include / at the end.
	% - MODE: 'startSimulation' or 'plot'

	default_arg('outputDir', 'SEAS_BP3QS_LU/');
	default_arg('MODE', 'startSimulation');
	default_arg('dipAngle', pi/2);
	default_arg('faultMotion', 'normal');

	% Start with options from input.SEAS_BP3QS_LU
	[domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = elastic.viscofriction.input.SEAS_BP3QS_LU(outputDir, MODE, dipAngle, faultMotion);
	Vp = 1e-9;


	% ------ Change BC related to tectonic plate movement.
	% Same BCs, but add sliding on side  boundaries -------
	def = domain.def;
	boundaryGroups = def.boundaryGroups;

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


	case 'normal'
		tectonic{1}.v = Vp;
	end

	% Switch sign compared to fault sliding.
	tectonic{2}.v = -tectonic{1}.v;
	%------------------------------------------------

end

