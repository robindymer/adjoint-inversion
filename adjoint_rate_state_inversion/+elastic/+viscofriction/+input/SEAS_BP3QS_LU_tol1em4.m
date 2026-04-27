function [domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = SEAS_BP3QS_LU_tol1em4(outputDir, MODE, dipAngle, faultMotion)
	% - outputDir: must include / at the end.
	% - MODE: 'startSimulation' or 'plot'

	default_arg('outputDir', 'SEAS_BP3QS_LU/');
	default_arg('MODE', 'startSimulation');
	default_arg('dipAngle', pi/2);
	default_arg('faultMotion', 'normal');

	% Start with options from input.SEAS_BP3QS_LU
	[domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = elastic.viscofriction.input.SEAS_BP3QS_LU(outputDir, MODE, dipAngle, faultMotion);

	% Change tolerance
	timeStepper.relTol = 1e-4;

end

