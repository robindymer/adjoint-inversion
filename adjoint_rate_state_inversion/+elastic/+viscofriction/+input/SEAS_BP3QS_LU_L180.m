function [domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = SEAS_BP3QS_LU_L180(outputDir, MODE, dipAngle, faultMotion)
	% - outputDir: must include / at the end.
	% - MODE: 'startSimulation' or 'plot'

	default_arg('outputDir', 'SEAS_BP3QS_LU/');
	default_arg('MODE', 'startSimulation');
	default_arg('dipAngle', pi/2);
	default_arg('faultMotion', 'normal');

	% Start with options from input.SEAS_BP3QS_LU
	[domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = elastic.viscofriction.input.SEAS_BP3QS_LU(outputDir, MODE, dipAngle, faultMotion);

	% Change domain size
	Wf = 40*1e3;
	Lx = 180*1e3;
	Ly = 180*1e3;

	stretchingDist = [];
	stretchingRatio = 20;

	def = domains.DippingFault(dipAngle, Wf, Lx, Ly, stretchingDist, stretchingDist, stretchingRatio, stretchingRatio);
	boundaryGroups = def.boundaryGroups;

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

end

