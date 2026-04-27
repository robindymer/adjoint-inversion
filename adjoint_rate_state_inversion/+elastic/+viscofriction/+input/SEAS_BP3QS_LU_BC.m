function [domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = SEAS_BP3QS_LU_BC(outputDir, MODE, dipAngle, faultMotion)
	% - outputDir: must include / at the end.
	% - MODE: 'startSimulation' or 'plot'

	default_arg('outputDir', 'SEAS_BP3QS_LU/');
	default_arg('MODE', 'startSimulation');
	default_arg('dipAngle', pi/2);
	default_arg('faultMotion', 'normal');

	% Start with options from input.SEAS_BP3QS_LU
	[domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output] = elastic.viscofriction.input.SEAS_BP3QS_LU(outputDir, MODE, dipAngle, faultMotion);
	Vp = 1e-9;

	%------- Domain change BC -----------
	def = domain.def;
	boundaryGroups = def.boundaryGroups;
	faultBoundaryGroups = domain.faultBoundaryGroups;

	% Traction BC (changed from parent input file) on west and east sides
	bcSidesN.boundary = boundaryGroups.sides;
	bcSidesN.type = {'n', 't'};

	bcSidesT.boundary = boundaryGroups.sides;
	bcSidesT.type = {'t', 't'};

	% Displacement BC (changed from parent input file) on bottom boundary
	bcBottomN.boundary = boundaryGroups.bottom;
	bcBottomN.type = {'n', 'd'};

	bcBottomT.boundary = boundaryGroups.bottom;
	bcBottomT.type = {'t', 'd'};

	% Free BC on surface
	bcSurfaceN.boundary = boundaryGroups.surface;
	bcSurfaceN.type = {'n', 't'};

	bcSurfaceT.boundary = boundaryGroups.surface;
	bcSurfaceT.type = {'t', 't'};

	bc = {bcSidesN, bcSidesT, bcBottomT, bcBottomN, bcSurfaceT, bcSurfaceN};

	domain = elastic.viscofriction.input.domainStruct(domain, def, bc, faultBoundaryGroups);

	% Add surface boundary groups (necessary for surface output)
	domain.surfaceBoundaryGroups = boundaryGroups.surface;

	% ------ Change BC related to tectonic plate movement -------
	% Bottom boundary moves with plate rate
	tectonic = cell(5, 1);
	for i = 1:numel(tectonic)
		tectonic{i} = struct;
	end

	% Fault sliding
	tectonic{1}.interface = struct;
	tectonic{1}.interface.boundaryGroups = boundaryGroups.tectonic;
	tectonic{1}.u0 = 0;

	% Bottom boundary, horizontal (sign for normal faulting)
	tectonic{2}.bc = struct;
	tectonic{2}.u0 = 0;
	tectonic{2}.bc.boundary = boundaryGroups.bottomRight;
	tectonic{2}.bc.type = {1, 'd'};
	tectonic{2}.v = cos(dipAngle)*Vp;

	tectonic{3}.bc = struct;
	tectonic{3}.u0 = 0;
	tectonic{3}.bc.boundary = boundaryGroups.bottomLeft;
	tectonic{3}.bc.type = {1, 'd'};
	tectonic{3}.v = -cos(dipAngle)*Vp;

	% Bottom boundary, vertical (sign for normal faulting)
	tectonic{4}.bc = struct;
	tectonic{4}.u0 = 0;
	tectonic{4}.bc.boundary = boundaryGroups.bottomRight;
	tectonic{4}.bc.type = {2, 'd'};
	tectonic{4}.v = -sin(dipAngle)*Vp;

	tectonic{5}.bc = struct;
	tectonic{5}.u0 = 0;
	tectonic{5}.bc.boundary = boundaryGroups.bottomLeft;
	tectonic{5}.bc.type = {2, 'd'};
	tectonic{5}.v = sin(dipAngle)*Vp;

	switch faultMotion
	case {'reverse', 'thrust'}
		tectonic{1}.v = -Vp;

		for i = 2:numel(tectonic)
			tectonic{i}.v = -tectonic{i}.v;
		end

	case 'normal'
		tectonic{1}.v = Vp;
	end


	%------------------------------------------------

end

