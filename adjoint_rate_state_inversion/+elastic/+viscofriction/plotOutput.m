function plotOutput(inputFile, outputDir)
	default_arg('outputDir', []);

	% MODE argument 'plot' avoids deleting files if simulation is currently running
	[domain, friction, preStress, material, mms, flowLaw, plotting, tectonic, timeStepper, output]...
	 = inputFile(outputDir, 'plot');

	yearScaling = elastic.helpers.secondsPerYear();

	for i = 1:numel(output.fields)

		% Initial figure setup
		fh = figure;
		xLabel = 'time (years)';

		% Read from file
		fieldName = output.fields{i};
		[t, field, s] = loadField(output, fieldName);

		% Divide stresses by 1e6 to get MPa
		if contains(fieldName, 'Stress')
			field = field/1e6;
		end

		% Plot scalars
		switch fieldName
		case {'slipVelocity','shearStressMax','shearStressMean','psiMin','psiMax'}
			semilogy(t/yearScaling, field);

			% Set appropriate y label
			switch fieldName
			case 'slipVelocity'
				ylabel('Max slip velocity (m/s)');

			case 'shearStressMax'
				ylabel('Max shear stress (MPa)');

			case 'shearStressMean'
				ylabel('Mean shear stress (MPa)');

			case 'psiMin'
				ylabel('Min psi');

			case 'psiMax'
				ylabel('Max psi');
			end
		end

		% Plot fault fields
		if contains(fieldName, 'Profile')

			% Setup spatial axis
			xFault = s.dat.xFault;
		    switch plotting.faultCoordinate
		    case 'x'
		        faultCoord = xFault(:,1)/1e3;
		        faultCoordLabel = 'x (km)';
		    case 'y'
		        faultCoord = xFault(:,2)/1e3;
		        faultCoordLabel = 'y (km)';
		    end

		    % Manipulate fields if required
		    switch fieldName
		    case 'slipVelocityProfile'
		    	field = log10(abs(field));
		    end

		    % Plot
			surf(t/yearScaling, faultCoord, field)
			ylabel(faultCoordLabel)
			shading interp
			view(0,90)
			colorbar

			% Set appropriate title
			switch fieldName
			case 'shearStressProfile'
				title('Shear stress (MPa)');
			case 'normalStressProfile'
				title('Normal stress (MPa)');
			case 'slipProfile'
				title('Slip (m)');
			case 'slipVelocityProfile'
				title('log10 of slip velocity in m/s');
			case 'psiProfile'
				title('State variable psi');
			end
		end

		% Finalize first figure
		figure(fh);
		xlabel(xLabel);

		% If we have the state variable psi, also plot state variable theta
		if strcmp(fieldName, 'psiProfile')
			b = evalOnLine(xFault, friction.b);
			V0 = friction.V0;
			f0 = friction.f0;
			L = friction.L;
			theta = 0*field;
			for j = 1:length(t)
				theta(:,j) = elastic.friction.psiToTheta(field(:,j), b, V0, f0, L);
			end
			figure;
			surf(t/yearScaling, faultCoord, theta)
			shading interp
			view(0,90)
			colorbar
			xlabel(xLabel);
			ylabel(faultCoordLabel);
			title('State variable theta');
		end

		% If we have the slip profile, also plot slip contours
		if strcmp(fieldName, 'slipProfile')
			plotSlipContours(output, field, t, faultCoord, faultCoordLabel);
		end

		% Plot surface fields
		if contains(fieldName, 'surface')

			% Setup spatial axis
			xSurface = s.dat.xSurface;
	        surfaceCoord = xSurface(:,1)/1e3;
	        surfaceCoordLabel = 'x (km)';

		    % Plot
			surf(t/yearScaling, surfaceCoord, field)
			xlabel('time (years)')
			ylabel(surfaceCoordLabel)
			shading interp
			view(0,90)
			colorbar

			% Set appropriate title
			switch fieldName
			case 'surfaceDisplacementX'
				title('Horizontal surface displacement (m)');
			case 'surfaceDisplacementY'
				title('Vertical surface displacement (m)');
			end
		end

	end

end

function fval = evalOnLine(X, f)
	if isa(f, 'function_handle')
		fval = f(X(:,1), X(:,2));
	else
		% Assume that f is a scalar
		fval = 0*X(:,1) + f;
	end
end

function plotSlipContours(output, slip, tSlip, faultCoord, faultCoordLabel)

	yearScaling = 365.25*24*3600;

	% load slipVelocity
	[t, slipVelocity] = loadField(output, 'slipVelocity');

	% Get unique time points
	[t, slipVelocity] = getUniqueTimePoints(t, slipVelocity);
	[tSlip, slip] = getUniqueTimePoints(tSlip, slip);

	% Set parameters
	cutOff = 1e-1;
	dt_slow = 20*yearScaling;
	dt_fast = 2;
	seismicColor = 'r';
	interseismicColor = 'b';

	% Determine if slip is seismic or aseismic
	quakeIndicator = double(slipVelocity > cutOff);
	tNew = newTimeAxis(quakeIndicator, t, dt_slow, dt_fast);
	% Interpolate quakeIndicator
	quakeIndicator = interp1(t, quakeIndicator, tNew);

	% Interpolate slip
	[m, ~] = size(slip);
	slipNew = [];
	for i = 1:m
		slipNew = [slipNew; interp1(tSlip, slip(i,:), tNew)];
	end

	% Plot curves
	figure;
	for i = 1:length(tNew)
		if quakeIndicator(i) > 1/2
			marker = seismicColor;
		else
			marker = interseismicColor;
		end
		plot(slipNew(:,i), faultCoord, marker);
		hold on;
	end
	xlabel('slip (m)')
	ylabel(faultCoordLabel)
end

function [t, field, s] = loadField(output, fieldName)

	fileName = [output.dir, fieldName, '.mat'];
	s = load(fileName);
	field = s.dat.field;
	t = s.dat.t;

end

function [t, field] = getUniqueTimePoints(t, field)
	[t, uniqueIndices] = unique(t);
	field = field(:, uniqueIndices);
end

function t = newTimeAxis(indicator, t, dt_1, dt_2)

	tTimes = transitionTimes(indicator, t);
	t = twoSpeedAxis(tTimes, dt_1, dt_2);

end

function transitionTimes = transitionTimes(indicator, t)

	transitionTimes = 0;
	for i = 2:length(t)
		if indicator(i) ~= indicator(i-1)
			transitionTimes = [transitionTimes, t(i)];
		end
	end

end

function t = twoSpeedAxis(transitionTimes, dt_1, dt_2)

	t = [];
	for i = 2:length(transitionTimes)
		if mod(i,2) == 0
			dt = dt_1;
		else
			dt = dt_2;
		end
		gap = transitionTimes(i) - transitionTimes(i-1);
		N = ceil(gap/dt);
		dt = gap/N;
		tLocal = (transitionTimes(i-1):dt:transitionTimes(i));
		t = [t, tLocal];
	end

end