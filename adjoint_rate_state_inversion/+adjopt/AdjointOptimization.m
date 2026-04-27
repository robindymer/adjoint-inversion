% Adjoint optimization

classdef AdjointOptimization < handle

	properties (Abstract)

		forwardDiscr
		adjointDiscr
		pars 			% Contains current parameter state

		m 				% Number of spatial grid points
		order			% Order of accuracy
		T 				% Final time
		H           	% Total space quadrature, for both components
		Ht 				% Time quadrature, a column vector of weights.
		dim				% Number of spatial dimensions

	end

	methods (Abstract)

		% PDE discretizations
		runForward(obj, T)
		runAdjoint(obj, T)
		updateForwardDiscr(obj, pars)
		updateAdjointDiscr(obj, pars)

		% Optimization-related computations
		M = computeMisfit(obj);
		grad = computeGradient(obj)
		grad = computeGradientFD(obj, deltaG)
		difference = compareGradients(obj, grad1, grad2)
		relnorm = gradientNorm(obj, grad) % Useful with e.g golden section search, may become
													 % obsolete.
		updateParameters(obj, pars, direction, steplength)

		% Sets the values of the inversion parameters to those prescribed by parsNative.
		setParameters(obj, parsNative)

		% For converting parameter sets from native to vector format and vice versa
		parsVec = parametersNativeToVector(obj, parsNative)
		parsNative = parametersVectorToNative(obj, parsVec)

		% For converting gradients from native to vector format and vice versa
		gradVec = gradientNativeToVector(obj, gradNative)
		gradNative = gradientVectorToNative(obj, gradVec)

		% Plotting
		[updatePar, updateRec, updateMis, figure_handle] = ...
		setupPlot(obj, trueParset, initialParset, M, M0, iter)

	end

end