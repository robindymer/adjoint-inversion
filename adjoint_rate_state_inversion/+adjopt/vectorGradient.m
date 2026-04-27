function [misfit, gradVec] = vectorGradient(parsVec, optObj)

	% Set new parameters
	parsNative = optObj.parametersVectorToNative(parsVec);
	optObj.setParameters(parsNative);
	optObj.updateForwardDiscr();

	% Compute gradient
	gradNative = optObj.computeGradient();
	gradVec = optObj.gradientNativeToVector(gradNative);

	% Compute misfit
	misfit = optObj.computeMisfit();

end