% Hessian-vector product for fmincon's HessianMultiplyFcn, HdpMat = H*dpMat.
% Each column of dpMat is a perturbation direction Delta p, and costs one Hessian-vector product.
function HdpMat = hessianVector(parsVec, dpMat, optObj)

	% Set new parameters
	parsNative = optObj.parametersVectorToNative(parsVec);
	optObj.setParameters(parsNative);
	optObj.updateForwardDiscr();

	HdpMat = zeros(size(dpMat));
	for k = 1:size(dpMat, 2)
		% Set perturbation direction and compute Hessian-vector product
		optObj.setPerturbationDirection(optObj.parametersVectorToNative(dpMat(:,k)));
		hessVecNative = optObj.computeHessianVector();
		HdpMat(:,k) = optObj.gradientNativeToVector(hessVecNative);
	end

end
