% Accepts scalars or function handles and returns struct with function handles
function material = materialStruct(rho, vs, vp)

	material = struct;

	if isa(rho, 'function_handle')
		material.rhoFun = rho;
	else
		material.rhoFun = @(x,y) 0*x + rho;
	end

	if isa(vs, 'function_handle')
		material.vsFun = vs;
	else
		material.vsFun = @(x,y) 0*x + vs;
	end

	if isa(vp, 'function_handle')
		material.vpFun = vp;
	else
		material.vpFun = @(x,y) 0*x + vp;
	end

	[lambdaFun, muFun] = elastic.speedsToModuli(material.vpFun, material.vsFun, material.rhoFun);
	material.lambdaFun = lambdaFun;
	material.muFun = muFun;

	material.CFun = elastic.isotropicStiffnessTensor(lambdaFun, muFun);

end