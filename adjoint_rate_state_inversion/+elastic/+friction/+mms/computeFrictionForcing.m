function F = computeFrictionForcing(f, V, psi, friction)

	a = friction.a;
	V0 = friction.V0;

	switch friction.coefficient
	case 'standard'
		F = f - a*asinh(abs(V)/(2*V0) * exp(psi/a));
	otherwise
		error('Not implemented');
	end

	F = matlabFunctionSizePreserving(F);

end