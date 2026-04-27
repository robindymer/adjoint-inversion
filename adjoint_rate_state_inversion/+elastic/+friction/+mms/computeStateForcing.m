function F = computeStateForcing(psi, V, f, friction, t)

	f0 = friction.f0;
	V0 = friction.V0;
	a = friction.a;
	b = friction.b;
	L = friction.L;

	psit = diff(psi, t);

	switch friction.law
	case 'slip'
		f_ss = f0 - (b-a)*log( (abs(V)+1e-14)/V0 );
		g = -(abs(V)/L)*(f - f_ss);
	case 'aging'
		g = (b*V0/L)*(exp((f0-psi)/b)-abs(V)/V0);
	end

	F = psit - g;

	F = matlabFunctionSizePreserving(F);

end