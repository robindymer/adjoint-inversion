function f = frictionCoefficientSteadyState(V, frictionPar)

	V = abs(V);
	a = frictionPar.a;
	b = frictionPar.b;
	V0 = frictionPar.V0;
	f0 = frictionPar.f0;

	switch frictionPar.coefficient
	case 'standard'
		f = f0 - (b-a).*log( max(V, 1e-14)/V0 );
	otherwise
		error('Not implemented');
	end

end