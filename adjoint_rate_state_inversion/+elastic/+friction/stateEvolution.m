function psit = stateEvolution(V, psi, f, frictionPar)

    V = abs(V);

    a = frictionPar.a;
    b = frictionPar.b;
    f0 = frictionPar.f0;
    V0 = frictionPar.V0;
    L = frictionPar.L;


	% This is using Psi (dimensionless, Psi = f0+b*log(V0*theta/L))
	switch frictionPar.law
	case 'slip'
		f_ss = elastic.friction.frictionCoefficientSteadyState(V, frictionPar);
		psit = -(V./L).*(f - f_ss);
	case 'aging'
		psit = (b.*V0./L).*(exp((f0-psi)./b)-V./V0);
	otherwise
		error('Not implemented');
	end

end