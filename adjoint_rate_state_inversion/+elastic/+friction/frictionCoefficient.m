function [f, frictionPar] = frictionCoefficient(V, psi, frictionPar)

	a = frictionPar.a;
	b = frictionPar.b;
	L = frictionPar.L;
	V0 = frictionPar.V0;
	f0 = frictionPar.f0;

	% Using Psi (dimensionless) as state variable
	switch frictionPar.coefficient
	case 'standard'
		f = a.*asinh(abs(V)/(2*V0).*exp(psi./a));

		% Return function handle to avoid repeating expensive switch statement
		frictionPar.coefficient = @(V, psi, fp) fp.a.*asinh(abs(V)/(2*fp.V0).*exp(psi./fp.a));
	otherwise
		error('Not implemented');
	end

end