function V = slipRate(psi, frictionPar)

	V = 0*psi;
	for i = 1:length(psi)
		res = @(V) elastic.friction.frictionLaw(V, psi, frictionPar);

		VMin = 0;
		VMax = ?;
		VBounds = sort([Vmin, Vmax]);

		V(i) = elastic.friction.bisection(res, VBounds);
	end

end