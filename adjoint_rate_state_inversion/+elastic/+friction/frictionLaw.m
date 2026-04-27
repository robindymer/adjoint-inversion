function res = frictionLaw(V, psi, shearStress, normalStress, frictionPar, frictionForcing)

	% Compute friction coefficient
	f = frictionPar.coefficient(V, psi, frictionPar);
	if ~isempty(frictionForcing)
		f = f + frictionForcing;
	end

	% --- Compute residual in force balance, including radiation damping ---
	% This uses normalStress>0 in compression convention, and
	% friction coefficient f always >=0, using sign(V) to capture
	% information about direction of sliding.
	% For V=0 it reduces to res = shearStress.
	res = shearStress - normalStress.*f.*sign(V) - frictionPar.eta.*V;

end