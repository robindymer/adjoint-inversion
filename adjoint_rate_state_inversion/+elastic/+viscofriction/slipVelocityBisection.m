function V = slipVelocityBisection(friction, shearStress, normalStress, psi, mms)
	% ----- Solve force balance for slip velocity ------------

	N = length(psi);
	V = zeros(size(psi));
	bracket = zeros(N, 2);

	tensile = normalStress <= 0;
	compressive = ~tensile;

	% If normal stress is tensile, set shear strength equal to zero and solve for V.
	V(tensile) = shearStress(tensile)./friction.eta(tensile);

	if ~isempty(mms)
		forcing = mms.frictionForcing(mms.t);
		frictionForcing = forcing(compressive);
	else
		frictionForcing = [];
	end

	% Bound solution, assuming maximum |V| is determined
	% by friction dropping to zero (zero shear strength)
	j = shearStress > 0;
	bracket(j, 1) = 0;
	bracket(j, 2) = shearStress(j)./friction.eta(j);
	bracket(~j, 1) = shearStress(~j)./friction.eta(~j);
	bracket(~j, 2) = 0;

	% Use bisection to solve for V where normal stress is compressive.
	tol = 1e-30;
	friction.eta = friction.eta(compressive);
	F = @(V) elastic.friction.frictionLaw(V, psi(compressive), shearStress(compressive), normalStress(compressive), friction, frictionForcing);
	V(compressive) = elastic.helpers.vectorBisection(F, bracket, tol);

end