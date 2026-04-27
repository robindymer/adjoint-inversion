function V = slipVelocity(friction, shearStress, normalStress, psi, mms)
	% ----- Solve force balance for slip velocity ------------

	frictionScalar = struct;
	frictionScalar.V0 = friction.V0;
	frictionScalar.f0 = friction.f0;
	frictionScalar.coefficient = friction.coefficient;
	frictionScalar.L = friction.L;

	V = zeros(size(psi));

	options = optimset('TolX',1e-30);

	for i = 1:length(V)
		frictionScalar.a = friction.a(i);
		frictionScalar.b = friction.b(i);
		frictionScalar.eta = friction.eta(i);

		if (normalStress(i) <= 0)
		    % Normal stress is tensile.
			% Set shear strength equal to zero and solve for V.
			V(i) = shearStress(i)/frictionScalar.eta;
		else
			if ~isempty(mms)
				forcing = mms.frictionForcing(mms.t);
				frictionForcing = forcing(i);
			else
				frictionForcing = [];
			end

			% Bound solution, assuming maximum |V| is determined
			% by friction dropping to zero (zero shear strength)
			if shearStress(i) > 0
				bracket = [0 shearStress(i)/frictionScalar.eta];
			else
				bracket = [shearStress(i)/frictionScalar.eta 0];
			end

			% Use previous solution as initial guess
			if i>1
				V0 = V(i-1);
			else
				V0 = bracket;
			end
			[V(i), ~, flag] = fzero( @(V) elastic.friction.frictionLaw(V, psi(i), shearStress(i), normalStress(i), frictionScalar, frictionForcing), V0, options);

			% If fzero found a solution outside the bracket, try again with bracket.
			if (V(i) <= bracket(1)) || (V(i) >= bracket(2))
				[V(i), ~, flag] = fzero( @(V) elastic.friction.frictionLaw(V, psi(i), shearStress(i), normalStress(i), frictionScalar, frictionForcing), bracket, options);
			end

			if flag~=1
				error('fzero failed');
			end

		end

	end

end