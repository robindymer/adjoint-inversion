function U_t = rhsFullyDynamicIntfForcing(t, U, ops, friction, preStress, mms, mexOps)

	% u: 	Displacements
	% v: 	Particle velocities, du/dt.
	% psi:  Fault state variable

	u = (ops.Eu)'*U;
	v = (ops.Ev)'*U;
	psi = (ops.Epsi)'*U;

	% Compute slip velocity on fault (note that the tangential directions are opposite)
	v_tangential_1 = (ops.eTangentFault1)'*v;
	v_tangential_2 = (ops.eTangentFault2)'*v;
	V = (v_tangential_1 + v_tangential_2);

	% Compute friction coefficient
	f = elastic.friction.frictionCoefficient(V, psi, friction);
	if ~isempty(mms)
		f = f + mms.frictionForcing(t);
	end


	% --- Compute frictional force ---

	% Perturbation in normal stress as average from two sides of fault;
	sigma = -1/2*(ops.tauNormalFault1' + ops.tauNormalFault2')*u;

	% Compressive normal stress
	sigmaTot = max(preStress.sigma0 + sigma, 0);
	shearStress = -sigmaTot.*f.*sign(V) - preStress.shear;
	%----------------------------------

	% Evolve displacements and velocities
	u_t = v;
	if ~isempty(mexOps)
		v_t = ops.D*u + elastic.mex.elasticOperatorMultiblock(u, mexOps.g, mexOps.RHOJi, mexOps.PHI, mexOps.ops);
	else
		v_t = ops.D*u;
	end
	v_t = v_t + (ops.SFault1+ops.SFault2)*shearStress;

	if ~isempty(mms)
		v_t = v_t + mms.forcing(t);
		v_t = v_t + ops.SFault1*mms.tractionForcing1(t) + ops.SFault2*mms.tractionForcing2(t);
		% v_t = v_t + ops.SNormalFault1*mms.normalTractionForcing1(t) + ops.SNormalFault2*mms.normalTractionForcing2(t);
	end

	% Evolve the state variable
	psi_t = elastic.friction.stateEvolution(V, psi, f, friction);
	if ~isempty(mms)
		psi_t = psi_t + mms.stateForcing(t);
	end

	U_t = ops.Eu*u_t + ops.Ev*v_t + ops.Epsi*psi_t;
end