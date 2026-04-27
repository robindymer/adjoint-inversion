function U0 = setupFullyDynamicInitialData(g, friction, ops, mms)

	if isempty(mms)
		zeroFun = @(x,y) [0*x; 0*x];
		u0 = grid.evalOn(g, zeroFun);
		v0 = grid.evalOn(g, zeroFun);
	else
		u0 = multiblock.evalOn(g, mms.u, 0);
		v0 = multiblock.evalOn(g, mms.ut, 0);
	end
	psi0 = elastic.helpers.evalOnLine(ops.xFault, friction.psi0);
	U0 = ops.Eu*u0 + ops.Ev*v0 + ops.Epsi*psi0;

end