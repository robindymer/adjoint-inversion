function e = computeError(u, t, discr, mms)

	u_exact = multiblock.evalOn(discr.grid, mms.u, t);
	e = discr.compareSolutions(u, u_exact);

end