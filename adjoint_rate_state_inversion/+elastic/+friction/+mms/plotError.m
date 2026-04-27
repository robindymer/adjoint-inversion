function plotError(u, psi, t, discr, ops, mms)

	u_exact = multiblock.evalOn(discr.grid, mms.u, t);

	figure;
	elastic.cmap;
	e = u-u_exact;
	e1 = e(1:2:end-1);
	s = multiblock.Surface(discr.grid, e1);
	colorbar;
	shading interp
	cmax = max(abs(e1));
    caxis([-cmax, cmax]);
    xlabel('x')
    ylabel('y')
    title('Error in u_1')

    figure;
    plot(ops.xFault(:,1), psi);
    hold on
    psi_exact = elastic.helpers.evalOnLine(ops.xFault, @(x,y) mms.psi(t,x,y));
    plot(ops.xFault(:,1), psi_exact);
    xlabel('x coordinate along fault')
    ylabel('\psi')


end