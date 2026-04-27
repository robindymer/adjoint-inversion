% Sets up a plot of the discretisation
% update is a function_handle accepting (t,u)
function [update, figure_handle] = stresses(discr, ops, U0, plotting, mms, tectonic)

    plotting_default = struct;
    plotting_default.xlims = [];
    default_struct('plotting', plotting_default);

    g = discr.grid;
    figure_handle = figure();
    elastic.cmap;

    gamma0 = (ops.Egamma)'*U0;
    delta0 = (ops.Edelta)'*U0;

    % Solve mechanical equilibrium eq to get u
    RHS = ops.Mech_gamma*gamma0 + ops.slipInserter_u*delta0;
    if ~isempty(mms)
        RHS = RHS + mms.mechForcing(0);
    end
    u = elastic.helpers.solveWithLU(ops.Mech_u_factorized, RHS);

    u1 = u(1:2:end-1);
    u2 = u(2:2:end);

    U_tangent = (ops.Edelta)'*U0;
    tau1 = -((ops.tauShearFault1_u)'*u + (ops.tauShearFault1_gamma)'*gamma0);
    tau2 = -((ops.tauShearFault2_u)'*u + (ops.tauShearFault2_gamma)'*gamma0);
    sigma1 = -((ops.tauNormalFault1_u)'*u + (ops.tauNormalFault1_gamma)'*gamma0);
    sigma2 = -((ops.tauNormalFault2_u)'*u + (ops.tauNormalFault2_gamma)'*gamma0);

    h1 = subplot(4,2,1);
    Sur1 = multiblock.Surface(g, u1);
    xlabel('x')
    ylabel('y')
    title('u_1')
    shading interp
    colorbar
    if ~isempty(plotting.xlims)
        xlim(plotting.xlims);
    end
    xlims1 = xlim;
    ylims1 = ylim;
    if ~isempty(plotting.xFault)
        drawLine(plotting.xFault);
    end
    if ~isempty(plotting.xSurface)
        drawLine(plotting.xSurface);
    end

    h2 = subplot(4,2,2);
    Sur2 = multiblock.Surface(g, u2);
    xlabel('x')
    ylabel('y')
    title('u_2')
    shading interp
    colorbar
    if ~isempty(plotting.xlims)
        xlim(plotting.xlims);
    end
    xlims2 = xlim;
    ylims2 = ylim;
    if ~isempty(plotting.xFault)
        drawLine(plotting.xFault);
    end
    if ~isempty(plotting.xSurface)
        drawLine(plotting.xSurface);
    end

    xFault = ops.xFault;
    switch plotting.faultCoordinate
    case 'x'
        faultCoord = xFault(:,1);
        faultCoordLabel = 'x (m)';
    case 'y'
        faultCoord = xFault(:,2);
        faultCoordLabel = 'y (m)';
    end

    h3 = subplot(4,2,3);
    L3 = plot(faultCoord, (ops.Epsi)'*U0);
    xlabel(faultCoordLabel)
    ylabel('State variable (psi)')
    if ~isempty(plotting.xlims)
        xlim(plotting.xlims);
    end
    timeTitleAxis = gca;

    h4 = subplot(4,2,4);
    L4 = plot(faultCoord, U_tangent);
    xlabel(faultCoordLabel)
    ylabel('Cumulative slip on fault (m)')
    if ~isempty(plotting.xlims)
        xlim(plotting.xlims);
    end

    h5 = subplot(4,2,5);
    L5 = plot(faultCoord, tau1/1e6);
    xlabel(faultCoordLabel)
    ylabel('Shear stress perturbation on side 1 (MPa)')
    if ~isempty(plotting.xlims)
        xlim(plotting.xlims);
    end

    h6 = subplot(4,2,6);
    L6 = plot(faultCoord, tau2/1e6);
    xlabel(faultCoordLabel)
    ylabel('Shear stress perturbation on side 2 (MPa)')
    if ~isempty(plotting.xlims)
        xlim(plotting.xlims);
    end

    h7 = subplot(4,2,7);
    L7 = plot(faultCoord, sigma1/1e6);
    xlabel(faultCoordLabel)
    ylabel('Compressive normal stress perturbation on side 1 (MPa)')
    if ~isempty(plotting.xlims)
        xlim(plotting.xlims);
    end

    h8 = subplot(4,2,8);
    L8 = plot(faultCoord, sigma2/1e6);
    xlabel(faultCoordLabel)
    ylabel('Compressive normal stress perturbation on side 2 (MPa)')
    if ~isempty(plotting.xlims)
        xlim(plotting.xlims);
    end

    function update_fun(t,U)
        gamma = (ops.Egamma)'*U;
        delta = (ops.Edelta)'*U;

        % Solve mechanical equilibrium eq to get u
        RHS = ops.Mech_gamma*gamma + ops.slipInserter_u*delta;
        if ~isempty(mms)
            RHS = RHS + mms.mechForcing(t);
        end
        % Add tectonic plate movement
        if ~isempty(tectonic)
            displ = tectonic.u0 + tectonic.v*t;
            RHS = RHS + ops.tectonicInserter_u * displ;
        end
        u = elastic.helpers.solveWithLU(ops.Mech_u_factorized, RHS);
        tau1 = -((ops.tauShearFault1_u)'*u + 0*(ops.tauShearFault1_gamma)'*gamma);
        tau2 = -((ops.tauShearFault2_u)'*u + 0*(ops.tauShearFault2_gamma)'*gamma);
        sigma1 = -((ops.tauNormalFault1_u)'*u + 0*(ops.tauNormalFault1_gamma)'*gamma);
        sigma2 = -((ops.tauNormalFault2_u)'*u + 0*(ops.tauNormalFault2_gamma)'*gamma);

        v1 = u(1:2:end-1);
        v2 = u(2:2:end);
        if ishandle(timeTitleAxis)
            title(timeTitleAxis, sprintf('T = %8.6f years',t/(365.25*24*3600)))
        end
        Sur1.ZData = v1;
        Sur1.CData = v1;
        Sur2.ZData = v2;
        Sur2.CData = v2;

        xlim(h1, xlims1);
        ylim(h1, ylims1);
        xlim(h2, xlims2);
        ylim(h2, ylims2);

        cmax = max(abs(Sur1.ZData));
        caxis(h1, [-cmax, cmax]);
        cmax = max(abs(Sur2.ZData));
        caxis(h2, [-cmax, cmax]);
        L3.YData = (ops.Epsi)'*U;
        L4.YData = delta;
        L5.YData = tau1/1e6;
        L6.YData = tau2/1e6;
        L7.YData = sigma1/1e6;
        L8.YData = sigma2/1e6;
        drawnow;
    end
    update = @(t,U) update_fun(t,U);
end

function drawLine(X)
    line(X(:,1), X(:,2), 0*X(:,1) + 1e10, 'color', 'black', 'linewidth', 1.5);
end