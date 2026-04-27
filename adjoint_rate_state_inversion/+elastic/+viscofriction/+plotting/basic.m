% Sets up a plot of the discretisation
% update is a function_handle accepting (t,u)
function [update, figure_handle] = basic(discr, ops, U0, plotting, mms, tectonic)

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

    h1 = subplot(4,1,1);
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

    h2 = subplot(4,1,2);
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

    a = gca;

    h3 = subplot(4,1,3);
    L3 = plot(ops.xFault(:,1), (ops.Epsi)'*U0);
    xlabel('x (m)')
    ylabel('State variable (psi)')
    if ~isempty(plotting.xlims)
        xlim(plotting.xlims);
    end

    h4 = subplot(4,1,4);
    L4 = plot(ops.xFault(:,1), U_tangent);
    xlabel('x (m)')
    ylabel('Cumulative slip on fault (m)')
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

        v1 = u(1:2:end-1);
        v2 = u(2:2:end);
        if ishandle(a)
            title(a,sprintf('u_2, T = %5.3f years',t/(365.25*24*3600)))
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
        drawnow;
    end
    update = @(t,U) update_fun(t,U);
end

function drawLine(X)
    line(X(:,1), X(:,2), 0*X(:,1) + 1e10, 'color', 'black', 'linewidth', 1.5);
end