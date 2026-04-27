% Sets up a plot of the discretisation
% update is a function_handle accepting (t,u)
function [update, figure_handle] = stateVariable(discr, ops, U0, plotting)

    plotting_default = struct;
    plotting_default.xlims = [];
    default_struct('plotting', plotting_default);

    g = discr.grid;
    figure_handle = figure();
    elastic.cmap;

    v0 = (ops.Ev)'*U0;
    v0_1 = discr.Ecomp{1}'*v0;
    v0_2 = discr.Ecomp{2}'*v0;

    U_tangent = ((ops.eTangentFault1)'+(ops.eTangentFault2)')*(ops.Eu)'*U0;

    h1 = subplot(4,1,1);
    Sur1 = multiblock.Surface(g, v0_1);
    xlabel('x')
    ylabel('y')
    title('v_x')
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
    Sur2 = multiblock.Surface(g, v0_2);
    xlabel('x')
    ylabel('y')
    title('v_y')
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
    ylabel('State variable psi)')
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

    function update_fun(t,U,E)
        v = (ops.Ev)'*U;
        v1 = E{1}'*v;
        v2 = E{2}'*v;
        U_tangent = ((ops.eTangentFault1)'+(ops.eTangentFault2)')*(ops.Eu)'*U;
        if ishandle(a)
            title(a,sprintf('v_y, T = %.3f',t))
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
        L4.YData = U_tangent;
        drawnow;
    end
    update = @(t,U) update_fun(t,U,discr.Ecomp);
end

function drawLine(X)
    line(X(:,1), X(:,2), 0*X(:,1) + 1e10, 'color', 'black', 'linewidth', 1.5);
end