% Sets up a plot of the discretisation
% update is a function_handle accepting (t,u)
function [update, figure_handle] = basic(discr, ops, U0)

    g = discr.grid;
    X = ops.xFault;
    figure_handle = figure();
    elastic.cmap;

    v0 = (ops.Ev)'*U0;
    v0_1 = discr.Ecomp{1}'*v0;
    v0_2 = discr.Ecomp{2}'*v0;

    U_tangent = ((ops.eTangentFault1)'+(ops.eTangentFault2)')*(ops.Eu)'*U0;

    h1 = subplot(3,1,1);
    Sur1 = multiblock.Surface(g, v0_1);
    xlabel('x')
    ylabel('y')
    title('v_x')
    shading interp
    colorbar
    xlims1 = xlim;
    ylims1 = ylim;
    line(X(:,1), X(:,2), 0*X(:,1) + 1e10, 'color', 'black', 'linewidth', 1.5);

    h2 = subplot(3,1,2);
    Sur2 = multiblock.Surface(g, v0_2);
    xlabel('x')
    ylabel('y')
    title('v_y')
    shading interp
    colorbar
    xlims2 = xlim;
    ylims2 = ylim;
    line(X(:,1), X(:,2), 0*X(:,1) + 1e10, 'color', 'black', 'linewidth', 1.5);

    a = gca;

    h3 = subplot(3,1,3);
    L3 = plot(X(:,1), U_tangent);
    xlabel('x (m)')
    ylabel('Cumulative slip on fault (m)')

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
        L3.YData = U_tangent;
        drawnow;
    end
    update = @(t,U) update_fun(t,U,discr.Ecomp);
end