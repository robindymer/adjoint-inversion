outDir = fileparts(mfilename('fullpath'));
[x, y] = domains.fractal_profile([-10 10], 1001);
[ymin, i] = min(y);
fprintf('min y = %.3f at x = %.3f\n', ymin, x(i));

f = figure('Visible','off','Position',[100 100 500 480]); hold on; box on
plot(x, y, 'k', 'LineWidth', 1.5)
axis equal; xlim([-10 10]); ylim([-10 10])
set(gca,'XTick',[-10 0 10],'YTick',[-10 0 10],'LineWidth',1.2,'FontSize',12)
text(0, 5, '$\Omega^+$', 'Interpreter','latex', 'FontSize', 16, 'HorizontalAlignment','center')
text(0, -5, '$\Omega^-$', 'Interpreter','latex', 'FontSize', 16, 'HorizontalAlignment','center')
text(0, 1.2, '$\Gamma$', 'Interpreter','latex', 'FontSize', 16, 'HorizontalAlignment','center')
xlabel('$x$ (km)', 'Interpreter','latex'); ylabel('$y$ (km)', 'Interpreter','latex')
exportgraphics(f, fullfile(outDir, 'fractal_fault_domain.png'), 'Resolution', 200);
exportgraphics(f, fullfile(outDir, 'fractal_fault_domain.pdf'), 'ContentType', 'vector');
