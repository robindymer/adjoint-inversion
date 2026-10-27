% Plots inversion results in the same way as direct_effect_inverse_crime.m, one figure per file.
% Handles both files saved by optimize_antiplaneshear (mat/inversion, variable history) and
% files saved by optimizer_benchmark/run_method.m (optimizer_benchmark/results, variable hist).
addSubpaths;
% resultsDir = fullfile('optimizer_benchmark', 'results');
resultsDir = fullfile('mat', 'inversion/inversion');
resultFiles = {'rough_lbfgs_maxiter200_m251_mp26_velocity_highres_vals_a00135.mat';
               'rough_newton-cg-ls_maxiter200_m251_mp26_velocity_highres_vals_a00135.mat'};
% resultFiles = {'res_newton_tr_repo_m251.mat'; 'res_fmincon_ip_lbfgs_nolb_m251.mat'};
iter = []; % Iteration to plot, [] for the final iteration of each file
assetsDir = 'assets';
if ~isfolder(assetsDir)
    mkdir(assetsDir);
end

for i = 1:numel(resultFiles)
    resultFile = resultFiles{i};
    filePath = fullfile(resultsDir, resultFile);
    if ~isfile(filePath)
        warning('Skipping %s: file not found', filePath);
        continue
    end

    vars = who('-file', filePath);
    if ismember('history', vars) % Saved by optimize_antiplaneshear
        S = load(filePath, 'history', 'truePars', 'initialPars');
        history = S.history;
        truePars = S.truePars;
        initialPars = S.initialPars;
    else % Saved by optimizer_benchmark/run_method.m
        S = load(filePath, 'hist');
        history = S.hist;
        m_p = size(history.x, 1);
        % These files do not store the true parameters. The true a on the parameter grid
        % does not depend on m, so it is taken from H0.mat (saved by time_evals.m at m = 36).
        H0 = load(fullfile('optimizer_benchmark', 'results', 'H0.mat'), 'xtrue');
        truePars.a = H0.xtrue;
        truePars.b = 0.011*ones(m_p, 1); % Constant b, see pars.rsFriction2DFractalFaultInversion
        initialPars.a = history.x(:,1);
    end
    m_p = size(history.x, 1); % Grid points on parameter grid

    if isempty(iter)
        it = size(history.x, 2) - 1; % Final iteration
    else
        it = iter;
    end

    x = linspace(-10,10,m_p); % Parameter grid
    if isempty(truePars) % Inversion of loaded highres data, where no true parameters are saved
        % True a used to generate the data, see pars.rsFriction2DFractalFaultForward. Evaluated on
        % the m_p = 1001 grid of the data, so that it is shown as a step as in the article.
        xTrue = linspace(-10,10,1001);
        if_velocity_weakening = (xTrue >= -5) & (xTrue <= 6);
        aTrue = 0.009*if_velocity_weakening + 0.013*(1-if_velocity_weakening);
        b = initialPars.b; % b is not inverted for
    else
        xTrue = x;
        aTrue = truePars.a;
        b = truePars.b;
    end

    % Hidden while drawing and exporting, since a tiling window manager otherwise resizes the figure
    fig = figure('Position', [100 100 900 500], 'Visible', 'off');
    % theme(fig, 'light');
    a_inverted = history.x(:,it+1); % Column 1 is iteration 0
    colors = get(groot, 'defaultAxesColorOrder');
    plot(xTrue,aTrue,'color',colors(1,:),'linewidth',2); hold on;
    plot(x,initialPars.a,':x','color',colors(2,:),'linewidth',2);
    plot(x,a_inverted,'k--x','linewidth',2);
    plot(x,b,'k:','linewidth',2);
    % y-limits from the reference curves, so that diverged values do not hide the rest of the plot.
    % Values outside the limits are marked at the edge of the axes, labeled with their actual value.
    yRef = [aTrue(:); initialPars.a(:); b(:)];
    yRange = max(yRef) - min(yRef);
    yAll = [yRef; a_inverted(:)];
    yl = [max(min(yAll), min(yRef) - yRange), min(max(yAll), max(yRef) + yRange)];
    above = a_inverted > yl(2);
    below = a_inverted < yl(1);
    plot(x(above), yl(2)*ones(nnz(above),1), 'r^', 'markersize', 10, 'markerfacecolor', 'r', 'HandleVisibility', 'off');
    plot(x(below), yl(1)*ones(nnz(below),1), 'rv', 'markersize', 10, 'markerfacecolor', 'r', 'HandleVisibility', 'off');
    text(x(above), yl(2)*ones(nnz(above),1), compose('  %.3g', a_inverted(above)), 'color', 'r', 'fontsize', 12, ...
        'VerticalAlignment', 'top');
    text(x(below), yl(1)*ones(nnz(below),1), compose('  %.3g', a_inverted(below)), 'color', 'r', 'fontsize', 12, ...
        'VerticalAlignment', 'bottom');
    hold off;
    xlabel('$x$ (km)','interpreter','latex');
    legendEntries = {'$\mathbf{a}$, true','$\mathbf{a}$, initial',sprintf('$\\mathbf{a}$, %d iterations',it),'$\mathbf{b}$'};
    if any(above | below) % Legend above the axes, so that it does not hide the markers at the edges
        ylim(yl);
        legend(legendEntries{:},'interpreter','latex','Location','northoutside','Orientation','horizontal');
    else
        ylim([yl(1), yl(2) + 0.35*diff(yl)]); % Headroom for the legend
        legend(legendEntries{:},'interpreter','latex','Location','northwest');
    end
    grid on;
    set(gca,'fontsize',16);

    exportgraphics(fig, fullfile(assetsDir, sprintf('%s_iter%d.png', erase(resultFile, '.mat'), it)), 'Resolution', 300);
    fig.Visible = 'on';
end
