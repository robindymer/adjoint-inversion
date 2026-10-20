% Plots inversion results in the same way as direct_effect_inverse_crime.m, one figure per file.
% Handles both files saved by optimize_antiplaneshear (mat/inversion, variable history) and
% files saved by optimizer_benchmark/run_method.m (optimizer_benchmark/results, variable hist).
addSubpaths;
% resultsDir = fullfile('optimizer_benchmark', 'results');
resultsDir = fullfile('mat', 'inversion');
resultFiles = {'rough_CG_maxiter200_m76_mp26_velocity_synth_vals_a00135.mat';
               'rough_lbfgs_maxiter200_m76_mp26_velocity_synth_vals_a00135.mat';
               'rough_newton-tr_maxiter200_m76_mp26_velocity_synth_vals_a00135.mat';
               'rough_newton-tr-noprec_maxiter200_m76_mp26_velocity_synth_vals_a00135.mat';
               'rough_fminunc-tr_maxiter200_m76_mp26_velocity_synth_vals_a00135.mat';
               'rough_newton-cg-ls_maxiter200_m76_mp26_velocity_synth_vals_a00135.mat'};
% resultFiles = {'res_newton_tr_repo_m251.mat'; 'res_fmincon_ip_lbfgs_nolb_m251.mat'};
iter = []; % Iteration to plot, [] for the final iteration of each file
assetsDir = fullfile('..', 'assets');

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

    fig = figure();
    theme(fig, 'light');
    x = linspace(-10,10,m_p); % Parameter grid
    a_inverted = history.x(:,it+1); % Column 1 is iteration 0
    plot(x,truePars.a, x,initialPars.a,':x', x,a_inverted,'k--x', x,truePars.b,'k:','linewidth',2);
    xlabel('$x$ (km)','interpreter','latex');
    legend('$\mathbf{a}$, true','$\mathbf{a}$, initial','$\mathbf{a}$, current','$\mathbf{b}$','interpreter','latex','Location','southwest');
    title(sprintf('%s, iter %d', strrep(erase(resultFile, '.mat'), '_', '\_'), it));
    grid on;
    set(gca,'fontsize',16);

    exportgraphics(fig, fullfile(assetsDir, sprintf('%s_iter%d.png', erase(resultFile, '.mat'), it)), 'Resolution', 300);
end
