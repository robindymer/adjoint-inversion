% Figures and LaTeX table for the high-resolution data inversions (notes/high_res_plan.md, steps 2-4).
% Compares lbfgs and newton-cg-ls by wall time to reach the minimizer of the objective that each run
% minimizes (J = F + R, with the smoothed R for TV). J* is the smallest J of the two runs with the
% same regularization, and the time to a level L is the first time with (J - J*)/(J0 - J*) <= L.
% The error in a is the relative error on the resolved nodes 6-25 (nodes 1-5 and 26 are not
% resolved by the receivers). Output in assets/high_res_*.pdf and assets/high_res_table.tex.
addSubpaths;
resultsDir = fullfile('mat', 'inversion', 'inversion'); % Results downloaded from the cluster
assetsDir = 'assets';
resolved = 6:25;
level = 1e-4; % a is within about 1% of its final value at this level
solvers = {'lbfgs', 'newton-cg-ls'};
solverNames = {'L-BFGS', 'Newton-CG'};
colors = get(groot, 'defaultAxesColorOrder');

% True a on the m_p = 26 parameter grid, from the inverse-crime run
S = load(fullfile('mat', 'inversion', 'rough_newton-cg-ls_maxiter200_m76_mp26_velocity_synth_vals_a00135.mat'), 'truePars');
aTrue = S.truePars.a(:);
xp = linspace(-10, 10, numel(aTrue)); % Parameter grid
xFine = linspace(-10, 10, 1001); % The true a as a step, on the grid of the data
aFine = 0.009*((xFine >= -5) & (xFine <= 6)) + 0.013*~((xFine >= -5) & (xFine <= 6));

run = @(solver, m, reg) loadRun(resultsDir, solver, m, reg, aTrue, resolved);

%% Figure 1: no regularization, convergence and final a
r = {run('lbfgs', 76, ''), run('newton-cg-ls', 76, '')};
fig = figure('Position', [100 100 1400 500], 'Visible', 'off');
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile;
for k = 1:2
    semilogy(r{k}.t, r{k}.F, '-', 'color', colors(k,:), 'linewidth', 2); hold on;
end
hold off;
xlabel('wall time (h)', 'interpreter', 'latex'); ylabel('misfit $F$', 'interpreter', 'latex');
legend(solverNames, 'interpreter', 'latex'); grid on; set(gca, 'fontsize', 16);
nexttile;
plot(xFine, aFine, 'k-', 'linewidth', 2); hold on;
for k = 1:2
    plot(xp, r{k}.a, '--x', 'color', colors(k,:), 'linewidth', 2);
end
hold off;
xlabel('$x$ (km)', 'interpreter', 'latex'); ylabel('$\mathbf{a}$', 'interpreter', 'latex');
yl = ylim; ylim([yl(1), yl(2) + 0.35*diff(yl)]); % Headroom for the legend
legend(['true', compose('%s (error %.2f)', string(solverNames'), [r{1}.errRes; r{2}.errRes])'], ...
       'interpreter', 'latex', 'Location', 'northwest', 'NumColumns', 3);
grid on; set(gca, 'fontsize', 16);
exportgraphics(fig, fullfile(assetsDir, 'high_res_noreg.pdf'));

%% Figure 2: convergence with regularization, relative objective against wall time
% The regularized cases where Newton-CG reaches the level first, used in Figs. 2 and 4 and the
% compact table
ncgCases = {76, 'tikhonov-grad_beta1e+03', 'tikhonov-grad, $\beta = 10^3$, $m = 76$';
            76, 'tikhonov-grad_beta1e+04', 'tikhonov-grad, $\beta = 10^4$, $m = 76$';
            76, 'tv_beta1e+02_eps1e-03', 'TV, $\beta = 100$, $\varepsilon = 10^{-3}$, $m = 76$';
            251, 'tv_beta1e+00', 'TV, $\beta = 1$, $\varepsilon = 10^{-4}$, $m = 251$'};
cases = ncgCases;
fig = figure('Position', [100 100 1400 900], 'Visible', 'off');
tiledlayout(2, 2, 'TileSpacing', 'compact');
for c = 1:size(cases, 1)
    r = {run('lbfgs', cases{c,1}, cases{c,2}), run('newton-cg-ls', cases{c,1}, cases{c,2})};
    Jstar = min(min(r{1}.J), min(r{2}.J));
    nexttile;
    for k = 1:2
        semilogy(r{k}.t, relObjective(r{k}, Jstar), '-', 'color', colors(k,:), 'linewidth', 2); hold on;
    end
    yline(level, 'k:', 'linewidth', 1.5, 'HandleVisibility', 'off');
    hold off;
    ylim([1e-5, 1]); yticks(10.^(-5:0));
    title(cases{c,3}, 'interpreter', 'latex');
    xlabel('wall time (h)', 'interpreter', 'latex');
    ylabel('$(J - J^*)/(J_0 - J^*)$', 'interpreter', 'latex');
    legend(solverNames, 'interpreter', 'latex'); grid on; set(gca, 'fontsize', 16);
end
exportgraphics(fig, fullfile(assetsDir, 'high_res_reg_convergence.pdf'));

%% Figure 3: error in a against beta
regs = {'tikhonov-grad', [1e3 1e4 1e5 1e6], 'tikhonov-grad';
        'tv', [1 10 100 1000], 'TV, $\varepsilon = 10^{-4}$'};
fig = figure('Position', [100 100 1400 500], 'Visible', 'off');
tiledlayout(1, 2, 'TileSpacing', 'compact');
for g = 1:size(regs, 1)
    nexttile;
    for k = 1:2
        err = nan(size(regs{g,2}));
        for j = 1:numel(regs{g,2})
            rr = run(solvers{k}, 76, sprintf('%s_beta%.0e', regs{g,1}, regs{g,2}(j)));
            if ~isempty(rr) && ~rr.checkpoint
                err(j) = rr.errRes;
            end
        end
        semilogx(regs{g,2}, err, '-o', 'color', colors(k,:), 'linewidth', 2); hold on;
    end
    r0 = {run('lbfgs', 76, ''), run('newton-cg-ls', 76, '')};
    for k = 1:2 % Unregularized error for reference
        yline(r0{k}.errRes, ':', 'color', colors(k,:), 'linewidth', 1.5, 'HandleVisibility', 'off');
    end
    hold off;
    title(regs{g,3}, 'interpreter', 'latex');
    xlabel('$\beta$', 'interpreter', 'latex'); ylabel('error in $\mathbf{a}$, nodes 6--25', 'interpreter', 'latex');
    legend(solverNames, 'interpreter', 'latex'); grid on; set(gca, 'fontsize', 16);
end
exportgraphics(fig, fullfile(assetsDir, 'high_res_beta_error.pdf'));

%% Figure 4: final a for the cases of Figure 2
fig = figure('Position', [100 100 1400 900], 'Visible', 'off');
tiledlayout(2, 2, 'TileSpacing', 'compact');
for c = 1:size(ncgCases, 1)
    r = {run('lbfgs', ncgCases{c,1}, ncgCases{c,2}), run('newton-cg-ls', ncgCases{c,1}, ncgCases{c,2})};
    nexttile;
    plot(xFine, aFine, 'k-', 'linewidth', 2); hold on;
    for k = 1:2
        plot(xp, r{k}.a, '--x', 'color', colors(k,:), 'linewidth', 2);
    end
    hold off;
    title(ncgCases{c,3}, 'interpreter', 'latex');
    xlabel('$x$ (km)', 'interpreter', 'latex'); ylabel('$\mathbf{a}$', 'interpreter', 'latex');
    yl = ylim; ylim([yl(1), yl(2) + 0.35*diff(yl)]); % Headroom for the legend
    legend(['true', compose('%s (error %.3f)', string(solverNames'), [r{1}.errRes; r{2}.errRes])'], ...
           'interpreter', 'latex', 'Location', 'northwest', 'NumColumns', 2);
    grid on; set(gca, 'fontsize', 16);
end
exportgraphics(fig, fullfile(assetsDir, 'high_res_ncg_a.pdf'));

%% Compact table for the cases of Figure 2: error, iterations, time to the level and speedup
fid = fopen(fullfile(assetsDir, 'high_res_table_ncg.tex'), 'w');
fprintf(fid, '%% Generated by high_res_report_figures.m\n');
fprintf(fid, '\\begin{tabular}{llrrrr}\n\\toprule\n');
fprintf(fid, 'Case & & error in $\\mathbf{a}$ & iter. & $t_{10^{%d}}$ (h) & speedup \\\\\n\\midrule\n', round(log10(level)));
for c = 1:size(ncgCases, 1)
    r = {run('lbfgs', ncgCases{c,1}, ncgCases{c,2}), run('newton-cg-ls', ncgCases{c,1}, ncgCases{c,2})};
    Jstar = min(min(r{1}.J), min(r{2}.J));
    tL = nan(1, 2);
    for k = 1:2
        j = find(relObjective(r{k}, Jstar) <= level, 1);
        if ~isempty(j)
            tL(k) = r{k}.t(j);
        end
    end
    if isnan(tL(1)) % L-BFGS does not reach the level, so the speedup is a lower bound
        tStr = {sprintf('$>$%.2f', r{1}.t(end)), sprintf('%.2f', tL(2))};
        speedup = sprintf('$>$%.1f', r{1}.t(end)/tL(2));
    else
        tStr = {sprintf('%.2f', tL(1)), sprintf('%.2f', tL(2))};
        speedup = sprintf('%.1f', tL(1)/tL(2));
    end
    caseStr = {ncgCases{c,3}, ''};
    speedStr = {'', speedup};
    for k = 1:2
        unfinished = '';
        if r{k}.checkpoint
            unfinished = '$^\dagger$';
        end
        fprintf(fid, '%s & %s & %.3f & %d%s & %s & %s \\\\\n', caseStr{k}, solverNames{k}, r{k}.errRes, ...
                r{k}.iter, unfinished, tStr{k}, speedStr{k});
    end
    if c < size(ncgCases, 1)
        fprintf(fid, '\\addlinespace\n');
    end
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);
type(fullfile(assetsDir, 'high_res_table_ncg.tex'));

%% Final a for the regularized cases where Newton-CG reaches the level first, one figure each, with
% the LaTeX figure environments (captions with the numbers) in assets/high_res_a_figures.tex
cases = {76, 'tikhonov-grad_beta1e+03', 'tikhonov-grad with $\beta = 10^3$';
         76, 'tikhonov-grad_beta1e+04', 'tikhonov-grad with $\beta = 10^4$';
         76, 'tv_beta1e+02_eps1e-03', 'TV with $\beta = 100$ and $\varepsilon = 10^{-3}$';
         251, 'tv_beta1e+00', 'TV with $\beta = 1$ and $\varepsilon = 10^{-4}$'};
fid = fopen(fullfile(assetsDir, 'high_res_a_figures.tex'), 'w');
fprintf(fid, '%% Generated by high_res_report_figures.m\n');
for c = 1:size(cases, 1)
    r = {run('lbfgs', cases{c,1}, cases{c,2}), run('newton-cg-ls', cases{c,1}, cases{c,2})};
    Jstar = min(min(r{1}.J), min(r{2}.J));
    tL = nan(1, 2);
    for k = 1:2
        j = find(relObjective(r{k}, Jstar) <= level, 1);
        if ~isempty(j)
            tL(k) = r{k}.t(j);
        end
    end
    if isnan(tL(1))
        timing = sprintf(['after %.2f~h, whereas L-BFGS does not reach it within its %d iterations ' ...
                          '(%.2f~h)'], tL(2), r{1}.iter, r{1}.t(end));
    else
        timing = sprintf('after %.2f~h against %.2f~h for L-BFGS (%.1f times faster)', tL(2), tL(1), tL(1)/tL(2));
    end
    fig = figure('Position', [100 100 900 500], 'Visible', 'off');
    plot(xFine, aFine, 'k-', 'linewidth', 2); hold on;
    for k = 1:2
        plot(xp, r{k}.a, '--x', 'color', colors(k,:), 'linewidth', 2);
    end
    hold off;
    xlabel('$x$ (km)', 'interpreter', 'latex'); ylabel('$\mathbf{a}$', 'interpreter', 'latex');
    yl = ylim; ylim([yl(1), yl(2) + 0.35*diff(yl)]); % Headroom for the legend
    legend(['true', solverNames], 'interpreter', 'latex', 'Location', 'northwest', 'NumColumns', 3);
    grid on; set(gca, 'fontsize', 16);
    name = sprintf('high_res_a_m%d_%s', cases{c,1}, cases{c,2});
    exportgraphics(fig, fullfile(assetsDir, [name, '.pdf']));

    unfinished = {'', ''};
    for k = 1:2
        if r{k}.checkpoint
            unfinished{k} = ', unfinished';
        end
    end
    fprintf(fid, ['\\begin{figure}\n  \\centering\n  \\includegraphics[width=0.8\\textwidth]{assets/compare/%s.pdf}\n', ...
                  '  \\caption{Final $\\mathbf{a}$, high-resolution data inverted on $m = %d$, $m_p = 26$, %s. ', ...
                  'Newton-CG reaches $(J - J^*)/(J_0 - J^*) = 10^{%d}$ %s. ', ...
                  'Relative error on nodes 6--25: %.3f (Newton-CG, %d iterations, %.2f~h%s) ', ...
                  'and %.3f (L-BFGS, %d iterations, %.2f~h%s).}\n  \\label{fig:highres-a-%d-%s}\n\\end{figure}\n\n'], ...
            name, cases{c,1}, cases{c,3}, round(log10(level)), timing, r{2}.errRes, r{2}.iter, r{2}.t(end), ...
            unfinished{2}, r{1}.errRes, r{1}.iter, r{1}.t(end), unfinished{1}, cases{c,1}, cases{c,2});
end
fclose(fid);

%% Figure 5: TV with beta = 100, effect of eps on Newton-CG. All runs on the exact objective F + beta*TV.
epsCases = {'tv_beta1e+02_eps1e-03', '$\varepsilon = 10^{-3}$';
            'tv_beta1e+02', '$\varepsilon = 10^{-4}$';
            'tv_beta1e+02_eps1e-05', '$\varepsilon = 10^{-5}$';
            'tv_beta1e+02_eps0e+00', '$\varepsilon = 0$'};
fig = figure('Position', [100 100 1400 500], 'Visible', 'off');
tiledlayout(1, 2, 'TileSpacing', 'compact');
for k = 1:2
    nexttile;
    names = {};
    for c = 1:size(epsCases, 1)
        rr = run(solvers{k}, 76, epsCases{c,1});
        if isempty(rr)
            continue
        end
        semilogy(rr.t, rr.Jexact, '-', 'color', colors(c,:), 'linewidth', 2); hold on;
        names{end+1} = epsCases{c,2}; %#ok<SAGROW>
        if rr.checkpoint
            names{end} = [names{end}, ' (unfinished)'];
        end
    end
    hold off;
    title(sprintf('%s, TV, $\\beta = 100$', solverNames{k}), 'interpreter', 'latex');
    xlabel('wall time (h)', 'interpreter', 'latex'); ylabel('$F + \beta\,\mathrm{TV}$ (exact)', 'interpreter', 'latex');
    legend(names, 'interpreter', 'latex'); grid on; set(gca, 'fontsize', 16);
end
exportgraphics(fig, fullfile(assetsDir, 'high_res_tv_eps.pdf'));

%% Table: all m = 76 runs, and the m = 251 runs used in the text
rows = {76, '', 'none', '--';
        76, 'tikhonov-grad_beta1e+03', 'tikhonov-grad', '$10^3$';
        76, 'tikhonov-grad_beta1e+04', 'tikhonov-grad', '$10^4$';
        76, 'tikhonov-grad_beta1e+05', 'tikhonov-grad', '$10^5$';
        76, 'tikhonov-grad_beta1e+06', 'tikhonov-grad', '$10^6$';
        76, 'tv_beta1e+00', 'TV', '$1$';
        76, 'tv_beta1e+01', 'TV', '$10$';
        76, 'tv_beta1e+02', 'TV', '$10^2$';
        76, 'tv_beta1e+02_eps1e-03', 'TV, $\varepsilon = 10^{-3}$', '$10^2$';
        76, 'tv_beta1e+03', 'TV', '$10^3$';
        251, '', 'none', '--';
        251, 'tv_beta1e+00', 'TV', '$1$'};
fid = fopen(fullfile(assetsDir, 'high_res_table.tex'), 'w');
fprintf(fid, '%% Generated by high_res_report_figures.m\n');
fprintf(fid, '\\begin{tabular}{rllrrrrrr}\n\\toprule\n');
fprintf(fid, '$m$ & $R$ & $\\beta$ & & $J$ & error in $\\mathbf{a}$ & iter. & time (h) & $t_{10^{%d}}$ (h) \\\\\n\\midrule\n', round(log10(level)));
for i = 1:size(rows, 1)
    r = {run('lbfgs', rows{i,1}, rows{i,2}), run('newton-cg-ls', rows{i,1}, rows{i,2})};
    Js = cellfun(@(x) min(x.J), r(~cellfun(@isempty, r)));
    Jstar = min(Js);
    tL = nan(1, 2);
    for k = 1:2
        if isempty(r{k})
            continue
        end
        relJ = relObjective(r{k}, Jstar);
        j = find(relJ <= level, 1);
        if ~isempty(j)
            tL(k) = r{k}.t(j);
        end
    end
    % The better solver of the pair in bold, for J, the error and the time to the level
    Jend = nan(1, 2);
    errRes = nan(1, 2);
    for k = 1:2
        if ~isempty(r{k})
            Jend(k) = r{k}.J(end);
            errRes(k) = r{k}.errRes;
        end
    end
    Jstr = boldSmaller(Jend, '%.3f');
    errStr = boldSmaller(errRes, '%.3f');
    tStr = boldSmaller(tL, '%.2f');
    for k = 1:2
        if k == 1
            lead = sprintf('%d & %s & %s', rows{i,1}, rows{i,3}, rows{i,4});
        else
            lead = ' & & ';
        end
        if isempty(r{k})
            fprintf(fid, '%s & %s & \\multicolumn{5}{c}{not run} \\\\\n', lead, solverNames{k});
            continue
        end
        unfinished = '';
        if r{k}.checkpoint
            unfinished = '$^\dagger$';
        end
        fprintf(fid, '%s & %s & %s & %s & %d%s & %.2f & %s \\\\\n', lead, solverNames{k}, Jstr{k}, ...
                errStr{k}, r{k}.iter, unfinished, r{k}.t(end), tStr{k});
    end
    if i < size(rows, 1)
        fprintf(fid, '\\addlinespace\n');
    end
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);
type(fullfile(assetsDir, 'high_res_table.tex'));

%% Figure 7: final a at m_p = 51 (step 5a)
% The error is the relative L2 error of the linearly interpolated a against the step on the
% resolved region x in [-6, 9.2], so that it can be compared with m_p = 26
win = (xFine >= -6) & (xFine <= 9.2);
cases51 = {'tikhonov-grad_beta1e+04', 'tikhonov-grad, $\beta = 10^4$';
           'tv_beta1e+01_eps1e-03', 'TV, $\beta = 10$, $\varepsilon = 10^{-3}$'};
fig = figure('Position', [100 100 1400 500], 'Visible', 'off');
tiledlayout(1, 2, 'TileSpacing', 'compact');
for c = 1:size(cases51, 1)
    nexttile;
    plot(xFine, aFine, 'k-', 'linewidth', 2); hold on;
    legendStr = {'true'};
    for k = 1:2
        S = load(fullfile(resultsDir, sprintf('rough_%s_maxiter200_m76_mp51_velocity_highres_vals_a00135_%s.mat', ...
                                              solvers{k}, cases51{c,1})), 'history');
        a = S.history.x(:,end);
        xp51 = linspace(-10, 10, numel(a));
        aInterp = interp1(xp51, a, xFine);
        err = norm(aInterp(win) - aFine(win))/norm(aFine(win));
        plot(xp51, a, '--x', 'color', colors(k,:), 'linewidth', 2);
        legendStr{end+1} = sprintf('%s (error %.3f, %d iter.)', solverNames{k}, err, numel(S.history.fval) - 1);
    end
    hold off;
    title(cases51{c,2}, 'interpreter', 'latex');
    xlabel('$x$ (km)', 'interpreter', 'latex'); ylabel('$\mathbf{a}$', 'interpreter', 'latex');
    yl = ylim; ylim([yl(1), yl(2) + 0.35*diff(yl)]); % Headroom for the legend
    legend(legendStr, 'interpreter', 'latex', 'Location', 'northwest');
    grid on; set(gca, 'fontsize', 16);
end
exportgraphics(fig, fullfile(assetsDir, 'high_res_mp51_a.pdf'));

% Loads one run; empty if the file does not exist. Unfinished runs are read from their checkpoint.
function r = loadRun(resultsDir, solver, m, reg, aTrue, resolved)
    base = sprintf('rough_%s_maxiter200_m%d_mp26_velocity_highres_vals_a00135', solver, m);
    if ~isempty(reg)
        base = [base, '_', reg];
    end
    r = [];
    path = fullfile(resultsDir, [base, '.mat']);
    r.checkpoint = ~isfile(path);
    if r.checkpoint
        path = fullfile(resultsDir, [base, '_checkpoint.mat']);
        if ~isfile(path)
            r = [];
            return
        end
    end
    S = load(path, 'history');
    h = S.history;
    if ~isfield(h, 'reg')
        h.reg = zeros(size(h.fval));
    end
    r.t = h.t/3600;
    r.J = h.fval;
    r.F = h.fval - h.reg;
    % Exact objective with eps = 0 for TV
    tok = regexp(reg, '^(tv|tikhonov-grad)_beta([^_]+)', 'tokens', 'once');
    if ~isempty(tok) && strcmp(tok{1}, 'tv')
        regExact = adjopt.regularization(struct('type', 'tv', 'beta', str2double(tok{2}), 'eps', 0), ...
                                         h.x(:,1), size(h.x, 1));
        r.Jexact = r.F + arrayfun(@(k) regExact.R(h.x(:,k)), (1:size(h.x, 2))');
    else
        r.Jexact = r.J;
    end
    r.a = h.x(:,end);
    r.iter = numel(h.fval) - 1;
    r.errRes = norm(r.a(resolved) - aTrue(resolved))/norm(aTrue(resolved));
end

function relJ = relObjective(r, Jstar)
    relJ = (r.J - Jstar)/(r.J(1) - Jstar);
end

% Formats the two values, with the smaller one in bold. No bold if they are equal as printed.
function s = boldSmaller(vals, fmt)
    s = cell(1, 2);
    for k = 1:2
        if isnan(vals(k))
            s{k} = '--';
        else
            s{k} = sprintf(fmt, vals(k));
        end
    end
    if strcmp(s{1}, s{2}) || all(isnan(vals))
        return
    end
    [~, k] = min(vals); % min ignores NaN, so a value beats a missing one
    s{k} = ['\textbf{', s{k}, '}'];
end
