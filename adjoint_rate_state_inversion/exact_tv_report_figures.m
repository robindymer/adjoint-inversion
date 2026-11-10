% Figures and table on TV with small or zero eps, all for m = 76 and m_p = 26:
% 1. high_res_tv_beta_eps.pdf: final a for the TV runs on the m = 1001 data (cluster runs, read from
%    mat/inversion/inversion/), for several beta and eps.
% 2. high_res_exact_tv.pdf: exact TV (eps = 0) on the m = 251 data, made locally with
%    scripts/high-res_runs_tv_exact.txt and scripts/high-res_runs_tv_exact_ncg.txt.
% 3. high_res_table_tv.tex: errors, final J and t_{10^-4} for both.
% Stopped runs are read from their checkpoints.
addSubpaths;
resultsDir = fullfile('mat', 'inversion');
clusterDir = fullfile('mat', 'inversion', 'inversion');
assetsDir = 'assets';
base = 'rough_%s_maxiter200_m76_mp26_velocity_highres_m251_mp251_vals_a00135_%s';
baseCluster = 'rough_%s_maxiter200_m76_mp26_velocity_highres_vals_a00135_%s';
level = 1e-4;
colors = get(groot, 'defaultAxesColorOrder');

S = load(fullfile(resultsDir, 'rough_newton-cg-ls_maxiter200_m76_mp26_velocity_synth_vals_a00135.mat'), 'truePars');
aTrue = S.truePars.a(:);
resolved = 6:25;
xp = linspace(-10, 10, numel(aTrue));
xFine = linspace(-10, 10, 1001);
aFine = 0.009*((xFine >= -5) & (xFine <= 6)) + 0.013*~((xFine >= -5) & (xFine <= 6));

load_run = @(solver, reg) loadHistory(fullfile(resultsDir, sprintf(base, solver, reg)), aTrue, resolved);
load_cluster = @(solver, reg) loadHistory(fullfile(clusterDir, sprintf(baseCluster, solver, reg)), aTrue, resolved);
solverNames2 = {'L-BFGS', 'Newton-CG'};

%% Figure 1: final a for TV on the m = 1001 data, over beta and eps
tvCases = {'tv_beta1e+00', 1, '10^{-4}'; 'tv_beta1e+01', 10, '10^{-4}'; 'tv_beta1e+02', 100, '10^{-4}';
           'tv_beta1e+03', 1000, '10^{-4}'; 'tv_beta1e+02_eps1e-03', 100, '10^{-3}';
           'tv_beta1e+02_eps1e-05', 100, '10^{-5}'; 'tv_beta1e+02_eps0e+00', 100, '0'};
fig = figure('Position', [100 100 1900 850], 'Visible', 'off');
tiledlayout(2, 4, 'TileSpacing', 'compact');
for c = 1:size(tvCases, 1)
    r = {load_cluster('lbfgs', tvCases{c,1}), load_cluster('newton-cg-ls', tvCases{c,1})};
    nexttile;
    plot(xFine, aFine, 'k-', 'linewidth', 2); hold on;
    names = {'true'};
    for k = 1:2
        plot(xp, r{k}.a, '--x', 'color', colors(k,:), 'linewidth', 2);
        names{end+1} = sprintf('%s (%.3f%s)', solverNames2{k}, r{k}.err, repmat(', stopped', 1, r{k}.stopped)); %#ok<SAGROW>
    end
    hold off;
    title(sprintf('$\\beta = %g$, $\\varepsilon = %s$', tvCases{c,2}, tvCases{c,3}), 'interpreter', 'latex');
    xlabel('$x$ (km)', 'interpreter', 'latex'); ylabel('$\mathbf{a}$', 'interpreter', 'latex');
    ylim([0.0075 0.0175]);
    legend(names, 'interpreter', 'latex', 'Location', 'north', 'fontsize', 12);
    grid on; set(gca, 'fontsize', 15);
end
exportgraphics(fig, fullfile(assetsDir, 'high_res_tv_beta_eps.pdf'));

%% Figure 2: exact TV on the m = 251 data

fig = figure('Position', [100 100 1900 500], 'Visible', 'off');
tiledlayout(1, 3, 'TileSpacing', 'compact');

% (a) L-BFGS, exact TV for three beta, and tikhonov-grad for reference
cases = {'tv_beta3e+00_eps0e+00', 'TV, $\beta = 3$';
         'tv_beta1e+01_eps0e+00', 'TV, $\beta = 10$';
         'tv_beta3e+01_eps0e+00', 'TV, $\beta = 30$';
         'tikhonov-grad_beta1e+04', 'tikhonov-grad, $\beta = 10^4$'};
nexttile;
plot(xFine, aFine, 'k-', 'linewidth', 2); hold on;
names = {'true'};
for c = 1:size(cases, 1)
    r = load_run('lbfgs', cases{c,1});
    plot(xp, r.a, '--x', 'color', colors(c,:), 'linewidth', 2);
    names{end+1} = sprintf('%s (%.3f%s)', cases{c,2}, r.err, repmat(sprintf(', %d it.', r.iter), 1, r.stopped)); %#ok<SAGROW>
    fprintf('L-BFGS %s: J = %.4f, error %.3f, median |Da| = %.1e, %d iterations, %.2f h\n', ...
            cases{c,1}, r.J(end), r.err, median(abs(diff(r.a))), r.iter, r.t(end));
end
hold off;
title('L-BFGS, exact TV ($\varepsilon = 0$)', 'interpreter', 'latex');
xlabel('$x$ (km)', 'interpreter', 'latex'); ylabel('$\mathbf{a}$', 'interpreter', 'latex');
yl = ylim; ylim([yl(1), yl(2) + 0.45*diff(yl)]); % Headroom for the legend
legend(names, 'interpreter', 'latex', 'Location', 'northwest', 'fontsize', 13);
grid on; set(gca, 'fontsize', 16);

% (b), (c) exact TV with beta = 10: L-BFGS, Newton-CG and Newton-CG from L-BFGS iteration 20
solvers = {'lbfgs', 'newton-cg-ls', 'lbfgs20-newton-cg-ls'};
solverNames = {'L-BFGS', 'Newton-CG', 'L-BFGS (20 it.) $\to$ Newton-CG'};
r = cellfun(@(s) load_run(s, 'tv_beta1e+01_eps0e+00'), solvers);
nexttile;
for k = 1:3
    i0 = 1 + 20*(k == 3); % The hybrid repeats the L-BFGS history up to iteration 20
    semilogy(r(k).t(i0:end), r(k).J(i0:end), '-', 'color', colors(k,:), 'linewidth', 2); hold on;
    fprintf('%s: J = %.4f, error %.3f, %d iterations, %.2f h%s\n', solvers{k}, r(k).J(end), r(k).err, ...
            r(k).iter, r(k).t(end), repmat(' (stopped)', 1, r(k).stopped));
end
hold off;
j = find(r(1).J <= r(2).J(end), 1);
fprintf('L-BFGS reaches the final Newton-CG J = %.3f at iteration %d, %.2f h\n', r(2).J(end), j - 1, r(1).t(j));
title('Exact TV, $\beta = 10$', 'interpreter', 'latex');
xlabel('wall time (h)', 'interpreter', 'latex'); ylabel('objective $J$', 'interpreter', 'latex');
legend(solverNames, 'interpreter', 'latex', 'fontsize', 13); grid on; set(gca, 'fontsize', 16);

nexttile;
plot(xFine, aFine, 'k-', 'linewidth', 2); hold on;
names = {'true'};
for k = 1:3
    plot(xp, r(k).a, '--x', 'color', colors(k,:), 'linewidth', 2);
    names{end+1} = sprintf('%s (%.3f%s)', solverNames{k}, r(k).err, repmat(', stopped', 1, r(k).stopped)); %#ok<SAGROW>
end
hold off;
title('Exact TV, $\beta = 10$', 'interpreter', 'latex');
xlabel('$x$ (km)', 'interpreter', 'latex'); ylabel('$\mathbf{a}$', 'interpreter', 'latex');
yl = ylim; ylim([yl(1), yl(2) + 0.35*diff(yl)]);
legend(names, 'interpreter', 'latex', 'Location', 'northwest', 'fontsize', 13);
grid on; set(gca, 'fontsize', 16);
exportgraphics(fig, fullfile(assetsDir, 'high_res_exact_tv.pdf'));

%% Table: error, final J and t_{10^-4} for L-BFGS and Newton-CG
% J* is the lowest J of the two solvers. A '>' time is the length of a run that did not reach the
% level, and a dagger marks a stopped run. t is left out when only L-BFGS was run.
rows = {'1001', 'tikhonov-grad_beta1e+04', 'tikhonov-grad', '$10^4$', '--';
        '1001', 'tv_beta1e+00', 'TV', '$1$', '$10^{-4}$';
        '1001', 'tv_beta1e+01', 'TV', '$10$', '$10^{-4}$';
        '1001', 'tv_beta1e+02_eps1e-03', 'TV', '$10^2$', '$10^{-3}$';
        '1001', 'tv_beta1e+02', 'TV', '$10^2$', '$10^{-4}$';
        '1001', 'tv_beta1e+02_eps1e-05', 'TV', '$10^2$', '$10^{-5}$';
        '1001', 'tv_beta1e+02_eps0e+00', 'TV', '$10^2$', '$0$';
        '1001', 'tv_beta1e+03', 'TV', '$10^3$', '$10^{-4}$';
        '251', 'tikhonov-grad_beta1e+04', 'tikhonov-grad', '$10^4$', '--';
        '251', 'tv_beta3e+00_eps0e+00', 'TV', '$3$', '$0$';
        '251', 'tv_beta1e+01_eps0e+00', 'TV', '$10$', '$0$';
        '251', 'tv_beta3e+01_eps0e+00', 'TV', '$30$', '$0$'};
fid = fopen(fullfile(assetsDir, 'high_res_table_tv.tex'), 'w');
fprintf(fid, '%% Generated by exact_tv_report_figures.m\n');
fprintf(fid, '\\begin{tabular}{rllrrrrrr}\n\\toprule\n');
fprintf(fid, '& & & \\multicolumn{2}{c}{error in $\\mathbf{a}$} & \\multicolumn{2}{c}{final $J$} & \\multicolumn{2}{c}{$t_{10^{-4}}$ (h)} \\\\\n');
fprintf(fid, '\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}\n');
fprintf(fid, 'data $m$ & $R$, $\\beta$ & $\\varepsilon$ & L-BFGS & Newton-CG & L-BFGS & Newton-CG & L-BFGS & Newton-CG \\\\\n\\midrule\n');
for i = 1:size(rows, 1)
    if i > 1 && ~strcmp(rows{i,1}, rows{i-1,1})
        fprintf(fid, '\\midrule\n');
    end
    if strcmp(rows{i,1}, '1001'), loader = load_cluster; else, loader = load_run; end
    r = {loader('lbfgs', rows{i,2}), []};
    try
        r{2} = loader('newton-cg-ls', rows{i,2});
    catch
    end
    Jstar = min(cellfun(@(rr) min(rr.J), r(~cellfun(@isempty, r))));
    [e, J, t] = deal(cell(1, 2));
    for k = 1:2
        if isempty(r{k})
            [e{k}, J{k}, t{k}] = deal('--');
            continue
        end
        dag = repmat('$^\dagger$', 1, r{k}.stopped);
        e{k} = sprintf('%.3f', r{k}.err);
        J{k} = sprintf('%.3f%s', r{k}.J(end), dag);
        j = find((r{k}.J - Jstar)/(r{k}.J(1) - Jstar) <= level, 1);
        if isempty(r{2})
            t{k} = '--'; % J* needs both solvers
        elseif isempty(j)
            t{k} = sprintf('$>$%.2f', r{k}.t(end));
        else
            t{k} = sprintf('%.2f', r{k}.t(j));
        end
    end
    fprintf(fid, '%s & %s, %s & %s & %s & %s & %s & %s & %s & %s \\\\\n', rows{i,1}, rows{i,3}, rows{i,4}, ...
            rows{i,5}, e{1}, e{2}, J{1}, J{2}, t{1}, t{2});
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);

function r = loadHistory(path, aTrue, resolved)
    r.stopped = ~isfile([path, '.mat']);
    if r.stopped
        path = [path, '_checkpoint'];
    end
    h = load([path, '.mat'], 'history').history;
    r.t = h.t/3600;
    r.J = h.fval;
    r.a = h.x(:,end);
    r.iter = numel(h.fval) - 1;
    r.err = norm(r.a(resolved) - aTrue(resolved))/norm(aTrue(resolved));
end
