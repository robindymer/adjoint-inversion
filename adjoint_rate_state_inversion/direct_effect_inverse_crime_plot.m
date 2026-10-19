% Plots an inversion result from optimizer_benchmark/results in the same way as
% direct_effect_inverse_crime.m.
addSubpaths;
resultsDir = fullfile('optimizer_benchmark', 'results');
resultFile = 'res_newton_tr_repo_m251.mat';
% resultFile = 'res_fmincon_ip_lbfgs_nolb_m251.mat';

S = load(fullfile(resultsDir, resultFile)); % Contains hist and method
history = S.hist;
m_p = size(history.x, 1); % Grid points on parameter grid

% The result files do not store the true parameters. The true a on the parameter grid
% does not depend on m, so it is taken from H0.mat (saved by time_evals.m at m = 36).
H0 = load(fullfile(resultsDir, 'H0.mat'), 'xtrue');
truePars.a = H0.xtrue;
truePars.b = 0.011*ones(m_p, 1); % Constant b, see pars.rsFriction2DFractalFaultInversion
initialPars.a = history.x(:,1);

fig = figure();
theme(fig, 'light');
x = linspace(-10,10,m_p); % Parameter grid
iter = size(history.x, 2) - 1; % Iteration to plot, final by default
% iter = 45;
a_inverted = history.x(:,iter+1); % Column 1 is iteration 0, see InvProblem.log
plot(x,truePars.a, x,initialPars.a,':x', x,a_inverted,'k--x', x,truePars.b,'k:','linewidth',2);
xlabel('$x$ (km)','interpreter','latex');
legend('$\mathbf{a}$, true','$\mathbf{a}$, initial','$\mathbf{a}$, current','$\mathbf{b}$','interpreter','latex','Location','southwest');
title(sprintf('%s, iter %d', strrep(erase(resultFile, '.mat'), '_', '\_'), iter));
grid on;
set(gca,'fontsize',16);

assetsDir = fullfile('..', 'assets');
exportgraphics(gcf, fullfile(assetsDir, sprintf('%s_iter%d.png', erase(resultFile, '.mat'), iter)), 'Resolution', 300);