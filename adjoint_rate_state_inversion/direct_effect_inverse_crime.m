addSubpaths;
maxIter = 200; % max number of optimization iterations
m = 76; % Grid points in x-direction
m_p = 26; % Grid points on parameter grid
scalingStudy = false; 
misfitType = 'velocity';
optPars = {'a'}; % parameter to invert for
val = 0.0135; % initial value
% lbfgs, newton-tr, CG, newton-tr-noprec, fminunc-tr, newton-cg-ls (see adjopt.optimize_general).
% Can be set before running the script, e.g. matlab -batch "opt_routine_str='newton-tr'; direct_effect_inverse_crime"
if ~exist('opt_routine_str', 'var')
    opt_routine_str = 'lbfgs';
end
% Regularization: 'none', 'tikhonov', 'tikhonov-grad', 'tv' or 'lasso' (see adjopt.regularization),
% with weight reg_beta. reg_eps smooths the absolute value in 'tv' and 'lasso'. E.g.
% matlab -batch "opt_routine_str='newton-cg-ls'; reg_type='tv'; reg_beta=1e-2; direct_effect_inverse_crime"
if ~exist('reg_type', 'var')
    reg_type = 'none';
end
if ~exist('reg_beta', 'var')
    reg_beta = 0;
end
regOpts = struct('type', reg_type, 'beta', reg_beta);
if exist('reg_eps', 'var')
    regOpts.eps = reg_eps;
end
% Hybrid: start from iteration x0_iter of the lbfgs run with the same settings, e.g. x0_iter=30
if ~exist('x0_iter', 'var')
    x0_iter = [];
end
fprintf('RUNNING: %s, regularization %s with beta = %g\n', opt_routine_str, reg_type, reg_beta);
[parsVec, misFit, flag, output, history, initialPars, truePars, parSet] = optimize_antiplaneshear(maxIter, m, m_p, misfitType, optPars, val, scalingStudy, [], [], [], opt_routine_str, 2, regOpts, x0_iter);

figure();
x = linspace(-10,10,m_p); % Parameter grid
a_inverted = history.x(:,end); % Inverted results at final iteration
plot(x,truePars.a, x,initialPars.a,':x', x,a_inverted,'k--x', x,truePars.b,'k:','linewidth',2);
xlabel('$x$ (km)','interpreter','latex');
legend('$\mathbf{a}$, true','$\mathbf{a}$, initial','$\mathbf{a}$, current','$\mathbf{b}$','interpreter','latex','Location','southwest');
grid on;
set(gca,'fontsize',16);
