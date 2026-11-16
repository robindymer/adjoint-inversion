% Verification of pcn.m on a posterior that is known exactly (script, Sec. 4.6, first test):
% the smoothness prior with a linear forward model y = G a + noise, so the posterior is Gaussian with
%
%   precision P_post = G'G/sigma^2 + P,   mean P_post \ (G'y/sigma^2 + P mu).
%
% No PDE solve. The sample mean and standard deviations should match the exact ones to within
% Monte Carlo error, which is large since the chain is strongly correlated: in a test run the mean
% was within 0.17 sd and the sds within 10 %.

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir, fullfile(fileparts(bayesDir), 'adjoint_rate_state_inversion', 'sbplib'));
m_p = 26;
prior = smoothness_prior(m_p);
rs = RandStream('mt19937ar', 'Seed', 3);
G = randn(rs, 10, m_p);
aTrue = prior.sample(1, rs);
sigma = 2e-3;
y = G*aTrue + sigma*randn(rs, 10, 1);

% Same fields as posterior.m that pcn uses
post.prior = prior;
post.prob.nFwd = 0;
post.Phi = @(a) 0.5*norm(G*a - y)^2/sigma^2;

% Exact posterior
Ppost = G'*G/sigma^2 + prior.P;
muPost = Ppost \ (G'*y/sigma^2 + prior.P*prior.mu);
sdPost = sqrt(diag(inv(Ppost)));

opts.nAdapt = 5000;
opts.printEvery = inf;
chain = pcn(post, prior.mu, 2e5, 0.1, opts);
A = chain.a(:, opts.nAdapt+1:end);
fprintf('s = %.3f, acceptance rate = %.2f\n', chain.s(end), mean(chain.accepted(opts.nAdapt+1:end)));
fprintf('max |mean - exact mean|/sd = %.3f, sd/exact sd in [%.3f, %.3f]\n', ...
    max(abs(mean(A, 2) - muPost)./sdPost), min(std(A, 0, 2)./sdPost), max(std(A, 0, 2)./sdPost));
