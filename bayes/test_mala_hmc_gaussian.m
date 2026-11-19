% Verification of mala.m and hmc.m on the exactly known Gaussian posterior of test_pcn_gaussian.m
% (smoothness prior, linear forward model y = G a + noise). No PDE solve. The sample mean and
% standard deviations should match the exact ones to within Monte Carlo error.

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir, fullfile(fileparts(bayesDir), 'adjoint_rate_state_inversion', 'sbplib'));
m_p = 26;
prior = smoothness_prior(m_p);
rs = RandStream('mt19937ar', 'Seed', 3);
G = randn(rs, 10, m_p);
aTrue = prior.sample(1, rs);
sigma = 2e-3;
y = G*aTrue + sigma*randn(rs, 10, 1);

post.prior = prior;
post.prob.nFwd = 0;
post.Phi = @(a) 0.5*norm(G*a - y)^2/sigma^2;
post.U = @(a) post.Phi(a) - prior.logp(a);
post.Ug = @(a) deal(post.U(a), G'*(G*a - y)/sigma^2 - prior.grad(a));

Ppost = G'*G/sigma^2 + prior.P;
muPost = Ppost \ (G'*y/sigma^2 + prior.P*prior.mu);
sdPost = sqrt(diag(inv(Ppost)));
report = @(name, A) fprintf('%-5s max |mean - exact mean|/sd = %.3f, sd/exact sd in [%.3f, %.3f]\n', ...
    name, max(abs(mean(A, 2) - muPost)./sdPost), min(std(A, 0, 2)./sdPost), max(std(A, 0, 2)./sdPost));

opts.nAdapt = 2000;
opts.printEvery = inf;
opts.rs = RandStream('mt19937ar', 'Seed', 1);
chain = mala(post, prior.mu, 1e5, 1e-3, opts);
kept = opts.nAdapt+1:numel(chain.U);
fprintf('MALA: h = %.2e, acceptance rate = %.2f\n', chain.h(end), mean(chain.accepted(kept)));
report('', chain.a(:, kept));

opts.nAdapt = 500;
opts.rs = RandStream('mt19937ar', 'Seed', 2);
chain = hmc(post, prior.mu, 1e4, 1e-2, 20, opts);
kept = opts.nAdapt+1:numel(chain.U);
fprintf('HMC:  eps = %.2e, nLeap = 20, acceptance rate = %.2f\n', chain.eps(end), mean(chain.accepted(kept)));
report('', chain.a(:, kept));
