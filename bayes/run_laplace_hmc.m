% Run HMC with mass matrix C^-1 (prior precision) on the noisy inverse-crime posterior (m = 36,
% m_p = 26, smoothness prior). Laplace approximation

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
m = 36; m_p = 26; eta = 0.02;
if ~exist('nIter', 'var'), nIter = 1500; end
if ~exist('nLeap', 'var'), nLeap = 10; end
if ~exist('nAdapt', 'var'), nAdapt = 100; end
if ~exist('eps0', 'var'), eps0 = 1e-5; end       % rough guess, 1/sqrt(largest Hessian eigenvalue); adapted
if ~exist('start', 'var'), start = 'prior'; end  % 'true': start at a_true, 'prior': at a prior draw
if ~exist('seed', 'var'), seed = 0; end

[prob, noise, aTrue] = setup_bayes_problem(m, m_p, eta, 1);
prior = smoothness_prior(m_p);
post = posterior(prob, noise, prior);

rs = RandStream('mt19937ar', 'Seed', seed);
switch start
    case 'true'
        a0 = aTrue;
    case 'prior'
        a0 = prior.sample(1, rs);
        while ~isfinite(post.Phi(a0))
            a0 = prior.sample(1, rs);
        end
end

resDir = fullfile(bayesDir, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
opts.nAdapt = nAdapt;
opts.rs = rs;
opts.saveFile = fullfile(resDir, sprintf('hmc_m%d_eta%g_L%d_%s_seed%d.mat', m, eta, nLeap, start, seed));
chain = hmc(post, a0, nIter, eps0, nLeap, opts);
save(opts.saveFile, 'chain', 'aTrue', 'noise', 'm', 'm_p', 'start', 'seed', 'nLeap');

kept = nAdapt+1:nIter;
fprintf('Final eps = %.2e, acceptance rate after burn-in = %.2f, %d forward + %d adjoint solves in %.0f s\n', ...
    chain.eps(end), mean(chain.accepted(kept)), chain.nFwd, chain.nAdj, chain.time);
