% Run prior-preconditioned MALA on the noisy inverse-crime posterior (m = 36, m_p = 26, smoothness
% prior). Laplace approximation

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
m = 36; m_p = 26; eta = 0.02;
if ~exist('nIter', 'var'), nIter = 15000; end
if ~exist('nAdapt', 'var'), nAdapt = 500; end
if ~exist('h0', 'var'), h0 = 1e-10; end          % rough guess, 1/(largest Hessian eigenvalue); adapted
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
opts.saveFile = fullfile(resDir, sprintf('mala_m%d_eta%g_%s_seed%d.mat', m, eta, start, seed));
chain = mala(post, a0, nIter, h0, opts);
save(opts.saveFile, 'chain', 'aTrue', 'noise', 'm', 'm_p', 'start', 'seed');

kept = nAdapt+1:nIter;
fprintf('Final h = %.2e, acceptance rate after burn-in = %.2f, %d forward + %d adjoint solves in %.0f s\n', ...
    chain.h(end), mean(chain.accepted(kept)), chain.nFwd, chain.nAdj, chain.time);
