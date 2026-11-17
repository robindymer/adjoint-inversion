% Run pCN on the noisy inverse-crime posterior (m = 36, m_p = 26, smoothness prior).
% The first nAdapt iterations tune s and are burn-in. The chain is saved to results/ as it runs.
% seed, start, nIter, nAdapt and s0 can be set before calling, e.g. for several chains in parallel
% (see run_chains.sh): matlab -batch "seed=2; start='prior'; run_pcn"

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
m = 36; m_p = 26; eta = 0.02;
if ~exist('nIter', 'var'), nIter = 2000; end
if ~exist('nAdapt', 'var'), nAdapt = 500; end
% TODO: Verify the s from Robbins-Monro manually
% Done by testing a few values of s until 20-30% acceptance rate is achieved
if ~exist('s0', 'var'), s0 = 1e-4; end         % rough guess from the posterior width estimate; adapted during burn-in
if ~exist('start', 'var'), start = 'true'; end % 'true': start at a_true, 'prior': at a prior draw
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
opts.saveFile = fullfile(resDir, sprintf('pcn_m%d_eta%g_%s_seed%d.mat', m, eta, start, seed));
chain = pcn(post, a0, nIter, s0, opts);
save(opts.saveFile, 'chain', 'aTrue', 'noise', 'm', 'm_p', 'start', 'seed');

kept = nAdapt+1:nIter;
fprintf('Final s = %.2e, acceptance rate after burn-in = %.2f, %d forward solves in %.0f s\n', ...
    chain.s(end), mean(chain.accepted(kept)), chain.nFwd, chain.time);
