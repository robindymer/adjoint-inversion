% Run pCN on the noisy inverse-crime posterior (m = 36, m_p = 26, smoothness prior).
% The first nAdapt iterations tune s and are burn-in. The chain is saved to results/ as it runs.

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
m = 36; m_p = 26; eta = 0.02;
nIter  = 2000;
nAdapt = 500;
% TODO: Verify the s from Robbins-Monro manually
% Done by testing a few values of s until 20-30% acceptance rate is achieved
s0     = 1e-4;     % rough guess from the posterior width estimate; adapted during burn-in
start  = 'true';   % 'true': start at a_true, 'prior': at a prior draw
seed   = 0;

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
