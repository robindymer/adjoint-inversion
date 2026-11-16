% Run H-pCN on the noisy inverse-crime posterior (m = 36, m_p = 26, smoothness prior): pCN with the
% Laplace approximation N(a_MAP, H^-1) as reference instead of the prior.
%
% The Hessian H is computed at a_true (26 Hv products) and cached in results/. The MAP estimate is
% approximated by one Newton step from a_true, which is close to the MAP for inverse-crime data.
% s is fixed (nAdapt = 0) and chosen by pilot runs for an acceptance rate of about 0.25. If the
% Laplace approximation were exact, every proposal would be accepted, even with s = 1.

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
m = 36; m_p = 26; eta = 0.02;
nIter  = 1000;
nAdapt = 0;
s      = 0.5;      % pilot value, tune for acceptance 0.2-0.3
seed   = 0;

[prob, noise, aTrue] = setup_bayes_problem(m, m_p, eta, 1);
prior = smoothness_prior(m_p);
post = posterior(prob, noise, prior);

resDir = fullfile(bayesDir, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
hessFile = fullfile(resDir, sprintf('hessian_m%d_eta%g.mat', m, eta));
if exist(hessFile, 'file')
    load(hessFile, 'H', 'lambda', 'aMap');
else
    [H, lambda] = posterior_hessian(post, aTrue);
    [uTrue, g] = post.Ug(aTrue);
    aMap = aTrue - H\g;   % one Newton step
    fprintf('Newton step: U(a_true) = %.4f, U(a_MAP) = %.4f\n', uTrue, post.U(aMap));
    save(hessFile, 'H', 'lambda', 'aMap', 'aTrue');
end
fprintf('Eigenvalues of the prior-preconditioned data Hessian: max %.2e, min %.2e, %d of %d above 1\n', ...
    lambda(1), lambda(end), sum(lambda > 1), m_p);

opts.ref = struct('mu', aMap, 'R', chol(H));
opts.nAdapt = nAdapt;
opts.rs = RandStream('mt19937ar', 'Seed', seed);
opts.saveFile = fullfile(resDir, sprintf('hpcn_m%d_eta%g_s%g_seed%d.mat', m, eta, s, seed));
chain = pcn(post, aMap, nIter, s, opts);
save(opts.saveFile, 'chain', 'aTrue', 'aMap', 'H', 'noise', 'm', 'm_p', 's', 'seed');

kept = nAdapt+1:nIter;
fprintf('s = %.2e, acceptance rate = %.2f, %d forward solves in %.0f s\n', ...
    s, mean(chain.accepted(kept)), chain.nFwd, chain.time);
