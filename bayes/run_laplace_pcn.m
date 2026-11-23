% Run H-pCN on the noisy inverse-crime posterior (m = 36, m_p = 26, smoothness prior): pCN with the
% Laplace approximation N(a_MAP, H^-1) as reference instead of the prior.
%
% The Hessian H is computed at a_true (26 Hv products) and cached in results/. The MAP estimate is
% approximated by one Newton step from a_true, which is close to the MAP for inverse-crime data.
% s is fixed (nAdapt = 0). If the Laplace approximation were exact, phi = U - |w|^2/2 would be
% constant and every proposal accepted, even with s = 1, which gives independent proposals. So s
% follows from the spread of phi in a pilot run: s = 1 if std(chain.phi) <~ 1, smaller otherwise.
% The 0.234 rule of the random walk does not apply.
% seed, start, nIter and s can be set before calling, e.g. for several chains in parallel
% (see run_chains.sh): matlab -batch "seed=2; start='laplace'; run_hpcn". With hessianOnly = true,
% only the Hessian is computed and cached, so that parallel chains do not all compute it.

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
m = 36; m_p = 26; eta = 0.02;
if ~exist('nIter', 'var'), nIter = 1000; end
nAdapt = 0;
if ~exist('s', 'var'), s = 0.5; end              % use s = 1 (independent proposals) if std(chain.phi) <~ 1
if ~exist('start', 'var'), start = 'map'; end    % 'map': at a_MAP, 'laplace': overdispersed Laplace draw
if ~exist('seed', 'var'), seed = 0; end
if ~exist('hessianOnly', 'var'), hessianOnly = false; end

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
    % Shortcut: one Newton step from a_true lands within 2.5e-4 of the newton-tr MAP (inverse crime
    % only; needs a_true). For high-res data use run_highres_map.m.
    aMap = aTrue - H\g;   % one Newton step
    fprintf('Newton step: U(a_true) = %.4f, U(a_MAP) = %.4f\n', uTrue, post.U(aMap));
    save(hessFile, 'H', 'lambda', 'aMap', 'aTrue');
end
fprintf('Eigenvalues of the prior-preconditioned data Hessian: max %.2e, min %.2e, %d of %d above 1\n', ...
    lambda(1), lambda(end), sum(lambda > 1), m_p);
if hessianOnly
    return
end

RH = chol(H);
rs = RandStream('mt19937ar', 'Seed', seed);
switch start
    case 'map'
        a0 = aMap;
    case 'laplace'
        % Draw from the Laplace approximation with twice its standard deviation, so that the
        % chains start dispersed around the posterior (script, Sec. 4.6)
        a0 = aMap + 2*(RH\randn(rs, m_p, 1));
        while ~isfinite(post.U(a0))
            a0 = aMap + 2*(RH\randn(rs, m_p, 1));
        end
end

opts.ref = struct('mu', aMap, 'R', RH);
opts.nAdapt = nAdapt;
opts.rs = rs;
opts.saveFile = fullfile(resDir, sprintf('hpcn_m%d_eta%g_s%g_%s_seed%d.mat', m, eta, s, start, seed));
chain = pcn(post, a0, nIter, s, opts);
save(opts.saveFile, 'chain', 'aTrue', 'aMap', 'H', 'noise', 'm', 'm_p', 's', 'start', 'seed');

kept = nAdapt+1:nIter;
fprintf('s = %.2e, acceptance rate = %.2f, %d forward solves in %.0f s\n', ...
    s, mean(chain.accepted(kept)), chain.nFwd, chain.time);
