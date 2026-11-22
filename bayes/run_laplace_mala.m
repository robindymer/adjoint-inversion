% Run MALA on the noisy inverse-crime posterior (m = 36, m_p = 26, smoothness prior) preconditioned
% with the Laplace approximation N(a_MAP, H^-1), i.e. preconditioner H^-1 instead of the prior
% covariance C.
%
% The Hessian and a_MAP are cached in results/ (shared with run_hpcn.m). In the whitened coordinates
% the posterior is close to N(0, I), so h ~ 1 scale (for a standard Gaussian in d = 26 the optimal
% h is about 1.65^2 d^(-1/3) ~ 0.9).

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
m = 36; m_p = 26; eta = 0.02;
if ~exist('nIter', 'var'), nIter = 15000; end
if ~exist('nAdapt', 'var'), nAdapt = 500; end
if ~exist('h0', 'var'), h0 = 0.5; end              % whitened coordinates; adapted
if ~exist('start', 'var'), start = 'laplace'; end  % 'map', 'laplace' (overdispersed draw), 'true', 'prior'
if ~exist('seed', 'var'), seed = 0; end

[prob, noise, aTrue] = setup_bayes_problem(m, m_p, eta, 1);
prior = smoothness_prior(m_p);
post = posterior(prob, noise, prior);

resDir = fullfile(bayesDir, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
hessFile = fullfile(resDir, sprintf('hessian_m%d_eta%g.mat', m, eta));
if exist(hessFile, 'file')
    load(hessFile, 'H', 'aMap');
else
    [H, lambda] = posterior_hessian(post, aTrue);
    [uTrue, g] = post.Ug(aTrue);
    aMap = aTrue - H\g;   % one Newton step
    fprintf('Newton step: U(a_true) = %.4f, U(a_MAP) = %.4f\n', uTrue, post.U(aMap));
    save(hessFile, 'H', 'lambda', 'aMap', 'aTrue');
end
RH = chol(H);

rs = RandStream('mt19937ar', 'Seed', seed);
switch start
    case 'map'
        a0 = aMap;
    case 'laplace'
        % Draw from the Laplace approximation with twice its standard deviation (script, Sec. 4.6)
        a0 = aMap + 2*(RH\randn(rs, m_p, 1));
        while ~isfinite(post.U(a0))
            a0 = aMap + 2*(RH\randn(rs, m_p, 1));
        end
    case 'true'
        a0 = aTrue;
    case 'prior'
        a0 = prior.sample(1, rs);
        while ~isfinite(post.Phi(a0))
            a0 = prior.sample(1, rs);
        end
end

opts.ref = struct('mu', aMap, 'R', RH);
opts.nAdapt = nAdapt;
opts.rs = rs;
opts.saveFile = fullfile(resDir, sprintf('mala_laplace_m%d_eta%g_%s_seed%d.mat', m, eta, start, seed));
chain = mala(post, a0, nIter, h0, opts);
save(opts.saveFile, 'chain', 'aTrue', 'aMap', 'H', 'noise', 'm', 'm_p', 'start', 'seed');

kept = nAdapt+1:nIter;
fprintf('Final h = %.2e, acceptance rate after burn-in = %.2f, %d forward + %d adjoint solves in %.0f s\n', ...
    chain.h(end), mean(chain.accepted(kept)), chain.nFwd, chain.nAdj, chain.time);
