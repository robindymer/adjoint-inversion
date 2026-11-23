% Run HMC on the noisy inverse-crime posterior (m = 36, m_p = 26, smoothness prior) with the Laplace
% approximation N(a_MAP, H^-1) as reference, i.e. mass matrix H (the posterior Hessian at the MAP)
% instead of the prior precision C^-1.
%
% The Hessian and a_MAP are cached in results/ (shared with run_hpcn.m). In the whitened coordinates
% the posterior is close to N(0, I), so eps ~ 1 scale: leapfrog is stable for eps < 2, and
% eps*nLeap ~ pi/2 (a quarter period) gives nearly independent proposals.

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
m = 36; m_p = 26; eta = 0.02;
if ~exist('nIter', 'var'), nIter = 1500; end
if ~exist('nLeap', 'var'), nLeap = 10; end
if ~exist('nAdapt', 'var'), nAdapt = 100; end
if ~exist('eps0', 'var'), eps0 = 0.15; end         % ~ (pi/2)/nLeap in whitened coordinates; adapted
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
    % Shortcut: one Newton step from a_true lands within 2.5e-4 of the newton-tr MAP (inverse crime
    % only; needs a_true). For high-res data use run_highres_map.m.
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
opts.saveFile = fullfile(resDir, sprintf('hmc_laplace_m%d_eta%g_L%d_%s_seed%d.mat', m, eta, nLeap, start, seed));
chain = hmc(post, a0, nIter, eps0, nLeap, opts);
save(opts.saveFile, 'chain', 'aTrue', 'aMap', 'H', 'noise', 'm', 'm_p', 'start', 'seed', 'nLeap');

kept = nAdapt+1:nIter;
fprintf('Final eps = %.2e, acceptance rate after burn-in = %.2f, %d forward + %d adjoint solves in %.0f s\n', ...
    chain.eps(end), mean(chain.accepted(kept)), chain.nFwd, chain.nAdj, chain.time);
