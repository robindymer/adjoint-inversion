% MAP and Laplace approximation for the high-res data with the hierarchical (GSBL) prior on the
% increments z = L a (port of python/hierarchical_prior.py):
%
%   a_1 ~ N(mu, gamma0^2),   z_i | theta_i ~ N(0, theta_i),   theta_i ~ InvGamma(beta, vartheta)
%
% A few theta_i can be large, so a(x) can jump at a few nodes and is nearly constant elsewhere,
% like the true a(x). The joint MAP of (a, theta) by the iterative alternating scheme (IAS):
%   theta-step  theta_i = (vartheta + z_i^2/2)/(beta + 3/2)   (closed form)
%   a-step      trust-region Newton-CG on U(a) with the Gaussian prior smoothness_prior(m_p, sqrt(theta)),
%               in the whitened coordinates w = R(a - mu) of that prior, since the small theta_i make
%               the prior precision badly conditioned
% interleaved with the noise inflation of run_highres_map.m (inflation = M(a)/M_noise). Then the
% Hessian of U(a | theta_MAP) at a_MAP, i.e. the Laplace approximation of a for fixed theta (it
% ignores the uncertainty in theta, so it is too narrow at the jumps).
%
% vartheta scales with the node spacing: 1e-8 at m_p = 26 (python), vartheta ~ dx otherwise.
% Saved to results/highres_map_m<m>_mp<m_p>_eta<eta>_gsbl<tag>.mat, in the format of
% run_highres_map.m. m, m_p, beta, vartheta, warmFile, maxOuter, innerIter, trOpts (trust region,
% defaults tuned on the m_p = 126 Gaussian MAP) can be set before calling. noHessian = true stops
% after the MAP (then lowrank_hessian.m, which uses the saved theta).

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
if ~exist('m', 'var'), m = 126; end
if ~exist('m_p', 'var'), m_p = 126; end
eta = 0.02;
if ~exist('beta', 'var'), beta = 0.5; end
if ~exist('vartheta', 'var'), vartheta = 1e-8*25/(m_p - 1); end
if ~exist('gamma0', 'var'), gamma0 = 2e-3; end
if ~exist('maxOuter', 'var'), maxOuter = 40; end
if ~exist('innerIter', 'var'), innerIter = 5; end
if ~exist('trOpts', 'var'), trOpts = struct('DeltaMax', 1e4, 'maxStep', 1e3, 'etaMax', 1e-3, 'maxCG', 60); end
if ~exist('tag', 'var'), tag = ''; end
if ~exist('noHessian', 'var'), noHessian = false; end
resDir = fullfile(bayesDir, 'results');
if ~exist('warmFile', 'var'), warmFile = fullfile(resDir, sprintf('highres_map_m%d_mp%d_eta%g.mat', m, m_p, eta)); end
outFile = fullfile(resDir, sprintf('highres_map_m%d_mp%d_eta%g_gsbl%s.mat', m, m_p, eta, tag));

[prob, noise, aTrue] = setup_highres_problem(m, m_p, eta, 1);
W = load(warmFile, 'aMap', 'inflation');
aMap = interp1(linspace(-10, 10, numel(W.aMap)), W.aMap, linspace(-10, 10, m_p)');
inflation = W.inflation;
fprintf('GSBL: beta = %g, vartheta = %.2e, warm start from %s (inflation %.3f)\n', beta, vartheta, ...
    warmFile, inflation);

L = diff(eye(m_p));
hist = [];
for outer = 1:maxOuter
    theta = (vartheta + (L*aMap).^2/2)/(beta + 1.5);
    prior = smoothness_prior(m_p, sqrt(theta), gamma0);
    noise.inflation = inflation;
    noise.sigma = sqrt(inflation)*noise.sigma0;
    post = posterior(prob, noise, prior);

    % a-step in whitened coordinates
    R = prior.R;
    toA = @(w) prior.mu + R\w;
    tropts = trOpts;
    tropts.maxIter = innerIter;
    tropts.precond = true;
    aOld = aMap;
    t0 = tic;
    [w, uMap, flag, output] = adjopt.trustRegionNewtonCG(@(w) fgWhite(post, R, toA(w)), ...
        @(w, V) R'\post.Uhv(toA(w), R\V), R*(aMap - prior.mu), [], tropts);
    aMap = toA(w);
    Mmap = prob.f(aMap);
    newInflation = max(1, Mmap/noise.Mnoise);
    da = norm(aMap - aOld)/norm(aOld);
    nJump = sum(abs(L*aMap) > 1e-3);
    hist(end+1, :) = [inflation, Mmap, newInflation, output.iterations, flag, da, nJump, toc(t0)]; %#ok<SAGROW>
    fprintf(['Outer %d: inflation %.3f, %d iterations (flag %d), M(a_MAP) = %.3f M_noise, |da|/|a| = %.2e, ' ...
        '%d increments > 1e-3, |a_MAP - a_true|/|a_true| = %.3e, solves %d/%d/%d\n'], outer, inflation, ...
        output.iterations, flag, Mmap/noise.Mnoise, da, nJump, norm(aMap - aTrue)/norm(aTrue), ...
        prob.nFwd, prob.nAdj, prob.nHv);
    save(outFile, 'aMap', 'aTrue', 'noise', 'inflation', 'theta', 'hist', 'beta', 'vartheta', 'gamma0', ...
        'm', 'm_p', 'eta');
    if da < 1e-3 && abs(newInflation - inflation) <= 0.05*inflation
        break
    end
    inflation = newInflation;
end
infHist = hist(:, 1:3);
if noHessian
    save(outFile, 'aMap', 'aTrue', 'noise', 'inflation', 'infHist', 'theta', 'hist', 'beta', ...
        'vartheta', 'gamma0', 'm', 'm_p', 'eta');
    fprintf('Saved %s (no Hessian)\n', outFile);
    return
end

% Laplace approximation of a | theta at the final (a, theta, inflation)
H = post.Uhv(aMap, eye(m_p));
H = (H + H')/2;
save(outFile, 'aMap', 'aTrue', 'noise', 'inflation', 'infHist', 'theta', 'hist', 'H', 'beta', ...
    'vartheta', 'gamma0', 'm', 'm_p', 'eta');
if min(eig(H)) <= 0
    error('The Hessian at a_MAP is not positive definite (min eigenvalue %.3e)', min(eig(H)));
end
fprintf('Solves: %d forward, %d adjoint, %d Hv\nSaved %s\n', prob.nFwd, prob.nAdj, prob.nHv, outFile);

function [u, g] = fgWhite(post, R, a)
    [u, g] = post.Ug(a);
    g = R'\g;
end
