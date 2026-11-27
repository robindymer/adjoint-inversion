% Metropolis-within-Gibbs for the high-res posterior with the hierarchical (GSBL) prior of
% gsbl_map.m, sampling (a, theta) jointly with the noise inflation fixed:
%
%   theta_i | a  ~ InvGamma(beta + 1/2, vartheta + z_i^2/2)                  exact draw
%   a | theta    one MALA step for U(a | theta) = Phi(a) - log p(a | theta), preconditioned by
%                G(theta) = H_data + P(theta), with H_data = H - P(theta_MAP) the data part of the
%                Hessian at the MAP (gsbl_map.m). G depends only on theta, so the step leaves
%                a | theta invariant.
%
% One forward and one adjoint solve per iteration (the data gradient of the current a is kept, the
% prior part is recomputed for the new theta). The step size h is adapted during the first nAdapt
% iterations towards 57.4 % acceptance (burn-in). Starts from an overdispersed draw around a_MAP.
% seed, nIter, nAdapt, h0, m, m_p, tag can be set before calling. Saved in the format of mala.m
% (chain.a, U, accepted, ...) plus chain.theta, to results/gsbl_highres_m<m>_mp<m_p>_eta<eta><tag>_seed<seed>.mat.

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
if ~exist('m', 'var'), m = 126; end
if ~exist('m_p', 'var'), m_p = 126; end
eta = 0.02;
if ~exist('nIter', 'var'), nIter = 600; end
if ~exist('nAdapt', 'var'), nAdapt = 100; end
if ~exist('h0', 'var'), h0 = 0.3; end
if ~exist('seed', 'var'), seed = 1; end
if ~exist('tag', 'var'), tag = ''; end
resDir = fullfile(bayesDir, 'results');
mapFile = fullfile(resDir, sprintf('highres_map_m%d_mp%d_eta%g_gsbl%s.mat', m, m_p, eta, tag));
saveFile = fullfile(resDir, sprintf('gsbl_highres_m%d_mp%d_eta%g%s_seed%d.mat', m, m_p, eta, tag, seed));

S = load(mapFile);
[prob, noise, aTrue] = setup_highres_problem(m, m_p, eta, S.inflation);
beta = S.beta; vartheta = S.vartheta; gamma0 = S.gamma0;
aMap = S.aMap;
L = diff(eye(m_p));
priorOf = @(theta) smoothness_prior(m_p, sqrt(theta), gamma0);
Hdata = S.H - priorOf(S.theta).P;
c = 1/(noise.sigma^2*noise.dt_d);

rs = RandStream('mt19937ar', 'Seed', seed);
% Overdispersed start: Laplace draw (fixed theta_MAP) with doubled spread
RH = chol(S.H);
a = aMap + 2*(RH\randn(rs, m_p, 1));
[phi, gphi] = dataTerm(prob, c, a);
while ~isfinite(phi)
    a = aMap + 2*(RH\randn(rs, m_p, 1));
    [phi, gphi] = dataTerm(prob, c, a);
end

h = h0; logH = log(h);
chain = struct('a', zeros(m_p, nIter), 'theta', zeros(m_p - 1, nIter), 'U', zeros(1, nIter), ...
    'Phi', zeros(1, nIter), 'alpha', zeros(1, nIter), 'accepted', false(1, nIter), 'h', zeros(1, nIter));
chain.opts = struct('nAdapt', nAdapt, 'beta', beta, 'vartheta', vartheta, 'mapFile', mapFile);
[nFwd0, nAdj0] = deal(prob.nFwd, prob.nAdj);
t0 = tic;
for k = 1:nIter
    % theta | a: InvGamma(beta + 1/2, vartheta + z^2/2) = scale/Gamma(beta + 1/2, 1)
    theta = (vartheta + (L*a).^2/2)./gammaDraw(rs, beta + 0.5, m_p - 1);
    prior = priorOf(theta);
    R = chol(Hdata + prior.P);

    % a | theta: preconditioned MALA, a' = a - h/2 G^-1 grad U + sqrt(h) R\xi
    u = phi - prior.logp(a);
    g = gphi - prior.grad(a);
    xi = randn(rs, m_p, 1);
    ap = a - 0.5*h*(R\(R'\g)) + sqrt(h)*(R\xi);
    [phip, gphip] = dataTerm(prob, c, ap);
    if isfinite(phip)
        up = phip - prior.logp(ap);
        gp = gphip - prior.grad(ap);
        logqFwd = -0.5*(xi'*xi);
        logqBwd = -sum((R*(a - ap + 0.5*h*(R\(R'\gp)))).^2)/(2*h);
        alpha = min(1, exp(u - up + logqBwd - logqFwd));
    else
        alpha = 0;
    end
    if rand(rs) < alpha
        a = ap; phi = phip; gphi = gphip; u = up;
        chain.accepted(k) = true;
    end
    chain.a(:, k) = a; chain.theta(:, k) = theta; chain.U(k) = u; chain.Phi(k) = phi;
    chain.alpha(k) = alpha; chain.h(k) = h;
    if k <= nAdapt
        logH = logH + k^(-0.6)*(alpha - 0.574);
        h = exp(logH);
    end
    if mod(k, 10) == 0
        fprintf('[%6.0fs] it %5d  Phi %10.4e  h %8.2e  acc. rate (last 10) %.2f  jumps > 1e-3: %d\n', ...
            toc(t0), k, phi, h, mean(chain.accepted(k-9:k)), sum(abs(L*a) > 1e-3));
        c2 = chain;
        for f = {'a', 'theta', 'U', 'Phi', 'alpha', 'accepted', 'h'}
            c2.(f{1}) = c2.(f{1})(:, 1:k);
        end
        c2.nFwd = prob.nFwd - nFwd0; c2.nAdj = prob.nAdj - nAdj0; c2.time = toc(t0);
        save(saveFile, '-struct', 'S', 'aTrue', 'aMap', 'H', 'noise', 'm', 'm_p');
        save(saveFile, '-append', 'seed');
        T.chain = c2;
        save(saveFile, '-struct', 'T', '-append');
    end
end
fprintf('Final h = %.2e, acceptance after burn-in = %.2f, %d forward solves in %.0f s\n', h, ...
    mean(chain.accepted(nAdapt+1:end)), prob.nFwd - nFwd0, toc(t0));

% Negative log likelihood Phi = c*M(a) and its gradient; inf if the solve fails or a <= 0
function [phi, g] = dataTerm(prob, c, a)
    if any(a <= 0)
        phi = inf; g = nan(size(a)); return
    end
    [f, gf] = prob.fg(a);
    phi = c*f; g = c*gf;
    if isnan(phi) || any(isnan(g))
        phi = inf;
    end
end

% n draws of Gamma(shape, 1) (Marsaglia & Tsang; shape < 1 by the boost U^(1/shape))
function x = gammaDraw(rs, shape, n)
    s = shape + (shape < 1);
    d = s - 1/3; cc = 1/sqrt(9*d);
    x = zeros(n, 1);
    for i = 1:n
        while true
            z = randn(rs); v = (1 + cc*z)^3;
            if v > 0 && log(rand(rs)) < 0.5*z^2 + d - d*v + d*log(v)
                x(i) = d*v; break
            end
        end
    end
    if shape < 1
        x = x.*rand(rs, n, 1).^(1/shape);
    end
end
