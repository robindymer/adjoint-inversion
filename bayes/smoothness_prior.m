% Gaussian smoothness prior on the increments of a (port of python/smoothness_prior.py):
%
%   a_1 ~ N(mu, gamma0^2),   a_{i+1} - a_i ~ N(0, gamma^2)
%
% i.e. p(a) ~ exp(-(a - mu)' P (a - mu) / 2) with P = L'L/gamma^2 + e_1 e_1'/gamma0^2, as
% gaussian_precision(n, gamma, 'differences', gamma0) in the course code.
%
% prior.logp(a), prior.grad(a)   log density (up to a constant) and its gradient
% prior.fromWhitened(w)          mu + R\w with P = R'R; w ~ N(0, I) gives a ~ prior (pCN works in w)
% prior.sample(n, rs)            n draws as columns, rs a RandStream (optional)
function prior = smoothness_prior(m_p, gamma, gamma0, mu)
    if nargin < 2, gamma = 1e-3; end
    if nargin < 3, gamma0 = 2e-3; end
    if nargin < 4, mu = 0.011; end   % b: neither VW nor VS is favoured

    L = diff(eye(m_p));         % (L a)_i = a_{i+1} - a_i
    P = L'*L/gamma^2;
    P(1,1) = P(1,1) + 1/gamma0^2;
    R = chol(P);                % upper triangular, P = R'R

    prior.m_p = m_p;
    prior.gamma = gamma;
    prior.gamma0 = gamma0;
    prior.mu = mu*ones(m_p, 1);
    prior.P = P;
    prior.R = R;
    prior.logp = @(a) -0.5*(a - prior.mu)'*P*(a - prior.mu);
    prior.grad = @(a) -P*(a - prior.mu);
    prior.fromWhitened = @(w) prior.mu + R\w;   % Cov = R^{-1} R^{-T} = P^{-1}
    prior.sample = @(n, rs) prior.fromWhitened(randnStream(rs, m_p, n));
end

function w = randnStream(rs, m, n)
    if nargin < 1 || isempty(rs)
        w = randn(m, n);
    else
        w = randn(rs, m, n);
    end
end
