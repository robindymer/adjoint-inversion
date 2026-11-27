% Smoothed total variation prior on the increments of a (cf. smoothness_prior.m):
%
%   a_1 ~ N(mu, gamma0^2),   p(a_{i+1} - a_i) ~ exp(-sqrt((a_{i+1} - a_i)^2 + eps^2)/b)
%
% i.e. a Laplace prior with scale b on the increments, smoothed at |Delta| < eps so that it is twice
% differentiable ('tv' in +adjopt/regularization.m). Favours piecewise constant a, but is not
% discretization invariant, so it is meant for MAP estimates rather than for sampling.
%
% prior.logp(a), prior.grad(a)   log density (up to a constant) and its gradient
% prior.hv(a, V)                 Hessian of -log p(a) applied to V (depends on a, unlike P)
% prior.mu, prior.R, prior.P     of the Gaussian smoothness prior with the same gamma0 and
%                                gamma = b*sqrt(pi/2), used only for whitening (run_highres_map.m)
%
% The default b = gamma*sqrt(2/pi) of the scaled Gaussian default gives the same mean |increment|.
function prior = tv_prior(m_p, b, eps, gamma0, mu)
    if nargin < 2 || isempty(b), b = 1e-3*sqrt(25/(m_p - 1))*sqrt(2/pi); end
    if nargin < 3 || isempty(eps), eps = 1e-5; end
    if nargin < 4, gamma0 = 2e-3; end
    if nargin < 5, mu = 0.011; end

    prior = smoothness_prior(m_p, b*sqrt(pi/2), gamma0, mu);
    prior.b = b;
    prior.eps = eps;
    L = sparse(diff(eye(m_p)));
    e1 = sparse(1, 1, 1, m_p, 1);
    mu1 = prior.mu(1);

    prior.logp = @(a) -sum(sqrt((L*a).^2 + eps^2) - eps)/b - (a(1) - mu1)^2/(2*gamma0^2);
    prior.grad = @(a) -L'*((L*a)./sqrt((L*a).^2 + eps^2))/b - e1*(a(1) - mu1)/gamma0^2;
    prior.hv = @(a, V) L'*((eps^2./((L*a).^2 + eps^2).^(3/2)).*(L*V))/b + e1*V(1,:)/gamma0^2;
    prior.sample = [];
end
