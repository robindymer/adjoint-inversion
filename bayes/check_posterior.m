% Sanity checks of the Bayesian setup before sampling (step 2 of the plan, first part):
%   1) likelihood scale: c*M(a_true) should be close to N/2 (measured ~0.84*N/2 earlier)
%   2) prior: squared Mahalanobis distance of a_true and prior draws
%   3) finite-difference check of the gradient of U

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
[prob, noise, aTrue] = setup_bayes_problem(36, 26, 0.02, 1);
prior = smoothness_prior(26);
post = posterior(prob, noise, prior);

% 1) Likelihood scale
Ut = post.c*prob.f(aTrue);
fprintf('c*M(a_true) = %.4e, N/2 = %.4e, ratio = %.3f\n', Ut, noise.N/2, Ut/(noise.N/2));

% 2) Prior
fprintf('-2 log p(a_true) = %.1f (draws: about m_p = %d)\n', -2*prior.logp(aTrue), prior.m_p);
draws = prior.sample(1e4, RandStream('mt19937ar', 'Seed', 0));
fprintf('Prior: P(a <= 0 somewhere) = %.3f\n', mean(any(draws <= 0, 1)));

% 3) Gradient check at a point away from the truth
rs = RandStream('mt19937ar', 'Seed', 2);
a0 = aTrue + 2e-4*randn(rs, size(aTrue));
dir = randn(rs, size(aTrue)); dir = dir/norm(dir);
[u0, g0] = post.Ug(a0);
fprintf('U(a0) = %.6e, directional derivative g''d = %.6e\n', u0, g0'*dir);
for h = [1e-5 1e-6 1e-7]
    fd = (post.U(a0 + h*dir) - post.U(a0 - h*dir))/(2*h);
    fprintf('  h = %.0e: central FD = %.6e, rel. error = %.2e\n', h, fd, abs(fd - g0'*dir)/abs(g0'*dir));
end
