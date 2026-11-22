% Deterministic regularized inversion (L-BFGS) on the same noisy inverse-crime data as the samplers
% (m = 36, m_p = 26, eta = 0.02, noise seed 1), for a fair comparison with the posterior.
%
% Objective (unscaled misfit, as in optimize_general):  M(a) + R(a)
% The posterior is U(a) = M(a)/(sigma^2 dt_d) - log p(a), so multiplying by sigma^2 dt_d:
%   'tikhonov-grad'  beta = sigma^2 dt_d / gamma^2 gives the smoothness prior on the increments, so
%                    the minimizer is (close to) the MAP estimate. It lacks the weak anchor of a_1
%                    at mu (gamma0), the only difference to the prior.
%   'tv'             beta = sigma^2 dt_d / b for a Laplace prior on the increments with scale b.
%                    Default b = gamma*sqrt(2/pi), the mean |increment| of the smoothness prior.
% betaScale multiplies the default beta, e.g. for a sweep. Choose beta by the discrepancy
% principle: the final misfit M should be close to the noise level M(a_true) (printed at the end).
%
% Set before calling, e.g. matlab -batch "regType='tv'; betaScale=3; run_deterministic_noisy"
% Results in adjoint_rate_state_inversion/mat/inversion/rough_<optRoutine>_..._synth_eta002_..._<reg>.mat

bayesDir = fileparts(mfilename('fullpath'));
codeDir = fullfile(fileparts(bayesDir), 'adjoint_rate_state_inversion');
m = 36; m_p = 26; eta = 0.02;
if ~exist('regType', 'var'), regType = 'tikhonov-grad'; end   % 'tikhonov-grad', 'tv' or 'none'
if ~exist('betaScale', 'var'), betaScale = 1; end
if ~exist('tvEps', 'var'), tvEps = 1e-4; end    % smoothing of |.| in TV, small against increments ~1e-3
if ~exist('maxIter', 'var'), maxIter = 200; end
if ~exist('a0', 'var'), a0 = 0.0135; end        % initial guess, as the noise-free run
if ~exist('optRoutine', 'var'), optRoutine = 'lbfgs'; end   % see optimize_antiplaneshear, e.g. 'newton-tr'
if ~exist('loadOptPath', 'var'), loadOptPath = ''; end      % continue from a saved run (absolute path);
                                                             % maxIter then counts the earlier iterations

% Noise and prior of the samplers (setup_bayes_problem.m, smoothness_prior.m)
addpath(bayesDir);
S = load(fullfile(bayesDir, 'results', sprintf('hessian_m%d_eta%g.mat', m, eta)), 'aMap', 'aTrue');
R = load(fullfile(bayesDir, 'results', sprintf('hpcn_m%d_eta%g_s1_laplace_seed1.mat', m, eta)), 'noise');
noise = R.noise;
prior = smoothness_prior(m_p);
s2dt = noise.sigma^2*noise.dt_d;

switch regType
    case 'tikhonov-grad'
        regOpts = struct('type', regType, 'beta', betaScale*s2dt/prior.gamma^2);
    case 'tv'
        regOpts = struct('type', regType, 'beta', betaScale*s2dt/(prior.gamma*sqrt(2/pi)), 'eps', tvEps);
    case 'none'
        regOpts = struct('type', 'none');
    otherwise
        error('Unknown regType %s', regType);
end
if isfield(regOpts, 'beta')
    fprintf('%s, beta = %.3e (betaScale %g), sigma^2 dt_d = %.3e\n', regType, regOpts.beta, betaScale, s2dt);
end

oldDir = cd(codeDir);
restoreDir = onCleanup(@() cd(oldDir));
[aOpt, fOpt, flag, output, history, initialPars, truePars] = optimize_antiplaneshear(maxIter, m, m_p, ...
    'velocity', {'a'}, a0, false, false, false, loadOptPath, optRoutine, 2, regOpts, eta);

% Same noise as the samplers?
assert(abs(history.noise.sigma - noise.sigma) < 1e-12*noise.sigma && history.noise.seed == noise.seed, ...
    'Noise differs from the Bayesian setup');

M = history.fval(end) - history.reg(end);
fprintf('Final misfit M = %.4e, regularization R = %.4e, %d iterations\n', M, history.reg(end), output.iterations);
fprintf('Noise level (discrepancy target) ~ 0.84 N sigma^2 dt_d / 2 = %.4e\n', 0.84*noise.N*s2dt/2);
fprintf('|a - a_true| / |a_true| = %.3e, |a - a_MAP| / |a_MAP| = %.3e\n', ...
    norm(aOpt - S.aTrue)/norm(S.aTrue), norm(aOpt - S.aMap)/norm(S.aMap));
