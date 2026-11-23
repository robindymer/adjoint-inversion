% MAP estimate, noise inflation and Laplace approximation for the high-res data
% (setup_highres_problem.m, default m = 151, m_p = 26, smoothness prior).
%
% 1) MAP of the posterior U(a) = M(a)/(sigma^2 dt_d) - log p(a) by trust-region Newton-CG with
%    exact Hessian-vector products (the full U, so the prior's anchor of a_1 is included).
% 2) Discrepancy principle for the model error: the misfit at the MAP should be the noise level
%    M_noise. If M(a_MAP) > M_noise, the model error is absorbed by inflating the noise variance,
%    inflation = M(a_MAP)/M_noise, and the MAP is recomputed (warm start), since the inflation
%    changes the weight of the prior. Repeated until the inflation changes by less than 5 %.
% 3) Hessian of U at the MAP from m_p Hessian-vector products (posterior_hessian.m).
%
% Saved to results/highres_map_m<m>_mp<m_p>_eta<eta>.mat (also after each round), which the
% chains of run_highres_hpcn.m load. m, m_p, maxIter can be set before calling. With resume = true,
% starts from the a_MAP and inflation saved there, e.g. if the job hit the time limit or the
% Hessian was not positive definite (MAP not converged): matlab -batch "resume=true; run_highres_map"

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
if ~exist('m', 'var'), m = 151; end
if ~exist('m_p', 'var'), m_p = 26; end
eta = 0.02;
if ~exist('maxIter', 'var'), maxIter = 60; end    % per round; inverse-crime runs needed 21
if ~exist('a0', 'var'), a0 = 0.0135; end          % constant initial guess, as before
if ~exist('resume', 'var'), resume = false; end
maxRounds = 4;

resDir = fullfile(bayesDir, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
outFile = fullfile(resDir, sprintf('highres_map_m%d_mp%d_eta%g.mat', m, m_p, eta));

[prob, noise, aTrue] = setup_highres_problem(m, m_p, eta, 1);
prior = smoothness_prior(m_p);

aMap = a0*ones(m_p, 1);
inflation = 1;
infHist = [];
if resume
    R = load(outFile, 'aMap', 'inflation', 'infHist');
    aMap = R.aMap; inflation = R.inflation; infHist = R.infHist;
    if ~isempty(infHist), inflation = max(1, infHist(end, 3)); end
    fprintf('Resuming from %s with inflation %.3f\n', outFile, inflation);
end
for round = 1:maxRounds
    noise.inflation = inflation;
    noise.sigma = sqrt(inflation)*noise.sigma0;
    post = posterior(prob, noise, prior);

    tropts.maxIter = maxIter;
    tropts.precond = true;
    tropts.OutputFcn = @(x, f, g, info) printIter(x, f, g, prob, post);
    t0 = tic;
    [aMap, uMap, flag, output] = adjopt.trustRegionNewtonCG(post.Ug, post.Uhv, aMap, zeros(m_p, 1), tropts);
    Mmap = prob.f(aMap);
    newInflation = max(1, Mmap/noise.Mnoise);
    infHist(end+1, :) = [inflation, Mmap, newInflation, output.iterations, flag, toc(t0)]; %#ok<SAGROW>
    fprintf(['Round %d: inflation %.3f, %d iterations (flag %d), M(a_MAP) = %.4e = %.3f M_noise, ' ...
        '|a_MAP - a_true|/|a_true| = %.3e\n'], round, inflation, output.iterations, flag, Mmap, ...
        Mmap/noise.Mnoise, norm(aMap - aTrue)/norm(aTrue));
    save(outFile, 'aMap', 'aTrue', 'noise', 'inflation', 'infHist', 'm', 'm_p', 'eta');
    if abs(newInflation - inflation) <= 0.05*inflation
        break
    end
    inflation = newInflation;
end
if round == maxRounds
    warning('Inflation not converged after %d rounds', maxRounds);
end

% Laplace approximation at the MAP, with the final inflation. H is saved before the check, so an
% indefinite H (MAP not converged) does not lose the solves.
H = post.Uhv(aMap, eye(m_p));
H = (H + H')/2;
save(outFile, 'aMap', 'aTrue', 'noise', 'inflation', 'infHist', 'H', 'm', 'm_p', 'eta');
if min(eig(H)) <= 0
    error(['The Hessian at a_MAP is not positive definite (min eigenvalue %.3e): the MAP is not ' ...
        'converged. Rerun with resume = true (and a larger maxIter).'], min(eig(H)));
end
R = prior.R;
lambda = sort(eig((R'\(H - prior.P))/R), 'descend');
fprintf('Eigenvalues of the prior-preconditioned data Hessian: max %.2e, %d of %d above 1\n', ...
    lambda(1), sum(lambda > 1), m_p);
fprintf('Solves: %d forward, %d adjoint, %d Hv\n', prob.nFwd, prob.nAdj, prob.nHv);
save(outFile, 'aMap', 'aTrue', 'noise', 'inflation', 'infHist', 'H', 'lambda', 'm', 'm_p', 'eta');
fprintf('Saved %s\n', outFile);

function stop = printIter(x, f, g, prob, post)
    fprintf('  U = %.6e, |g| = %.3e, M = %.4e, solves %d/%d/%d\n', f, norm(g), (f + post.prior.logp(x))/post.c, ...
        prob.nFwd, prob.nAdj, prob.nHv);
    stop = false;
end
