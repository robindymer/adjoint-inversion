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
% Saved to results/highres_map_m<m>_mp<m_p>_eta<eta><tag>.mat (also after each round), which the
% chains of run_highres_hpcn.m load. m, m_p, maxIter, gammaFactor (prior, see highres_prior.m) can
% be set before calling. warmFile: start from the a_MAP of another MAP file, interpolated to m_p.
% whiten = true runs Newton-CG in the whitened prior coordinates w = R(a - mu) (needed at large m_p,
% where the prior precision is badly conditioned and steps in a stay tiny). a_MAP is checkpointed
% at every iteration. trOpts: extra options of the trust region (e.g. DeltaMax, maxStep). uTol > 0 ends
% a round when U has decreased by less than uTol over the last 5 iterations (U is a negative log
% density, so a change below 1 is negligible); 0 keeps only the gradient test.
% inflation0: initial inflation. With resume = true, starts from the a_MAP and inflation saved there
% (or inflation0 if set), e.g. if the job hit the time limit or the Hessian was not positive
% definite (MAP not converged):
% matlab -batch "resume=true; run_highres_map"
% priorType = 'tv' uses the smoothed TV prior (highres_prior.m). noHessian = true stops after the MAP
% (the Hessian then needs lowrank_hessian.m or hessian_columns.m).

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
if ~exist('m', 'var'), m = 151; end
if ~exist('m_p', 'var'), m_p = 26; end
eta = 0.02;
if ~exist('maxIter', 'var'), maxIter = 150; end   % per round; 60 was hit in every round
if ~exist('a0', 'var'), a0 = 0.0135; end          % constant initial guess, as before
if ~exist('resume', 'var'), resume = false; end
if ~exist('maxRounds', 'var'), maxRounds = 8; end
if ~exist('gammaFactor', 'var'), gammaFactor = 1; end
if ~exist('warmFile', 'var'), warmFile = ''; end
if ~exist('whiten', 'var'), whiten = false; end
if ~exist('trOpts', 'var'), trOpts = struct(); end
if ~exist('uTol', 'var'), uTol = 0; end
if ~exist('priorType', 'var'), priorType = 'gaussian'; end
if ~exist('noHessian', 'var'), noHessian = false; end

resDir = fullfile(bayesDir, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
[prior, tag] = highres_prior(m_p, gammaFactor, priorType);
outFile = fullfile(resDir, sprintf('highres_map_m%d_mp%d_eta%g%s.mat', m, m_p, eta, tag));

[prob, noise, aTrue] = setup_highres_problem(m, m_p, eta, 1);

aMap = a0*ones(m_p, 1);
if ~isempty(warmFile)
    W = load(warmFile, 'aMap');
    aMap = interp1(linspace(-10, 10, numel(W.aMap)), W.aMap, linspace(-10, 10, m_p)');
    fprintf('Warm start from %s\n', warmFile);
end
inflation = 1;
if exist('inflation0', 'var'), inflation = inflation0; end   % initial guess, saves a round
infHist = [];
if ~resume
    save(outFile, 'aMap', 'aTrue', 'noise', 'inflation', 'infHist', 'm', 'm_p', 'eta', 'gammaFactor');
end
if resume
    R = load(outFile, 'aMap', 'inflation', 'infHist');
    aMap = R.aMap; inflation = R.inflation; infHist = R.infHist;
    if ~isempty(infHist), inflation = max(1, infHist(end, 3)); end
    if exist('inflation0', 'var'), inflation = inflation0; end
    fprintf('Resuming from %s with inflation %.3f\n', outFile, inflation);
end
for round = 1:maxRounds
    noise.inflation = inflation;
    noise.sigma = sqrt(inflation)*noise.sigma0;
    post = posterior(prob, noise, prior);

    tropts = trOpts;
    tropts.maxIter = maxIter;
    tropts.precond = true;
    t0 = tic;
    if whiten
        R = prior.R;
        toA = @(w) prior.mu + R\w;
        tropts.OutputFcn = @(w, f, g, info) printIter(toA(w), f, g, info, prob, post, outFile, inflation, uTol);
        [w, uMap, flag, output] = adjopt.trustRegionNewtonCG(@(w) fgWhite(post, R, toA(w)), ...
            @(w, V) R'\post.Uhv(toA(w), R\V), R*(aMap - prior.mu), [], tropts);
        aMap = toA(w);
    else
        tropts.OutputFcn = @(x, f, g, info) printIter(x, f, g, info, prob, post, outFile, inflation, uTol);
        [aMap, uMap, flag, output] = adjopt.trustRegionNewtonCG(post.Ug, post.Uhv, aMap, zeros(m_p, 1), tropts);
    end
    Mmap = prob.f(aMap);
    newInflation = max(1, Mmap/noise.Mnoise);
    infHist(end+1, :) = [inflation, Mmap, newInflation, output.iterations, flag, toc(t0)]; %#ok<SAGROW>
    fprintf(['Round %d: inflation %.3f, %d iterations (flag %d), M(a_MAP) = %.4e = %.3f M_noise, ' ...
        '|a_MAP - a_true|/|a_true| = %.3e\n'], round, inflation, output.iterations, flag, Mmap, ...
        Mmap/noise.Mnoise, norm(aMap - aTrue)/norm(aTrue));
    save(outFile, 'aMap', 'aTrue', 'noise', 'inflation', 'infHist', 'm', 'm_p', 'eta', 'gammaFactor');
    converged = abs(newInflation - inflation) <= 0.05*inflation;
    if converged
        break
    end
    inflation = newInflation;
end
if ~converged
    warning('Inflation not converged after %d rounds', maxRounds);
end
inflation = noise.inflation;   % the inflation a_MAP and H belong to, which the chains load
if noHessian
    save(outFile, 'aMap', 'aTrue', 'noise', 'inflation', 'infHist', 'm', 'm_p', 'eta', 'gammaFactor');
    fprintf('Saved %s (no Hessian)\n', outFile);
    return
end

% Laplace approximation at the MAP, with the final inflation. H is saved before the check, so an
% indefinite H (MAP not converged) does not lose the solves.
H = post.Uhv(aMap, eye(m_p));
H = (H + H')/2;
save(outFile, 'aMap', 'aTrue', 'noise', 'inflation', 'infHist', 'H', 'm', 'm_p', 'eta', 'gammaFactor');
if min(eig(H)) <= 0
    error(['The Hessian at a_MAP is not positive definite (min eigenvalue %.3e): the MAP is not ' ...
        'converged. Rerun with resume = true (and a larger maxIter).'], min(eig(H)));
end
R = prior.R;
lambda = sort(eig((R'\(H - prior.P))/R), 'descend');
fprintf('Eigenvalues of the prior-preconditioned data Hessian: max %.2e, %d of %d above 1\n', ...
    lambda(1), sum(lambda > 1), m_p);
fprintf('Solves: %d forward, %d adjoint, %d Hv\n', prob.nFwd, prob.nAdj, prob.nHv);
save(outFile, 'aMap', 'aTrue', 'noise', 'inflation', 'infHist', 'H', 'lambda', 'm', 'm_p', 'eta', 'gammaFactor');
fprintf('Saved %s\n', outFile);

function stop = printIter(x, f, g, info, prob, post, outFile, inflation, uTol)
    persistent uHist
    if strcmp(info, 'init')
        uHist = [];
    end
    uHist(end+1) = f;
    fprintf('  U = %.6e, |g| = %.3e, M = %.4e, solves %d/%d/%d  %s\n', f, norm(g), ...
        (f + post.prior.logp(x))/post.c, prob.nFwd, prob.nAdj, prob.nHv, info);
    if exist(outFile, 'file')   % checkpoint for resume = true
        aMap = x; %#ok<NASGU>
        save(outFile, 'aMap', 'inflation', '-append');
    end
    stop = uTol > 0 && numel(uHist) > 5 && uHist(end-5) - uHist(end) < uTol;
    if stop
        fprintf('  U decreased by %.2e over 5 iterations: stop\n', uHist(end-5) - uHist(end));
    end
end

function [u, g] = fgWhite(post, R, a)
    [u, g] = post.Ug(a);
    g = R'\g;
end
