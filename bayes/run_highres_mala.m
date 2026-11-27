% H-MALA for the high-res data: MALA preconditioned with the Laplace approximation N(a_MAP, H^-1)
% from run_highres_map.m, with the inflated noise found there (cf. run_laplace_mala.m).
% Starts from an overdispersed Laplace draw. Default nIter = 500 (1 forward + 1 adjoint solve per
% iteration) matches the budget of the 1000 H-pCN iterations.
% seed, nIter, nAdapt, h0, m, m_p can be set before calling: matlab -batch "seed=2; run_highres_mala"

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
if ~exist('m', 'var'), m = 151; end
if ~exist('m_p', 'var'), m_p = 26; end
eta = 0.02;
if ~exist('nIter', 'var'), nIter = 500; end
if ~exist('nAdapt', 'var'), nAdapt = 50; end
if ~exist('h0', 'var'), h0 = 0.5; end    % whitened coordinates; adapted
if ~exist('seed', 'var'), seed = 1; end
if ~exist('gammaFactor', 'var'), gammaFactor = 1; end   % prior, see highres_prior.m

resDir = fullfile(bayesDir, 'results');
[prior, tag] = highres_prior(m_p, gammaFactor);
L = load(fullfile(resDir, sprintf('highres_map_m%d_mp%d_eta%g%s.mat', m, m_p, eta, tag)), ...
    'aMap', 'H', 'inflation');
aMap = L.aMap; H = L.H;
[prob, noise, aTrue] = setup_highres_problem(m, m_p, eta, L.inflation);
post = posterior(prob, noise, prior);

RH = chol(H);
rs = RandStream('mt19937ar', 'Seed', seed);
a0 = aMap + 2*(RH\randn(rs, m_p, 1));
while ~isfinite(post.U(a0))
    a0 = aMap + 2*(RH\randn(rs, m_p, 1));
end

opts.ref = struct('mu', aMap, 'R', RH);
opts.nAdapt = nAdapt;
opts.rs = rs;
opts.saveEvery = 10;
opts.saveFile = fullfile(resDir, sprintf('mala_highres_m%d_mp%d_eta%g%s_seed%d.mat', m, m_p, eta, tag, seed));
start = 'laplace';
chain = mala(post, a0, nIter, h0, opts);
save(opts.saveFile, 'chain', 'aTrue', 'aMap', 'H', 'noise', 'm', 'm_p', 'start', 'seed');

kept = nAdapt+1:nIter;
fprintf('Final h = %.2e, acceptance rate after burn-in = %.2f, %d forward + %d adjoint solves in %.0f s\n', ...
    chain.h(end), mean(chain.accepted(kept)), chain.nFwd, chain.nAdj, chain.time);
