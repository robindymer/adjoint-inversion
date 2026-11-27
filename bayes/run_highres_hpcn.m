% H-pCN for the high-res data: pCN with the Laplace approximation N(a_MAP, H^-1) from
% run_highres_map.m as reference, and the inflated noise found there.
%
% Starts from a draw of the Laplace approximation with doubled standard deviation (overdispersed,
% script Sec. 4.6). s = 1 gives independent proposals and is right if the posterior is close to
% the Laplace approximation, i.e. if std(chain.phi) <~ 1; otherwise lower s (e.g. 0.5).
% seed, nIter, s, m, m_p can be set before calling: matlab -batch "seed=2; run_highres_hpcn"
% The chain is saved as it runs, in the format of run_laplace_pcn.m, so python/diagnostic.py works:
%   python diagnostic.py "../results/hpcn_highres_m151_mp26_eta0.02_s1_seed*.mat"

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
if ~exist('m', 'var'), m = 151; end
if ~exist('m_p', 'var'), m_p = 26; end
eta = 0.02;
if ~exist('nIter', 'var'), nIter = 1000; end
if ~exist('s', 'var'), s = 1; end
if ~exist('seed', 'var'), seed = 1; end
if ~exist('gammaFactor', 'var'), gammaFactor = 1; end   % prior, see highres_prior.m
if ~exist('nAdapt', 'var'), nAdapt = 0; end   % > 0: adapt s towards 25 % acceptance (burn-in)

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
opts.saveFile = fullfile(resDir, sprintf('hpcn_highres_m%d_mp%d_eta%g%s_s%g_seed%d.mat', m, m_p, eta, tag, s, seed));
start = 'laplace';
chain = pcn(post, a0, nIter, s, opts);
save(opts.saveFile, 'chain', 'aTrue', 'aMap', 'H', 'noise', 'm', 'm_p', 's', 'start', 'seed');

fprintf('s = %g, acceptance rate = %.2f, std(phi) = %.2f, %d forward solves in %.0f s\n', ...
    s, mean(chain.accepted), std(chain.phi), chain.nFwd, chain.time);
