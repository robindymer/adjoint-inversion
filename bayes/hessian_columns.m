% Columns cols of the Hessian of U at the a_MAP of a MAP file (run_highres_map.m), so that the m_p
% Hessian-vector products can run in several processes; combine with combine_hessian.m.
% With dirFile (a .mat with directions D, m_p x r), H*D(:, cols) instead of columns of the identity,
% for the low-rank Hessian of lowrank_hessian.m.
% Set mapFile, cols, gammaFactor (default 1) before calling:
%   matlab -batch "mapFile='results/...mat'; cols=1:63; hessian_columns"
% Saved to <mapFile without .mat>_Hcols<first>-<last>.mat.

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
if ~exist('gammaFactor', 'var'), gammaFactor = 1; end
S = load(mapFile, 'aMap', 'inflation', 'm', 'm_p', 'eta');
[prob, noise] = setup_highres_problem(S.m, S.m_p, S.eta, S.inflation);
if ismember('theta', who('-file', mapFile))   % GSBL MAP (gsbl_map.m): a | theta
    T = load(mapFile, 'theta', 'gamma0');
    prior = smoothness_prior(S.m_p, sqrt(T.theta), T.gamma0);
else
    prior = highres_prior(S.m_p, gammaFactor);
end
post = posterior(prob, noise, prior);
if exist('dirFile', 'var')
    E = getfield(load(dirFile, 'D'), 'D'); %#ok<GFLD>
    outBase = dirFile(1:end-4);
else
    E = eye(S.m_p);
    outBase = mapFile(1:end-4);
end
post.Ug(S.aMap);   % forward and adjoint solve at a_MAP, needed by the Hessian-vector products
t0 = tic;
Hcols = post.Uhv(S.aMap, E(:, cols));
fprintf('%d Hessian-vector products in %.0f s\n', numel(cols), toc(t0));
save(sprintf('%s_Hcols%d-%d.mat', outBase, cols(1), cols(end)), 'Hcols', 'cols');
