% Combine the Hessian columns of hessian_columns.m into H in the MAP file, as run_highres_map.m
% does after the MAP: H, lambda (eigenvalues of the prior-preconditioned data Hessian), and the final
% misfit in infHist. Set mapFile (and gammaFactor) before calling.

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
if ~exist('gammaFactor', 'var'), gammaFactor = 1; end
S = load(mapFile);
files = dir(sprintf('%s_Hcols*.mat', mapFile(1:end-4)));
H = nan(S.m_p);
for k = 1:numel(files)
    C = load(fullfile(files(k).folder, files(k).name));
    H(:, C.cols) = C.Hcols;
end
if any(isnan(H(:)))
    error('Hessian columns missing');
end
H = (H + H')/2;
prior = highres_prior(S.m_p, gammaFactor);
R = prior.R;
lambda = sort(eig((R'\(H - prior.P))/R), 'descend');
fprintf('min eig(H) = %.3e; prior-preconditioned data Hessian: max %.2e, %d of %d above 1\n', ...
    min(eig(H)), lambda(1), sum(lambda > 1), S.m_p);
S.H = H; S.lambda = lambda;
save(mapFile, '-struct', 'S');
fprintf('Saved H to %s\n', mapFile);
