% Low-rank Hessian at a_MAP by a single-pass randomized eigendecomposition of the
% prior-preconditioned data Hessian (Halko, Martinsson & Tropp 2011, Sec. 5.5):
%
%   H = P + R' (R'^-1 H_data R^-1) R ~ P + R' V Lambda V' R,   P = R'R the prior precision
%
% Directions with eigenvalue << 1 are informed by the prior only and are dropped. Two steps:
%   step = 'directions'  draw Omega (m_p x r Gaussian) and save D = R\Omega to <mapFile>_dirs.mat;
%                        then run hessian_columns with dirFile on chunks of 1:r (several processes)
%   step = 'combine'     Y = R'^-1 H_data D, Q = orth(Y), B = (Q'Y)(Q'Omega)^-1, B = V Lambda V',
%                        and save H, lambda (all r eigenvalues) to the MAP file
% Set mapFile, step, r (default 100) and gammaFactor (default 1) before calling.

bayesDir = fileparts(mfilename('fullpath'));
addpath(bayesDir);
if ~exist('gammaFactor', 'var'), gammaFactor = 1; end
if ~exist('r', 'var'), r = 100; end
S = load(mapFile);
if ismember('theta', who('-file', mapFile))   % GSBL MAP (gsbl_map.m): a | theta
    T = load(mapFile, 'theta', 'gamma0');
    prior = smoothness_prior(S.m_p, sqrt(T.theta), T.gamma0);
else
    prior = highres_prior(S.m_p, gammaFactor);
end
R = prior.R;
dirFile = sprintf('%s_dirs.mat', mapFile(1:end-4));
switch step
    case 'directions'
        Omega = randn(RandStream('mt19937ar', 'Seed', 0), S.m_p, r);
        D = R\Omega;
        save(dirFile, 'Omega', 'D');
        fprintf('Saved %d directions to %s\n', r, dirFile);
    case 'combine'
        W = load(dirFile);
        r = size(W.D, 2);
        HD = nan(S.m_p, r);
        files = dir(sprintf('%s_Hcols*.mat', dirFile(1:end-4)));
        for k = 1:numel(files)
            C = load(fullfile(files(k).folder, files(k).name));
            HD(:, C.cols) = C.Hcols;
        end
        if any(isnan(HD(:)))
            error('Hessian-vector products missing');
        end
        Y = R'\(HD - prior.P*W.D);   % prior-preconditioned data Hessian times Omega
        [Q, ~] = qr(Y, 0);
        B = (Q'*Y)/(Q'*W.Omega);
        B = (B + B')/2;
        [U, Lam] = eig(B);
        [lambda, idx] = sort(diag(Lam), 'descend');
        V = Q*U(:, idx);
        keep = lambda > 0;
        H = prior.P + R'*(V(:, keep)*diag(lambda(keep))*V(:, keep)')*R;
        H = (H + H')/2;
        fprintf(['Low rank r = %d: lambda max %.2e, %d above 1, smallest kept %.2e ' ...
            '(should be << 1); min eig(H) = %.3e\n'], r, lambda(1), sum(lambda > 1), ...
            min(lambda(keep)), min(eig(H)));
        S.H = H; S.lambda = lambda; S.lowRank = r;
        save(mapFile, '-struct', 'S');
        fprintf('Saved H to %s\n', mapFile);
end
