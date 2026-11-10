% Adds noise to high-res data
% Assume data_noisy = data + N(0, sigma^2), iid at uniform sample times t_n
%
% The high-res recordings are stored as spline function handles
% d_r(t) = ppval(pp_r, t) (see setSyntheticReceiverData). We therefore
%   1) sample each clean spline at t_n = 0:dt_d:T,
%   2) add iid noise with sigma = eta * max|v| (max over all receivers),
%      where eta is typically a few percent
%   3) fit a new spline through the noisy samples and store it in the same
%      receiverData format, so optimize_general loads it unchanged.
% The clean and noisy samples are also stored as plain arrays,
% can be useful when splines is not needed (or used by Python)

dataDir = '../adjoint_rate_state_inversion/mat/fractal_fault_m1001_mp1001';
eta     = 0.02;   % relative noise level
dt_d    = 0.01;   % data sampling interval [s]: 100 Hz seismometer. Resolves
                  % the signal (99.9% of energy below 21.5 Hz, m=251 data).
                  % Only sigma^2*dt_d enters the likelihood.
seed    = 1;

S = load(fullfile(dataDir, 'receiverData.mat'));
receiverData = S.receiverData;
nRec = numel(receiverData.recordings);

% Recover the time span from the stored spline (its breakpoints are the
% solver's output times)
pp0 = functions(receiverData.recordings{1}).workspace{1}.data_pp;
t = pp0.breaks(1):dt_d:pp0.breaks(end);
nT = numel(t);

% 1) Sample clean data
d_clean = zeros(nRec, nT);
for r = 1:nRec
    d_clean(r,:) = receiverData.recordings{r}(t);
end

% 2) Add iid Gaussian noise
sigma = eta * max(abs(d_clean(:)));
rng(seed);
d_noisy = d_clean + sigma * randn(nRec, nT);

% 3) Back to spline handles
for r = 1:nRec
    pp = spline(t, d_noisy(r,:));
    receiverData.recordings{r} = @(tq) ppval(pp, tq);
end

% Plain-array copies and noise metadata (ignored by optimize_general)
receiverData.noise = struct('sigma', sigma, 'eta', eta, 'dt_d', dt_d, ...
    'seed', seed, 't', t, 'd_clean', d_clean, 'd_noisy', d_noisy);

outFile = fullfile(dataDir, sprintf('receiverData_noisy_eta%g.mat', eta));
save(outFile, 'receiverData', '-v7.3', '-nocompression');
fprintf('sigma = %.3e (eta = %g, max|v| = %.3e), %d receivers x %d samples\n', ...
    sigma, eta, max(abs(d_clean(:))), nRec, nT);
fprintf('Saved %s\n', outFile);
