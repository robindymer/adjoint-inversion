% Discretization (model) error of the coarse forward model against the high-resolution data
% (m = 1001, mat/fractal_fault_m1001_mp1001), to choose the smallest m for the high-res run.
%
% For each m, the misfit M(a_true) against the CLEAN high-res recordings is pure model error
% (discretization of the PDE + representation of a on the m_p = 26 node grid). It is compared
% with the noise level of the noisy data (eta = 0.02, dt_d = 0.01 s, add_noise.m),
%   M_noise = 1/2 sum_r int (d_noisy - d_clean)^2 dt  ~  sigma^2 dt_d N / 2.
% Model error is negligible if M_model << M_noise.
% Set ms and m_p before calling, e.g. matlab -batch "ms=[36 51]; m_p=101; discretization_error"

bayesDir = fileparts(mfilename('fullpath'));
codeDir = fullfile(fileparts(bayesDir), 'adjoint_rate_state_inversion');
oldDir = cd(codeDir);
restoreDir = onCleanup(@() cd(oldDir));
addSubpaths;
if ~exist('ms', 'var'), ms = [36 51 76 101]; end
if ~exist('m_p', 'var'), m_p = 26; end
order = 8;
dataDir = fullfile(codeDir, 'mat', 'fractal_fault_m1001_mp1001');

S = load(fullfile(dataDir, 'receiverData.mat'));
clean = S.receiverData;
Nz = load(fullfile(dataDir, 'receiverData_noisy_eta0.02.mat'));
nz = Nz.receiverData.noise;
% Noise level over the receivers used in the inversion (receiverSpacing 2), not all in the file
opts0 = struct('m', 36, 'm_p', 26, 'misfitType', 'velocity', 'inversionParameters', {{'a'}}, ...
    'initialGuessScalings', [], 'initialGuessValues', 0.011, 'receiverSpacing', 2);
parset0 = pars.rsFriction2DFractalFaultInversion(opts0);
used = cellfun(@(x) find(cellfun(@(p) isequal(p, x), Nz.receiverData.positions), 1), parset0.receivers.x);
nRec = numel(used);
Nused = nRec*numel(nz.t);
Mnoise = 0.5*sum(trapz(nz.t, (nz.d_noisy(used,:) - nz.d_clean(used,:)).^2, 2));
fprintf('Noise: sigma = %.3e, dt_d = %g, %d of %d receivers, N = %d, sigma^2 dt_d N/2 = %.4e, M_noise (trapz) = %.4e\n', ...
    nz.sigma, nz.dt_d, nRec, size(nz.d_clean, 1), Nused, nz.sigma^2*nz.dt_d*Nused/2, Mnoise);

res = zeros(numel(ms), 3);
for i = 1:numel(ms)
    m = ms(i);
    opts = struct('m', m, 'm_p', m_p, 'misfitType', 'velocity', 'inversionParameters', {{'a'}}, ...
        'initialGuessScalings', [], 'initialGuessValues', 0.011, 'receiverSpacing', 2);
    [parset, trueParset] = pars.rsFriction2DFractalFaultInversion(opts);
    tmp = AntiplaneShear2DRSFrictionOpt(trueParset, {'a'}, order);
    aTrue = tmp.pars.a;
    clear tmp
    optObj = AntiplaneShear2DRSFrictionOpt(parset, {'a'}, order);
    rec = cell(numel(parset.receivers.x), 1);
    for r = 1:numel(parset.receivers.x)
        j = find(cellfun(@(p) isequal(p, parset.receivers.x{r}), clean.positions), 1);
        rec{r} = clean.recordings{j};
    end
    optObj.receiverData = rec;
    optObj.setParameters(optObj.parametersVectorToNative(aTrue));
    optObj.updateForwardDiscr();
    t0 = tic;
    optObj.runForward(false);
    tFwd = toc(t0);
    M = optObj.computeMisfit();
    res(i,:) = [M, M/Mnoise, tFwd];
    fprintf('m = %4d: M_model = %.4e, M_model/M_noise = %.3f, forward solve %.0f s\n', ...
        m, M, M/Mnoise, tFwd);
end
save(fullfile(bayesDir, 'results', sprintf('discretization_error_mp%d.mat', m_p)), 'ms', 'res', 'Mnoise', 'nz');
