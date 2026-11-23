% Bayesian inversion of a(x) with the noisy high-resolution data (no inverse crime).
%
% The data are the m = 1001 recordings of mat/fractal_fault_m1001_mp1001 with iid Gaussian noise
% N(0, sigma0^2), sigma0 = eta*max|v|, sampled at dt_d = 0.01 s (add_noise.m). The forward model is
% solved on the coarser grid m with a(x) on m_p nodes, so the misfit also contains model error
% (discretization of the PDE and of a; see discretization_error.m). This is accounted for by
% inflating the noise variance (simple approach of the course script):
%
%   sigma^2 = inflation*sigma0^2,   -log L(a) = M(a)/(sigma^2*dt_d)
%
% with the inflation factor chosen by the discrepancy principle in run_highres_map.m.
%
% prob       InvProblem: prob.f(a) misfit, prob.fg(a) misfit and gradient, prob.Hv(a,V)
% noise      struct with eta, sigma0, inflation, sigma (inflated, used by posterior.m), dt_d,
%            N (data used: receivers x samples), Mnoise (misfit of the pure noise with the
%            solver's time quadrature, i.e. the discrepancy target without model error)
% aTrue      true a sampled at the m_p nodes (for plots only; the data were made with a sharp step)
function [prob, noise, aTrue] = setup_highres_problem(m, m_p, eta, inflation)
    bayesDir = fileparts(mfilename('fullpath'));
    codeDir = fullfile(fileparts(bayesDir), 'adjoint_rate_state_inversion');
    addpath(codeDir, fullfile(codeDir, 'optimizer_benchmark'));
    addSubpaths;   % sbplib, with default_arg

    default_arg('m', 151);
    default_arg('m_p', 26);
    default_arg('eta', 0.02);
    default_arg('inflation', 1);

    dataDir = fullfile(codeDir, 'mat', 'fractal_fault_m1001_mp1001');
    parSetOpts.initialGuessScalings = [];
    parSetOpts.initialGuessValues = 0.011;    % not used
    parSetOpts.misfitType = 'velocity';
    parSetOpts.inversionParameters = {'a'};
    parSetOpts.m = m;
    parSetOpts.m_p = m_p;
    parSetOpts.receiverSpacing = 2;           % 88 receivers
    order = 8;
    [parset, trueParset] = pars.rsFriction2DFractalFaultInversion(parSetOpts);

    tmp = AntiplaneShear2DRSFrictionOpt(trueParset, {'a'}, order);
    aTrue = tmp.pars.a;
    clear tmp

    % Noisy and clean high-res recordings at the receivers of parset
    Nz = load(fullfile(dataDir, sprintf('receiverData_noisy_eta%g.mat', eta)));
    Cl = load(fullfile(dataDir, 'receiverData.mat'));
    nRec = numel(parset.receivers.x);
    dNoisy = cell(nRec, 1);
    dClean = cell(nRec, 1);
    used = zeros(nRec, 1);
    for r = 1:nRec
        used(r) = find(cellfun(@(p) isequal(p, parset.receivers.x{r}), Nz.receiverData.positions), 1);
        dNoisy{r} = Nz.receiverData.recordings{used(r)};
        j = find(cellfun(@(p) isequal(p, parset.receivers.x{r}), Cl.receiverData.positions), 1);
        dClean{r} = Cl.receiverData.recordings{j};
    end
    nzInfo = Nz.receiverData.noise;
    clear Nz Cl

    optObj = AntiplaneShear2DRSFrictionOptHessian(parset, {'a'}, order);
    optObj.receiverData = dNoisy;
    prob = InvProblem(optObj, aTrue);

    noise.eta = eta;
    noise.seed = nzInfo.seed;
    noise.sigma0 = nzInfo.sigma;
    noise.inflation = inflation;
    noise.sigma = sqrt(inflation)*nzInfo.sigma;
    noise.dt_d = nzInfo.dt_d;
    noise.N = nRec*numel(nzInfo.t);

    % Misfit of the pure noise with the solver's quadrature (needs the time grid of one forward solve)
    Mtrue = prob.f(aTrue);
    T = optObj.forwardTimeIntegrationData.T;
    noise.Mnoise = 0;
    for r = 1:nRec
        noise.Mnoise = noise.Mnoise + 0.5*((dNoisy{r}(T) - dClean{r}(T)).^2)*optObj.Ht;
    end
    noise.MnoiseFormula = noise.sigma0^2*noise.dt_d*noise.N/2;

    fprintf(['High-res noisy data: m = %d, m_p = %d, sigma0 = %.3e (eta = %g), dt_d = %g, N = %d\n' ...
        '  M_noise = %.4e (sigma0^2 dt_d N/2 = %.4e), M(a_true nodes) = %.4e, inflation = %.3f\n'], ...
        m, m_p, noise.sigma0, eta, noise.dt_d, noise.N, noise.Mnoise, noise.MnoiseFormula, Mtrue, inflation);
end
