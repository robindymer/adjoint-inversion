% Inverse-crime problem with noisy data for Bayesian inversion of a(x).
%
% As optimizer_benchmark/setup_problem.m, but iid Gaussian noise N(0, sigma^2) with
% sigma = eta*max|v| is added to the synthetic recordings at the solver's output times,
% so dt_d is the solver's time step. The likelihood is -log L(a) = M(a)/(sigma^2*dt_d).
%
% prob    InvProblem: prob.f(a) misfit, prob.fg(a) misfit and gradient, prob.Hv(a,V)
% noise   struct with eta, seed, sigma, dt_d, N (number of data = receivers x times)
% aTrue   true a at the m_p parameter nodes
function [prob, noise, aTrue] = setup_bayes_problem(m, m_p, eta, seed)
    bayesDir = fileparts(mfilename('fullpath'));
    codeDir = fullfile(fileparts(bayesDir), 'adjoint_rate_state_inversion');
    addpath(codeDir, fullfile(codeDir, 'optimizer_benchmark'));
    addSubpaths;   % sbplib, with default_arg

    default_arg('m', 36);
    default_arg('m_p', 26);
    default_arg('eta', 0.02);
    default_arg('seed', 1);

    parSetOpts.initialGuessScalings = [];
    parSetOpts.initialGuessValues = 0.011;    % not used by the samplers
    parSetOpts.misfitType = 'velocity';
    parSetOpts.inversionParameters = {'a'};
    parSetOpts.m = m;
    parSetOpts.m_p = m_p;
    parSetOpts.receiverSpacing = 2;           % 88 receivers
    order = 8;
    [parset, trueParset] = pars.rsFriction2DFractalFaultInversion(parSetOpts);

    % Clean synthetic data
    trueOptObj = AntiplaneShear2DRSFrictionOptHessian(trueParset, {'a'}, order);
    trueOptObj.runForward(false);
    data = trueOptObj.forwardReceiverRecordings;
    timeData = trueOptObj.forwardTimeIntegrationData;
    k = trueOptObj.k;
    aTrue = trueOptObj.pars.a;
    clear trueOptObj;

    % Add noise, with its own stream so the global RNG is untouched
    optObj = AntiplaneShear2DRSFrictionOptHessian(parset, {'a'}, order, k);
    tData = optObj.removeStagedData(timeData.T, timeData.nStages);
    noise.eta = eta;
    noise.seed = seed;
    noise.sigma = eta*max(abs(data(:)));
    noise.dt_d = mean(diff(tData));
    noise.N = size(data, 1)*numel(tData);
    rs = RandStream('mt19937ar', 'Seed', seed);
    data = data + noise.sigma*randn(rs, size(data));
    optObj.setSyntheticReceiverData(data, timeData);

    prob = InvProblem(optObj, aTrue);
    fprintf('Noisy inverse-crime data: m = %d, sigma = %.3e (eta = %g), dt_d = %.3e, N = %d\n', ...
        m, noise.sigma, eta, noise.dt_d, noise.N);
end
