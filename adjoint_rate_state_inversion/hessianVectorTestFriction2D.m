clear all;
close all;
addSubpaths;

QUICK_TEST = true;
DO_VALIDATION_SWEEP = false;
DO_EPSILON_CONSISTENCY_CHECK = true;

order = 8;
inversionPars = {'a'};
initialGuessScaling = 1.1;

parsetFun = @pars.rsFriction2DFractalFaultVerification;
parsetOpts.inversionParameters = inversionPars;
parsetOpts.initialGuessScalings = initialGuessScaling;
parsetOpts.misfitType = 'velocity';

if QUICK_TEST
    % Lower-cost setup for fast sanity checks.
    order = 4;
    parsetOpts.m = 41;
    % Keep inversion grid equal to simulation grid to avoid GlueOps projection path.
    parsetOpts.m_p = parsetOpts.m;
    parsetOpts.k = 0.02;
    parsetOpts.receiverSpacing = 2;
end

[parset, trueParset] = parsetFun(parsetOpts);

% Generate synthetic data at true model.
trueOpt = AntiplaneShear2DRSFrictionOpt(trueParset, inversionPars, order);
trueOpt.runForward(false, [], [], false);
syntheticData = {trueOpt.forwardReceiverRecordings, trueOpt.forwardTimeIntegrationData};

% Build inversion object.
opt = AntiplaneShear2DRSFrictionOpt(parset, inversionPars, order);
opt.setSyntheticReceiverData(syntheticData{:});

% Direction for Hessian-vector product.
direction = struct;
direction.a = ones(size(opt.pars.a));
direction.a = direction.a / norm(direction.a);

% Compute Hessian-vector product.
epsHv = 1e-6;
hv = opt.computeHessianVector(direction, epsHv, false, false);

fprintf('Computed Hessian-vector product with epsilon = %.2e\n', epsHv);
fprintf('||H*dp||_2 for a: %e\n', norm(hv.a));

if DO_EPSILON_CONSISTENCY_CHECK
    if QUICK_TEST
        epsCheck = [1e-5, 1e-6, 1e-7];
    else
        epsCheck = [1e-4, 1e-5, 1e-6];
    end

    hvVecs = cell(numel(epsCheck),1);
    hvNorms = zeros(numel(epsCheck),1);
    for i = 1:numel(epsCheck)
        hv_i = opt.computeHessianVector(direction, epsCheck(i), false, false);
        hvVecs{i} = opt.gradientNativeToVector(hv_i);
        hvNorms(i) = norm(hvVecs{i});
    end

    fprintf('\nEpsilon consistency check for H*dp:\n');
    for i = 1:numel(epsCheck)
        fprintf('  eps = %.1e, ||H*dp||_2 = %e\n', epsCheck(i), hvNorms(i));
    end

    for i = 1:numel(epsCheck)-1
        relDiff = norm(hvVecs{i} - hvVecs{i+1})/(norm(hvVecs{i+1}) + eps);
        fprintf('  rel diff [%.1e vs %.1e] = %e\n', epsCheck(i), epsCheck(i+1), relDiff);
    end
end

if DO_VALIDATION_SWEEP
    % Keep this short in quick mode.
    if QUICK_TEST
        epsList = 10.^(-4:-1:-6);
    else
        epsList = 10.^(-3:-1:-7);
    end
    [epsListOut, relErr] = opt.validateHessianVectorFD(direction, epsList, false);

    fprintf('\nValidation against a tighter finite-difference reference:\n');
    for i = 1:numel(epsListOut)
        fprintf('  eps = %.1e, relative error = %e\n', epsListOut(i), relErr(i));
    end

    figure;
    loglog(epsListOut, relErr, '-o', 'LineWidth', 2);
    grid on;
    xlabel('epsilon');
    ylabel('relative error in H*dp');
    title('Hessian-vector finite-difference validation');
end
