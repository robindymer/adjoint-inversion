function [optObj, x0, truePars, xtrue] = setup_problem(m, m_p, val)
    % Mirrors adjopt.optimize_general up to the optimizer call, using the Hessian class.
    parSetOpts.initialGuessScalings = [];
    parSetOpts.initialGuessValues = val;
    parSetOpts.misfitType = 'velocity';
    parSetOpts.inversionParameters = {'a'};
    parSetOpts.m = m;
    parSetOpts.m_p = m_p;
    order = 8;
    [parset, trueParset] = pars.rsFriction2DFractalFaultInversion(parSetOpts);
    trueOptObj = AntiplaneShear2DRSFrictionOptHessian(trueParset, {'a'}, order);
    trueOptObj.runForward(false);
    trueReceiverData = trueOptObj.forwardReceiverRecordings;
    trueTimeIntegrationData = trueOptObj.forwardTimeIntegrationData;
    k = trueOptObj.k;
    truePars = trueOptObj.pars;
    clear trueOptObj;
    optObj = AntiplaneShear2DRSFrictionOptHessian(parset, {'a'}, order, k);
    optObj.setSyntheticReceiverData(trueReceiverData, trueTimeIntegrationData);
    x0 = optObj.parametersNativeToVector(optObj.pars);
    xtrue = truePars.a;
end
