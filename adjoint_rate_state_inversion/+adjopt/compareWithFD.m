function compareWithFD(optObjFunc, parSetFunc, invPars, parSetOpts, deltaParExp, order, dim)

default_arg('dim', 2);
default_arg('order', 2);
default_arg('deltaParExp', -6);
default_arg('parSetOpts', struct);
default_arg('invPars', {'rho'});
default_arg('parSetFunc', @pars.rhoInversionParsReceiverLine);
default_arg('optObjFunc', @Elastic2DOpt);

[parset, trueParset] = parSetFunc(parSetOpts);

% Run with true parameters and save data
trueOptObj = optObjFunc(trueParset, invPars, order);
trueOptObj.runForward(false);
trueReceiverData = trueOptObj.forwardReceiverRecordings;
k = trueOptObj.k;
clear trueOptObj;

% Create optimization object
adjOpt = optObjFunc(parset, invPars, order, k);
adjOpt.receiverData = trueReceiverData;

% Compute gradient with adjoint method
grad = adjOpt.computeGradient();

% Compute gradient with brute force FD
dg = 10^deltaParExp;

tic
gradFD = adjOpt.computeGradientFD(dg);

% Compute errors
err = adjOpt.compareGradients(gradFD, grad);
% for j = 1:numel(invPars)
% 	err{j} = abs(gradFD{j}-grad{j});
% end

figure
adjOpt.plotGradient(gradFD)
title('FD')

figure
adjOpt.plotGradient(grad)
title('Adjoint')

figure
adjOpt.plotGradient(err);
title('Absolute error')

toc

