function runFDConvergence(optObjFunc, parSetFunc, invPars, parSetOpts, deltaParExps, order, dim)

default_arg('dim', 2);
default_arg('order', 2);
default_arg('deltaParExps', -1:-1/2:-12);
default_arg('parSetOpts', struct);
default_arg('invPars', {'rho'});
default_arg('parSetFunc', @pars.rhoInversionParsReceiverLine);
default_arg('optObjFunc', @Elastic2DOpt);

linewidth = 2;
legSize = 13;
labelSize = 16;
tickSize = 14;

[parset, trueParset] = parSetFunc(parSetOpts);

% Run with true parameters and save data
trueOptObj = optObjFunc(trueParset, invPars, order);
trueOptObj.runForward(false);
trueReceiverData = trueOptObj.forwardReceiverRecordings;
trueTimeIntegrationData = trueOptObj.forwardTimeIntegrationData;
clear trueOptObj;

% Create optimization object
adjOpt = optObjFunc(parset, invPars, order);
adjOpt.setSyntheticReceiverData(trueReceiverData,trueTimeIntegrationData);

% Compute gradient with adjoint method
grad = adjOpt.computeGradient();

% Compute gradient with brute force FD
DG = 10.^deltaParExps;
err = zeros(size(DG));

% If adaptive timestepping, change to non-adaptive and 
% use timesteps from forward solve
if adjOpt.tsOpts.forwardMethod.adaptive
	adjOpt.tsOpts.forwardMethod.adaptive = false;
	adjOpt.tsOpts.k = adjOpt.forwardTimeIntegrationData.k;
end

tic
str = util.replace_string('','   %d %%',0);
for i = 1:length(DG)
	dg = DG(i);
	gradFD = adjOpt.computeGradientFD(dg);

	difference = adjOpt.compareGradients(grad, gradFD);
	err(i) = adjOpt.gradientNorm(difference)./adjOpt.gradientNorm(grad);

	str = util.replace_string(str,'   %.2f %%',i/length(DG)*100);
end
str = util.replace_string(str,'');
toc

figure
gca.FontSize = tickSize;
loglog(DG, err, '-', 'linewidth', linewidth);
hold on
ref = 2*err(1)/DG(1) * DG;
loglog(DG, ref, '--k', 'linewidth', linewidth);
hold off
xlabel('\Delta parameter', 'fontsize', labelSize)
ylabel('Normalized relative l^2-error', 'fontsize', labelSize)
lgd = legend('Naive gradient - adjoint gradient','1st order reference');
lgd.FontSize = legSize;
lgd.Location = 'north';

