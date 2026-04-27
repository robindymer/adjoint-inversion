function marmousiTests(downSampling, plotFlag)

default_arg('downSampling',40);
default_arg('plotFlag',false);

% ------- Setup ----------------------
order = 4;
if downSampling > 20
	order = 2;
end

opts = struct;
opts.nRecEl = 1;
opts.nRecAc = 0;
opts.nS = 2;
opts.T = 4;
invPars = {'source'};

opts.nRec = opts.nRecEl;

parSetFuncGood = @pars.marmousi2;
optObjFuncGood = @ElasticAcoustic2DMultiblockOpt;

parSetFuncBad = @pars.marmousi2elastic;
optObjFuncBad = @Elastic2DMultiblockOpt;

% For saving runs
frameSpacing = 20;
fNameAdjGood = sprintf('goodAdjDS%d.mat',downSampling);
fNameFwdGood = sprintf('goodFwdDS%d.mat',downSampling);
fNameAdjBad = sprintf('badAdjDS%d.mat',downSampling);

saveDataGood.frameSpacing = frameSpacing;
saveDataBad.frameSpacing = frameSpacing;

saveDataGood.filename = fNameAdjGood;
saveDataBad.filename = fNameAdjBad;

saveDataGoodForward.frameSpacing = frameSpacing;
saveDataGoodForward.filename = fNameFwdGood;
%---------------------------------------------



[parset, trueParset] = parSetFuncGood(opts, downSampling);
disp('Parset created');

% Run with true parameters and save data
trueOptObj = optObjFuncGood(trueParset, invPars, order);
disp('Optobj created');
disp('Running true forward ...');
tic
trueOptObj.runForward(plotFlag);
t = toc;
disp(['... done! ' num2str(t) ' s.']);
trueReceiverData = trueOptObj.forwardReceiverRecordings;
k = trueOptObj.k;
trueParset = trueOptObj.pars;
clear trueOptObj;

% ----- Create good optimization object -----------
optGood = optObjFuncGood(parset, invPars, order, k);
optGood.receiverData = trueReceiverData;

disp('Running forward ...');
tic
optGood.runForward(plotFlag, [], saveDataGoodForward);
t = toc;
disp(['... done! ' num2str(t) ' s.']);
disp('Updating adjoint discr ...');
tic
optGood.updateAdjointDiscr();
t = toc;
disp(['... done! ' num2str(t) ' s.']);
disp('Running adjoint ...');
tic
optGood.runAdjoint(plotFlag, [], saveDataGood);
t = toc;
disp(['... done! ' num2str(t) ' s.']);
forwardReceiverRecordings = optGood.forwardReceiverRecordings;
clear optGood
%-------------------------------------------------

% ----- Create bad optimization object -----------
parsetBad = parSetFuncBad(opts, downSampling);
optBad = optObjFuncBad(parsetBad, invPars, order, k);
optBad.receiverData = trueReceiverData.elastic;

% Take forward receiver recordings from elastic-acoustic simulation
optBad.forwardReceiverRecordings = forwardReceiverRecordings.elastic;

disp('Updating bad adjoint discr ...');
tic
optBad.updateAdjointDiscr();
t = toc;
disp(['... done! ' num2str(t) ' s.']);
disp('Running bad adjoint ...');
tic
optBad.runAdjoint(plotFlag, [], saveDataBad);
t = toc;
disp(['... done! ' num2str(t) ' s.']);
clear optBad;
%-------------------------------------------------
