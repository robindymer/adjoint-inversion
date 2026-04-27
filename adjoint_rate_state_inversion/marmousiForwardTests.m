function marmousiForwardTests(downSampling, plotFlag, filename, parSetFunc, sourceLocation)

default_arg('downSampling',40);
default_arg('plotFlag',false);
default_arg('filename', 'elAc');
default_arg('parSetFunc', @pars.marmousi2ElasticSources)
default_arg('sourceLocation', 'bottom')

% ------- Setup ----------------------
progressBar = true;
order = 4;
if downSampling > 20
	order = 2;
end

opts = struct;
opts.nS = 1;
opts.T = 3;
invPars = {'source'};

opts.location = sourceLocation;
optObjFunc = @ElasticAcoustic2DMultiblockOpt;

% For saving runs
frameSpacing = 20;
if downSampling == 4
	frameSpacing = 50;
end
fName = [filename 'DS' num2str(downSampling) '.mat'];

saveData.frameSpacing = frameSpacing;
saveData.filename = fName;
%---------------------------------------------

[parset, trueParset] = parSetFunc(opts, downSampling);
disp('Parset created');

% Run with true parameters and save data
trueOptObj = optObjFunc(trueParset, invPars, order);
disp('Optobj created');
disp('Running forward ...');
tic
trueOptObj.runForward(plotFlag, [], saveData, progressBar);
t = toc;
disp(['... done! ' num2str(t) ' s.']);

