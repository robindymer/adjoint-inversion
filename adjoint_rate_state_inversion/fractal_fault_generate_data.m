%clear all;
addSubpaths;
order = 8;

parset = @pars.rsFriction2DFractalFaultForward;
paropts.misfitType = 'velocity';
[~, parset_true] = parset(paropts);

plotFlag = false;
progressBar = true;
saveOpts.frameSpacing = 20;
%saveOpts = [];
saveDir = sprintf('mat/fractal_fault_m%d_mp%d',parset_true.m,parset_true.m_p);
success = mkdir(saveDir);
if ~success
    error('Failed to create directory %s',saveDir);
end
saveOpts.saveDir = saveDir;

% Create and save synthetic receiver data
data_opt = AntiplaneShear2DRSFrictionOpt(parset_true, [], order);
data_opt.runForward(plotFlag,parset_true.T,saveOpts,progressBar)
data_opt.setSyntheticReceiverData(data_opt.forwardReceiverRecordings,data_opt.forwardTimeIntegrationData);

% 
receiverData.misfitType = paropts.misfitType;
receiverData.recordings = data_opt.receiverData;
receiverData.positions = parset_true.receivers.x;
receiverData.blockIds = parset_true.receivers.blockIds;
filename = sprintf('%s/receiverData.mat',saveDir);
save(filename, 'receiverData', '-v7.3','-nocompression');
