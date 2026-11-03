%clear all;
addSubpaths;
order = 8;

parset = @pars.rsFriction2DFractalFaultForward;
paropts.misfitType = 'velocity';
% Grid and parameter grid size, m = m_p = 1001 by default. Can be set before running the script,
% e.g. matlab -batch "data_m=251; fractal_fault_generate_data". Smaller grids save no frames.
if exist('data_m', 'var')
    paropts.m = data_m;
    paropts.m_p = data_m;
end
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
if exist('data_m', 'var')
    saveOpts = [];
end

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
