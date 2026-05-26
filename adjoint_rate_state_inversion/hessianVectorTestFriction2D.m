
clear all;
addSubpaths;
order = 8;
inversionPars = {'a'};
initialGuessScaling = 1.1;
T = 1;
m = 51;

parset = @pars.rsFriction2DFractalFaultVerification; delta_p = 1e-9;

parsetOpts.inversionParameters = inversionPars;
parsetOpts.initialGuessScalings = initialGuessScaling;
parsetOpts.T = T;
parsetOpts.m = m;
[parset, parset_true] = parset(parsetOpts);

plotFlag = true;
progressBar = true;

% Get synthetic data
% TODO: Add option for loading data from mat files
data_opt = AntiplaneShear2DRSFrictionOptHessian(parset_true, [], order);
data_opt.runForward(false,[],[], progressBar);
syntheticData = {data_opt.forwardReceiverRecordings,data_opt.forwardTimeIntegrationData};

% [delta_grad_disp, grad_adj_disp, grad_fd_disp] = compare_gradients(parset, order, syntheticData, 'displacement', inversionPars, delta_p, plotFlag, progressBar);
tic
[delta_hessVec, hessVec_vel, hessVec_fd_vel]    = compare_gradients(parset, order, syntheticData, 'velocity', inversionPars, delta_p, plotFlag, progressBar);
toc

function [delta_hessVec, hessVec, hessVec_fd] = compare_gradients(parset, order, syntheticData, misfitType, inversionPars, delta_p, plotFlag, progressBar)
    parset.misfitType = misfitType;
    adj_opt = AntiplaneShear2DRSFrictionOptHessian(parset, inversionPars, order);
    adj_opt.setSyntheticReceiverData(syntheticData{:});
    hessVec = adj_opt.computeHessian(plotFlag, progressBar);

    % Dummy
    delta_hessVec = 0;
    hessVec_fd = 0;

    % grad_adj = adj_opt.computeGradient(plotFlag,progressBar);
    % if adj_opt.tsOpts.forwardMethod.adaptive
    %     adj_opt.tsOpts.forwardMethod.adaptive = false;
    %     adj_opt.tsOpts.k = adj_opt.forwardTimeIntegrationData.k;
    % end
        
    % grad_fd = adj_opt.computeGradientFD(delta_p);
    % delta_grad = adj_opt.compareGradients(grad_adj,grad_fd);     
    % figure();
    % fprintf('%s Misfit = %f\n',misfitType,adj_opt.computeMisfit());
    % for i = 1:numel(inversionPars)
    %     parName = inversionPars{i};
    %     plot(delta_grad.(parName),'LineWidth',2);
    %     hold on;
    %     fprintf('Gradient L2-norm absolute error parameter %s: %e\n', parName, norm(delta_grad.(parName)));
    %     fprintf('Gradient L2-norm relative error parameter %s: %e\n', parName, norm(delta_grad.(parName))/norm(grad_adj.(parName)));
    % end
    % xlabel('Fault index');
    % ylabel('Absolute error')
    % legend(inversionPars{:});
    % title(sprintf('Error using %s misfit', misfitType));
    % hold off;
end