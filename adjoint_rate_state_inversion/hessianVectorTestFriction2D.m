
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
parsetOpts.delta_p = delta_p;
[parset, parset_true] = parset(parsetOpts);

plotFlag = true;
progressBar = true;

% Get synthetic data
% TODO: Add option for loading data from mat files
data_opt = AntiplaneShear2DRSFrictionOptHessian(parset_true, [], order);
data_opt.runForward(false,[],[], progressBar);
syntheticData = {data_opt.forwardReceiverRecordings,data_opt.forwardTimeIntegrationData};

% [delta_grad_disp, grad_adj_disp, grad_fd_disp] = compare_gradients(parset, order, syntheticData, 'displacement', inversionPars, delta_p, plotFlag, progressBar);
[delta_hessVec, hessVec_vel, hessVec_fd_vel]    = compare_gradients(parset, order, syntheticData, 'velocity', inversionPars, delta_p, plotFlag, progressBar);

function [delta_hessVec, hessVec, hessVec_fd] = compare_gradients(parset, order, syntheticData, misfitType, inversionPars, delta_p, plotFlag, progressBar)
    parset.misfitType = misfitType;
    adj_opt = AntiplaneShear2DRSFrictionOptHessian(parset, inversionPars, order);
    adj_opt.setSyntheticReceiverData(syntheticData{:});

    tic
    hessVec = adj_opt.computeHessianVector(plotFlag, progressBar);
    toc
    
    tic
    hessVec_fd = adj_opt.computeHessianVectorFD(delta_p);
    toc

    % Trickery!
    hessVec_struct = struct();
    hessVec_struct.a = hessVec;
    hessVec = hessVec_struct;

    delta_hessVec = adj_opt.compareGradients(hessVec,hessVec_fd);     
    figure();
    fprintf('%s Misfit = %f\n',misfitType,adj_opt.computeMisfit());
    for i = 1:numel(inversionPars)
        parName = inversionPars{i};
        plot(delta_hessVec.(parName),'LineWidth',2);
        hold on;
        fprintf('Hessian vector L2-norm absolute error parameter %s: %e\n', parName, norm(delta_hessVec.(parName)));
        fprintf('Hessian vector L2-norm relative error parameter %s: %e\n', parName, norm(delta_hessVec.(parName))/norm(hessVec.(parName)));
    end
    xlabel('Fault index');
    ylabel('Absolute error')
    legend(inversionPars{:});
    title(sprintf('Error using %s misfit', misfitType));
    hold off;
end