
clear all;
addSubpaths;
order = 8;
inversionPars = {'a'};
initialGuessScaling = 1.1;
T = 6;
m = 51;
m_p = 51;

parset = @pars.rsFriction2DFractalFaultVerification; 
delta_p = 1e-9; % For FD (also second pertubation?)

parsetOpts.inversionParameters = inversionPars;
parsetOpts.initialGuessScalings = initialGuessScaling;
parsetOpts.T = T;
parsetOpts.m = m;
parsetOpts.m_p = m_p;
parsetOpts.k_idx = 10;
parsetOpts.delta_p = delta_p;
[parset, parset_true] = parset(parsetOpts);

plotFlag = false;
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
    hessVec_fd = adj_opt.computeHessianVectorFD(1e-5);
    toc

    % Trickery!
    hessVec_struct = struct();
    if ~isstruct(hessVec)
        hessVec_struct.a = hessVec;
        hessVec = hessVec_struct;
    end

    delta_hessVec = adj_opt.compareGradients(hessVec,hessVec_fd);     
    figure();
    fprintf('%s Misfit = %f\n',misfitType,adj_opt.computeMisfit());
    for i = 1:numel(inversionPars)
        parName = inversionPars{i};
        plot(delta_hessVec.(parName),'LineWidth',2);
        hold on;
        % fprintf('Norm analytical: %e\n', norm(hessVec.a));
        % fprintf('Norm FD: %e\n', norm(hessVec_fd.a));
        % disp(hessVec.a(1:11));
        % disp(hessVec_fd.a(1:11));
        % disp(hessVec.a(1:11) ./ hessVec_fd.a(1:11));
        fprintf('Hessian vector L2-norm absolute error parameter %s: %e\n', parName, norm(delta_hessVec.(parName)));
        fprintf('Hessian vector L2-norm relative error parameter %s: %e\n', parName, norm(delta_hessVec.(parName))/norm(hessVec.(parName)));
    end
    xlabel('Fault index');
    ylabel('Absolute error')
    legend(inversionPars{:});
    title(sprintf('Error using %s misfit', misfitType));
    hold off;

    % % Test different FD step sizes
    % fd_steps = logspace(-5, -10, 6); 

    % for step = fd_steps
    %     tic
    %     hessVec_fd = adj_opt.computeHessianVectorFD(step);
    %     toc
        
    %     delta_hessVec = adj_opt.compareGradients(hessVec, hessVec_fd); 
        
    %     fprintf('\n--- FD Step Size: %e ---\n', step);
    %     for i = 1:numel(inversionPars)
    %         parName = inversionPars{i};
    %         rel_error = norm(delta_hessVec.(parName)) / norm(hessVec.(parName));
    %         fprintf('Parameter %s relative error: %e\n', parName, rel_error);
    %     end
    % end
end