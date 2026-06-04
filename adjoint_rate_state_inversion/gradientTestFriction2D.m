
clear all;
addSubpaths;
order = 8;
inversionPars = {'a'};
initialGuessScaling = 1.1;
T = 5;
m = 51;
m_p = 11;

parset = @pars.rsFriction2DFractalFaultVerification; delta_p = 1e-9;
parsetOpts.delta_p = delta_p;
parsetOpts.T = T;
parsetOpts.m = m;
parsetOpts.m_p = m_p;
parsetOpts.k_idx = 10;

parsetOpts.inversionParameters = inversionPars;
parsetOpts.initialGuessScalings = initialGuessScaling;
[parset, parset_true] = parset(parsetOpts);

plotFlag = false;
progressBar = true;

% Get synthetic data
% TODO: Add option for loading data from mat files
data_opt = AntiplaneShear2DRSFrictionOpt(parset_true, [], order);
data_opt.runForward(false,[],[], progressBar);
syntheticData = {data_opt.forwardReceiverRecordings,data_opt.forwardTimeIntegrationData};

% [delta_grad_disp, grad_adj_disp, grad_fd_disp] = compare_gradients(parset, order, syntheticData, 'displacement', inversionPars, delta_p, plotFlag, progressBar);
[delta_grad_vel, grad_adj_vel, grad_fd_vel]    = compare_gradients(parset, order, syntheticData, 'velocity', inversionPars, delta_p, plotFlag, progressBar);


function [delta_grad, grad_adj, grad_fd] = compare_gradients(parset, order, syntheticData, misfitType, inversionPars, delta_p, plotFlag, progressBar)
    parset.misfitType = misfitType;
    adj_opt = AntiplaneShear2DRSFrictionOpt(parset, inversionPars, order);
    adj_opt.setSyntheticReceiverData(syntheticData{:});
    grad_adj = adj_opt.computeGradient(plotFlag,progressBar);
    if adj_opt.tsOpts.forwardMethod.adaptive
        adj_opt.tsOpts.forwardMethod.adaptive = false;
        adj_opt.tsOpts.k = adj_opt.forwardTimeIntegrationData.k;
    end
        
    grad_fd = adj_opt.computeGradientFD(delta_p);
    delta_grad = adj_opt.compareGradients(grad_adj,grad_fd);     
    figure();
    fprintf('%s Misfit = %f\n',misfitType,adj_opt.computeMisfit());
    for i = 1:numel(inversionPars)
        parName = inversionPars{i};
        plot(delta_grad.(parName),'LineWidth',2);
        hold on;
        fprintf('Gradient L2-norm absolute error parameter %s: %e\n', parName, norm(delta_grad.(parName)));
        fprintf('Gradient L2-norm relative error parameter %s: %e\n', parName, norm(delta_grad.(parName))/norm(grad_adj.(parName)));
    end
    xlabel('Fault index');
    ylabel('Absolute error')
    legend(inversionPars{:});
    title(sprintf('Error using %s misfit', misfitType));
    hold off;

    % Plot grad_fd and grad_adj for visual comparison
    % figure();
    % for i = 1:numel(inversionPars)
    %     parName = inversionPars{i};
    %     subplot(numel(inversionPars),1,i);
    %     plot(grad_adj.(parName),'LineWidth',2);
    %     hold on;
    %     plot(grad_fd.(parName),'--','LineWidth',2);
    %     xlabel('Fault index');
    %     ylabel('Gradient');
    %     legend('Adjoint Gradient','Finite Difference Gradient');
    %     title(sprintf('Gradient comparison for parameter %s using %s misfit', parName, misfitType));
    % end

    % Test different FD step sizes and plot them on the same figure for comparison
    figure();
    % First plot adjoint
    hold on;
    for i = 1:numel(inversionPars)
        parName = inversionPars{i};
        plot(grad_adj.(parName),'LineWidth',2);
    end
    fd_steps = logspace(-7, -11, 5); 

    for step = fd_steps
        tic
        grad_fd = adj_opt.computeGradientFD(step);
        toc
        % Then plot FD
        for i = 1:numel(inversionPars)
            parName = inversionPars{i};
            plot(grad_fd.(parName),'--','LineWidth',2);
        end

        delta_grad = adj_opt.compareGradients(grad_adj, grad_fd);
        
        fprintf('\n--- FD Step Size: %e ---\n', step);
        for i = 1:numel(inversionPars)
            parName = inversionPars{i};
            rel_error = norm(delta_grad.(parName)) / norm(grad_adj.(parName));
            fprintf('Parameter %s relative error: %e\n', parName, rel_error);
        end
    end
    hold off;
end