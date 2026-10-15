
clear all;
addSubpaths;
order = 6;
inversionPars = {'a'};
initialGuessScaling = 1.1;
T = 4;
m = 26;
m_p = 26;

parset = @pars.rsFriction2DFractalFaultVerification; 
delta_p = 1e-9; % For FD (also second pertubation?)

parsetOpts.inversionParameters = inversionPars;
parsetOpts.initialGuessScalings = initialGuessScaling;
parsetOpts.T = T;
parsetOpts.m = m;
parsetOpts.m_p = m_p;
parsetOpts.k_idx = 17;
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
    hessVec_fd = adj_opt.computeHessianVectorFD(delta_p);
    toc

    % Trickery!
    hessVec_struct = struct();
    if ~isstruct(hessVec)
        hessVec_struct.a = hessVec;
        hessVec = hessVec_struct;
    end

    delta_hessVec = adj_opt.compareGradients(hessVec,hessVec_fd);
    disp("Delta Hessian-Vector:")
    delta_hessVec.a
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

    figure();
    % Plot hessvec and hessVec_fd for visual comparison
    for i = 1:numel(inversionPars)
        parName = inversionPars{i};
        subplot(numel(inversionPars),1,i);
        plot(hessVec.(parName), 'b-', 'LineWidth', 2);
        hold on;
        plot(hessVec_fd.(parName), 'r--', 'LineWidth', 2);
        legend('Adjoint Hessian-Vector', 'FD Hessian-Vector');
        title(sprintf('Hessian-Vector for parameter %s', parName));
        xlabel('Fault index');
        ylabel('Hessian-Vector value');
        hold off;
    end

    % Relative error of the FD Hessian-vector product vs FD step size
    % Test different FD step sizes
    fd_steps = logspace(-4, -10, 7);
    rel_err = zeros(numel(inversionPars), numel(fd_steps));
    for j = 1:numel(fd_steps)
        tic
        hessVec_fd = adj_opt.computeHessianVectorFD(fd_steps(j));
        toc
        delta_hessVec = adj_opt.compareGradients(hessVec, hessVec_fd);
        for i = 1:numel(inversionPars)
            parName = inversionPars{i};
            rel_err(i,j) = norm(delta_hessVec.(parName)) / norm(full(hessVec.(parName)));
            fprintf('step %.1e, %s: relative error %e\n', fd_steps(j), parName, rel_err(i,j));
        end
    end

    figure();
    loglog(fd_steps, rel_err, 'o-', 'LineWidth', 2);
    hold on;
    % First-order reference slope (forward difference), anchored at the largest step
    loglog(fd_steps, rel_err(1,1)*fd_steps/fd_steps(1), 'k--', 'LineWidth', 1.5);
    hold off;
    grid on;
    xlabel('FD step \delta_p');
    ylabel('Relative error');
    legend([inversionPars, {'O(\delta_p)'}], 'Location', 'best');
    title(sprintf('Hessian-vector FD error vs step size (%s misfit)', misfitType));
end