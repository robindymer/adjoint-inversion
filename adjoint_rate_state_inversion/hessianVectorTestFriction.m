clear all;
addSubpaths;
order = 4;
%parset = @pars.rsFrictionErickson;
%parset = @pars.rsFrictionLinear;
%parset = @pars.rsFrictionTimeDepNonlinear;
parset = @pars.rsFrictionOffFaultTransportHessian;

[parset, parset_true] = parset(); % Specify misfit later.

% Get synthetic data
% TODO: Add option for loading data from mat files
data_opt = AntiplaneShearRSFrictionOptHessian(parset_true, [], order);
data_opt.runForward(false);
syntheticData = {data_opt.forwardReceiverRecordings,data_opt.forwardTimeIntegrationData};

%[delta_grad_disp, grad_adj_disp, grad_fd_disp] = compare_gradients(parset, order, syntheticData, 'displacement');
[delta_hessVec, hessianVector, hessianVector_fd]    = compare_hessians(parset, order, syntheticData, 'velocity');

% Prints "Failed try" - this is due to adaptive time stepping, 
% RK4 checks if tol is met at each step, and if not, reduces the time step and tries again.
function [delta_hessVec, hessianVector, hessianVector_fd] = compare_hessians(parset, order, syntheticData, misfitType)
    parset.misfitType = misfitType;
    % TODO: Find better name than adj_opt
    adj_opt = AntiplaneShearRSFrictionOptHessian(parset, [], order);
    adj_opt.setSyntheticReceiverData(syntheticData{:});

    tic
    hessianVector = adj_opt.computeHessianVector();
    toc
    disp(hessianVector);
    tic
    hessianVector_fd = adj_opt.computeHessianVectorFD(1e-6);
    toc
    disp(hessianVector_fd);

    delta_hessVec = abs(hessianVector-hessianVector_fd)./(abs(hessianVector));
    fprintf('%s misfit:\n',misfitType);
    fprintf('HessVec: %e, \t HessVec FD: %e\n',full(hessianVector), hessianVector_fd);
    fprintf('Pointwise relative error: %e\n',delta_hessVec);
    % Makes sense that it is larger since hessian value is larger
    fprintf('Hessian vector norm absolute error: %e\n', adj_opt.gradientNorm(hessianVector-hessianVector_fd));
    fprintf('Hessian vector norm relative error: %e\n', adj_opt.gradientNorm(hessianVector-hessianVector_fd)/adj_opt.gradientNorm(hessianVector));
    
    % Dummy assigning
    % delta_grad = 0;
    % grad_adj = 0;
    % grad_fd = 0;
    grad_adj = adj_opt.computeGradient();
    if adj_opt.tsOpts.forwardMethod.adaptive
        adj_opt.tsOpts.forwardMethod.adaptive = false;
        adj_opt.tsOpts.k = adj_opt.forwardTimeIntegrationData.k;
    end

    grad_fd = adj_opt.computeGradientFD(1e-6);

    delta_grad = abs(grad_adj-grad_fd)./(abs(grad_adj));
    fprintf('%s misfit:\n',misfitType);
    fprintf('Grad Adj: %e, \t Grad FD: %e\n',full(grad_adj), grad_fd);
    fprintf('Pointwise relative error: %e\n',delta_grad);
    fprintf('Gradient norm absolute error: %e\n', adj_opt.gradientNorm(grad_adj-grad_fd));
    fprintf('Gradient norm relative error: %e\n', adj_opt.gradientNorm(grad_adj-grad_fd)/adj_opt.gradientNorm(grad_adj));
end