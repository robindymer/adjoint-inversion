clear all;
order = 4;
%parset = @pars.rsFrictionErickson;
%parset = @pars.rsFrictionLinear;
%parset = @pars.rsFrictionTimeDepNonlinear;
parset = @pars.rsFrictionOffFaultTransport;

[parset, parset_true] = parset(); % Specify misfit later.

% Get synthetic data
% TODO: Add option for loading data from mat files
data_opt = AntiplaneShearRSFrictionOpt(parset_true, [], order);
data_opt.runForward(true);
syntheticData = {data_opt.forwardReceiverRecordings,data_opt.forwardTimeIntegrationData};

%[delta_grad_disp, grad_adj_disp, grad_fd_disp] = compare_gradients(parset, order, syntheticData, 'displacement');
[delta_grad_vel, grad_adj_vel, grad_fd_vel]    = compare_gradients(parset, order, syntheticData, 'velocity');


function [delta_grad, grad_adj, grad_fd] = compare_gradients(parset, order, syntheticData, misfitType)
    parset.misfitType = misfitType;
    adj_opt = AntiplaneShearRSFrictionOpt(parset, [], order);
    adj_opt.setSyntheticReceiverData(syntheticData{:});
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