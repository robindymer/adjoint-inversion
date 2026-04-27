clear all;
order = 4;

parset = @pars.rsFriction2DSinusoidal;
paropts.misfitType = 'velocity';
[parset, parset_true] = parset(paropts);

inversionPars = {'a'};

plotFlag = true;

% Create synthetic data
% TODO: Add option for loading data from mat files
data_opt = AntiplaneShear2DRSFrictionOpt(parset_true, [], order);
data_opt.runForward();
syntheticData = {data_opt.forwardReceiverRecordings,data_opt.forwardTimeIntegrationData};

adj_opt = AntiplaneShear2DRSFrictionOpt(parset, inversionPars, order);
adj_opt.setSyntheticReceiverData(syntheticData{:});

% Can use
% grad = adj_opt.computeGradient(plotFlag);
% or manually:
adj_opt.runForward(plotFlag);
adj_opt.updateAdjointDiscr();
adj_opt.runAdjoint(plotFlag);
grad = adj_opt.gradientFormula();

%% Plotting
nStages = adj_opt.forwardTimeIntegrationData.nStages;
t = adj_opt.removeStagedData(adj_opt.forwardTimeIntegrationData.T,nStages);

% Receiver recordings
figure();
nrec = length(adj_opt.forwardReceiverRecordings);
legend_str = cell(nrec,1);
for i = 1:nrec
    rec = adj_opt.removeStagedData(adj_opt.forwardReceiverRecordings{i},nStages);
    plot(t,rec,'LineWidth',2);
    hold on;
    legend_str{i} = sprintf('Receiver %d',i);
end
hold off;
xlabel('$t$','interpreter','latex')
ylabel('Recorded data','interpreter','latex');
legend(legend_str{:});

% Gradient
figure();
plot(grad.a,'LineWidth',2);
xlabel('Fault index','interpreter','latex');
ylabel('$\mathcal{F}_a$','Interpreter','latex');

% Forward and adjoint slip velocity (V*) and state evolution
[nf, ~] = size(adj_opt.forwardFaultVariables.V);
figure();
subplot(2,2,1)
V = adj_opt.removeStagedData(adj_opt.forwardFaultVariables.V,nStages);
surface(t,1:nf,V);
xlabel('$t$','interpreter','latex')
ylabel('Fault index','interpreter','latex')
title('$V^*$','interpreter','latex')
shading flat;
axis tight
colorbar;
subplot(2,2,2)
Psi = adj_opt.removeStagedData(adj_opt.forwardFaultVariables.Psi,nStages);
surface(t,1:nf,Psi);
xlabel('$t$','interpreter','latex')
ylabel('Fault index','interpreter','latex')
title('$\Psi$','interpreter','latex')
shading flat;
axis tight
colorbar;
subplot(2,2,3)
V_adj = fliplr(adj_opt.removeStagedData(adj_opt.adjointFaultVariables.V,nStages)); % Flip from reverse to forward time
surface(t,1:nf,V_adj);
xlabel('$t$','interpreter','latex')
ylabel('Fault index','interpreter','latex')
title('$V^{\dagger*}$','interpreter','latex')
shading flat;
axis tight
colorbar;
subplot(2,2,4)
Psi_adj = fliplr(adj_opt.removeStagedData(adj_opt.adjointFaultVariables.V,nStages)); % Flip from reverse to forward time
surface(t,1:nf,Psi_adj);
xlabel('$t$','interpreter','latex')
ylabel('Fault index','interpreter','latex')
title('$\Psi^\dagger$','interpreter','latex')
shading flat;
axis tight
colorbar;
drawnow

% Derivatives of F and G w.r.t V and Psi
figure();
subplot(2,2,1)
F_V = fliplr(adj_opt.removeStagedData(adj_opt.adjointDiscr.friction.data.F_V,nStages));  % Flip from reverse to forward time
surface(t,1:nf,F_V);
xlabel('$t$','interpreter','latex')
ylabel('Fault index')
title('$F_V$','interpreter','latex')
shading flat;
axis tight
colorbar;
subplot(2,2,2)
G_V = fliplr(adj_opt.removeStagedData(adj_opt.adjointDiscr.friction.data.G_V,nStages));  % Flip from reverse to forward time
surface(t,1:nf,G_V);
xlabel('$t$','interpreter','latex')
ylabel('Fault index','interpreter','latex')
title('$G_V$','interpreter','latex')
shading flat;
axis tight
colorbar;
subplot(2,2,3)
F_Psi = fliplr(adj_opt.removeStagedData(adj_opt.adjointDiscr.friction.data.F_Psi,nStages)); % Flip from reverse to forward time
surface(t,1:nf,F_Psi);
xlabel('$t$','interpreter','latex')
ylabel('Fault index','interpreter','latex')
title('$F_\Psi$','interpreter','latex')
shading flat;
axis tight
colorbar;
subplot(2,2,4)
G_psi = fliplr(adj_opt.removeStagedData(adj_opt.adjointDiscr.friction.data.G_Psi,nStages)); % Flip from reverse to forward time
surface(t,1:nf,G_psi);
xlabel('$t$','interpreter','latex')
ylabel('Fault index','interpreter','latex')
title('$G_\Psi$','interpreter','latex')
shading flat;
axis tight
colorbar;
drawnow