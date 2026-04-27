%clear all;
order = 8;
% 
parset = @pars.rsFriction2DFractalFaultInversion;
inversionPars = {'a'};
paropts.misfitType = 'velocity';
paropts.inversionParameters = inversionPars;
paropts.initialGuessValues = 0.0135;
paropts.receiverSpacing = 2;
paropts.doFilter = false;
[parset, parset_true] = parset(paropts);

plotFlag = true;
progressBar = true;

% Generate this data by running fractal_fault_generate_data.m with the 
% pars.rsFriction2DFractalFaultForward parameter set.
loadData = false;
loadPath = 'mat/fractal_fault_m1001_mp1001/receiverData.mat'; 

% Create optimization object
adj_opt = AntiplaneShear2DRSFrictionOpt(parset, inversionPars, order);

% Data
if loadData
    fprintf('Loading data from %s\n',loadPath);
    load(loadPath);
    nr = numel(parset_true.receivers.x);
    adj_opt.receiverData = cell(nr,1);
    for i = 1:nr
        for j = 1:numel(receiverData.positions)
            if receiverData.positions{j} == parset_true.receivers.x{i}
                break
            end
        end
        adj_opt.receiverData{i} = receiverData.recordings{j};
    end
else
    disp('Generating synthetic data');
    data_opt = AntiplaneShear2DRSFrictionOpt(parset_true, inversionPars, order);
    data_opt.runForward(plotFlag, parset_true.T, [], progressBar);
    adj_opt.setSyntheticReceiverData(data_opt.forwardReceiverRecordings, data_opt.forwardTimeIntegrationData);
end




% Perform forward and backward solve, computing the gradient.
disp('Computing gradient');
adj_opt.runForward(plotFlag, parset.T, [], progressBar);
adj_opt.updateAdjointDiscr();
adj_opt.runAdjoint(plotFlag, parset.T, [], progressBar);
grad = adj_opt.gradientFormula();

%% Plotting
domain = parset.domain;
fault_id = domain.boundaryGroups.fault_minus;
nStages = adj_opt.forwardTimeIntegrationData.nStages;
t = adj_opt.removeStagedData(adj_opt.forwardTimeIntegrationData.T,nStages);
X_fault = adj_opt.forwardDiscr.grid.getBoundary(fault_id);
x = X_fault(:,1);



% Gradient
figure();
plot(grad.a,'LineWidth',2);
xlabel('Parameter index','interpreter','latex');
ylabel('$\frac{\partial \mathcal{F}}{\partial p}$','interpreter','latex');

% Forward and adjoint slip velocity (V*) and state evolution
figure();
Vfwd = adj_opt.removeStagedData(adj_opt.forwardFaultVariables.V,nStages);
surface(x,t,Vfwd');
xlabel('$x$ (km)','interpreter','latex')
ylabel('$t$ (s)','interpreter','latex')
title('$V^*$','interpreter','latex')
shading flat;
axis tight
colorbar;
figure();
Psifwd = adj_opt.removeStagedData(adj_opt.forwardFaultVariables.Psi,nStages);
surface(x,t,Psifwd');
xlabel('$x$ (km)','interpreter','latex')
ylabel('$t$ (s)','interpreter','latex')
title('$\Psi$','interpreter','latex')
shading flat;
axis tight
colorbar;

figure();
Vadj = fliplr(adj_opt.removeStagedData(adj_opt.adjointFaultVariables.V,nStages));
surface(x,t,Vadj');
xlabel('$x$ (km)','interpreter','latex')
ylabel('$t$ (s)','interpreter','latex')
title('$V^{\dagger*}$','interpreter','latex')
shading flat;
axis tight
colorbar;
figure();
Psiadj = fliplr(adj_opt.removeStagedData(adj_opt.adjointFaultVariables.Psi,nStages));
surface(x,t,Psiadj');
xlabel('$x$ (km)','interpreter','latex')
ylabel('$t$ (s)','interpreter','latex')
title('$\Psi^\dagger$','interpreter','latex')
shading flat;
axis tight
colorbar;

% Derivatives of F and G w.r.t V and Psi
figure();
F_V = fliplr(adj_opt.removeStagedData(adj_opt.adjointDiscr.friction.data.F_V,nStages));  % Flip from reverse to forward time
surface(x,t,F_V');
xlabel('$x$ (km)','interpreter','latex')
ylabel('$t$ (s)','interpreter','latex')
title('$F_V$','interpreter','latex')
shading flat;
axis tight
colorbar;
figure();
G_V = fliplr(adj_opt.removeStagedData(adj_opt.adjointDiscr.friction.data.G_V,nStages));  % Flip from reverse to forward time
surface(x,t,G_V');
xlabel('$x$ (km)','interpreter','latex')
ylabel('$t$ (s)','interpreter','latex')
title('$G_V$','interpreter','latex')
shading flat;
axis tight
colorbar;
figure();
F_Psi = fliplr(adj_opt.removeStagedData(adj_opt.adjointDiscr.friction.data.F_Psi,nStages)); % Flip from reverse to forward time
surface(x,t,F_Psi');
xlabel('$x$ (km)','interpreter','latex')
ylabel('$t$ (s)','interpreter','latex')
title('$F_\Psi$','interpreter','latex')
shading flat;
axis tight
colorbar;
figure();
G_Psi = fliplr(adj_opt.removeStagedData(adj_opt.adjointDiscr.friction.data.G_Psi,nStages)); % Flip from reverse to forward time
surface(x,t,G_Psi');
xlabel('$x$ (km)','interpreter','latex')
ylabel('$t$ (s)','interpreter','latex')
title('$G_\Psi$','interpreter','latex')
shading flat;
axis tight
colorbar;
drawnow