close all;
clear all;
addSubpaths;
order = 8;
inversionPars = {'a'};
initialGuessScaling = 1.1;
deltaParExps = [-5:-1:-8,-8.3,-8.5:-0.5:-9.5,-10:-1:-12];

parset = @pars.rsFriction2DFractalFaultVerification;
parSetOpts.inversionParameters = inversionPars;
parSetOpts.initialGuessScalings = initialGuessScaling;

parSetOpts.misfitType = 'displacement';
adjopt.runFDConvergence(@AntiplaneShear2DRSFrictionOpt,parset,inversionPars,parSetOpts,deltaParExps,order,2);
drawnow;
figure(1);
ax1 = gca;

parSetOpts.misfitType = 'velocity';
adjopt.runFDConvergence(@AntiplaneShear2DRSFrictionOpt,parset,inversionPars,parSetOpts,deltaParExps,order,2);
drawnow;
figure(2);
ax2 = gca;
% Fix plots
copyobj(ax2.Children(2), ax1);
uistack(ax1.Children(2),'up',1);

ax1.Children(2).Color = [0.8500 0.3250 0.0980]; % Red/Orange
ax1.Children(2).Marker = '+';
ax1.Children(2).MarkerFaceColor = [0.8500 0.3250 0.0980]; % Red/Orange

ax1.Children(3).Color = [0 0.4470 0.7410]; % Blue
ax1.Children(3).Marker = 'o';
ax1.Children(3).MarkerFaceColor = [0 0.4470 0.7410]; % Blue

grid(ax1,'on')
legend(ax1,'Displacement misfit', 'Velocity misfit', '1st-order reference', 'interpreter','latex');
xlabel(ax1,'$\Delta a$','Interpreter','latex');
ylabel(ax1,'$e(\Delta a)$','Interpreter','latex');
axis(ax1,[1e-12, 1e-5, 1e-6, 1e-1]);
set(ax1,'fontsize',16);
close(figure(2));







