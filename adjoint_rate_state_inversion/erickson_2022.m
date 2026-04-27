clear all;
% close all;
compute_semianalytic = true;

order = 4;

%% Parameters
[~,parameters] = pars.rsFrictionErickson();

% Unpack parameter struct
opset = parameters.opset;
m = parameters.m;
T = parameters.T;
xlims = parameters.xlims;
material = parameters.material;
bc = parameters.bc;
friction = parameters.friction;
sources = parameters.sources;
initialconditions = parameters.initialconditions;
interp_data = parameters.interpolate_data;
ts_method = parameters.tsOpts.forwardMethod;

%% Create discretization
discr = elastic.AntiplaneShearRSFrictionFwdDiscr(opset, m, xlims, order, material, bc, ...
                    friction, sources, initialconditions);
[ts,N] = discr.getTimestepper(ts_method,T);

plot_opts.plot_variables = 'u';
plot_opts.axis_u = [-1,1,0,1];
upd = discr.setupPlot(plot_opts);
repr = discr.getTimeSnapshot(0);
upd(repr);
gca;
title(sprintf('t = %f',repr.t));
legend('approx');
pause();
repr = discr.getTimeSnapshot(ts);
if ~ts_method.adaptive
    for i = 1:N
        if mod(i,20) == 0
            repr = discr.getTimeSnapshot(ts);
            upd(repr); 
            gca;
            title(sprintf('t = %f',repr.t));
            legend('approx');
            drawnow;
        end
        ts.step();
    end
else
    ts = time.runAdaptiveTS(ts, T);   
end
repr = discr.getTimeSnapshot(ts);
upd(repr); 
gca;
title(sprintf('t = %f',repr.t));
legend('approx');
drawnow;
hold on;

if compute_semianalytic
    a = friction.params.a;
    x_s = initialconditions.x_s;
    sigma = initialconditions.sigma;
    tau = @(V) friction.funs.tau(0, V, 0, a);
    [u_m, u_p, V] = semianalytic_solution_erickson_2022(tau,x_s,sigma);
    x_m = discr.grid.grids{1}.points();
    x_p = discr.grid.grids{2}.points();
    plot(x_m,u_m(x_m,T),x_p,u_p(x_p,T),'linewidth',2);
    legend('approx','exact');
end