clear all;

order = 4;
ic_method = 'erickson2022';

%% Parameters
opts = struct;
opts.ic_method = ic_method;

[~, params] = pars.rsFrictionOffFaultTransport(opts);
%[params, ~] = pars.rsFrictionOffFaultTransport();

% Initial conditions
opset = params.opset;
m = params.m;
xlims = params.xlims;
bc = params.bc;
friction = params.friction;
sources = [];
material = params.material;
initialconditions = params.initialconditions;
interpolate_data = params.interpolate_data;
T = params.T;

%% Create discretization

discr = elastic.AntiplaneShearRSFrictionFwdDiscr(opset, m, xlims, order, material, bc, friction, sources, initialconditions);
method = params.tsOpts.forwardMethod;
k = params.tsOpts.k;
[ts,N] = discr.getTimestepper(method,T,k);

plot_opts.plot_variables = 'trajectory';
plot_opts.axis_u = [xlims(1),xlims(3),-0.001,0.001];
plot_opts.axis_v = [xlims(1),xlims(3),-0.01,0.01];
upd = discr.setupPlot(plot_opts);
repr = discr.getTimeSnapshot(0);
upd(repr);
pause();
delta_t_plot = 0.0001*T/(m-1);
tic;
if method.adaptive
    ts = time.runAdaptiveTS(ts, T, @(next_plot_time)plot_callback(discr, ts, upd, next_plot_time, delta_t_plot));
    
else
    next_plot_time = delta_t_plot;
    for i = 1:N
        ts.step();
        next_plot_time = plot_callback(discr,ts,upd,next_plot_time,delta_t_plot);
    end
end
toc;

function next_plot_time = plot_callback(discr, ts, upd, next_plot_time, delta_t_plot)
    if ts.t < next_plot_time
        return;
    else
        repr = discr.getTimeSnapshot(ts);
        upd(repr);
        drawnow;
        next_plot_time = ts.t + delta_t_plot;
    end
end