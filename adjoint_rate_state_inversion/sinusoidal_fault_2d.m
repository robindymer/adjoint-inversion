order = 4;
m = 401;

plot_opts = struct;
plot_opts.plot_variables = 'all';

%% Parameters
[~, params] = pars.rsFriction2DSinusoidal();

% Initial conditions
opset = params.opset;
bc = params.bc;
friction = params.friction;
sources = [];
material = params.material;
initialconditions = params.initialconditions;
T = params.T;
domain = params.domain;

%% Create discretization
discr = elastic.AntiplaneShear2DRSFrictionFwdDiscr(opset, domain, m, order, material, bc, friction, sources, initialconditions);

method = params.tsOpts.forwardMethod;
k = params.tsOpts.k;
[ts,N] = discr.getTimestepper(method,T,k);

upd = discr.setupPlot(plot_opts);
repr = discr.getTimeSnapshot(0);
upd(repr);
% pause();
delta_t_plot = 1*T/(m-1);
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
w = ts.getV();
u = discr.E.u*w;

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