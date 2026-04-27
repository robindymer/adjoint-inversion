addSubpaths;
order = 8;

plot_opts = struct;
plot_opts.plot_variables = 'v';
plot_opts.axlims = [-10 10 -10 10 -10 10];

%% Parameters
% High-resolution setup using true parameter value
[~, params] = pars.rsFriction2DFractalFaultForward();
% Lower-resolution inversion setup using initial parameter value
%opts.initialGuessValues = 0.0135;
%params = pars.rsFriction2DFractalFaultInversion(opts);

% Initial conditions
opset = params.opset;
bc = params.bc;
friction = params.friction;
sources = [];
material = params.material;
initialconditions = params.initialconditions;
T = params.T;
domain = params.domain;
m = params.m;

%% Create discretization
discr = elastic.AntiplaneShear2DRSFrictionFwdDiscr(opset, domain, m, order, material, bc, friction, sources, initialconditions);

method = params.tsOpts.forwardMethod;
k = params.tsOpts.k;
[ts,N] = discr.getTimestepper(method,T,k);

upd = discr.setupPlot(plot_opts);
repr = discr.getTimeSnapshot(0);
upd(repr);
tic;
time_stamps = [2,3,5];
%time_stamps = [];
if ~isempty(time_stamps)
    for it = 1:length(time_stamps)
        t = time_stamps(it);
        ts.evolve(t, true);
        repr = discr.getTimeSnapshot(ts);
        upd(repr);
        filestr = sprintf("fractal_fault_m%d_t%d.fig",m,t); 
        savefig(filestr);
    end
else
    for i = 1:N
        ts.step();
        repr = discr.getTimeSnapshot(ts);
        upd(repr);
        drawnow;
    end
    repr = discr.getTimeSnapshot(ts);
    upd(repr);
end
toc;
w = ts.getV();
u = discr.E.u*w;