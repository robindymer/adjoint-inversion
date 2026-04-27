clear all;
close all;

% Discretization
ms = 17*2.^[1:3];
orders = 4:2:8;

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
ts_method = parameters.tsOpts.forwardMethod;


%% Analytic solution
a = friction.params.a;
x_s = initialconditions.x_s;
sigma = initialconditions.sigma;
tau = @(V) friction.funs.tau(0, V, 0, a);
[u_m, u_p] = semianalytic_solution_erickson_2022(tau,x_s,sigma);
u_analytic = @(x_m,x_p,t) [u_m(x_m,t);u_p(x_p,t)];

%% Create discretization
no = length(orders);
nm = length(ms);
err = cell(no,1);
q = cell(no,1);
for i = 1:no
    o = orders(i);
    eH = zeros(nm,1);
    i
    for j = 1:nm
        j
        m = ms(j);
        discr = elastic.AntiplaneShearRSFrictionFwdDiscr(opset, m, xlims, o, material, bc, friction, sources, initialconditions);
        [ts,N] = discr.getTimestepper(ts_method,T);
        if ts_method.adaptive
            time.runAdaptiveTS(ts, T);
        else
            ts.stepN(N,true);
        end
        w = ts.getV();
        u = discr.E.u*w;
        x_m = discr.grid.grids{1}.points();
        x_p = discr.grid.grids{2}.points();
        e = u-u_analytic(x_m,x_p,T);
        H = discr.H;
        eH(j) = sqrt(e'*H*e);
    end
    err{i} = eH;
    q{i} = log(eH(1:end-1)./eH(2:end))./log((ms(2:end)-1)./(ms(1:end-1)-1));
end       
convergenceTable(sprintf('L2, t = %f', T), num2cell(orders), ms, err, q,'tex');