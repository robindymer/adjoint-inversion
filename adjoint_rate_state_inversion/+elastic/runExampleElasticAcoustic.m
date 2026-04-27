function runExampleElasticAcoustic(m,T,order,point_sources,check_ev_flag)

default_arg('m', [51, 51]);
default_arg('T', 10);
default_arg('order', 4);
default_arg('check_ev_flag',false)
default_arg('point_sources',[])

% ---- Forcing ------------------
f1 = @(x,y,t) 0*x;
f2 = @(x,y,t) 0*x;
f = @(t,x,y) 0*x;
F.elastic = {{f1, f2}, {f1, f2}};
F.acoustic = f;
%--------------------------------

% --- Material coefficients -----
l1 = @(x,y) 0*x + 1;
l2 = @(x,y) 0*x + 1;
m1 = @(x,y) 0*x + 1;
m2 = @(x,y) 0*x + 1;

rho_acoustic = @(x,y) 0*x + 1;
c = @(x,y) 0*x + 1;

lambda = {l1, l2};
mu = {m1, m2};
rho = {[],[]};

parameters.elastic.lambda = lambda;
parameters.elastic.mu = mu;
parameters.elastic.rho = rho;

parameters.acoustic.rho = rho_acoustic;
parameters.acoustic.c = c;
% -------------------------------

% ---- Domain ----------------
xl = -0.6; xm = 0; xr = 0.6;
yl = -1; ym = 0; yr = 1;
xlim = [xl, xr];
ylim = [yr, yl];
ms = {m, m};

%-- Elastic domain south of acoustic --
xlim = [xl, xm, xr];
ylim = [yr, ym];
domain.acoustic.xlim = xlim;
domain.acoustic.ylim = ylim;
domain.acoustic.m = ms;
domain.acoustic.def = multiblock.domain.Rectangle(xlim, ylim);
domain.acoustic.interfaceGroup = 'S';

xlim = [xl, xm, xr];
ylim = [ym, yl];
domain.elastic.xlim = xlim;
domain.elastic.ylim = ylim;
domain.elastic.m = ms;
domain.elastic.def = multiblock.domain.Rectangle(xlim, ylim);
domain.elastic.interfaceGroup = 'N';
% ------------------------------

% ---- Boundary conditions ----
% Acoustic
bcw.boundary = domain.acoustic.def.boundaryGroups.W;
bcw.type = 'd';

bcn.boundary = domain.acoustic.def.boundaryGroups.N;
bcn.type = 'n';

bce.boundary = domain.acoustic.def.boundaryGroups.E;
bce.type = 'd';
bc.acoustic = {bcw, bce, bcn};

% Elastic
bcws1.boundary = domain.elastic.def.boundaryGroups.WS;
bcws1.type = {1, 'd'};

bcws2.boundary = domain.elastic.def.boundaryGroups.WS;
bcws2.type = {2, 'd'};

bce1.boundary = domain.elastic.def.boundaryGroups.E;
bce1.type = {1, 't'};

bce2.boundary = domain.elastic.def.boundaryGroups.E;
bce2.type = {2, 't'};

bc.elastic = {bcws1, bcws2, bce1, bce2};
% -----------------------------

% --- Point sources -----------
pointSources = [];
% -----------------------------

% Create discretization object
discr = elastic.elasticAcousticDiscr(domain, order, parameters, bc, F, pointSources);
if check_ev_flag
	testElasticAcousticMatrix(discr);
end

% Initial data functions
x0 = 0;
y0 = -0.5;
u0 = @(t,x,y) 3*[exp(-((x-x0)/0.2).^2 - ((y-y0)/0.2).^2) ; exp(-((x-x0)/0.2).^2 - ((y-y0)/0.2).^2)];
ut0 = @(t,x,y) [0*x; 0*x];
phi0 = @(t,x,y) 0*x;
phit0 = @(t,x,y) 0*x;

% Initial data
u0 = grid.evalOn(discr.g.elastic, @(x,y) u0(0,x,y));
ut0 = grid.evalOn(discr.g.elastic, @(x,y) ut0(0,x,y));
phi0 = grid.evalOn(discr.g.acoustic, @(x,y) phi0(0,x,y));
phit0 = grid.evalOn(discr.g.acoustic, @(x,y) phit0(0,x,y));

v0 = [u0; phi0];
vt0 = [ut0; phit0];
discr.v0 = v0;
discr.v0t = vt0;


% Time-step
noname.animate(discr, [], T);

end