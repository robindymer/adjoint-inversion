% A line of nRec equispaced receivers from point x0 to point x1

function [pars, truePars] = muElAc2Blocks(opts)

% Default values
default_arg('opts', struct);

opts = defaultField(opts, 'm', {[6, 7],[6,7]});
% opts = defaultField(opts, 'x0', [-0.6, -0.5]);
% opts = defaultField(opts, 'x1', [0.6, -0.5]);
opts = defaultField(opts, 'nRec', 10);

ms = opts.m;
% x0 = opts.x0;
% x1 = opts.x1;
nRec = opts.nRec;

% ---- Domain ----------------
xl = -1; xm = 0; xr = 1;
yl = -1; ym = 0; yr = 1;
T = 5;
dim = 2;

%-- Elastic domain south of acoustic --
xlim = [xl, xm, xr];
% xlim = [xl, xr];
ylim = [yr, ym];
domain.acoustic.xlim = xlim;
domain.acoustic.ylim = ylim;
domain.acoustic.m = ms;
domain.acoustic.def = multiblock.domain.Rectangle(xlim, ylim);
domain.acoustic.interfaceGroup = 'S';

xlim = [xl, xm, xr];
% xlim = [xl, xr];
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
bcw1.boundary = domain.elastic.def.boundaryGroups.W;
bcw1.type = {1, 'd'};

bcw2.boundary = domain.elastic.def.boundaryGroups.W;
bcw2.type = {2, 'd'};

bcs1.boundary = domain.elastic.def.boundaryGroups.S;
bcs1.type = {1, 'd'};

bcs2.boundary = domain.elastic.def.boundaryGroups.S;
bcs2.type = {2, 'd'};

bce1.boundary = domain.elastic.def.boundaryGroups.E;
bce1.type = {1, 't'};

bce2.boundary = domain.elastic.def.boundaryGroups.E;
bce2.type = {2, 't'};

bc.elastic = {bcw1, bcw2, bcs1, bcs2, bce1, bce2};
% -----------------------------

% Material parameters
% lambda = @(x,y) 0*x +1;
% mu = @(x,y) 0*x +1;
% rho = @(x,y) 0*x +1;

lambda = @(x,y) 2 + sin(x+y);
mu = @(x,y) 2 + cos(x+y);
rho = @(x,y) 2 + sin(x-y);

rho_acoustic = @(x,y) 0*x + 1;
c = @(x,y) 0*x + 1;

%--- Elastic sources ---
A = 30;
sigma = 0.2;
t0 = 1;
g = @(t) A*exp(-(t-t0).^2 / sigma^2 );
sources = struct;
sources.elastic.x = {[0.1, -0.5]};
sources.elastic.g = {{g, g}};
sources.elastic.x = [];
%-----------------------

%--- Elastic receivers ---
% Receiver line
receivers = struct;
x = linspace(xl,xr,nRec);
y = linspace(ym,ym,nRec);
receivers.elastic.x = cell(1, nRec);
for i = 1:nRec
	receivers.elastic.x{i} = [x(i), y(i)];
end
% receivers.elastic.x = [];
%----------------------

%--- Acoustic sources ---
A = 1e0;
sigma = 0.2;
t0 = 1;
g = @(t) A*exp(-(t-t0).^2 / sigma^2 );

nS = 10;
x = linspace(xl,xr,nS);
y = linspace(yr,yr,nS);
sources.acoustic.x = cell(1, nS);
for i = 1:nS
	sources.acoustic.x{i} = [x(i), y(i)];
	sources.acoustic.g{i} = g;
end
% sources.acoustic.x = {[0.5, 0.75]};
% sources.acoustic.g = {g};
%-----------------------

%--- Acoustic receivers ---
x = linspace(xl,xr,nRec);
y = linspace(0.5,0.5,nRec);
receivers.acoustic.x = cell(1, nRec);
for i = 1:nRec
	receivers.acoustic.x{i} = [x(i), y(i)];
end
% receivers.acoustic.x = [];
% receivers.acoustic.x = {[0.5, 0.25]};
%----------------------

% Perturb mu
perturbedMu = {@(x,y) 2*mu(x,y), @(x,y) 1.5*mu(x,y)} ;
% perturbedMu = {@(x,y) 2*mu(x,y)} ;

% Put other parameters in cell arrays
lambda = {lambda, lambda};
mu = {mu, mu};
rho = {rho, rho};
% lambda = {lambda};
% mu = {mu};
% rho = {rho};



% Return struct
pars = struct;
pars.domain = domain;
pars.elastic.lambda = lambda;
pars.elastic.mu = perturbedMu;
pars.elastic.rho = rho;

pars.acoustic.c = c;
pars.acoustic.rho = rho_acoustic;

pars.sources = sources;
pars.receivers = receivers;
pars.bc = bc;
pars.T = T;
pars.m = ms;
pars.dim = dim;

truePars = pars;
truePars.elastic.mu = mu;

end
