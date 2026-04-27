% A line of nRec equispaced receivers at ocean bottom.
function [pars, truePars] = marmousi2ElasticSources(opts, downSampling)

% Default values
default_arg('downSampling', 20);
default_arg('opts', struct);

opts = defaultField(opts, 'nS', 1);
opts = defaultField(opts, 'T', 1);

nS = opts.nS;
T = opts.T;
T = 3;

dim = 2;

% Read from .mat file
% filename = sprintf('+pars/marmousi2padded_downsampled_%d.mat',downSampling);
% filename = sprintf('+pars/marmousi2_downsampled_%d.mat',downSampling);
% filename = sprintf('+pars/marmousi2smoothed_downsampled_%d.mat',downSampling);
% filename = sprintf('+pars/marmousi2const_downsampled_%d.mat',downSampling);
filename = sprintf('+pars/marmousi2crop_downsampled_%d.mat',downSampling);
load(filename);


%-- Elastic domain south of acoustic --

% Acoustic
xmin_w = min(min(X_w));
xmax_w = max(max(X_w));
ymin_w = min(min(Z_w));
ymax_w = max(max(Z_w));
[Ny, Nx] = size(X_w);

xlim = [xmin_w, xmax_w];
ylim = [ymax_w, ymin_w];

domain.acoustic.xlim = xlim;
domain.acoustic.ylim = ylim;
domain.acoustic.m = {[Nx, Ny]};
domain.acoustic.def = multiblock.domain.Rectangle(xlim, ylim);
domain.acoustic.interfaceGroup = 'S';

% Elastic
xmin_e = min(min(X_e));
xmax_e = max(max(X_e));
ymin_e = min(min(Z_e));
ymax_e = max(max(Z_e));
[Ny, Nx] = size(X_e);

xlim = [xmin_e, xmax_e];
ylim = [ymax_e, ymin_e];

domain.elastic.xlim = xlim;
domain.elastic.ylim = ylim;
domain.elastic.m = {[Nx, Ny]};
domain.elastic.def = multiblock.domain.Rectangle(xlim, ylim);
domain.elastic.interfaceGroup = 'N';
% ------------------------------


% ---- Boundary conditions ----
% Acoustic
bcw.boundary = domain.acoustic.def.boundaryGroups.W;
bcw.type = 'n';

bcn.boundary = domain.acoustic.def.boundaryGroups.N;
bcn.type = 'd';

bce.boundary = domain.acoustic.def.boundaryGroups.E;
bce.type = 'n';
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
bce1.type = {1, 'd'};

bce2.boundary = domain.elastic.def.boundaryGroups.E;
bce2.type = {2, 'd'};

bc.elastic = {bcw1, bcw2, bcs1, bcs2, bce1, bce2};
% -----------------------------

% Material parameters
mu = rho_e.*vs_e.^2;
lambda = rho_e.*vp_e.^2 - 2*mu;
lambda = lambda(:);
mu = mu(:);
rho = rho_e(:);

rho_acoustic = rho_w(:);
c = vp_w(:);

% ===== Sources ==========
sources = struct;

xmin = xlim(1);
xmax = xlim(2);
L = xmax_e - xmin_e;
xl = xmin_e + 3/8*L;
xr = xmin_e + 5/8*L;

%--- Elastic sources ---
A = 5e-1;

% Ricker wavelet
t0 = 1.2;
sigma = 1;
g = @(t) A*rickerWavelet(t, t0, sigma);
g0 = @(t) 0*t;

x = 10700;
D = 0;
y = linspace(ymax_e-D,ymax_e-D,nS);
sources.elastic.x = cell(1, nS);
for i = 1:nS
	sources.elastic.x{i} = [x(i), y(i)];
	sources.elastic.g{i} = {g, g};
	% sources.elastic.g{i} = {g0, g};
end
%-----------------------

%--- Acoustic sources ---
sources.acoustic.x = [];
sources.acoustic.g = [];
% =============================

% ===== Receivers ==========
receivers = struct;

%--- Elastic receivers ---
receivers.elastic.x = [];
%----------------------

%--- Acoustic receivers ---
receivers.acoustic.x = [];
%----------------------


% Perturb rho
% perturbedMu = {@(x,y) 2*mu(x,y), @(x,y) 1.5*mu(x,y)} ;
perturbedRho = {1.1*rho};

% Put other parameters in cell arrays
lambda = {lambda};
mu = {mu};
rho = {rho};

% Return struct
pars = struct;
pars.domain = domain;
pars.elastic.lambda = lambda;
pars.elastic.mu = mu;
pars.elastic.rho = perturbedRho;

pars.acoustic.c = c;
pars.acoustic.rho = rho_acoustic;

pars.sources = sources;
pars.receivers = receivers;
pars.bc = bc;
pars.T = T;
% pars.m = ms;
pars.dim = dim;

truePars = pars;
truePars.elastic.rho = rho;

end
