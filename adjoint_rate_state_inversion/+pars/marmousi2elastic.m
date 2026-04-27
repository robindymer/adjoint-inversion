% A line of nRec equispaced receivers at ocean bottom.
function [pars, truePars] = marmousi2elastic(opts, downSampling)

% Default values
default_arg('downSampling', 20);
default_arg('opts', struct);

opts = defaultField(opts, 'nRec', 10);
opts = defaultField(opts, 'nS', 2);
opts = defaultField(opts, 'T', 1);

nRec = opts.nRec;
nS = opts.nS;
T = opts.T;

dim = 2;

% Read from .mat file
filename = sprintf('+pars/marmousi2padded_downsampled_%d.mat',downSampling);
% filename = sprintf('+pars/marmousi2_downsampled_%d.mat',downSampling);
load(filename);


%-- Domain --
xmin = min(min(X));
xmax = max(max(X));
ymin = min(min(Z));
ymax = max(max(Z));
[Ny, Nx] = size(X);
m = {[Nx, Ny]};

xlim = [xmin, xmax];
ylim = [ymax, ymin];

domain = multiblock.domain.Rectangle(xlim, ylim);
% ------------------------------

% ---- Boundary conditions ----
bc = struct;
W = domain.boundaryGroups.W;
E = domain.boundaryGroups.E;
S = domain.boundaryGroups.S;
N = domain.boundaryGroups.N;

BCW1 = bcStruct(W, {1, 'd'}, []);
BCW2 = bcStruct(W, {2, 'd'}, []);
BCE1 = bcStruct(E, {1, 'd'}, []);
BCE2 = bcStruct(E, {2, 'd'}, []);

% Free surface
BCN1 = bcStruct(N, {1, 't'}, []);
BCN2 = bcStruct(N, {2, 't'}, []);

BCS1 = bcStruct(S, {1, 'd'}, []);
BCS2 = bcStruct(S, {2, 'd'}, []);

bc = {BCW1, BCW2, BCE1, BCE2, BCS1, BCS2, BCN1, BCN2};
% -----------------------------

% Material parameters
mu = rho.*vs.^2;
lambda = rho.*vp.^2 - 2*mu;
lambda = lambda(:);
mu = mu(:);
rho = rho(:);

% rho_acoustic = rho_w(:);
% c = vp_w(:);

% ===== Sources ==========
sources = struct;

xmin = xlim(1);
xmax = xlim(2);
L = xmax - xmin;
xl = xmin + 3/8*L;
xr = xmin + 5/8*L;

% A = 0;
% sigma = 0.2;
% t0 = 1;
% g = @(t) A*exp(-(t-t0).^2 / sigma^2 );
% sources.x = {[0.1, -0.5]};
% sources.g = {{g, g}};
sources.x = [];
sources.g = [];

% =============================

% ===== Receivers ==========
receivers = struct;

ymax_e = max(max(Z_e));

x = linspace(xl,xr,nRec);
if nRec == 1
	x = (xl+xr)/2;
end
y = linspace(ymax_e,ymax_e,nRec);
receivers.x = cell(1, nRec);
for i = 1:nRec
	receivers.x{i} = [x(i), y(i)];
end
% receivers.x = [];
%----------------------


% Perturb rho (Perturb exact rho to get marmousi rho in forward/adjoint)
% perturbedMu = {@(x,y) 2*mu(x,y), @(x,y) 1.5*mu(x,y)} ;
perturbedRho = {1.1*rho};

% Put other parameters in cell arrays
% lambda = {lambda, lambda};
% mu = {mu, mu};
% rho = {rho, rho};
lambda = {lambda};
mu = {mu};
rho = {rho};

% Return struct
pars = struct;
pars.xlim = xlim;
pars.ylim = ylim;
pars.lambda = lambda;
pars.mu = mu;
pars.rho = rho;

pars.sources = sources;
pars.receivers = receivers;
pars.bc = bc;
pars.T = T;
pars.m = m;
pars.dim = dim;

truePars = pars;
truePars.rho = perturbedRho;

end

function s = bcStruct(name, type, data);
	s = struct;
	s.boundary = name;
	s.type = type;
	s.data = data;
end

function I = myIntegral(f, a, b)
	I = zeros(size(b));
	for i = 1:length(b)
		I(i) = integral(f,a,b(i));
	end
end
