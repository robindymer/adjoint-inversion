% A line of nRec equispaced receivers from point x0 to point x1

function [pars, truePars] = muInversionParsReceiverLine(opts)

% Default values
default_arg('opts', struct);

opts = defaultField(opts, 'm', [7, 7]);
opts = defaultField(opts, 'x0', [0, 0.5]);
opts = defaultField(opts, 'x1', [1, 0.5]);
opts = defaultField(opts, 'nRec', 100);

m = opts.m;
x0 = opts.x0;
x1 = opts.x1;
nRec = opts.nRec;

m = opts.m;
x0 = opts.x0;
x1 = opts.x1;
nRec = opts.nRec;

% Domain
dim = 2;
xlim = {0, 1};
ylim = {0, 1};
T = 3;

% Boundary conditions
bc = struct;
bc.names = {'w','e','n','s','w','e','n','s'};
bc.types = { {1, 't'}, {1, 't'}, {1, 't'}, {1, 't'},...
			 {2, 't'}, {2, 't'}, {2, 't'}, {2, 't'} };
bc.funcs = {[], [], [], [], [], [], [], []};

% Material parameters
lambda = @(x,y) 0*x + 1 + 1/4*cos(x).*sin(y);
mu = @(x,y) 0*x +1 + 1/4*cos(2*x).*sin(y);
rho = @(x,y) 0*x +1 + 1/4*cos(x).*sin(2*y);

% Sources
A = 30;
sigma = 0.2;
t0 = 1;
g = @(t) A*exp(-(t-t0).^2 / sigma^2 );
sources = struct;
sources.x = {[0.1, 0.1]};
sources.g = {{g, g}};

% Perturb mu
perturbedMu = @(x,y) 2*mu(x,y);

% Receiver line
receivers = struct;
x = linspace(x0(1),x1(1),nRec);
y = linspace(x0(2),x1(2),nRec);
receivers.x = cell(1, nRec);
for i = 1:nRec
	receivers.x{i} = [x(i), y(i)];
end

% Return struct
pars = struct;
pars.xlim = xlim;
pars.ylim = ylim;
pars.lambda = lambda;
pars.mu = perturbedMu;
pars.rho = rho;
pars.sources = sources;
pars.receivers = receivers;
pars.bc = bc;
pars.T = T;
pars.m = m;
pars.dim = dim;

truePars = pars;
truePars.mu = mu;