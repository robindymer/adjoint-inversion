% A line of nRec equispaced receivers from point x0 to point x1

function [pars, truePars] = rho2Blocks(opts)

% Default values
default_arg('opts', struct);

opts = defaultField(opts, 'm', {[7, 7],[7,7]});
opts = defaultField(opts, 'x0', [-1, 0.5]);
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
xlim = [-1, 0, 1];
ylim = [1, 0];
T = 3;

domain = multiblock.domain.Rectangle(xlim, ylim);


% Boundary conditions
bc = struct;
W = domain.boundaryGroups.W;
E = domain.boundaryGroups.E;
S = domain.boundaryGroups.S;
N = domain.boundaryGroups.N;

BCW1 = bcStruct(W, {1, 't'}, []);
BCW2 = bcStruct(W, {2, 'd'}, []);
BCE1 = bcStruct(E, {1, 't'}, []);
BCE2 = bcStruct(E, {2, 'd'}, []);
BCN1 = bcStruct(N, {1, 't'}, []);
BCN2 = bcStruct(N, {2, 't'}, []);
BCS1 = bcStruct(S, {1, 'd'}, []);
BCS2 = bcStruct(S, {2, 'd'}, []);

bc = {BCW1, BCW2, BCE1, BCE2, BCS1, BCS2, BCN1, BCN2};

% Material parameters
% lambda = @(x,y) 0*x +1;
% mu = @(x,y) 0*x +1;
% rho = @(x,y) 0*x +1;

lambda = @(x,y) 2 + sin(x+y);
mu = @(x,y) 2 + cos(x+y);
rho = @(x,y) 2 + 1/2*sin(x-y);

% Sources
A = 30;
sigma = 0.2;
t0 = 1;
g = @(t) A*exp(-(t-t0).^2 / sigma^2 );
sources = struct;
sources.x = {[0.1, 0.1]};
sources.g = {{g, g}};

% Perturb rho
perturbedRho = {@(x,y) 2*rho(x,y), @(x,y) 1.5*rho(x,y)} ;

% Put other parameters in cell arrays
lambda = {lambda, lambda};
mu = {mu, mu};
rho = {rho, rho};

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
pars.mu = mu;
pars.rho = perturbedRho;
pars.sources = sources;
pars.receivers = receivers;
pars.bc = bc;
pars.T = T;
pars.m = m;
pars.dim = dim;

truePars = pars;
truePars.rho = rho;

end

function s = bcStruct(name, type, data);
	s = struct;
	s.boundary = name;
	s.type = type;
	s.data = data;
end