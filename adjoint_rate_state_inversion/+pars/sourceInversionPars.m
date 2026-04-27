function [pars, truePars] = sourceInversionPars(opts)

% Default values
default_arg('opts', struct);

opts = defaultField(opts, 'm', [7, 7]);

m = opts.m;

% Domain
dim = 2;
xlim = {0, 1};
ylim = {0, 1};
T = 3;

% Boundary conditions
bc = struct;
bc.names = {'w','e','n','s','w','e','n','s'};
bc.types = { {1, 'd'}, {1, 'd'}, {1, 'd'}, {1, 'd'},...
			 {2, 'd'}, {2, 'd'}, {2, 'd'}, {2, 'd'} };
bc.funcs = {[], [], [], [], [], [], [], []};

% Material parameters
lambda = @(x,y) 0*x +1;
mu = @(x,y) 0*x +1;
rho = @(x,y) 0*x +1;

% Sources
A = 30;
sigma = 0.2;
t0 = 1;
g = @(t) A*exp(-(t-t0).^2 / sigma^2 );
sources = struct;
sources.x = {[0.1, 0.1]};
sources.g = {{g, g}};

% Perturb sources
perturbedSources = sources;
G = sources.g;
for i = 1:numel(G)
	for d = 1:dim
		% perturbedSources.g{i}{d} = @(t) G{i}{d}(t) + A/4*sin(d*t);
		perturbedSources.g{i}{d} = @(t) 2*G{i}{d}(t-1);
		% perturbedSources.g{i}{d} = @(t) 0.9*G{i}{d}(t);
	end
end

% Receivers
receivers = struct;
receivers.x = {[0.8, 0.8]};
% receivers.x = {[0.8, 0.8],[0.1,0.9], [0.2,0.2], [0.9,0.1]};
% receivers.x = {[0.8, 0.8],[0.1,0.9], [0.2,0.2], [0.9,0.1], [0.5, 0.5],...
% 			 	[0.4,0.6], [0.6,0.4], [0, 0], [0.25, 1], [1, 0.25]};

% Return struct
pars = struct;
pars.xlim = xlim;
pars.ylim = ylim;
pars.lambda = lambda;
pars.mu = mu;
pars.rho = rho;
pars.sources = perturbedSources;
pars.receivers = receivers;
pars.bc = bc;
pars.T = T;
pars.m = m;
pars.dim = dim;

truePars = pars;
truePars.sources= sources;

