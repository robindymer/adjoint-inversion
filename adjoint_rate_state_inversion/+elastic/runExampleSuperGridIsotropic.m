function [e] = runExampleSuperGridIsotropic(m,T,order,plot_flag,mSG,stretchingRatio,stencil,checkEvFlag)

default_arg('stencil', 'narrow');
default_arg('plot_flag',true)
default_arg('checkEvFlag', false);
default_arg('order',6)
default_arg('T',10)
default_arg('m',21);
default_arg('mSG', 101);
default_arg('stretchingRatio', 2^7);

% Material parameters
rho = @(x,y) 0*x+1;
lambda = @(x,y) 0*x + 1;
mu = @(x,y) 0*x + 1;

dim = 2;
C = cell(dim,dim,dim,dim);
d = @kroneckerDelta;
for i = 1:dim
	for j = 1:dim
		for k = 1:dim
			for l = 1:dim
					C{i,j,k,l} = @(x,y) lambda(x,y)*d(i,j)*d(k,l) + ...
										mu(x,y)*(d(i,k)*d(j,l) + d(i,l)*d(j,k)) ;
			end
		end
	end
end

% Forcing
F = [];

% Domain
L = 1;
h = L/(m-1);

SG = struct;
SG.mDOI = m;
SG.mSG = mSG;
SG.ratio = stretchingRatio;

def = elastic.SuperGridSquare(L, SG);
% figure
% def.show([],true,m);
% pause
switch stencil
case 'staggered'
	g = def.getLebedevGrid(m);
otherwise
	g = def.getGrid(m);
end
DOI_IDs = def.DOI_IDs;

nDOIGrids = numel(DOI_IDs);
if nDOIGrids > 1
	g_DOI = g.grids(DOI_IDs);
	conn = g.connections(DOI_IDs, DOI_IDs);
	g_DOI = multiblock.Grid(g_DOI, conn);
else
	g_DOI = g.grids{DOI_IDs};
end

bc1.boundary = def.boundaryGroups.all;
bc1.type = {1,'t'};
bc1.data = [];

bc2.boundary = def.boundaryGroups.all;
bc2.type = {2,'t'};
bc2.data = [];

bc = {bc1, bc2};

% SuperGrid.gamma = {0, 5e-2, 0};
SuperGrid.gamma = [];
SuperGrid.DOI_IDs = DOI_IDs;
SuperGrid.W_IDs = def.W_IDs;
SuperGrid.E_IDs = def.E_IDs;
SuperGrid.S_IDs = def.S_IDs;
SuperGrid.N_IDs = def.N_IDs;

% Point source
% pointSource.x = {[0, 0]};
pointSource.x = {[0.5, 0.47]};
pointSource.blockIds = [5];

f0 = 4;
g0 = @(t) 1e1*rickerWavelet(t, 1/f0, f0);
pointSource.g = {{g0, @(t) 2*g0(t)}};

% Create discretization object
switch stencil
case {'narrow', 'upwind'}
	discr = elastic.elasticAnisotropicSupergridDiscr(g, order, C, rho, F, bc, pointSource, [], SuperGrid, [], stencil);
case 'staggered'
	discr = elastic.elasticAnisotropicStaggeredSupergridDiscr(g, order, C, rho, F, bc, pointSource, SuperGrid);
end

if checkEvFlag
	testDiscrMatrix([],[],[],discr);
end

% Initial data
A = 3*0;
sigma = 0.1*L;

f_gauss = @(x,y) A*exp( -1/sigma^2*((x-0).^2 + (y-0).^2) );
f_zero  = @(x,y) 0*x;

u0 = @(x,y) [f_gauss(x,y); f_gauss(x,y)];
u0t = @(x,y) [f_zero(x,y); f_zero(x,y)];

switch stencil
case 'staggered'
	discr.v0 = grid.evalOnStaggered(discr.grid, u0);
	discr.v0t = grid.evalOnStaggered(discr.grid, u0t);
otherwise
	discr.v0 = grid.evalOn(discr.grid, u0);
	discr.v0t = grid.evalOn(discr.grid, u0t);
end

% Solve
if plot_flag
	noname.animate(discr, [], T);
	e = [];
else
	[ts,N] = discr.getTimestepper([],T);
	fprintf('Computing o = %d, m = %d.  \n',order,m(1));
	ts.stepN(N,true);
	[u, ~] = ts.getV;
	e = [];

	% Compare
	% e = discr.compareSolutionsAnalytical(....
	% 			 u, u_exact, T );
end

% Return meshsize
% hx = (xr-xl)/(m(1)-1);
% hy = (yr-yl)/(m(2)-1);
% h = min(hx,hy);

end