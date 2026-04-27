function [e] = runExampleSuperGridAnisotropic(m,T,order,plot_flag,mSG,stretchingRatio)

default_arg('plot_flag',true)
default_arg('order',6)
default_arg('T',10)
default_arg('m',21);
default_arg('mSG', 101);
default_arg('stretchingRatio', 2^7);

% Material parameters
rho = @(x,y) 0*x+1;

c_diag = @(x,y,i,j,k,l) sin(i*x+j*y) + 1/2*sin(k*x-l*y) + 8;
c_off = @(x,y,i,j,k,l) 1/4*sin(i*x+j*y) + 1/8*sin(k*x-l*y) + 0;

dim = 2;
C = cell(dim,dim,dim,dim);
d = @kroneckerDelta;
for i = 1:dim
	for j = 1:dim
		for k = 1:dim
			for l = 1:dim
				if i==k && j==l
					C{i,j,k,l} = @(x,y) c_diag(x,y,i,j,k,l);
				else
					C{i,j,k,l} = @(x,y) c_off(x,y,i,j,k,l);
				end
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
g = def.getGrid(m);
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

SuperGrid.gamma = {0, 1e-1, 0};
SuperGrid.DOI_IDs = DOI_IDs;
SuperGrid.W_IDs = def.W_IDs;
SuperGrid.E_IDs = def.E_IDs;
SuperGrid.S_IDs = def.S_IDs;
SuperGrid.N_IDs = def.N_IDs;

% --- Point sources -----------
pointSources = struct;
t0 = 1;
s = 0.2;
g0 = @(t) 1*exp( -1/s^2 * (t-t0).^2 );
pointSources.g = {{g0, @(t)2*g0(t)}};
pointSources.x = {[0,0]};
pointSources.blockIds = [5];
% -----------------------------

% Create discretization object
discr = elastic.elasticAnisotropicSupergridDiscr(g, order, C, rho, F, bc, pointSources, [], SuperGrid);

% Initial data
A = 0;
sigma = 0.1*L;

f_gauss = @(x,y) A*exp( -1/sigma^2*((x-0).^2 + (y-0).^2) );
f_zero  = @(x,y) 0*x;

u0 = @(x,y) [f_gauss(x,y); f_gauss(x,y)];
u0t = @(x,y) [f_zero(x,y); f_zero(x,y)];

discr.v0 = grid.evalOn(discr.grid, u0);
discr.v0t = grid.evalOn(discr.grid, u0t);

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