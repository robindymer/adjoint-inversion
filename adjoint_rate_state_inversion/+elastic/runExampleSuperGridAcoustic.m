function runExampleSuperGridAcoustic(m,T,order,plot_flag,mSG,stretchingRatio)

default_arg('plot_flag',true)
default_arg('order',6)
default_arg('T',10)
default_arg('m',21);
default_arg('mSG', 101);
default_arg('stretchingRatio', 2^7);

% Material parameters
a = 1;
b = 1^2;
c = sqrt(a*b);

% Forcing
F = [];

% Domain
L = 1;
h = L/(m-1);

SG = struct;
SG.mDOI = m;
SG.mSG = mSG;
SG.ratio = stretchingRatio;

def = elastic.SuperGridChannel(L, SG);
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

% Boundary data, 'n' for Neumann, 'd' for Dirichlet.
bc.boundary = def.boundaryGroups.all;

bc.type = 'd';
bc.data = [];
bc = {bc};

% --- Point sources -----------
pointSources = struct;
t0 = 1;
s = 0.2;
g_ac = @(t) 1*exp( -1/s^2 * (t-t0).^2 );
pointSources.g = {g_ac};
pointSources.x = {[0,0]};
pointSources.blockIds = [2];
% -----------------------------

SuperGrid.gamma = {0, 1e-1, 0};
SuperGrid.DOI_IDs = DOI_IDs;
SuperGrid.W_IDs = def.W_IDs;
SuperGrid.E_IDs = def.E_IDs;
SuperGrid.S_IDs = def.S_IDs;
SuperGrid.N_IDs = def.N_IDs;

% Create discretization object
discr = elastic.acousticDiscrCurveSupergrid(g, order, a, b, F, bc, pointSources, [], SuperGrid);

% Initial data
A = 0;
sigma = 0.1*L;
u = @(x,y) A*exp( -1/sigma^2*((x-0).^2 + (y-0).^2) );
ut = @(x,y) 0*x;
discr.v0 = grid.evalOn(discr.grid, u);
discr.v0t = grid.evalOn(discr.grid, ut);

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
end

end