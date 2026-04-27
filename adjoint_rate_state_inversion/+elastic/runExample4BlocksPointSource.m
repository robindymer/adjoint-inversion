function runExample4BlocksPointSource(m,T,order,point_sources,check_ev_flag)

default_arg('m', [51, 51]);
default_arg('T', 10);
default_arg('order', 4);
default_arg('check_ev_flag',false)
default_arg('point_sources',[])

% Domain
xlim = [-8,0,8];
ylim = [3,0,-3];
domain = multiblock.domain.Rectangle(xlim, ylim);
ms = {m, m, m, m};

% ----- Boundary conditions ------
f0 = @(t,x,y)0*x;

domain = multiblock.domain.Rectangle(xlim, ylim);

bc = struct;
W = domain.boundaryGroups.W;
E = domain.boundaryGroups.E;
S = domain.boundaryGroups.S;
N = domain.boundaryGroups.N;

BCW1 = bcStruct(W, {1, 'd'}, f0);
BCW2 = bcStruct(W, {2, 'd'}, f0);
BCE1 = bcStruct(E, {1, 'd'}, f0);
BCE2 = bcStruct(E, {2, 'd'}, f0);
BCN1 = bcStruct(N, {1, 'd'}, f0);
BCN2 = bcStruct(N, {2, 'd'}, f0);
BCS1 = bcStruct(S, {1, 'd'}, f0);
BCS2 = bcStruct(S, {2, 'd'}, f0);

bc = {BCW1, BCW2, BCE1, BCE2, BCS1, BCS2, BCN1, BCN2};
% -------------------------------

% ---- Forcing ------------------
F0 = @(x,y,t) 0*x;
F1 = {F0, F0};
F2 = {F0, F0};
F = {F1, F2, F1, F2};
%--------------------------------

% --- Material coefficients -----
l1 = @(x,y) 0*x + 1;
l2 = @(x,y) 0*x + 1;
m1 = @(x,y) 0*x + 1;
m2 = @(x,y) 0*x + 1;

lambda = {l1, l2, l1, l2};
mu = {m1, m2, m1, m2};
rho = {[],[],[],[]};
% -------------------------------

% Zero initial data
u1 = @(x,y,t) 0*x;
u2 = u1;
u1t = u1;
u2t = u1;

% Initial data functions
u_exact = @(x,y,t) [u1(x,y,t); u2(x,y,t)];
ut_exact = @(x,y,t) [u1t(x,y,t); u2t(x,y,t)];

% Point source
if isempty(point_sources)
	point_sources = struct;
	point_sources.x = {[0.8*xlim(end), 0.5], [0, ylim(1)], [-7, 1], [-6, -1], [1,-1] };
	t0 = 1;
	sigma = 0.1;
	A = 40;
	g = @(t) A*exp( -(t-t0).^2/sigma^2 );
	gx = g;
	gy = @(t) 0.5*g(t);
	point_sources.g = {{gx, gy}, {@(t) gx(t), gy}, {gx, gy}, {gx, gy}, {gx, gy}};
end

% Create discretization object
discr = elastic.elasticDiscrInterface(m, order, xlim, ylim, lambda, mu, rho, F, bc,...
												 point_sources);

if check_ev_flag
	testDiscrMatrix([],[],[],discr);
end

% Initial data
discr.v0 = grid.evalOn(discr.grid, @(x,y) u_exact(x,y,0));
discr.v0t = grid.evalOn(discr.grid, @(x,y) ut_exact(x,y,0));

% Time-step
noname.animate(discr, [], T);

end

function s = bcStruct(name, type, data);
	s = struct;
	s.boundary = name;
	s.type = type;
	s.data = data;
end