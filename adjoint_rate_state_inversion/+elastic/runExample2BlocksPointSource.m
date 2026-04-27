function runExample2BlocksPointSource(m,T,order,point_sources,check_ev_flag)

default_arg('m', [51, 51]);
default_arg('T', 10);
default_arg('order', 4);
default_arg('check_ev_flag',false)
default_arg('point_sources',[])

% ----- Boundary conditions ------
f0 = @(t,x,y)0*x;

% First component, both blocks
bcw1_1.boundary = {1 ,'w'};
bcw1_1.type = {1 ,'d'};
bcw1_1.data = f0;

bcn1_1.boundary = {1 ,'n'};
bcn1_1.type = {1 ,'t'};
bcn1_1.data = f0;

bcs1_1.boundary = {1 ,'s'};
bcs1_1.type = {1, 'd'};
bcs1_1.data = f0;

bce2_1.boundary = {2 ,'e'};
bce2_1.type = {1 ,'d'};
bce2_1.data = f0;

bcs2_1.boundary = {2 ,'s'};
bcs2_1.type = {1, 'd'};
bcs2_1.data = f0;

bcn2_1.boundary = {2 ,'n'};
bcn2_1.type = {1,'t'};
bcn2_1.data = f0;

% Second component, both blocks
bcw1_2.boundary = {1 ,'w'};
bcw1_2.type = {2 ,'d'};
bcw1_2.data = f0;

bcn1_2.boundary = {1 ,'n'};
bcn1_2.type = {2 ,'t'};
bcn1_2.data = f0;

bcs1_2.boundary = {1 ,'s'};
bcs1_2.type = {2, 'd'};
bcs1_2.data = f0;

bce2_2.boundary = {2 ,'e'};
bce2_2.type = {2 ,'d'};
bce2_2.data = f0;

bcs2_2.boundary = {2 ,'s'};
bcs2_2.type = {2, 'd'};
bcs2_2.data = f0;

bcn2_2.boundary = {2 ,'n'};
bcn2_2.type = {2,'t'};
bcn2_2.data = f0;

bc = {bcw1_1, bcn1_1, bcs1_1, bce2_1, bcs2_1, bcn2_1,...
	  bcw1_2, bcn1_2, bcs1_2, bce2_2, bcs2_2, bcn2_2};
% -------------------------------

% ---- Forcing ------------------
F0 = @(x,y,t) 0*x;
F1 = {F0, F0};
F2 = {F0, F0};
F = {F1, F2};
%--------------------------------

% --- Material coefficients -----
l1 = @(x,y) 0*x + 1;
l2 = @(x,y) 0*x + 1;
m1 = @(x,y) 0*x + 1;
m2 = @(x,y) 0*x + 1;

lambda = {l1, l2};
mu = {m1, m2};
rho = {[],[]};
% -------------------------------

% Zero initial data
u1 = @(x,y,t) 0*x;
u2 = u1;
u1t = u1;
u2t = u1;

% Initial data functions
u_exact = @(x,y,t) [u1(x,y,t); u2(x,y,t)];
ut_exact = @(x,y,t) [u1t(x,y,t); u2t(x,y,t)];

xlim = [-8,0,8];
ylim = [3,-3];
ms = {m, m};

% Point source
if isempty(point_sources)
	point_sources = struct;
	point_sources.x = {[0.8*xlim(end), 0.5], [0, ylim(1)], [-7, 0]};
	t0 = 1;
	sigma = 0.1;
	A = 40;
	g = @(t) A*exp( -(t-t0).^2/sigma^2 );
	gx = g;
	gy = @(t) 0.5*g(t);
	point_sources.g = {{gx, gy}, {@(t) -gx(t), gy}, {gx, gy}};
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