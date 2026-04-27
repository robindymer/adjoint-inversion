function runExample2Blocks(m,T,order,check_ev_flag)

default_arg('m', [51, 51]);
default_arg('T', 10);
default_arg('order', 4);
default_arg('check_ev_flag',false)

% ----- Boundary conditions ------
f0 = @(t,x,y)0*x;

% First component, both blocks
bcw1_1.boundary = {1 ,'w'};
bcw1_1.type = {1 ,'d'};
bcw1_1.data = f0;

bcn1_1.boundary = {1 ,'n'};
bcn1_1.type = {1 ,'d'};
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
bcn2_1.type = {1,'d'};
bcn2_1.data = f0;

% Second component, both blocks
bcw1_2.boundary = {1 ,'w'};
bcw1_2.type = {2 ,'d'};
bcw1_2.data = f0;

bcn1_2.boundary = {1 ,'n'};
bcn1_2.type = {2 ,'d'};
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
bcn2_2.type = {2,'d'};
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

xlim = [-1,0,1];
ylim = [1,0];
ms = {m, m};

% Create discretization object
discr = elastic.elasticDiscrInterface(m, order, xlim, ylim, lambda, mu, rho, F, bc);

if check_ev_flag
	testDiscrMatrix([],[],[],discr);
end

% Initial data
v0 = grid.evalOn(discr.grid, @v0_fun);
v0t = grid.evalOn(discr.grid, @v0t_fun);
discr.v0 = v0;
discr.v0t = v0t;

% Time-step
noname.animate(discr, [], T);

end

function v0 = v0_fun(x,y)
x0 = 0.5;
y0 = 0.5;
sigma = 0.05;
v0 = 5*[(exp(-((x-x0).^2+(y-y0).^2)/sigma^2)); ...
      0*x];
end

function v0t = v0t_fun(x,y)
v0t = [0*x; 0*x];
end

