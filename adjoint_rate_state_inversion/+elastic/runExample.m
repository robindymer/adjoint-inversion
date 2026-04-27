function runExample(m,T,order)

default_arg('m', [51, 51]);
default_arg('T', 10);
default_arg('order', 4);

% Create discretization object
BD.names = {'w','e','n','s','w','e','n','s'};
BD.types = { {1, 'd'}, {1, 't'}, {1, 'd'}, {1, 't'},...
			 {2, 'd'}, {2, 't'}, {2, 't'}, {2, 'd'} };
BD.funcs = {[], [], [], [], [], [], [], []};
discr = elastic.elasticDiscr(m, order, [], [], [], [], [], [], BD);

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

