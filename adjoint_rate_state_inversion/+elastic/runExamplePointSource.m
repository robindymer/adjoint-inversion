function [e, h] = runExamplePointSource(m,T,order,point_sources,plot_flag,check_ev_flag,pars,domain,saveData)
default_arg('saveData',[]);
default_arg('domain',[]);
default_arg('pars',[]);
default_arg('check_ev_flag',false);
default_arg('plot_flag',true);
default_arg('point_sources',[]);
default_arg('order',4);
default_arg('T',10);
default_arg('m',[81,31]);

% Equation parameters
if isempty(pars)
	lambda = @(x,y) 2 + cos((x+y)/4);
	mu = @(x,y) 3 + sin((x+y)/8);
	rho = @(x,y) 1 + 1/2*cos((x+y)/4);
else
	lambda = pars.lambda;
	mu = pars.mu;
	rho = pars.rho;
end

% Zero initial data
u1 = @(x,y,t) 0*x;
u2 = u1;
u1t = u1;
u2t = u1;

% Forcing (not including point sources)
F = [];

% Domain
if isempty(domain)
	xl = -8;
	xr = 8;

	yl = -3;
	yr = 3;
else
	xl = domain.xl;
	xr = domain.xr;
	yl = domain.yl;
	yr = domain.yr;
end
xlim = {xl, xr};
ylim = {yl, yr};

% Point source
if isempty(point_sources)
	point_sources = struct;
	point_sources.x = {[0, -2], [0, yr]};
	t0 = 1;
	sigma = 0.1;
	A = 40;
	g = @(t) A*exp( -(t-t0).^2/sigma^2 );
	gx = g;
	gy = @(t) 0.5*g(t);
	point_sources.g = {{gx, gy}, {@(t) -gx(t), gy}};
end

% Boundary data, Mixed homogeneous BC. 't' for traction, 'd' for Dirichlet.
BD = struct;
BD.names = {'w','e','n','s','w','e','n','s'};
BD.types = { {1, 'd'}, {1, 'd'}, {1, 't'}, {1, 'd'},...
			 {2, 'd'}, {2, 'd'}, {2, 't'}, {2, 'd'} };
BD.funcs = {[], [], [], [], [], [], [], []};

% Create discretization object
discr = elastic.elasticDiscr(m, order, xlim, ylim, lambda, mu, rho, F, BD, point_sources);
if check_ev_flag
	testDiscrMatrix([],[],[],discr);
end

% Initial data functions
u_exact = @(x,y,t) [u1(x,y,t); u2(x,y,t)];
ut_exact = @(x,y,t) [u1t(x,y,t); u2t(x,y,t)];

% Initial data
discr.v0 = grid.evalOn(discr.grid, @(x,y) u_exact(x,y,0));
discr.v0t = grid.evalOn(discr.grid, @(x,y) ut_exact(x,y,0));

% Solve
if isempty(saveData)
	if plot_flag
		noname.animate(discr, [], T);
	else
		[ts,N] = discr.getTimestepper([],T);
		fprintf('Computing o = %d, m = %d.  \n',order,m(1));
		ts.stepN(N,true);
	end
else
	[ts,N] = discr.getTimestepper([],T);
	di = saveData.timeStepsPerFrame;
	j = 1;
	for i = 0:di:N
		[u, t] = ts.stepN(di,false);
		U(:,j+1) = u;
		tvec(j+1) = t;
		j = j+1;
	end

	filename = saveData.filename;
	save(filename,'U','tvec','discr');

end