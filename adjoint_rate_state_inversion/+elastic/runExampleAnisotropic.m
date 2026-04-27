function runExampleAnisotropic(m,T,order,point_sources,check_ev_flag)

default_arg('m', [51, 51]);
default_arg('T', 10);
default_arg('order', 4);
default_arg('check_ev_flag',false)
default_arg('point_sources',[])

% Domain
% xlim = [-8,0,8];
% ylim = [3,0,-3];
xlim = [-1,1];
ylim = [1,0];
domain = multiblock.domain.Rectangle(xlim, ylim);
% ms = {m, m, m, m};
ms = {m};

% ----- Boundary conditions ------
f0 = @(t,x,y)0*x;
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
F0 = @(t,x,y) 0*x;
F1 = {F0, F0};
F2 = {F0, F0};
F = {F1, F2, F1, F2};
%--------------------------------

% --- Elastic stiffness tensor -----
c_diag = @(x,y) 0*x + 3;
c_off = @(x,y) 0*x + 1;

isotropic = true;
lambda = 1;
mu = 1;
d = @kroneckerDelta;

dim = 2;
C = cell(dim,dim,dim,dim);
for i = 1:dim
	for j = 1:dim
		for k = 1:dim
			for l = 1:dim
				if i==j && i==k && i==l
					C{i,j,k,l} = c_diag;
				else
					C{i,j,k,l} = c_off;
				end

				if isotropic
					C{i,j,k,l} = @(x,y) 0*x + lambda*d(i,j)*d(k,l) + ...
								+ mu*(d(i,k)*d(j,l) + d(i,l)*d(j,k));
				end

			end
		end
	end
end
% rho = {[],[],[],[]};
rho = @(x,y) 0*x + 1;
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
	point_sources.x = {[0.8*xlim(end), 0.5], [0, ylim(1)], [-0.5, 0.2], [-0.5, 0.5], [0,0.5] };
	t0 = 1;
	sigma = 0.2;
	A = 1;
	g = @(t) A*exp( -(t-t0).^2/sigma^2 );
	gx = g;
	gy = @(t) 0.5*g(t);
	point_sources.g = {{gx, gy}, {@(t) gx(t), gy}, {gx, gy}, {gx, gy}, {gx, gy}};
end

% Create discretization object
discr = elastic.elasticAnisotropicDiscrInterface(ms, order, xlim, ylim, C, rho, F, bc,...
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

function s = bcStruct(name, type, data)
	s = struct;
	s.boundary = name;
	s.type = type;
	s.data = data;
end