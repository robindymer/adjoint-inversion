function runExampleElasticAcousticCurve(m,T,order,check_ev_flag)

default_arg('check_ev_flag', false)
default_arg('order', 4)
default_arg('T', 10)
default_arg('m', 41);


% ---- Domain ----------------
% Acoustic disk inside elastic annulus
R0 = 0.5;
R1 = 1.5;
R2 = 4;
bc = struct;

% Acoustic
domain.acoustic.m = m;

% Circle
% domain.acoustic.def = multiblock.domain.Circle(R1);
% domain.acoustic.interfaceGroup = domain.acoustic.def.boundaryGroups.all;
% bc.acoustic = [];

% Annulus
domain.acoustic.def = Annulus(R0,R1);
domain.acoustic.interfaceGroup = domain.acoustic.def.boundaryGroups.outer;
bc_ac.boundary = domain.acoustic.def.boundaryGroups.inner;
bc_ac.type = 'd';
bc_ac.data = [];
bc.acoustic = {bc_ac};

% Elastic
domain.elastic.m = m;
domain.elastic.def = Annulus(R1,R2);
domain.elastic.interfaceGroup = domain.elastic.def.boundaryGroups.inner;

bc1.boundary = domain.elastic.def.boundaryGroups.outer;
bc1.type = {1, 'd'};
bc1.data = [];

bc2.boundary = domain.elastic.def.boundaryGroups.outer;
bc2.type = {2, 'd'};
bc2.data = [];

bc.elastic = {bc1, bc2};

g = struct;
g.acoustic = domain.acoustic.def.getGrid(m);
g.elastic = domain.elastic.def.getGrid(m);

domain.g = g;
% ------------------------------


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
rho = @(x,y) 0*x + 1;
% -------------------------------


% --- Acoustic material -------
rho_acoustic = @(x,y) 0*x + 1;
c = @(x,y) 0*x + 1;
%-------------------------------

% ---- Forcing ------------------
F.elastic = [];
F.acoustic = [];
%--------------------------------

% --- Material coefficients -----
parameters.elastic.C = C;
parameters.elastic.rho = rho;

parameters.acoustic.rho = rho_acoustic;
parameters.acoustic.c = c;
% -------------------------------

% --- Point sources -----------
pointSources = struct;
pointSources.acoustic = struct;
t0 = 1;
s = 0.2;
g_ac = @(t) 10*exp( -1/s^2 * (t-t0).^2 );
pointSources.acoustic.g = {g_ac, g_ac};
pointSources.acoustic.x = {[1, 0.5], [-1, -0.5]};
pointSources.acoustic.blockIds = [2, 4];
% -----------------------------

% Create discretization object
discr = elastic.elasticAcousticCurveDiscr(domain, order, parameters, bc, F, pointSources);
if check_ev_flag
	testElasticAcousticMatrix(discr);
end

% Initial data
u0 = @(x,y) [0*x; 0*x];
ut0 = @(x,y) [0*x; 0*x];
s = 0.1;
phi0 = @(x,y) 0*10*exp(-1/s^2*( (x-1).^2 + (y-0).^2 ));
phit0 = @(x,y) 0*x;

u0 = grid.evalOn(g.elastic, u0);
ut0 = grid.evalOn(g.elastic, ut0);
phi0 = grid.evalOn(g.acoustic, phi0);
phit0 = grid.evalOn(g.acoustic, phit0);

v0 = [u0; phi0];
vt0 = [ut0; phit0];
discr.v0 = v0;
discr.v0t = vt0;

% Solve
noname.animate(discr, [], T);