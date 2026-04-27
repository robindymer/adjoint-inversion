function [ops, plotting] = setupOps(discr, faultBoundaryGroups, plotting)
	plotting_default = struct;
	plotting_default.faultBoundaryGroup = [];
	plotting_default.surfaceBoundaryGroup = [];
	default_struct('plotting', plotting_default);

	diffOp = discr.diffOp;

	faultType = {'t', 't'};
	[~, faultS1] = diffOp.boundary_condition(faultBoundaryGroups{1}, faultType);
	[~, faultS2] = diffOp.boundary_condition(faultBoundaryGroups{2}, faultType);

	% For mms data on normal components of traction, the following may be used:
	faultType = {'n', 't'};
	[~, faultNormalS1] = diffOp.boundary_condition(faultBoundaryGroups{1}, faultType);
	[~, faultNormalS2] = diffOp.boundary_condition(faultBoundaryGroups{2}, faultType);

	eTangentFault1 = diffOp.getBoundaryOperator('et', faultBoundaryGroups{1});
	eTangentFault2 = diffOp.getBoundaryOperator('et', faultBoundaryGroups{2});

	eNormalFault1 = diffOp.getBoundaryOperator('en', faultBoundaryGroups{1});
	eNormalFault2 = diffOp.getBoundaryOperator('en', faultBoundaryGroups{2});

	tau_n_1 = diffOp.getBoundaryOperator('tau_n', faultBoundaryGroups{1});
	tau_n_2 = diffOp.getBoundaryOperator('tau_n', faultBoundaryGroups{2});

	m = length(discr.D);
	[mFault, ~] = size(eTangentFault1');

	eu = speye(m, m);
	ev = speye(m, m);
	Eu = cell2mat({eu, 0*ev, sparse(m, mFault)})';
	Ev = cell2mat({0*eu, ev, sparse(m, mFault)})';
	Epsi = cell2mat({sparse(mFault, m), sparse(mFault, m), speye(mFault, mFault)})';

	ops = struct;
	ops.D = discr.D;

	ops.Eu = Eu;
	ops.Ev = Ev;
	ops.Epsi = Epsi;

	ops.SFault1 = faultS1;
	ops.SFault2 = faultS2;

	ops.SNormalFault1 = 1/2*faultNormalS1;
	ops.SNormalFault2 = 1/2*faultNormalS2;

	ops.eTangentFault1 = eTangentFault1;
	ops.eTangentFault2 = eTangentFault2;

	ops.eNormalFault1 = eNormalFault1;
	ops.eNormalFault2 = eNormalFault2;

	ops.tauNormalFault1 = tau_n_1;
	ops.tauNormalFault2 = tau_n_2;

	ops.m = m;
	ops.mFault = mFault;

	ops.xFault = discr.grid.getBoundary(faultBoundaryGroups{1});

	if ~isempty(plotting.surfaceBoundaryGroup)
		plotting.xSurface = discr.grid.getBoundary(plotting.surfaceBoundaryGroup);
	else
		plotting.xSurface = [];
	end
	if ~isempty(plotting.faultBoundaryGroup)
		plotting.xFault = discr.grid.getBoundary(plotting.faultBoundaryGroup);
	else
		plotting.xFault = ops.xFault;
	end

end