% Creates a 2D staggered grid of Lebedev (checkerboard) type
% Primal grid: equidistant with m points.
% Dual grid: m + 1 points, h/2 spacing first point.
% First grid line is "primal", 2nd is "dual", etc.
%
% Examples
%   g = grid.Lebedev2d(m, xlims, ylims)
%   g = grid.Lebedev2d([21, 31], {0,2}, {0,3})
function g = lebedev2d(m, xlims, ylims, opSet)

    default_arg('opSet', @(m,lim) sbp.D1StaggeredUpwind(m,lim,2));

    if ~iscell(xlims) || numel(xlims) ~= 2
        error('grid:lebedev2D:InvalidLimits','The limits should be cell arrays with 2 elements.');
    end

    if ~iscell(ylims) || numel(ylims) ~= 2
        error('grid:lebedev2D:InvalidLimits','The limits should be cell arrays with 2 elements.');
    end

    if xlims{1} > xlims{2}
        error('grid:lebedev2D:InvalidLimits','The elements of the limit must be increasing.');
    end

    if ylims{1} > ylims{2}
        error('grid:lebedev2D:InvalidLimits','The elements of the limit must be increasing.');
    end

    opsX = opSet(m(1), xlims);
    xp = opsX.x_primal;
    xd = opsX.x_dual;

    opsY = opSet(m(2), ylims);
    yp = opsY.x_primal;
    yd = opsY.x_dual;

    % 4 Cartesian grids with spacing h
    % 2 grids for displacements (u)
    % 2 grids for stresses (sigma)
    % Density needs to be evaluated on the u grids
    % The stiffness tensor is evaluated on the sigma grids

    gu1 = grid.Cartesian(xp, yp);
    gu2 = grid.Cartesian(xd, yd);
    gs1 = grid.Cartesian(xd, yp);
    gs2 = grid.Cartesian(xp, yd);

    gu = {gu1, gu2};
    gs = {gs1, gs2};

    dim = 2;
    g = grid.Staggered(dim, gu, gs);

end



