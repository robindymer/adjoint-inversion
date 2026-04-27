% Creates a 1D staggered grid of dimension length(m).
% over the interval xlims
% Primal grid: equidistant with m points.
% Dual grid: m + 1 points, h/2 spacing first point.
% Examples
%   g = grid.primal_dual_1D(m, xlim)
%   g = grid.primal_dual_1D(11, {0,1})
function [g_primal, g_dual] = primalDual1D(m, xlims)

    if ~iscell(xlims) || numel(xlims) ~= 2
        error('grid:primalDual1D:InvalidLimits','The limits should be cell arrays with 2 elements.');
    end

    if xlims{1} > xlims{2}
        error('grid:primalDual1D:InvalidLimits','The elements of the limit must be increasing.');
    end

    xl = xlims{1};
    xr = xlims{2};
    h = (xr-xl)/(m-1);

    % Primal grid
    g_primal = grid.equidistant(m, xlims);
    g_primal.h = h;

    % Dual grid
    x = [xl; linspace(xl+h/2, xr-h/2, m-1)'; xr];
    g_dual = grid.Cartesian(x);
    g_dual.h = h;
end