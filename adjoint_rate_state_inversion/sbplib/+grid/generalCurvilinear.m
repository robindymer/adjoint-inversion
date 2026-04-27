% Creates a curvilinear grid of dimension length(m).
% over the logical domain xi_lim, eta_lim, using a provided
% coordinate mapping. The grid point distribution is
% specified is 
% Examples:
%   g = grid.generalCurvilinear(mapping, [mx, my], xlim, ylim, {'equidist'}, {'boundaryopt', order, 'accurate'})
function g = generalCurvilinear(mapping, m, varargin)
    n = length(m);

    % Check that parameters matches dimensions
    matchingParams = false;
    if length(varargin) == 2*n
            matchingParams = iscell([varargin{1:2*n}]);
    end
    assert(matchingParams,'grid:generalCurvilinear:NonMatchingParameters','The number of parameters per dimensions do not match.');
    
    X = [];
    h = [];
    inds_periodic = [];
    for i = 1:n
        lim = varargin{i};
        opts = varargin{i+n};
        gridtype = opts{1};
        switch gridtype
        case 'equidist'
            gridgenerator = @()util.get_grid(lim{1},lim{2},m(i));
        case 'boundaryopt'
            order = opts{2};
            stencil_type = opts{3};
            switch stencil_type
            case {'Accurate','accurate','A','acc'}
                gridgenerator = @()sbp.grid.accurateBoundaryOptimizedGrid(lim,m(i),order);
            case {'Minimal','minimal','M','min'}
                gridgenerator = @()sbp.grid.minimalBoundaryOptimizedGrid(lim,m(i),order);
            case {'accurate-upwind'}
                gridgenerator = @()sbp.grid.upwindAccurateBoundaryOptimizedGrid(lim,m(i),order);
            end
        case 'periodic'
            gridgenerator = @()util.get_periodic_grid(lim{1},lim{2},m(i));
            inds_periodic = [inds_periodic, i];
        otherwise
            error("grid type %s not supported. Must be one of 'equidist', 'boundaryopt', 'periodic'",gridtype);
        end
        try
            [X{i},h(i)] = gridgenerator();
        catch exception % Propagate any errors in the grid generation functions.
            msgText = getReport(exception);
            error('grid:boundaryOptimizedCurvilinear:InvalidParameter',msgText)
        end
    end
    g = grid.Curvilinear(mapping, X{:});
    g.logic.h = h;
    for i = inds_periodic
        g.logic.lim{i}{2} = g.logic.lim{i}{2}+h(i);
    end
end